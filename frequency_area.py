# Script for calculating joint frequency-area distribution from a VCF
# Cian Williams cdw45@cam.ac.uk

import numpy as np
import pandas as pd
import allel
import zarr
import argparse
import re
from scipy.spatial import ConvexHull, distance
from multiprocessing import Pool
from pathlib import Path
from pyproj import Transformer


"""
Compute frequency-area distribution for SNPs in a VCF file.
 
Input:
    --vcf       path to input VCF file
    --samples   path to sample metadata CSV (columns: sample_id, longitude, latitude, cohort, ploidy)
    --out       path to output TSV
    --zarr      path to zarr store (will be created if it doesn't exist)
    --chrom     optional chromosome name to filter to
    --chunk     chunk size for processing (default: 10000)
    --threads   number of threads for parallel processing (default: 1)
 
Output:
    TSV with columns: position, ref, alt, mean_daf, area
"""

# ── CLI ──────────────────────────────────────────────────────────────────────
 
def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute SNP frequency-area distributions from a VCF."
    )
    parser.add_argument("--vcf",     required=True,  help="Path to input VCF file")
    parser.add_argument("--samples", required=True,  help="Path to sample metadata CSV")
    parser.add_argument("--out",     required=True,  help="Path to output TSV")
    parser.add_argument("--zarr",    required=True,  help="Path to zarr store (created if absent)")
    parser.add_argument("--chrom",   default=None,   help="Chromosome to filter to (optional)")
    parser.add_argument("--chunk",   type=int, default=10_000, help="Chunk size (default: 10000)")
    parser.add_argument("--threads", type=int, default=1,      help="Parallel threads (default: 1)")
    return parser.parse_args()

# ── VCF → ZARR ───────────────────────────────────────────────────────────────
 
def load_or_convert_zarr(vcf_path: str, zarr_path: str) -> zarr.hierarchy.Group:
    """Convert VCF to Zarr if not already done, then open it."""
    if not Path(zarr_path).exists():                                        # If it doesn't already exist
        print(f"Converting VCF to Zarr at {zarr_path} ...")
        allel.vcf_to_zarr(
            vcf_path, zarr_path,
            fields=["calldata/GT", "variants/POS", "variants/REF",
                    "variants/ALT", "variants/CHROM", "variants/DA"],       # DA is the derived allele
            overwrite=True                                                  # Probably redundant but a good safety check
        )
        print("Conversion complete.")
    else:                                                                    # If it does already exist
        print(f"Using existing Zarr store at {zarr_path}")
    return zarr.open(zarr_path, mode="r")                                    # Return the opened zarr store


# ── ALLELE FREQUENCIES ────────────────────────────────────────────────────────
 
def cohort_mean_af(genotype_chunk: allel.GenotypeArray,
                   samples: pd.DataFrame,
                   da: np.ndarray) -> np.ndarray:
    """
    Compute mean derived allele frequency across cohorts,
    accounting for mixed ploidy on sex chromosomes.
    Samples w/ ploidy=1 are pseudo-diploid encoded so only the first
    haplotype is used to avoid double-counting.

    Returns array of shape [n_variants, n_derived_alleles] with mean DAF per derived allele.
    To account for uneven sampling between populations (cohorts) we calculate frequency for each population separately
    and then take the mean
    """
    cohorts = samples["cohort"].unique()
    cohort_freqs = []

    for cohort in cohorts:                                                                  # Loop through each cohort
        cohort_samples = samples[samples["cohort"] == cohort]                               # Select only samples in that cohort

        diploid_idx = cohort_samples[cohort_samples["ploidy"] == 2].index.to_numpy()        # Get np arrays with indices of diploid individuals
        haploid_idx = cohort_samples[cohort_samples["ploidy"] == 1].index.to_numpy()        # Haploids

        # guard: skip cohorts with no samples in this chunk's sample list (eg if there are samples that aren't actually in our VCF)
        # Or if we somehow chunked the VCF by samples (don't do this)
        diploid_idx = diploid_idx[diploid_idx < genotype_chunk.shape[1]]
        haploid_idx = haploid_idx[haploid_idx < genotype_chunk.shape[1]]

        if len(diploid_idx) == 0 and len(haploid_idx) == 0:
            continue                                                                        # Skip this cohort

        # Calculate allele frequencies
        ac_diploid = (genotype_chunk[:, diploid_idx, :].count_alleles(max_allele=3)         # All variants, diploid individuals, all ploidy positions
                      if len(diploid_idx) > 0                                               # count_alleles will return an array of shape n_variants, 4
                      else np.zeros((genotype_chunk.shape[0], 4), dtype=int))

        ac_haploid = (allel.GenotypeArray(genotype_chunk[:, haploid_idx, 0:1]).count_alleles(max_allele=3) # the 0:1 is critical here - we keep only the first haplotype (which is identical to the second haplotype)
                      if len(haploid_idx) > 0
                      else np.zeros((genotype_chunk.shape[0], 4), dtype=int))

        ac_combined = ac_diploid + ac_haploid                                               # Total counts for all alleles (n_variants, 4)
        total = ac_combined.sum(axis=1, keepdims=True).astype(float)                        # Total alleles observed at each site (n_variants, 1)

        # index into ac using DA to get derived allele count for each variant
        derived_counts = ac_combined[np.arange(len(da)), da].reshape(-1, 1)                 # Grab counts at indices i, da[i] and reshape to (n_variants, 1) (if your site is biallelic)
                                                                                            # I flippin love numpy. so elegant
        with np.errstate(invalid="ignore"):                                                 # Ignore divide by 0 warnings (since np.where will compute derived/total before evaluating the condition)
            freq = np.where(total > 0, derived_counts / total, np.nan)
        cohort_freqs.append(freq)

    if not cohort_freqs:
        n_variants = genotype_chunk.shape[0]
        return np.full((n_variants, 1), np.nan)                                              # If no valid cohorts were found, return an array of nan (n_variants, 1)

    return np.nanmean(np.stack(cohort_freqs, axis=0), axis=0)                                # Get mean among cohorts

# ── SPATIAL SPREAD ────────────────────────────────────────────────────────────
 
def compute_area(carrier_locs: np.ndarray) -> float:                                         # Takes a (n_carriers, 2) array of locations passed from process_chunk()
    """
    Compute spatial spread of allele carriers.
    - 1 location  → NaN
    - 2 locations → Euclidean distance
    - 3+ locations → Convex hull area
    """
    locs = np.unique(carrier_locs, axis=0)                                                      # Get unique locations
    n = len(locs)
    if n > 2:
        try:
            transformer = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)        # Reproject from lon/lat to equal area (in m). always_xy means that the order stays longitude, latitude
            x, y = transformer.transform(locs[:, 0], locs[:, 1]) 
            area = ConvexHull(np.column_stack([x, y])).volume                                   # Area in m^2

            return area
        except Exception:                                                                       # If the points are all colinear
            return np.nan
    elif n == 2:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:6933") # WGS84 → equal-area
        x, y = transformer.transform(locs[:, 1], locs[:, 0])  # note: lat, lon order
        locs_adjusted = np.column_stack([x, y])
        #return distance.euclidean(locs_adjusted[0], locs_adjusted[1])                            # For 2 points Euclidean distance also carries info but I'm not sure how to integrate it with area so have left it for now
        return np.nan
    return np.nan

# ── CHUNK PROCESSING ──────────────────────────────────────────────────────────
 
def process_chunk(args_tuple):
    """
    Process one chunk of variants.
    Returns a list of dicts, one per (variant, alt allele) pair.
    """
    gt_chunk, pos_chunk, ref_chunk, alt_chunk, da_chunk, samples, mean_af = args_tuple
 
    coords = samples[["longitude", "latitude"]].to_numpy()
    records = []
 
    # precompute carrier masks for all alt alleles at once:
    # shape [n_alt_alleles, n_variants, n_samples]
    n_variants, n_samples, _ = gt_chunk.shape                               # gt_chunk has shape (n_variants, n_samples, ploidy)
    n_alts = alt_chunk.shape[1]                                             # This will just be 3. We don't use it though
    
    ploidy = samples["ploidy"].to_numpy()
    haploid_mask = ploidy == 1                                              # Identify haploid individuals

    for i in range(n_variants):                                             # For each variant
        gt = gt_chunk[i]                                                    # shape (n_samples, ploidy)
 
        af = mean_af[i, 0]                                                  # Arbitrarily pick the first variant (just a safeguard as we should have only biallelic SNPs)
        if np.isnan(af) or af == 0:                                         # Skip the variant if there are no derived alleles
            continue

        allele_code = int(da_chunk[i])                                      # Get index (0 or 1) of the derived allele for this variant

        # Work out which samples actually carry the derived allele
        carrier_mask = (gt == allele_code).any(axis=1)
        carrier_mask[haploid_mask] = (gt[haploid_mask, 0] == allele_code)   # Not strictly necessary with pseudo-diploids but doesn't hurt
        carriers = np.where(carrier_mask)[0]                                # [0] is needed since np.where returns a tuple

        area = np.nan
        if len(carriers) > 0:                                               # If there are carriers of the derived allele
            area = compute_area(coords[carriers])                           # Get the minimum convex hull for carriers

        records.append({
            "position": int(pos_chunk[i]),
            "ref":      ref_chunk[i],
            "alt":      alt_chunk[i],
            "derived":  allele_code,
            "mean_daf": float(af),
            "area":     area,
        })
 
    return records

# ── Convert derived allele to integer index ──────────────────────────────────────────────────────────────────────
def da_base_to_index(da_bases, ref, alt):
    """
    Convert derived allele bases to allele indices.
    e.g. if ref='A', alt=['G', '', ''] and da='G', returns 1
    e.g. if ref='A', alt=['G', '', ''] and da='A', returns 0
    """
    indices = np.zeros(len(da_bases), dtype=int)                                # This will store the index of the derived allele
    for i, (base, r, alts) in enumerate(zip(da_bases, ref, alt)):               # Takes the ith element of da_bases, ref, and alt respectively
        if isinstance(base, bytes):                                             # Decode zarr bytes into python string
            base = base.decode()
        if isinstance(r, bytes):
            r = r.decode()
        alts = [a.decode() if isinstance(a, bytes) else a for a in alts]

        if base == r:                                                           # If derived = reference
            indices[i] = 0
        else:
            for j, a in enumerate(alts):
                if base == a:
                    indices[i] = j + 1                                          # +1 is necessary due to python's 0-indexing
                    break                                                       # Stop at the first match (not strictly necessary since this script expects biallelic sites)
    return indices

# ── MAIN ──────────────────────────────────────────────────────────────────────
 
def main():
    args = parse_args()
 
    # load sample metadata
    samples = pd.read_csv(args.samples)
    required_cols = {"sample_id", "longitude", "latitude", "cohort", "sex"}             # Check all the columns we need are there
    if not required_cols.issubset(samples.columns):
        raise ValueError(f"Metadata CSV must contain columns: {required_cols}")
    samples = samples.reset_index(drop=True)                                            # Reset all row indices (in case of messy/filtered df)
 
    # load / convert zarr
    callset = load_or_convert_zarr(args.vcf, args.zarr)                                 # Convert VCF to zarr
 
    positions = callset["variants/POS"][:]
    ref       = callset["variants/REF"][:]
    alt       = callset["variants/ALT"][:]
    chroms    = callset["variants/CHROM"][:]
    da_bases        = callset["variants/DA"][:]
    da = da_base_to_index(da_bases, ref, alt)                                           # Convert derived allele strings to indices (0=ref 1=alt)
    genotype  = allel.GenotypeArray(callset["calldata/GT"])
 
    # optional chromosome filter
    if args.chrom:
        mask = chroms == args.chrom.encode() if isinstance(chroms[0], bytes) else chroms == args.chrom      # Get a mask for all the positions that are on our selected chromosome
        positions = positions[mask]
        ref       = ref[mask]
        alt       = alt[mask]
        genotype  = genotype[mask]
        da        = da[mask]
        print(f"Filtered to chromosome {args.chrom}: {len(positions)} variants")
    else:
        print(f"Processing all chromosomes: {len(positions)} variants")
 
    # build chunks for parallelisation
    chunk_size = args.chunk
    n_variants = len(positions)
    chunk_args = []
 
    for start in range(0, n_variants, chunk_size):                                                  # Looping over variants in steps of chunk_size
        end = min(start + chunk_size, n_variants)
        gt_chunk  = genotype[start:end]
        pos_chunk = positions[start:end]
        ref_chunk = ref[start:end]
        alt_chunk = alt[start:end]
        da_chunk  = da[start:end]
        maf       = cohort_mean_af(gt_chunk, samples, da_chunk)                                     # Compute the derived allele frequencies for that chunk (chunk_size, 1)
        chunk_args.append((gt_chunk, pos_chunk, ref_chunk, alt_chunk, da_chunk, samples, maf))      # List of tuples of args for each chunk
 
    print(f"Processing {len(chunk_args)} chunks of up to {chunk_size} variants ...")
 
    # run chunks (parallel or serial)
    if args.threads > 1:
        with Pool(processes=args.threads) as pool:
            results = pool.map(process_chunk, chunk_args)
    else:
        results = [process_chunk(c) for c in chunk_args]
 
    # flatten and save
    all_records = [rec for chunk in results for rec in chunk]                                               # Convoluted but necessary since results is a list of lists (each entry is the output from 1 chunk)
    df = pd.DataFrame(all_records, columns=["position", "ref", "alt", "derived", "mean_daf", "area"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Done. {len(df)} records written to {args.out}")
 
 
if __name__ == "__main__":
    main()