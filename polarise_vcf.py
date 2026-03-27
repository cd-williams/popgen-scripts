# Script to polarise a VCF containing samples from one species (A) and 2 outgroups (B and C)
# Cian Williams cdw45@cam.ac.uk

import argparse
import cyvcf2
from collections import defaultdict

"""
Polarise SNPs in a VCF file for one species (A) and 2 outgroups (B and C).
NOTE this method is quite naive, it just finds sites that are fixed in B and C and polymorphic in A
 
Input:
    --vcf       path to input VCF file
    --samples   path to sample metadata TSV (columns: sample_id, species (A, B or C))
    --out       path to output VCF file
    --thresh    proportion of individuals in B or C that can have missing calls
 
Output:
    VCF with 2 new fields in the INFO column (AA and DA)
"""

# ── CLI ──────────────────────────────────────────────────────────────────────
 
def parse_args():
    parser = argparse.ArgumentParser(
        description="Polarise alleles using parsimony"
    )
    parser.add_argument("--vcf",     required=True,  help="Path to input VCF file")
    parser.add_argument("--samples", required=True,  help="Path to sample metadata TSV")
    parser.add_argument("--out",     required=True,  help="Path to output VCF file")
    parser.add_argument("--thresh",     required=True,  help="Proportion of individuals that can have missing calls in B or C")
    return parser.parse_args()


# ── Function to get alleles that are fixed in B and C ──────────────────────────────────────────────────────────────────────
def get_species_allele(variant, alleles, indices, missing_thresh):
    """
    Returns the fixed allele for a species if fixed, else None.
    missing_thresh: maximum tolerated proportion of missing genotypes.
    """
    observed = []
    for idx in indices:
        gt = variant.genotypes[idx] # returns a list of [allele1, allele2, phased (1/0)]
        for a in gt[:-1]:           # strip phasing flag
            if a != -1:             # -1 denotes a missing call ('.' in the original VCF)
                observed.append(a)

    if not observed:
        return None                 # all missing

    # check missing data threshold
    n_possible = len(indices) * 2   # assumes diploidy
    if len(observed) / n_possible < (1 - missing_thresh):
        return None

    allele_set = set(observed)
    if len(allele_set) == 1:        # If the species is fixed for a single allele
        return alleles[list(allele_set)[0]]
    
    return None  # polymorphic within species

# ── Main ──────────────────────────────────────────────────────────────────────
args = parse_args()

thresh = float(args.thresh)

# Assign samples to species
sample_map = defaultdict(list)

with open(args.samples) as f:
    for line in f:
        sample, species = line.strip().split()
        sample_map[species].append(sample)


# Extract sample indices from VCF
vcf = cyvcf2.VCF(args.vcf)
samples = vcf.samples

b_indices = [samples.index(s) for s in sample_map["B"]]
c_indices = [samples.index(s) for s in sample_map["C"]]

# Add fields to VCF header
vcf.add_info_to_header({'ID': 'AA', 'Number': '1', 'Type': 'String',
                        'Description': 'Ancestral allele (B+C consensus)'})
vcf.add_info_to_header({'ID': 'DA', 'Number': '1', 'Type': 'String',
                        'Description': 'Derived allele (B+C consensus)'})


# Assign ancestral/derived in a new VCF
w = cyvcf2.Writer(args.out, vcf, mode="wz")
# Keep track of polarisation to print to the screen
polarised = 0
unpolarised = 0

for variant in vcf:
    alleles = [variant.REF] + variant.ALT              # Adding 2 lists makes a list [REF, ALT1, ALT2]

    b_allele = get_species_allele(variant, alleles, b_indices, thresh)
    c_allele = get_species_allele(variant, alleles, c_indices, thresh)

    if b_allele and c_allele and b_allele == c_allele: # both species are fixed for an identical allele
        aa = b_allele
        da = [a for a in alleles if a != aa]
        variant.INFO["AA"] = aa
        variant.INFO["DA"] = da[0] if da else "."      # Arbitrarily pick the first one (doesn't matter since VCF should be biallelic only)
        polarised += 1
    else:
        variant.INFO["AA"] = "."
        variant.INFO["DA"] = "."
        unpolarised += 1


    w.write_record(variant)
    
    # Print progress every 10k sites
    if (polarised+unpolarised) % 10000 == 0:
        print(f"Processed {polarised+unpolarised} sites, {polarised} polarised")

w.close()

# Report how many variants were polarised
print(f"Was able to polarise {polarised} out of {polarised+unpolarised} sites")