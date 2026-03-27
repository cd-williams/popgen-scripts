# popgen-scripts
Reusable scripts for miscellaneous popgen analyses

### polarise_vcf.py
Script that takes a biallelic SNP vcf with individuals from 3 species: A (focal), B and C (outgroups).

Performs parsimony-based polarisation of ancestral/derived alleles in A. I.e. if an allele is fixed in B and C 
but polymorphic in A then the B-C allele is assigned as derived. The script adds 2 new INFO fields to each site in
the VCF, AA and DA.

Arguments:
- `--vcf` path to input VCF file
- `--samples` TSV where column 1 is sample names and column 2 is the species (A, B or C)
- `--out` path to output VCF file
- `--thresh` proportion of missing calls to tolerate in species B or C 


### freqyency_area.py
Script that takes a polarised biallelic SNP vcf for georeferenced samples and computes the joint frequency-area spectrum as described in [Rehmann et al 2025](https://academic.oup.com/mbe/article/42/6/msaf141/8169206). 

Arguments:
- `--vcf` path to input VCF file
- `--samples` sample metadata CSV. Columns must include sample_id, longitude, latitude, cohort, ploidy. "Cohort" is required if you want to account for uneven spatial sampling. If you're not bothered just fill this column with a single value
- `--out` path to output CSV file
- `--zarr` path to zarr store (will be created if it doesn't exist). Zarr is a super efficient structure for storing blocked arrays (ideal for genomic data) and facilitates parallelisation
- `--chrom` (optional) chromosome name to filter to
- `--chunk` parallelisation involves splitting the VCF into chunks, this controls the chunk size (default 10000)
- `--threads` number of threads for parallel processing (default 1)
