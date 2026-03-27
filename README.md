# popgen-scripts
Reusable scripts for miscellaneous popgen analyses

### polarise_vcf.py
Script that takes a biallelic SNP vcf with individuals from 3 species: A (focal), B and C (outgroups).

Performs parsimony-based polarisation of ancestral/derived alleles in A (ie if an allele is fixed in B and C) 
but polymorphic in A then the B-C allele is assigned as derived. The script adds 2 new INFO fields to each site in
the VCF, AA and DA.

Arguments:
- `--vcf` path to input VCF file
- `--samples` TSV where column 1 is sample names and column 2 is the species (A, B or C)
- `--out` path to output VCF file
- `--thresh` proportion of missing calls to tolerate in species B or C 
