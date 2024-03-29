### v1.2 - 05/11/2022
    Fixed installation issue by adapting scripts to Python 3.9 and Biopython 1.79
### v1.1.1 - 04/25/2021
    Updated find_differential.py script to speed up the finding of differential regions by integrating bedtools
### v1.1 - 04/13/2021
    Fixed the installation error caused by outdated logomaker package dependency (https://pypi.org/project/logomaker/)
    Fixed bugs in concatenating files and maggie_vcf_input.py
### v1.0.1 - 11/22/2020
    Fixed bugs that will not terminate parallele processes when computing motif scores
### v1.0 - 10/15/2020
    Official release of MAGGIE v1.0!
### v0.3.4 - 06/12/2020
    Add script to find differential regions based on two BED files
### v0.3.3 - 05/10/2020
    Add function to split the input variants given in VCF file based on genomic annotations in GFF3 file
### v0.3.2 - 05/03/2020
    Fix bugs for VCF input; add a script to generate and conduct analysis on simulated sequences; upload histone QTL files
### v0.3.1 - 04/26/2020
    Update modules to extract sequences from reference genomes; add modules to generate simulated sequences
### v0.3.0 - 03/21/2020
    Now allow VCF file as input besides FASTA files!
### v0.2.2 - 01/12/2020
    Add information showing in the outputs, including numbers of motif mutation, summary and distribution of motif score differences
### v0.2.1 - 01/11/2020
    Incorporate redundant motifs from JASPAR 2020 database; \
    Show logo of the top motif of the merged list in HTML output; \
    Upload allele-specific CTCF, strain-specific PU.1, and QTLs datasets; \
    Update environment configuration; \
    Add progress bar.
### v0.2.0 - 09/30/2019
    Update computation pipeline and algorithm. Now accept FASTA files as inputs
### v0.1.0 - 01/08/2019
    Initial release.
