### Data description
QTLs of lymphoblastoid cell lines were downloaded from [Degner et al., 2012](https://doi:10.1038/nature10808) for DNase and from [Grubert et al., 2015](http://dx.doi.org/10.1016/j.cell.2015.07.048) for several different histone modifications. 

For each type of QTLs, we provide the following files:

\[QTL\]\.vcf: information about QTL in VCF file, including chromosome, position, variant ID, reference and alternative allele, and effect size. DNase QTLs are based on the hg18 genome, while histone QTLs are based on the hg19 genome. 

\[QTL\]\_high.fa: 100-bp sequences centering on the alleles associated with higher traits

\[QTL\]\_low.fa: 100-bp sequences centering on the alleles associated with lower traits

The example command line to process these files is:
```bash
python ./bin/maggie_vcf_input.py \
./data/QTL/DNase/dsQTL.vcf \
hg18 \
-e 6 \
-o ./data/QTL/DNase/maggie_output/
```
or
```bash
python ./bin/maggie_fasta_input.py \
./data/QTL/DNase/dsQTL_high.fa \
./data/QTL/DNase/dsQTL_low.fa \
-o ./data/QTL/DNase/maggie_output/
```
