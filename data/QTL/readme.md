### Data description
DNase QTLs were downloaded from [Degner et al., 2012](https://doi:10.1038/nature10808)). Histone QTLs were downloaded from [Grubert et al., 2015](http://dx.doi.org/10.1016/j.cell.2015.07.048). 

100-bp sequences were extracted centering on the QTL SNPs using bedtools "getfasta". DNase QTLs are based on the hg18 genome, while histone QTLs are based on the hg19 genome. 

\[QTL\]\_high.fa: alleles associated with higher traits

\[QTL\]\_low.fa: the other alleles associated with lower traits

The example command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
./data/QTL/DNase/dsQTL_high.fa \
./data/QTL/DNase/dsQTL_low.fa \
-o ./data/PU1/
```
