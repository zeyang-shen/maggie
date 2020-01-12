Data description: CTCF allele-specific binding sites were downloaded from [Shi et al., 2016](https://doi.org/10.1093/nar/gkw691) which were identified for GM12878 cells. 

CTCF_binding_alleles.fa: sequences associated with CTCF binding

CTCF_nonbinding_alleles.fa: the other alleles not associated with CTCF binding; come as pairs with sequences in "CTCF_binding_alleles.fa"

Sequences were extracted from GRCh37 genome using bedtools "getfasta".

The command line to process these files is:
```bash
python ./bin/maggie_computation_fasta_input.py \
./data/ASB/CTCF_binding_alleles.fa \
./data/ASB/CTCF_nonbinding_alleles.fa \
-o ./data/ASB/
```
