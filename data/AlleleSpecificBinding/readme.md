### Data description
CTCF allele-specific binding sites identified in GM12878 cells were downloaded from [Shi et al., 2016](https://doi.org/10.1093/nar/gkw691) and were then extracted from GRCh37 genome with a fixed size of 100 base pairs. 

CTCF_binding_alleles.fa: sequences associated with CTCF binding

CTCF_nonbinding_alleles.fa: the other alleles not associated with CTCF binding; come as pairs with sequences in "CTCF_binding_alleles.fa"

The command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
./data/AlleleSpecificBinding/CTCF_binding_alleles.fa \
./data/AlleleSpecificBinding/CTCF_nonbinding_alleles.fa \
-o ./data/AlleleSpecificBinding/maggie_output/ \
-p 8
```
