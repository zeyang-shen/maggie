### Data description
Strain-specific activated or repressed regulatory elements in response to KLA were identified based on ATAC-seq and H3K27ac ChIP-seq data of macrophages for the four genetically diverse strains of mice: C57BL/6J (C57), NOD/ShiLtJ (NOD), PWK/PhJ (PWK), and SPRET/EiJ (SPRET) ([Link et al., 2018](https://doi.org/10.1016/j.cell.2018.04.018)). 

./Activated/input_seq/{X}_EnhAct_spec_{Y}_\[ref/mut\].fa: activated elements specific to X when comparing X and Y; "ref" indicates sequences from strain X, and "mut" from strain Y

./Repressed/input_seq/{X}_EnhRep_spec_{Y}_\[ref/mut\].fa: repressed elements specific to X when comparing X and Y; "ref" indicates sequences from the genome of strain X, and "mut" from strain Y

Sequences were extracted using MMARGE "extract_sequences" ([Link et al., 2018](https://doi.org/10.1093/nar/gky491)).

The command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
$(ls -1 ./data/KLA_response/Activated/input_seq/*_ref.fa | paste -sd ",") \
$(ls -1 ./data/KLA_response/Activated/input_seq/*_mut.fa | paste -sd ",") \
-o ./data/KLA_response/Activated/
```

It will concatenate all the positive sequences ('ref') and all the negative sequences ('mut') before doing MAGGIE analysis. 