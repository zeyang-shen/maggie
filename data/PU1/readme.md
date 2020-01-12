### Data description
PU.1 strain-specific binding sites were identified based on PU.1 ChIP-seq data of macrophages for the two genetically diverse strains of mice: C57BL/6J (C57) and BALB/cJ (BALB) ([Link et al., 2018](https://doi.org/10.1016/j.cell.2018.04.018)). 

BALB_PU1_SpecificBinding_BALB_seqs.fa: BALB-specific PU.1 binding sites, sequences from BALB genome

BALB_PU1_SpecificBinding_C57_seqs.fa: pair with "BALB_PU1_SpecificBinding_BALB_seqs.fa", sequences from C57 genome

C57_PU1_SpecificBinding_C57_seqs.fa: C57-specific PU.1 binding sites, sequences from C57 genome

C57_PU1_SpecificBinding_BALB_seqs.fa: pair with "C57_PU1_SpecificBinding_C57_seqs.fa", sequences from BALB genome

Sequences were extracted using MMARGE "extract_sequences" ([Link et al., 2018](https://doi.org/10.1093/nar/gky491)).

The command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
./data/PU1/BALB_PU1_SpecificBinding_BALB_seqs.fa,./data/PU1/C57_PU1_SpecificBinding_C57_seqs.fa \
./data/PU1/BALB_PU1_SpecificBinding_C57_seqs.fa,./data/PU1/C57_PU1_SpecificBinding_BALB_seqs.fa \
-o ./data/PU1/
```

It will concatenate both BALB-specific and C57-specific sequences before doing MAGGIE analysis. 