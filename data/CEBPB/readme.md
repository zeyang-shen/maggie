### Data description
C/EBPb strain-specific binding sites were identified based on C/EBPb ChIP-seq data of macrophages for the two genetically diverse strains of mice: C57BL/6J (C57) and BALB/cJ (BALB) ([Link et al., 2018](https://doi.org/10.1016/j.cell.2018.04.018)). 

BALBC_notx_CEBPB_spec_C57Bl6_ref.fa: BALB-specific C/EBPb binding sites, sequences from BALB genome

BALBC_notx_CEBPB_spec_C57Bl6_mut.fa: pair with "BALBC_notx_CEBPB_spec_C57Bl6_ref.fa", sequences from C57 genome

C57Bl6_notx_CEBPB_spec_BALBC_ref.fa: C57-specific C/EBPb binding sites, sequences from C57 genome

C57Bl6_notx_CEBPB_spec_BALBC_mut.fa: pair with "C57Bl6_notx_CEBPB_spec_BALBC_ref.fa", sequences from BALB genome

Sequences were extracted using MMARGE "extract_sequences" ([Link et al., 2018](https://doi.org/10.1093/nar/gky491)).

The command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
./data/CEBPB/BALBC_notx_CEBPB_spec_C57Bl6_ref.fa,./data/CEBPB/C57Bl6_notx_CEBPB_spec_BALBC_ref.fa \
./data/CEBPB/BALBC_notx_CEBPB_spec_C57Bl6_mut.fa,./data/CEBPB/C57Bl6_notx_CEBPB_spec_BALBC_mut.fa \
-o ./data/CEBPB/
```

It will concatenate both BALB-specific and C57-specific sequences before doing MAGGIE analysis. 