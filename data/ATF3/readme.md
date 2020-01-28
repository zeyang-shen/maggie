### Data description
ATF3 strain-specific binding sites were identified based on ATF3 ChIP-seq data of macrophages for the two genetically diverse strains of mice: C57BL/6J (C57) and BALB/cJ (BALB) ([Fonseca et al., 2019](https://doi.org/10.1038/s41467-018-08236-0)). 

balbc_veh_atf3_spec_c57bl6_ref.fa: BALB-specific ATF3 binding sites, sequences from BALB genome

balbc_veh_atf3_spec_c57bl6_mut.fa: pair with "balbc_veh_atf3_spec_c57bl6_ref.fa", sequences from C57 genome

c57bl6_veh_atf3_spec_balbc_ref.fa: C57-specific ATF3 binding sites, sequences from C57 genome

c57bl6_veh_atf3_spec_balbc_mut.fa: pair with "c57bl6_veh_atf3_spec_balbc_ref.fa", sequences from BALB genome

Sequences were extracted using MMARGE "extract_sequences" ([Link et al., 2018](https://doi.org/10.1093/nar/gky491)).

The command line to process these files is:
```bash
python ./bin/maggie_fasta_input.py \
./data/ATF3/balbc_veh_atf3_spec_c57bl6_ref.fa,./data/ATF3/c57bl6_veh_atf3_spec_balbc_ref.fa \
./data/ATF3/balbc_veh_atf3_spec_c57bl6_mut.fa,./data/ATF3/c57bl6_veh_atf3_spec_balbc_mut.fa \
-o ./data/ATF3/
```

It will concatenate both BALB-specific and C57-specific sequences before doing MAGGIE analysis. 