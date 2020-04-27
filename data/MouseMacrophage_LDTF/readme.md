### Data description
PU.1 and C/EBPb strain-specific binding sites were identified based on ChIP-seq data of macrophages for the two genetically diverse strains of mice: C57BL/6J (C57) and BALB/cJ (BALB) ([Link et al., 2018](https://doi.org/10.1016/j.cell.2018.04.018)). ATF3 strain-specific binding sites were also from C57 and BALB ([Fonseca et al., 2019](https://doi.org/10.1038/s41467-018-08236-0)). Sequences were extracted from their respective genomes using MMARGE "extract_sequences" ([Link et al., 2018](https://doi.org/10.1093/nar/gky491)).

BALB_\[TF\]\_C57_ref.fa: BALB-specific TF binding sites, sequences from BALB genome

BALB_\[TF\]\_C57_mut.fa: BALB-specific TF binding sites, sequences from C57 genome

C57_\[TF\]\_BALB_ref.fa: C57-specific TF binding sites, sequences from C57 genome

C57_\[TF\]\_BALB_mut.fa: C57-specific TF binding sites, sequences from BALB genome

The example command line for concatenating both BALB-specific and C57-specific sequences and running MAGGIE analysis on these binding sites is:
```bash
python ./bin/maggie_fasta_input.py \
./data/MouseMacrophage_LDTF/PU1/BALB_PU1_spec_C57_ref.fa,./data/MouseMacrophage_LDTF/PU1/C57_PU1_spec_BALB_ref.fa \
./data/MouseMacrophage_LDTF/PU1/BALB_PU1_spec_C57_mut.fa,./data/MouseMacrophage_LDTF/PU1/C57_PU1_spec_BALB_mut.fa \
-o ./data/MouseMacrophage_LDTF/PU1/maggie_output/
```