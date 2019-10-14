Data description: binding sites called from [ENCODE K562 CTCF ChIP-seq data](https://www.encodeproject.org/experiments/ENCSR000DWE/) which have at least 0.4 fold change difference between two alleles of K562 cells. 

CTCF_orig.fa: sequences with stronger binding signals of CTCF

CTCF_mut.fa: sequences with weaker binding signals of CTCF

Genomes of two alleles were first generated using [MMARGE](https://academic.oup.com/nar/article/46/14/7006/5035172), and the ChIP-seq reads were mapped to both alleles and assigned to one of the alleles with a higher mapping score.
