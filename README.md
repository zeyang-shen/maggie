[![python-version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/zeyang-shen/maggie/issues)
[![DOI](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btaa476.svg)](https://doi.org/10.1093/bioinformatics/btaa476)

# MAGGIE
MAGGIE provides a framework for identifying DNA sequence motifs mediating transcription factor binding and function. By leveraging measurements and genetic variation information from different genotypes (human individuals, animal strains, or alleles), MAGGIE associates the mutation of DNA sequence motif with various types of epigenomic features, including but not limited to transcription factor binding, open chromatin, histone modification, and stimulus response of regulatory elements. 

Here is the overview of the method:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/method.png" width="900" height="280">
</p>

## Installing MAGGIE
First, copy the github folder and go into the "maggie" folder:
```bash
git clone https://github.com/zeyang-shen/maggie.git
cd maggie
```
Next, configure an environment where MAGGIE can work. Anaconda is required for environment setup and package installation (https://www.anaconda.com/download/#macos). 

After installing the anaconda, run the following command to automatically create an environment named "maggie" with required dependencies:
```bash
conda env create --file environment.yml
```
After setting up the environment, activate it:
```bash
conda activate maggie
```
Now you are ready to run your own analysis under the "maggie" folder!

## Quick Usage
All the executable scripts are stored in the `bin/` directory. Here is the usage of MAGGIE on a toy example of CTCF allele-specific binding sites stored in FASTA files.

Let's first go into the cloned folder:
```bash
cd maggie
```
Then you can run the script for FASTA inputs as below:
```bash
python ./bin/maggie_fasta_input.py \
./data/AlleleSpecificBinding/CTCF_binding_alleles.fa \
./data/AlleleSpecificBinding/CTCF_nonbinding_alleles.fa \
-o ./data/AlleleSpecificBinding/maggie_output/ \
-p 8
```
After the job is done, open the "mergedSignificant.html" file at "data/AlleleSpecificBinding/maggie_output/" with your web browser and take a look at the significant motifs. 

Alternatively, you can add the `bin/` directory to your `PATH` in order to execute those scripts from anywhere:
```bash
export PATH=/path/to/your/cloned/maggie/bin:$PATH
```
Then you can execute the previous script by `maggie_fasta_input.py` directly.

Go to our [tutorials](https://github.com/zeyang-shen/maggie/wiki/Tutorial) for usage of MAGGIE in other cases. 

## Example output
MAGGIE will display significant motifs in the HTML format. Here is an example for CTCF allele-specific binding sites:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/html_example.png" width="900" height="200">
</p>

Header: total number of samples

Column 1: ranking based on absolute value of -log10(p-value)

Column 2: merged motifs based on a high correlation among changes of their motif scores

Column 3: PWM logo for the motif with lowest p-value

Column 4: signed -log10(p-value) and 90% confidence interval

Column 5: # and percentage of sequences with motif mutations

Column 6: # sequences with higher motif scores in the positive set and its fraction of Column 5 

Column 7: # sequences with higher motif scores in the negative set and its fraction of Column 5 

Column 8: median value of non-zero motif score differences

Column 9: mean value of non-zero motif score differences

Column 10: distribution of non-zero motif score differences

## Documentation
Please go to our [wiki page](https://github.com/zeyang-shen/maggie/wiki) for more detailed usage of MAGGIE.

## Citation
If you use our findings, the software, or the NF-kb ChIP-seq data at [GEO:GSE144070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144070), please cite

[Shen, et al. MAGGIE: leveraging genetic variation to identify DNA sequence motifs mediating transcription factor binding and function. Bioinformatics, 2020](https://doi.org/10.1093/bioinformatics/btaa476)

## Contact
If you enconter a problem when using the software, you can
1. check the [FAQ](https://github.com/zeyang-shen/maggie/wiki/FAQ) page
2. post an issue on [Issue](https://github.com/zeyang-shen/maggie/issues) section
3. or email Zeyang Shen by zes017@ucsd.edu

## License
[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/maggie/blob/master/LICENSE)

## Contributors
MAGGIE was developed primarily by Zeyang Shen, with contributions and suggestions by Marten Hoeksema and Zhengyu Ouyang. Supervision for the project was provided by Christopher Glass and Christopher Benner. 
