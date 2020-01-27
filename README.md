# MAGGIE
MAGGIE provides a framework for analyzing transcription factor motifs within DNA regulatory elements and for detecting effects of mutations on motifs when studying different conditions, individuals, genotypes, etc. 

The overview of the method is as below:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/method.png" width="800" height="300">
</p>

## Environment setup and installation
Anaconda is required for environment setup and package installation (https://www.anaconda.com/download/#macos). After installing the anaconda, run the following commands to automatically create an environment named "maggie" with required dependencies so that your home environment will not be affected:
```bash
conda env create --file environment.yml
```
After setting up a new environment, activate the environment:
```bash
conda activate maggie
```

Finally, copy the github folder and you are ready to run your own analysis:
```bash
git clone https://github.com/zeyang-shen/maggie.git
```

## Quick Usage
### Python package
```python
import maggie
```

### Command lines
All the command line tools are stored in the 'bin' directory. Below is an example of using MAGGIE with FASTA files as inputs:
```bash
cd maggie
python ./bin/maggie_fasta_input.py \
./data/ASB/CTCF_binding_alleles.fa \
./data/ASB/CTCF_nonbinding_alleles.fa \
-o ./data/ASB/ \
-p 1
```

## Interpretation of outputs
An example of significant hits displayed in the HTML format for CTCF allele-specific binding sites looks like this:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/html_example.png" width="900" height="300">
</p>

Column 1: ranking of significant hits based on absolute value of -log10(p-value) \\
Column 2: motif names within the same cluster based on a high correlation among changes of their motif scores \\
Column 3: PWM logo for the motif with lowest p-value in the cluster; the rest motifs usually look quite similar \\
Column 4: signed -log10(p-value) based on the original inputs and its 90% confidence interval after doing bootstrapping \\
Column 5: number of input sequences with mutations on the specific motif and its percentage in the total sequences \\
Column 6: number of input sequences with higher motif scores in the positive sequences and its percentage in the mutated sequences \\
Column 7: number of input sequences with higher motif scores in the negative sequences and its percentage in the mutated sequences \\
Column 8: the median value of all the non-zero motif score differences \\
Column 9: the mean value of all the non-zero motif score differences \\
Column 10: the distribution of all the non-zero motif score differences \\

## Contact
If you enconter a problem when using the software, you can
1. post an issue on [Issue](https://github.com/zeyang-shen/maggie/issues) section
2. or email Zeyang Shen by zes017@ucsd.edu

## License

[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/maggie/blob/master/LICENSE)

## Contributors
MAGGIE was developed primarily by Zeyang Shen, with contributions and suggestions by Marten Hoeksema and Zhengyu Ouyang. Supervision for the project was provided by Christopher K. Glass and Chris Benner. 
