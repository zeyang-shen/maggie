# MAGGIE
MAGGIE provides a framework for identifying DNA sequence motifs mediating transcription factor binding and function. By leveraging measurements and genetic variation information from different genotypes (human individuals, animal strains, or alleles), MAGGIE associates the mutation of DNA sequence motif with various types of epigenomic features, including but not limited to transcription factor binding, open chromatin, histone modification, and stimulus response of regulatory elements. 

The overview of the method is as below:


<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/method.png" width="800" height="250">
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

Finally, copy the github folder and you are ready to run your own analysis under the "maggie" folder:
```bash
git clone https://github.com/zeyang-shen/maggie.git
cd maggie
```

## Quick Usage
### Python package
```python
from maggie import score, utils

# read positive and negative sequences
pos_seq_file = './data/ASB/CTCF_binding_alleles.fa'
neg_seq_file = './data/ASB/CTCF_nonbinding_alleles.fa'
pos_seq_dict = utils.read_fasta(pos_seq_file)
neg_seq_dict = utils.read_fasta(neg_seq_file)

# load motif PWMs
motif_dict = score.load_motifs('./data/JASPAR2020_CORE_vertebrates_motifs/')

# test for one motif, and store p-values and score differences
name, ID, _, pvals, score_diffs = score.test_one_motif(motif_dict['CTCF$MA0139.1'], pos_seq_dict, neg_seq_dict)

# print results
print('Signed -log10(p-value) for %s-%s is: %.2f' % (name, ID, pvals[0]))
```

### Command lines
All the command line tools are stored in the 'bin' directory. Below is an example of using MAGGIE with FASTA files as inputs:
```bash
python ./bin/maggie_fasta_input.py \
./data/ASB/CTCF_binding_alleles.fa \ #positive sequences
./data/ASB/CTCF_nonbinding_alleles.fa \ #negative sequences
-o ./data/ASB/ \ #output path; results will be stored in a folder named "maggie_output"
-p 1 #number of processors to use, default: 1
```

## Example outputs
An example of significant hits displayed in the HTML format for CTCF allele-specific binding sites looks like this:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/html_example.png" width="900" height="200">
</p>

Column 1: ranking of significant hits based on absolute value of -log10(p-value)

Column 2: motif names within the same cluster based on a high correlation among changes of their motif scores

Column 3: PWM logo for the motif with lowest p-value in the cluster; the rest motifs usually look quite similar

Column 4: signed -log10(p-value) based on the original inputs and its 90% confidence interval after doing bootstrapping

Column 5: number of input sequences with mutations on the specific motif and its percentage in the total sequences

Column 6: number of input sequences with higher motif scores in the positive sequences and its percentage in the mutated sequences 

Column 7: number of input sequences with higher motif scores in the negative sequences and its percentage in the mutated sequences 

Column 8: the median value of all the non-zero motif score differences

Column 9: the mean value of all the non-zero motif score differences

Column 10: the distribution of all the non-zero motif score differences


## Contact
If you enconter a problem when using the software, you can
1. post an issue on [Issue](https://github.com/zeyang-shen/maggie/issues) section
2. or email Zeyang Shen by zes017@ucsd.edu

## License

[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/maggie/blob/master/LICENSE)

## Contributors
MAGGIE was developed primarily by Zeyang Shen, with contributions and suggestions by Marten Hoeksema and Zhengyu Ouyang. Supervision for the project was provided by Christopher Glass and Christopher Benner. 
