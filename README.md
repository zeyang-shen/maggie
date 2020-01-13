# MAGGIE <img src="https://github.com/zeyang-shen/maggie/blob/master/image/Maggie_half.png" width="120" height="100">
MAGGIE provides a framework for analyzing transcription factor motifs within DNA regulatory elements and for detecting effects of mutations on motifs when studying different conditions, individuals, genotypes, etc. 

The overview of the method is as below:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/method.png" width="400" height="597">
</p>

## Environment setup and installation
Anaconda is required for environment setup and package installation (https://www.anaconda.com/download/#macos). After installing the anaconda, run the following commands to automatically create an environment named "maggie" with required dependencies so that your home environment will not be affected:
```bash
conda env create --file environment.yml
```
After setting up a new environment, activate the environment and you are ready to run your own analysis:
```bash
conda activate maggie
```

Finally, clone the github folder:
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

## Contact
If you enconter a problem when using the software, you can
1. post an issue on [Issue](https://github.com/zeyang-shen/maggie/issues) section
2. or email Zeyang Shen by zes017@ucsd.edu

## License

[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/maggie/blob/master/LICENSE)

## Contributors
MAGGIE was developed primarily by Zeyang Shen, with contributions and suggestions by Marten Hoeksema and Zhengyu Ouyang. Supervision for the project was provided by Christopher K. Glass and Chris Benner. 
