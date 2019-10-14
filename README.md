# MAGGIE <img src="https://github.com/zeyang-shen/maggie/blob/master/image/Maggie_half.png" width="120" height="100">
MAGGIE provides a framework for analyzing transcription factor motifs within DNA regulatory elements and for detecting effects of mutations on motifs when studying different conditions, individuals, genotypes, etc. 

The overview of the method is as below:

<p align="center">
<img src="https://github.com/zeyang-shen/maggie/blob/master/image/method.png" width="400" height="597">
</p>

## Environment setup and installation
Anaconda is required for environment setup and package installation (https://www.anaconda.com/download/#macos). After installing the anaconda, run the following commands to automatically create an environment named "maggie" with required dependencies so that your home environment will not be affected:
```
conda env create --file environment.yml
```
And then use pip to install required packages listed in requirements.txt:
```
pip install -r requirements.txt
```

After setting up a new environment, activate the environment and you are ready to run your own analysis:
```
conda activate maggie
```

## Quick Usage
### Python package
```
import maggie
```

### Command lines
All the command line tools are stored in the 'bin' directory. Below is an example of using MAGGIE with FASTA files as inputs:
```
cd maggie
python ./bin/maggie_computation_fasta_input.py ./examples/data/ENCODE_K562_CTCF/CTCF_orig.fa ./examples/data/ENCODE_K562_CTCF/CTCF_mut.fa -o ./examples/data/ENCODE_K562_CTCF/ -p 1
```


## Citation


## Contact
If you enconter a problem when using the software, you can
1. search your question in [FAQ](https://github.com/...)
2. or post an issue on [Issue](https://github.com/zeyang-shen/maggie/issues) section
3. or email Zeyang Shen by zes017@ucsd.edu

## License

[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/maggie/blob/master/LICENSE)

## Contributors
MAGGIE was developed primarily by Zeyang Shen, with contributions and suggestions by Marten Hoeksema and Zhengyu Ouyang. Supervision for the project was provided by Christopher K. Glass and Chris Benner. 
