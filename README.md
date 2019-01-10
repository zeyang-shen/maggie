![maggie_simpson](https://user-images.githubusercontent.com/33430859/50864683-abfbd000-1357-11e9-8973-2bf7e3ca6dcd.png)

# MAGGIE
MAGGIE provides a framework for analyzing transcription factor motifs within DNA regulatory elements and for detecting effects of mutations on motifs when studying different conditions, individuals, genotypes, etc. 

## Environment setup and installation
Anaconda is required for environment setup and package installation (https://www.anaconda.com/download/#macos). After installing the anaconda, run the following commands to create an environment with required dependencies so that your home environment will not be affected:
```
conda env create --file environment.yml
```

Alternatively, you may also set up an environment using requirements.txt file:
```
conda env create maggie
pip install -r requirements.txt
```

After setting up a new environment, activate the environment and you are ready to use all the fantastic tools:
```
source activate maggie
```

## Quick Usage
### Python package
```
import maggie
```

### Command lines
All the command line tools are stored in the 'bin' directory. Below is an example of using one of the commands:
```
cd maggie
python ./bin/compute_score_diff.py -h
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
MAGGIE was developed primarily by Zeyang Shen, with significant contributions and suggestions by Marten Hoeksema. Supervision for the project was provided by Professor Christopher K. Glass. 
