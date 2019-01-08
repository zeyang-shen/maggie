from setuptools import setup

with open("README.txt", 'r') as f:
    long_description = f.read()

setup(name='VarAnalysis',
      version='0.1.0',
      description='A module to analyze DNA mutations',
      license='GNU GPLv3',
      long_description=long_description,
      author='Zeyang Shen',
      author_email='zes017@ucsd.edu',
      #url="http://www.foopackage.com/",
      packages=['VarAnalysis'],  #same as name
      install_requires=[
          'numpy', 
          'pandas',
          'biopython==1.70', 
          'scipy==1.1.0',
          'matplotlib==2.0.2',
          'seaborn==0.7.1'
      ],
      zip_safe=False,
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Programming Language :: Python :: 3',
          'Framework :: Jupyter',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ]
)

