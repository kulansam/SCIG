#!/usr/bin/env python

"""Description
Setup script for SCIG: Specifying the Cell Identity Genes in Single-cells 
@author:  Kulandaisamy Arulsamy
@contact: KulandaiSamy.Arulsamy@childrens.harvard.edu
"""

import sys
from setuptools import setup, find_packages


def main():
    
    setup(name="SCIG",
          version="1.0",
          description="SCIG: Specifying the Cell Identity Genes in Single-cells",
          author='Kulandaisamy Arulsamy',
          author_email='KulandaiSamy.Arulsamy@childrens.harvard.edu',
          packages=find_packages(),
          url='https://github.com/kulansam/SCIG',
          scripts=['src/SCIG.py',
                   ],
          include_package_data=True,
          package_data={
              '': ['data/*.pkl', 'test/*.txt'],
          },
          license='MIT',
          install_requires=[
              'numpy','pandas','qnorm','regex','scikit-learn','rpy2','bioinfokit','anndata','scanpy','argparse','scipy']
          )

if __name__ == '__main__':
    main()