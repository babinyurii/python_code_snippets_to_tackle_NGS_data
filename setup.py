# this import is a quick fix of warning 'install_requires"
# https://stackoverflow.com/questions/8295644/pypi-userwarning-unknown-distribution-option-install-requires

from distutils.core import setup
import setuptools

setup(
  name = 'sequman',
  packages = ['sequman'],   
  version = '0.0.2',      
  license='MIT',        
  description = 'functions for routine work with sequence data',   
  author = 'Yuriy Babin',                  
  author_email = 'babin.yurii@gmail.com',      
  url = 'https://github.com/babinyurii/sequman', 
  download_url = 'https://github.com/babinyurii/sequman/archive/v_0.0.2.tar.gz',
  keywords = ['bioinformatics', 'sequencing'],   
  install_requires=[            
          'pandas',
          'biopython',
          'pyvcf',
          'matplotlib',
          'numpy',
          'seaborn'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
