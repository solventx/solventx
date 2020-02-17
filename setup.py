from setuptools import setup

setup(name='solventx',
      version='0.0.1',
      packages=['solventx',],
      description='Solvent Extraction configuration design model',
      author = 'Nwike I',
      author_email='nwike.iloeje@gmail.com',
      install_requires=['scipy>=1.2.0','numpy>=1.15.4', cantera>=2.4.0\
            'pandas>=0.24.2', 'seaborn>=0.9.0','pytest>=5.0.1'],	  
      )
