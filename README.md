
**Status:** Expect regular updates and bug fixes.
# Solvent Extraction process design and simulation model

Solvent Extraction process simulator

## Links
* Source code repository: https://github.com/solventx/solventx

## Installation
You can install the module directly from github with following commands:
```
git clone https://github.com/solventx/solventx.git
cd solventx
pip install -e .
```

### Dependencies  
#### Graphviz
[Graphviz](https://graphviz.gitlab.io/download/). Please add Graphviz dot executable and its directory to your system's PATH variable

Example (for Windows): 
`C:\Program Files (x86)\Graphviz2.38\bin` and 
`C:\Program Files (x86)\Graphviz2.38\bin\dot.exe`

#### Cantera
Cantera is not available with pip and need to be installed within a Conda environment. You can either install in an existing environment or within a freshly created environment. Detailed instructions can be found [here.](https://cantera.org/install/conda-install.html)

Additionally, add the file **elementz.xml** located in **/dependencies/** to the Cantera sys directory. 

Example (for Windows): `C:\Users\username\AppData\Local\Continuum\anaconda3\envs\myenv\Lib\site-packages\cantera\data`

#### Other Python Modules
SciPy, Numpy, Matlplotlib, Pandas, Seaborn

## Using the module
The module can be imported as a normal python module: 

```python
import solventx
```

Try out the Jupyter notebooks with a demo [here.]


## Issues
Please feel free to raise an issue for bugs or feature requests.

## Who is responsible?
Core developer:  
- Nwike Iloeje ciloeje@anl.gov  

Contributor:  
- Siby Jose Plathottam splathottam@anl.gov 
- Blake Richey blake.e.richey@gmail.com

## Acknowledgement  

## Citation
If you use this code please cite it as:

```
@misc{solventx,
  title = {{solventx}: Solvent extraction process simulator.},
  author = "{Nwike Iloeje}",
  howpublished = {\url{https://github.com/solventx/solventx}},
  url = "https://github.com/solventx/solventx",
  year = 2019,
  note = "[Online; accessed 28-August-2019]"
}
```

## Copyright and License  
Copyright © 2019, UChicago Argonne, LLC

SOLVENT EXTRACTION PROCESS SIMULATOR (SolventX) is distributed under the terms of [BSD-3 OSS License.](LICENSE.md)
