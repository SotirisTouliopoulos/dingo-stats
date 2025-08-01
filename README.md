
# Welcome to dingo-stats, a python package for statistics on flux sampling datasets!


## Installation

To run the package locally we recommend to create a virtual conda environment.

Make sure you have `conda` installed. To create the environment type:

```
conda create -n dingo-stats python=3.10.18
```

And then type `y` when this message appears : "The following NEW packages will be INSTALLED: ..."

To activate the environment type:

```
conda activate dingo-stats
```


`dingo-stats` relies on some functions implemented in the `dingo` package.

To install `dingo` type the following inside the activated environment:

```
CFLAGS="-I/usr/include/suitesparse" pip install --force-reinstall dingo-walk==0.1.6
```


And then install the rest dependencies of the package with:

```
python -m pip install --force -r requirements.txt 
```


## Structure of the `dingo-stats` repository

- `src`: folder contains all the source code of the package
- `tests`: provides tutorials on how to use the functions found in `src`
- `analysis`: contains extended tutorials aiming to produce biologically significant knowledge
- `ext_data`: contains all external data (e.g. models, maps) used in `tests`


## Suggested order to explore implemented functions

To better understand how to use the functions from `dingo-stats` you may examine the jupyter notebooks under the `tests` folder in the following order:

- `tests/load_modify_sample.ipynb` alongside `tests/pathways.ipynb`
- `tests/loopless.ipynb`
- `tests/escher_maps.ipynb`
- `tests/distributions_comparison.ipynb`
- `tests/copulas.ipynb`
- `tests/correlations.ipynb`
- `tests/clustering.ipynb`
- `tests/graphs.ipynb`

Parts of the functions are used with the same way multiple times across the notebooks. 
To see a complete analysis from loading a model to performing all-possible statistics check the:

- `analysis/complete_workflow.ipynb`


## How to build the ReadTheDocs locally

```
cd docs
make clean
make html
```