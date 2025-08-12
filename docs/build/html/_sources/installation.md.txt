
# Installation

## Installation of dingo

if you want to import `dingo` inside the `dingo-stats` package you can skip this part and go directly to the `Installation of dingo-stats` section.

To install dingo you can follow instructions on the [github](https://github.com/GeomScale/dingo) page

Or install it through `pip`:

```bash
pip install dingo-walk
```

## Installation of dingo-stats

To run the package locally we recommend to create a virtual conda environment.

Make sure you have `conda` installed. To create the environment type:

```bash
conda create -n dingo-stats python=3.10.18
```

And then type `y` when this message appears : "The following NEW packages will be INSTALLED: ..."

To activate the environment type:

```bash
conda activate dingo-stats
```


`dingo-stats` relies on some functions implemented in the `dingo` package.

To install `dingo` type the following inside the activated environment:

```bash
CFLAGS="-I/usr/include/suitesparse" pip install --force-reinstall dingo-walk==0.1.6
```


And then install the rest dependencies of the package with:

```bash
python -m pip install --force -r requirements.txt
```


## How to build the ReadTheDocs of `dingo-stats` locally

```
cd docs
make clean
make html
```