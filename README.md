# StormExplorers

## Introduction
Extends [ASVLite](https://github.com/resilient-swarms/ASVLite) to simulate a tropical storm and wave glider dynamics in the storm.

## Dependency
Some of the dependencies, such as cdsapi, are available only in `conda-forge`
```
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Install packages.
```
sudo apt install build-essential swig python3-dev
conda install cartopy netcdf4 bs4 requests toml tqdm pandas cdsapi pyproj shapely adjusttext scikit-learn xgboost seaborn cython jupyterlab gxx_linux-64 gcc_linux-64 swig auto-sklearn
```

Follow the steps defined [here](https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key) to install the CDS API key.

Although not a dependency for ASVLite, having the following tools will be useful to explore netCDF (.nc) files.

```
sudo apt install ncview netcdf-bin
```

Useful commands:
`ncdump` - view information and data in the `nc` file.
`ncview` - visualisation of data in the `nc` file. 

## Build instruction

Compile ASVLite Python wrapper. 
``` 
git clone --recurse-submodules https://github.com/resilient-swarms/StormExplorers.git
cd StormExplorers/dependency/ASVLite/wrapper/cython
python setup.py build_ext -i
```

Compile StormExplorers Cython source.
```
cd source
python setup.py build_ext -i
```