from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
asvlite_wrapper_source_dir = asvlite_wrapper_dir.joinpath("source")
sys.path.insert(0, str(asvlite_wrapper_dir))
sys.path.insert(0, str(asvlite_wrapper_source_dir))

extensions = []
extensions.append(Extension("cds", sources=["cds.pyx"]))
extensions.append(Extension("epsg", sources=["epsg.pyx"]))
extensions.append(Extension("netcdf", sources=["netcdf.pyx"]))
extensions.append(Extension("thrust_calibrator", sources=["thrust_calibrator.pyx"]))

asvlite_include_dir = root_dir.joinpath("dependency", "ASVLite", "include")
setup(ext_modules=cythonize(extensions, language_level = "3"), include_dirs=[np.get_include(), asvlite_include_dir])