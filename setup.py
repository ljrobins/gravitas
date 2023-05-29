from setuptools import setup, find_packages, Extension
import os
import numpy as np

_SOURCES = [os.path.join('gravitas', x) for x in os.listdir('gravitas') if '.c' == x[-2:]]
_INCDIR = ['gravitas', np.get_include()]
# _LIB_DIR
print(_SOURCES)
setup(
    name='gravitas',
    version='0.1.0',
    packages=find_packages(),
    license='GPL-2',
    requires=[
        'numpy',
    ],
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'gravitas._grav',
            # the files to compile into our module relative to ``setup.py``
            sources=_SOURCES,
            include_dirs=_INCDIR
        ),
    ],
    include_package_data=True,  # Include non-Python files in packages
    zip_safe=False,  # Allow the package to be unzipped without modification
)