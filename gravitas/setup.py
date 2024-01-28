from distutils.core import setup
from pybind11.setup_helpers import Pybind11Extension
import eigency
import platform

std_arg = '-std=c++17'
opt_arg = '-O3'
if platform.system() == 'Windows':
    std_arg = '/std:c++17'
    opt_arg = '/O2'

ext_modules = [
    Pybind11Extension(
        name="gravitas",
        sources=["python_bindings.cpp"],
        include_dirs=['.', *tuple(eigency.get_includes())],
        extra_compile_args=[std_arg, opt_arg],
    ),
]

setup(name = 'gravitas',
    ext_modules=ext_modules)