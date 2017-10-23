from setuptools import setup
from Cython.Build import cythonize
import pysam

setup(
    name='count_splits',
    ext_modules=cythonize('helpers.pyx'),
    include_dirs=pysam.get_include()
)
