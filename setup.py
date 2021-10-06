import setuptools
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

with open('README.md', 'r') as fh:
  long_description = fh.read()

setuptools.setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=[Extension("uio_tools/eosinterx",
                           ["uio_tools/eosinterx.pyx"],
                           include_dirs=[np.get_include()],
                           extra_compile_args=[
                               "-march=native", "-fopenmp", "-O3"],
                           extra_link_args=["-fopenmp"]
                           )],
    name="uio_tools",
    version="0.0.1",
    author="Siddhant A. Deshmukh",
    author_email="siddhant593@gmail.com",
    description="A collection of tools to read in and analyse UIO files from CO5BOLD model atmospheres.",
    long_description=long_description,
    url="https://github.com/SiddhantDeshmukh/uio_tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6"
)
