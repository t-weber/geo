# geo
Basic geometry algorithms and calculations.

[![DOI: 10.5281/zenodo.4297475](https://zenodo.org/badge/DOI/10.5281/zenodo.4297475.svg)](https://doi.org/10.5281/zenodo.4297475)


## Compile
    ./get_ext.sh
    mkdir build  &&  pushd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2 -march=native -s" -DCMAKE_VERBOSE_MAKEFILE=True ..  &&  make -j4
