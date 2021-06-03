# geo
Basic geometry algorithms and calculations.

## Compile
    ./get_ext.sh
    mkdir build  &&  pushd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2 -march=native -s" -DCMAKE_VERBOSE_MAKEFILE=True ..  &&  make -j4
