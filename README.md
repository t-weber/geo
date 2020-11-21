# geo
Basic geometry algorithms.


## Compile
    mkdir build  &&  pushd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2 -march=native -s" -DCMAKE_VERBOSE_MAKEFILE=True ..  &&  make -j4
