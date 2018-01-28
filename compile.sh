#!/bin/sh

export CMAKE_BUILD_TYPE=Debug 
export BUILD_SHARED_LIBS=ON 
export CROSS=0 
export EXTRA_CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../../install-master"
rm -rf build
mkdir build
cd build
cmake -DBUILD_EXAMPLES=ON -DBUILD_SHARED_LIBS=$BUILD_SHARED_LIBS -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE $EXTRA_CMAKE_FLAGS ..
make VERBOSE=1
make test VERBOSE=1
make install VERBOSE=1
