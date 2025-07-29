#!/usr/bin/env bash

git clone https://bitbucket.org/blaze-lib/blaze.git blaze.git
pushd blaze.git
cmake . -DCMAKE_INSTALL_PREFIX=$(readlink -f ..)/blaze.installed
make install
popd

mkdir build
cd build

cmake ..                                                                       \
  -DCMAKE_CXX_COMPILER=${BUILD_PREFIX}/bin/$(basename ${CXX})                  \
  -DCMAKE_INSTALL_PREFIX=${PREFIX}                                             \
  -DCMAKE_BUILD_TYPE=Release                                                   \
  -Dblaze_ROOT=$(readlink -f ../blaze.installed)                               \
  -DEnableMPI=OFF                                                              \
  -DTests=ON                                                                   \
  -DDocumentation=OFF

make -j2 VERBOSE=1
ctest --output-on-failure
make install
