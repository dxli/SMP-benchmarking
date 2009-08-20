#!/bin/bash

aclocal
libtoolize --copy --force --automake
autoconf
aclocal
autoheader
automake --add-missing --foreign
CFLAGS="-O2 -march=native -fomit-frame-pointer -pipe"
CXXFLAGS="${CFLAGS}"
./configure && make
