#!/bin/bash -l

# TODO: Do this for all the available coefficient files.
g++ -std=c++11 -mpclmul -O4 -DCOEFFSFILE='"coeffs1.h"' -DN=11 count_bigq.cpp -o a1.out
