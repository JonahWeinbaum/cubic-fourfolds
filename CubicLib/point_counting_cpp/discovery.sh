#!/bin/bash

run_discovery_job(){
    # Compile
    g++ -DN=5 -std=c++11 -DARM -O3 -DCOEFFSFILE="\"coeffs$1.h\"" count_bigq.cpp -o a.$1.out

    # Run
    ./a.$1.out
}

# For loop 5 times
for i in {1..1}
do
    run_discovery_job $i
done

#seq 1 196 | parallel

