#include <assert.h>
#include "Fq.h"
#include <iostream>
#include <ctime>


ff2k_t test_ff2k_mult(ff2k_t x, ff2k_t y) {
  // Stock multiplication algorithm.
  unsigned a = x, b = y, ab = 0;

  while (b != 0) {
    // b = const + higher-deg part
    // if const = 1 then ret += a
    if (b & 1) ab ^= a;
    // kill const, divide higher-deg part by t
    b >>= 1;
    // multiply a by t
    a <<= 1;
    if (a & q) a ^= p;
  }

  return ab;
}


int main() {

  if (initialized_ff_bitsize() == 11) {
    // Tests only work for q=2^11.

  
    // Test Fq multiplication.
    assert (ff2k_mult(3, 3) == 5);
    assert (ff2k_mult(3, 0) == 0);
    assert (ff2k_mult(0, 3) == 0);
    assert (ff2k_mult(1, 3) == 3);
    assert (ff2k_mult(3, 1) == 3);

    // Test Fq division
    assert (ff2k_divi(ff2k_mult(3,5), 5) == 3);

    std::clock_t cputime = std::clock();
  
    // Test gcd_degree
    unsigned f[] = {9, 3, 0, 0, 1}; // t^4 + (k.1+1) * t + (k.1^3 + 1);
    unsigned h[] = {211, 1, 1798, 0}; // k.1^1704*t^2 + t + k.1^1832

    assert (gcd_degree(f, h, 4, 3) == 2);
    
    // Test monic root counting.
    unsigned f2[] = {9, 3, 0, 0, 1}; // t^4 + (k.1+1) * t + (k.1^3 + 1);
    assert (count_poly_roots(f2, 4) == 2);

    // Test with a nontrivial leading coefficient.
    unsigned f3[] = {9, 3, 0, 0, 2}; // k.1 * t^4 + (k.1+1) * t + (k.1^3 + 1);
    assert (count_poly_roots(f3, 4) == 2);

    std::cout << "Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC)
              << std::endl;

  } else {
    // Test Fq multiplication.
    for (unsigned i=0; i < (1 << 17); i++) {
      for (unsigned j=0; j<10000; j++) {
        assert (ff2k_mult(i, j) == test_ff2k_mult(i, j));
      }
    }

    std::cout << "Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC)
              << std::endl;

  }
  
   
  return 0;
}
