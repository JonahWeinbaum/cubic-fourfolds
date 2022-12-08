#include <assert.h>
#include "Fq.h"
#include <iostream>
#include <ctime>

int main() {

  // Tests only work for q=2^11.
  assert (initialized_ff_bitsize() == 11);
  
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

  std::cout << "Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC) << std::endl;

  return 0;
}
