#include <assert.h>
#include "Fq.h"
#include <iostream>
#include <ctime>


ff2k_t test_ff2k_mult(ff2k_t x, ff2k_t y) {
  // Stock multiplication algorithm.
  long a = (long)x, b = (long)y, ab = 0;

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

  return (ff2k_t)ab;
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
    std::clock_t cputime = std::clock();

    // Test Fq multiplication.
    for (unsigned i=0; i < (1 << 17); i++) {
      for (unsigned j=0; j<10000; j++) {
        assert (ff2k_mult(i, j) == test_ff2k_mult(i, j));
      }
    }
    
    std::cout << "Multiplication: Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC)
              << std::endl;

    for (unsigned i=0; i<100; i++)
      for (unsigned j=0; j<100; j++)
        for (unsigned a=1; a<50; a++) {
          
          unsigned c = ff2k_mult(ff2k_mult(i,j), a);
          unsigned b = ff2k_mult(i ^ j, a);
        
          unsigned g[] = {c, b, a};
          unsigned rts[] = {0,0,0};

          quadratic_roots(g, rts);
          if (i == j) {
            assert(count_quadratic_roots(g) == 1);
            g[0] = c;
            g[1] = b;
            g[2] = a;

            assert(rts[0] == 1 && rts[1] == i);
          } else {
            assert(count_quadratic_roots(g) == 2);
            g[0] = c;
            g[1] = b;
            g[2] = a;
          
            assert((rts[0] == 2 && rts[1] == i && rts[2] == j) ||
                   (rts[0] == 2 && rts[1] == j && rts[2] == i)
                   );
          }
        }


    // In the event that N is odd, we can easily check if certain polynomials
    // don't have a root
    for (unsigned i=0; i<100; i=i+2) // Even numbers
      for (unsigned a=1; a<50; a++) {
        if (initialized_ff_bitsize() % 2 == 1) {
          unsigned* g = new unsigned[3];
          
          g[0] = ff2k_mult(a, i ^ 1);
          g[1] = a;
          g[2] = a;         
          assert(count_quadratic_roots(g) == 0);

          g[0] = ff2k_mult(a, i ^ 1); // Ensure the last element is not a trace.
          g[1] = a;
          g[2] = a;
          unsigned rts[3] = {0,0,0};
          
          quadratic_roots(g, rts);
          assert(rts[0] == 0);
            
        }
        
      }

    std::cout << "Quadratic roots: Total time: "
              << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC)
              << std::endl;

        
  }

  // Quick check about integer sizes.
  __int128 Q = (1 << 22);
  std::string magmaQstring = "309485083608321363567181825";
  assert(int128_to_string(Q*Q*Q*Q + Q*Q*Q + Q + 1).compare(magmaQstring) == 0);
   
  return 0;
}
