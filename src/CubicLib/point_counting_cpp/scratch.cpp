// #include <wmmintrin.h>
#include <ctype.h>
#include <iostream>
#include <ctime>
#include <bit>
#include <bitset>
// #include "tableio.h"

const unsigned n = 15;
const unsigned FINITEFIELDBITSIZE = n;

// our irreducible polynomials
const unsigned polynomials[] = { 1, // placeholder
  1<<1, // for F_2, take x
  (1<<2) + (1<<1) + 1, // for F_{2^2}, take x^2 + x + 1
  (1<<3) + (1<<1) + 1, // for F_{2^3}, take x^3 + x + 1
  (1<<4) + (1<<1) + 1, // etc.
  (1<<5) + (1<<2) + 1,
  (1<<6) + (1<<1) + 1,
  (1<<7) + (1<<1) + 1,
  (1<<8) + (1<<4) + (1<<3) + (1<<1) + 1,
  (1<<9) + (1<<1) + 1,
  (1<<10) + (1<<3) + 1,
  (1<<11) + (1<<2) + 1,
  (1<<12) + (1<<3) + 1,
  (1<<13) + (1<<4) + (1<<3) + (1<<1) + 1,
  (1<<14) + (1<<5) + 1,
  (1<<15) + (1<<1) + 1,
  (1<<16) + (1<<5) + (1<<3) + (1<<1) + 1,
  (1<<17) + (1<<3) + 1,
  (1<<18) + (1<<3) + 1,
  (1<<19) + (1<<5) + (1<<2) + (1<<1) + 1,
  (1<<20) + (1<<3) + 1,
  (1<<21) + (1<<2) + 1,
  (1<<22) + (1<<1) + 1 };


// #ifdef WITHCACHE

// unsigned mult(unsigned a, unsigned b) {
//   return mult_table[a][b];
// }

// # elif defined FINITEFIELDBITSIZE


const unsigned p = polynomials[FINITEFIELDBITSIZE];

using ff2k_t = unsigned;

static inline ff2k_t _ff2k_pclmul (ff2k_t a, ff2k_t b)
{
  uint32_t c = 0;
  for (int i = 0; i < 15; i++) {
    c ^= ((a >> i) & 1) ? b : 0;
    c = c >> 1 | c << (15 - 1);
    c &= 0x7FFF;
    c ^= ((c >> 15) & 1) ? ((1<<15) + (1<<1) + 1) : 0;
  }
  //c ^= ((1<<22) + (1<<1) + 1)*((c >> 22) & 1);
  return c;
}

// static inline uint32_t _ff2k_pclmul(unsigned a, unsigned b)
// {
//   const unsigned n = FINITEFIELDBITSIZE;
//   register __m128i A, B, C;
//   A[0]=a; B[0] = b;
//   C = _mm_clmulepi64_si128 (A,B,0);
//   unsigned ab = (unsigned)C[0];
  
//   // Reduce modulo the polynomial.
//   unsigned pxi = p << n-2;
//   for (int i = 2*n-2; i >= n; i--) {
//     // If the i-th bit is 1, reduce by the polynomial.
//     if (ab & (1<<i)) ab ^= pxi;
//     pxi = pxi >> 1;
//   }     
//   return ab; 
// } 

// #endif

// Test the code
int main() {
  std::clock_t cputime = std::clock();
  
  // To benchmark, we are going to use GF(2^13) arithmetic.
  const unsigned q = 1<<n;

  // unsigned x=0;
  // for (int i=0; i < q; i++) {
  //   for (int j=0; j < q; j++) {
  //     x += _ff2k_pclmul(i, j);
  //   }
  // }

  std::cout << std::popcount((unsigned)1123) << std::endl;
  
  //std::cout << "Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC) << std::endl;
  //std::cout << "Dummy value: " << x << std::endl;
  
  return 0;

  // unsigned q = 1<<n;
  // unsigned** mult = read_table(q, q, "mult_" + std::to_string(q));
  // unsigned** mult_new = new unsigned*[q];
  
  // for (int i=0; i<q; i++) {
  //   mult_new[i] = new unsigned[q];
  //   for (int j=0; j<q; j++) {
  // //     uint32_t ret = _ff2k_pclmul(i, j);
  // //     // Reduce modulo the polynomial.
  // //     unsigned pxi = p << n-2;
  // //     for (int i = 2*n-2; i >= n; i--) {
  // //       // If the i-th bit is 1, reduce by the polynomial.
  // //       if (ret & (1<<i)) ret ^= pxi;
  // //       pxi = pxi >> 1;
  // //     }     
  //     mult_new[i][j] =  _ff2k_pclmul(i, j);
  //   }
  // }

 // compare_2_table(mult, mult_new, q, q);
}
