#include <wmmintrin.h>
#include <ctype.h>
#include <iostream>
#include "tableio.h"

const unsigned n = 11;
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
  const unsigned n = FINITEFIELDBITSIZE;
  register __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  ff2k_t ab = (ff2k_t)C[0];
  
  // Reduce modulo the polynomial.
  unsigned pxi = p << n-2;
  for (int i = 2*n-2; i >= n; i--) {
    // If the i-th bit is 1, reduce by the polynomial.
    if (ab & (1<<i)) ab ^= pxi;
    pxi = pxi >> 1;
  }     
  return ab;
}

// #endif

// Test the code
int main() {
    unsigned q = 1<<n;
    unsigned** mult = read_table(q, q, "mult_" + std::to_string(q));
    unsigned** mult_new = new unsigned*[q];
    
    for (int i=0; i<q; i++) {
      mult_new[i] = new unsigned[q];
      for (int j=0; j<q; j++) {
	mult_new[i][j] = _ff2k_pclmul(i, j);
      }
    }

    compare_2_table(mult, mult_new, q, q);
    return 0;
}
