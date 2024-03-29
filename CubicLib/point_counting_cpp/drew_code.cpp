#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <wmmintrin.h>
#include <ctype.h>

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

const unsigned p = polynomials[n];

static inline uint32_t _ff2k_pclmul(unsigned a, unsigned b)
{
  const unsigned n = FINITEFIELDBITSIZE;
  register __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  unsigned ab = (unsigned)C[0];
  
  // Reduce modulo the polynomial.
  unsigned pxi = p << n-2;
  for (int i = 2*n-2; i >= n; i--) {
    // If the i-th bit is 1, reduce by the polynomial.
    if (ab & (1<<i)) ab ^= pxi;
    pxi = pxi >> 1;
  }     
  return ab; 
} 

static inline uint64_t _ff2k_pclmul_delred(unsigned a, unsigned b)
{
  register __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  return C[0];  
} 

unsigned reduce_cubic(int64_t ab) {

  // Reduce modulo the polynomial.
  int64_t pxi = (int64_t)p << n-2;
  for (int i = 3*n-3; i >= n; i--) {
    // If the i-th bit is 1, reduce by the polynomial.
    if (ab & (1<<i)) ab ^= pxi;
    pxi = pxi >> 1;
  }     
  return (unsigned)ab; 

}

int main() {
  std::clock_t cputime = std::clock();
  
  // To benchmark, we are going to use GF(2^13) arithmetic.
  const unsigned q = 1<<n;

  unsigned x=0;
  for (int i=0; i < q; i++) {
    for (int j=0; j < q; j++) {
      x += _ff2k_pclmul(i, j);
    }
  }

  std::cout << "Total time: " << (std::clock() - cputime) * (1.0/CLOCKS_PER_SEC) << std::endl;
  std::cout << "Dummy value: " << x << std::endl;
  
  return 0;
}
