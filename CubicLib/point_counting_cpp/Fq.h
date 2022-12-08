#include <ctype.h>
#include <assert.h>
#include <iostream>
#include "constants.h"


int initialized_ff_bitsize() {
  return FINITEFIELDBITSIZE;
}

#ifndef ARM
#include <wmmintrin.h>

inline ff2k_t ff2k_mult(ff2k_t a, ff2k_t b)
{
  const unsigned n = FINITEFIELDBITSIZE;
  register __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  ff2k_t ab = (ff2k_t)C[0];
  
  // Reduce modulo the polynomial.
  ff2k_t pxi = p << n-2;
  for (int i = 2*n-2; i >= n; i--) {
    // If the i-th bit is 1, reduce by the polynomial.
    if (ab & (1<<i)) ab ^= pxi;
    pxi = pxi >> 1;
  }     
  return ab; 
} 

// For evaluating complicated polynomials it can help to delay reduction.
static inline uint64_t _ff2k_pclmul_delred(unsigned a, unsigned b)
{
  register __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  return C[0];  
} 

#else

ff2k_t ff2k_mult(ff2k_t x, ff2k_t y) {
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
#endif

inline ff2k_t ff2k_square(ff2k_t a) {
  return ff2k_mult(a, a);
}

// NOTE: No division by zero check.
inline ff2k_t ff2k_inv(ff2k_t a) {
  const int n = FINITEFIELDBITSIZE;

  // Use Lagrange to invert. that is, compute a^(2^n-2).
  // Note the binary expansion of 2^n-2 = 111...110.
  ff2k_t inva = 1;
  ff2k_t sqac = ff2k_square(a);
  for (int i=1; i<n; i++) {
    inva = ff2k_mult(inva, sqac);
    sqac = ff2k_square(sqac);
  }

  
  return inva;
}

ff2k_t ff2k_divi(ff2k_t a, ff2k_t b) {
  return ff2k_mult(a, ff2k_inv(b));
}


/*
unsigned reduce_cubic(uint64_t ab) {

  // Reduce modulo the polynomial.
  uint64_t pxi = (uint64_t)p << n-2;
  for (int i = 3*n-3; i >= n; i--) {
    // If the i-th bit is 1, reduce by the polynomial.
    if (ab & (1<<i)) ab ^= pxi;
    pxi = pxi >> 1;
  }     
  return (unsigned)ab; 
}
*/

////////////////////////////////////////////////////////////////////////////////
//
// Polynomial root counting.
//
////////////////////////////////////////////////////////////////////////////////

void make_monic_nodiv(unsigned* f, int d) {
  // Use a division-free algorithm to make f monic via the transform
  // f(x) -> a^(d-1) f(x/a). It is assumed the leading coefficient of f is non-zero and
  // that d = degree(f).

  unsigned powa = 1;
  for (int i=1; i<=d; i++) {
    f[d-i] = ff2k_mult(f[d-i], powa);
    powa = ff2k_mult(powa, f[d]);
  }
  f[d] = 1;
  return;
}

int gcd_degree(unsigned* A, unsigned* B, int dA, int dB) {

  int degA = dA; // Note: degA must be accurate. degB is just a bound.
  int degB = dB;

  // Possibly correct the degree of B.
  while (degB >= 0) {
    if (B[degB] != 0) break;
    degB--;
  }
  
  /////////////////
  // Do Euclid's algorithm.
  
  while (degB > 0) {

    // Main Euclid loop.
    for (int i=0; i <=degA - degB; i++) {
      int lt = A[degA-i];

      // Scale up A. (Instead of making B monic.)
      // Note we deliberately avoid touching the leading coefficient, since
      // we still need this value for the reductions.
      for (int j=0; j<degA-i; j++) {
        A[j] = ff2k_mult(A[j], B[degB]);
      }

      // Perform the reduction.
      for (int j=1; j <= degB; j++) {
        A[degA-i-j] ^= ff2k_mult(B[degB-j], lt);
      }
      A[degA-i] = 0; // Discard the leading coefficient, which is theoretically killed.
    }
    
    // Swap
    unsigned* C = A;
    A = B;
    B = C;
    degA = degB;
    degB = degB - 1;

    while (degB >= 0) {
      if ( B[degB] != 0) break;
      degB--;
    }
  }

  // Final count.
  return (degB < 0) ? degA : 0;
}


int count_poly_roots(unsigned* poly, int deg) {
  // Probably should put this in the header file, so that the compiler can
  // see that the degree is bounded by 4. Maybe even template this...

  if (deg==1) {
    if ((poly[deg] == 0) and (poly[0] == 0)) return q; // Identically zero
    return (poly[deg] == 0);
  }

  const unsigned n = FINITEFIELDBITSIZE;
  
  // The key trick is to quickly compute the GCD with x^q-x, then check
  // the degree.

  if (poly[deg] == 0) return count_poly_roots(poly, deg-1);

  // Declare f. Note this mutates poly.
  unsigned* f = poly;
  make_monic_nodiv(f, deg);
  
  // q = 2^n.
  // Compute x^q mod f(x). 
  unsigned powx[deg];
  unsigned sqpowx[2*deg-1];
 
  // Initialize
  for (int ll=0; ll <= deg-1; ll++) {
    powx[ll] = 0;
  }
  for (int ll=0; ll <= 2*deg-2; ll++) {
    sqpowx[ll] = 0;
  }          
  powx[1] = 1;

  // Square-Reduce loop
  for (int i=0; i < n; i++) {
    // Square
    for (int j=0; j<=deg-1; j++) {
      sqpowx[2*j] = ff2k_square(powx[j]);
    }

    // Reduce
    for (int j = deg-2; j>=0; j--) {
      if (sqpowx[j+deg] != 0) {
        for (int ll=0; ll < deg; ll++) {
          sqpowx[j+ll] ^= ff2k_mult(sqpowx[j+deg], f[ll]);
        }
      }
    }

    // Repeat
    for (int ll=0; ll <= deg-1; ll++) {
      powx[ll] = sqpowx[ll];
    }
    for (int ll=0; ll <= 2*deg-2; ll++) {
      sqpowx[ll] = 0;
    }          
  }

  // std::cout << powx[0] << " " << powx[1] << " " << powx[2] << " " << powx[3] << std::endl;
  
  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  return gcd_degree(f, powx, deg, deg-1);
}
