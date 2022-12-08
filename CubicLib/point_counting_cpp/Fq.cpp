#include <ctype.h>
#include <assert.h>
#include <iostream>

// Setup aliasing for the value of q.
#if defined N
const unsigned FINITEFIELDBITSIZE = N;
#elif defined n
const unsigned FINITEFIELDBITSIZE = n;
#else
static_assert(false, "No value of N was specified during compilation.");
#endif

//const unsigned n = FINITEFIELDBITSIZE;

// REMARK: These are not Conway polynomials!
// Defining polynomials for Fq.
const unsigned polynomials_addington[] = { 1, // placeholder
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
  (1<<15) + (1<<1) + 1, };
/* with n > 15 you'll just run out of memory
  (1<<16) + (1<<5) + (1<<3) + (1<<1) + 1,
  (1<<17) + (1<<3) + 1,
  (1<<18) + (1<<3) + 1,
  (1<<19) + (1<<5) + (1<<2) + (1<<1) + 1,
  (1<<20) + (1<<3) + 1,
  (1<<21) + (1<<2) + 1,
  (1<<22) + (1<<1) + 1 }; */

const unsigned t = 1;

// These are Conway polynomials.
const unsigned polynomials[] = {1, // placeholder
  (t << 1)  + 1,
  (t << 2)  + (t << 1) + 1,
  (t << 3)  + (t << 1) + 1,
  (t << 4)  + (t << 1) + 1,
  (t << 5)  + (t << 2) + 1,
  (t << 6)  + (t << 4) + (t << 3) + (t << 1) + 1,
  (t << 7)  + (t << 1) + 1,
  (t << 8)  + (t << 4) + (t << 3) + (t << 2) + 1,
  (t << 9)  + (t << 4) + 1,
  (t << 10) + (t << 6) + (t << 5) + (t << 3) + (t << 2) + (t << 1) + 1,
  (t << 11) + (t << 2) + 1,
  (t << 12) + (t << 7) + (t << 6) + (t << 5) + (t << 3) + (t << 1) + 1,
  (t << 13) + (t << 4) + (t << 3) + (t << 1) + 1,
  (t << 14) + (t << 7) + (t << 5) + (t << 3) + 1,
  (t << 15) + (t << 5) + (t << 4) + (t << 2) + 1,
  (t << 16) + (t << 5) + (t << 3) + (t << 2) + 1,
  (t << 17) + (t << 3) + 1,
  (t << 18) + (t << 12) + (t << 10) + (t << 1) + 1,
  (t << 19) + (t << 5)  + (t << 2) + (t << 1) + 1,
  (t << 20) + (t << 10) + (t << 9) + (t << 7) + (t << 6) + (t << 5) + (t << 4) + (t << 1) + 1,
  (t << 21) + (t << 6)  + (t << 5) + (t << 2) + 1,
  (t << 22) + (t << 12) + (t << 11) + (t << 10) + (t << 9) + (t << 8) + (t << 6) + (t << 5) + 1
};


const unsigned q = 1 << FINITEFIELDBITSIZE;
const unsigned p = polynomials[FINITEFIELDBITSIZE];
using ff2k_t = unsigned;

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

/*
int count_linear_roots(unsigned a, unsigned b) {
  if ((a==0) and (b==0)) return q; // Identically zero
  return (a==0);
}

int count_quadratic_roots(unsigned a, unsigned b, unsigned c) {

  if (a == 0) return count_linear_roots(b, c);

  // It's likely that by transforming to a standard form
  //
  //    x^2 + x + A
  //
  // and checking whether A lies in the space of traces, we can do things
  // much faster. For the sake of testing, we hold off on this approach.

  // Replace f 
  unsigned f[] = {c,b,a};
  make_monic_nodiv(f, 2);
  
  // q = 2^n.
  // Compute x^q mod f(x).
  unsigned powx[]   = {0,1};
  unsigned sqpowx[] = {0,0,0};
  
  for (int i=0; i < n; i++) {
    // Square
    for (int j=0; j<=1; j++) {
      sqpowx[2*j] = ff2k_square(powx[j]);
    }

    // Reduce
    for (int j=1; j>=0; j--) {
      if (sqpowx[j+3] != 0) {
        sqpowx[j+1] ^= ff2k_mult(sqpowx[j+3], f[1]);
        sqpowx[j+0] ^= ff2k_mult(sqpowx[j+3], f[0]);
      }
    }

    // Repeat
    powx[0] = sqpowx[0];
    powx[1] = sqpowx[1];

    sqpowx[0] = 0;
    sqpowx[1] = 0;
    sqpowx[2] = 0;
  }

  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  unsigned coeffs[] = {c,b,a};
  return gcd_degree(coeffs, powx, 2, 1);
}

int count_cubic_roots(unsigned a, unsigned b, unsigned c, unsigned d) {

  // The key trick is to quickly compute the GCD with x^q-x, then check
  // the degree.

  if (a == 0) return count_quadratic_roots(b, c, d);

  // Replace f 
  unsigned f[] = {d,c,b,a};
  make_monic_nodiv(f, 3);

  // q = 2^n.
  // Compute x^q mod f(x). Note that we scale up x rather than scale down f.
  unsigned powx[]   = {0,1,0};
  unsigned sqpowx[] = {0,0,0,0,0};
  
  for (int i=0; i < n; i++) {
    // Square
    for (int j=0; j<=2; j++) {
      sqpowx[2*j] = ff2k_square(powx[j]);
    }

    // Reduce
    for (int j=2; j>=0; j--) {
      if (sqpowx[j+3] != 0) {
        sqpowx[j+2] ^= ff2k_mult(sqpowx[j+3], f[2]);
        sqpowx[j+1] ^= ff2k_mult(sqpowx[j+3], f[1]);
        sqpowx[j+0] ^= ff2k_mult(sqpowx[j+3], f[0]);
      }
    }

    // Repeat
    powx[0] = sqpowx[0];
    powx[1] = sqpowx[1];
    powx[2] = sqpowx[2];

    sqpowx[0] = 0;
    sqpowx[1] = 0;
    sqpowx[2] = 0;
    sqpowx[3] = 0;
    sqpowx[4] = 0;
  }

  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  unsigned coeffs[] = {d,c,b,a};
  return gcd_degree(coeffs, powx, 3, 2);
}



int count_quartic_roots(unsigned a, unsigned b, unsigned c, unsigned d, unsigned e) {

  // The key trick is to quickly compute the GCD with x^q-x, then check
  // the degree.

  if (a == 0) return count_cubic_roots(b, c, d, e);

  // Replace f 
  unsigned f[] = {e,d,c,b,a};
  make_monic_nodiv(f, 4);
  
  // q = 2^n.
  // Compute x^q mod f(x). 
  unsigned powx[]   = {0,1,0,0};
  unsigned sqpowx[] = {0,0,0,0,0,0,0};
  
  for (int i=0; i < n; i++) {
    // Square
    for (int j=0; j<=3; j++) {
      sqpowx[2*j] = ff2k_square(powx[j]);
    }

    // Reduce
    for (int j=2; j>=0; j--) {
      if (sqpowx[j+4] != 0) {
        sqpowx[j+3] ^= ff2k_mult(sqpowx[j+4], f[3]);
        sqpowx[j+2] ^= ff2k_mult(sqpowx[j+4], f[2]);
        sqpowx[j+1] ^= ff2k_mult(sqpowx[j+4], f[1]);
        sqpowx[j+0] ^= ff2k_mult(sqpowx[j+4], f[0]);
      }
    }

    // Repeat
    powx[0] = sqpowx[0];
    powx[1] = sqpowx[1];
    powx[2] = sqpowx[2];
    powx[3] = sqpowx[3];

    sqpowx[0] = 0;
    sqpowx[1] = 0;
    sqpowx[2] = 0;
    sqpowx[3] = 0;
    sqpowx[4] = 0;
    sqpowx[5] = 0;
    sqpowx[6] = 0;

  }

  // std::cout << powx[0] << " " << powx[1] << " " << powx[2] << " " << powx[3] << std::endl;
  
  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  return gcd_degree(f, powx, 4, 3);
}

int count_quadratic_roots(unsigned* f) {
  return count_quadratic_roots(f[2], f[1], f[0]);
}

int count_cubic_roots(unsigned* f) {
  return count_cubic_roots(f[3], f[2], f[1], f[0]);
}

int count_quartic_roots(unsigned* f) {
  return count_quartic_roots(f[4], f[3], f[2], f[1], f[0]);
}
*/

int count_poly_roots(unsigned* poly, int deg) {

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
