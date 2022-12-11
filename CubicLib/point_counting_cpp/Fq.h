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
  __m128i A, B, C;
  A[0]=a; B[0] = b;
  C = _mm_clmulepi64_si128 (A,B,0);
  ff2k_t ab = (ff2k_t)C[0];
  
  // Reduce modulo the polynomial.
  ff2k_t pxi = n > 1 ? (p << (n-2)) : p; // Should compile away.
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
  __m128i A, B, C;
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

ff2k_t ff2k_sqrt(ff2k_t a) {
  const int n = FINITEFIELDBITSIZE;
  ff2k_t s = a;
  for (int i=1; i<=n-1; i++) {
    s = ff2k_square(s); // Compute a^(2^(n-1)).
  }
  return s;
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
// Arf invariants
//
////////////////////////////////////////////////////////////////////////////////

// In the special case b=1, we save a few operations.
int Arf_invariant_b_equals_1(unsigned a, unsigned c) {
  const unsigned  n = FINITEFIELDBITSIZE;
  const unsigned* trace_basis = trace_bases[n];
  const unsigned* pretrace_basis = pretrace_bases[n];

  unsigned acbb  = ff2k_mult(a, c);
  
  for (int i = 0; i < n-1; i++) {
    if (acbb & trace_basis[i]) {
       acbb ^= trace_basis[i];
    }
  }

  // The trace is zero if and only if what remains is zero.
  return (int)(acbb != 0);
}

// Arf_invariant(unsigned a, unsigned b, unsigned c)
// 
// Given a quadratic polynomial ax^2 + bx + c with b != 0, determine the
// Arf invariant (Equal to 0 or 1).
//
int Arf_invariant(unsigned a, unsigned b, unsigned c) {
  
  unsigned binv2 = ff2k_square(ff2k_inv(b));
  unsigned acbb  = ff2k_mult(ff2k_mult(a, c), binv2);

  return Arf_invariant_b_equals_1(1, acbb);
}

int Arf_invariant_mu2_b_equals_1(unsigned a, unsigned c) {
  return Arf_invariant_b_equals_1(a, c) == 1 ? -1 : 1;
}

int Arf_invariant_mu2(unsigned X, unsigned Y, unsigned Z) {
  return Arf_invariant(X, Y, Z) == 1 ? -1 : 1;
}


////////////////////////////////////////////////////////////////////////////////
//
// Polynomial utilities.
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


////////////////////////////////////////////////////////////////////////////////
//
// Polynomial root counting.
//
////////////////////////////////////////////////////////////////////////////////

int count_poly_roots(unsigned* f, int deg);

int count_quadratic_roots(unsigned* f) {
  // Return an array {n, r1, r2}. The first number indicates the number
  // of distinct roots. The remaining values are the actual roots.
  //
  // It is assumed that f is genuinely a quadratic.
  const int n = FINITEFIELDBITSIZE;
  
  // ax^2 + bx + c.
  unsigned c = f[0];
  unsigned b = f[1];
  unsigned a = f[2];

  if (a==0)
    throw std::invalid_argument("Leading coefficient must be non-zero.");
  
  if (b==0) return 1;  // The quadratic has a double root.

  int arf = Arf_invariant(f[2], f[1], f[0]);

  if (!(arf == 0 || arf == 1)) {
    std::cerr << arf << std::endl;
    std::cerr << (int)(2 != 0) << std::endl;
    assert(arf == 0 || arf == 1);
  }
  return (arf ^ 1) << 1;
}


int count_quartic_roots(ff2k_t* poly) {
  // The polynomial needs to be a genuine quartic.

  const unsigned n = FINITEFIELDBITSIZE;
  const int deg = 4;
  
  // Declare f. Note this mutates poly.
  ff2k_t* f = poly;

  // First we use a divisionless algorithm to make f depressed.
  ff2k_t b = f[3];
  ff2k_t d = f[1];

  if (b != 0) {
    // Multiply by b^4 and absorb.
    ff2k_t bb   = ff2k_square(b);
    ff2k_t bbb  = ff2k_mult(bb, b);
    ff2k_t bbbb = ff2k_square(bb);
    
    f[3] = bb;
    f[2] = ff2k_mult(f[2], bb);
    f[1] = ff2k_mult(f[1], bbb);
    f[0] = ff2k_mult(f[0], bbbb);

    // Now shift by alpha = sqrt(b'/d') = sqrt(db).
    ff2k_t alpha = ff2k_sqrt(ff2k_mult(b, d));

    ff2k_t A0 = f[0];
    ff2k_t powal = alpha;
    for (int i=1; i<=4; i++) {
      A0 ^= ff2k_mult(f[i], powal);
      powal = ff2k_mult(powal, alpha);
    }
    ff2k_t A2 = f[2] ^ ff2k_mult(f[3], alpha);

    // Compute x^4 * f(1/x + alpha).
    f[0] = f[4];
    f[1] = f[3];
    f[2] = A2;
    f[3] = 0;
    f[4] = A0;

    if (A0 == 0) {
      // The polynomial A0 + A2 x^2 + f[3] x^3 + f[4] x^4 has alpha as a root.
      //
      // We count out the multiplicity and use generic functionality on what remains.
      if (A2 == 0)
        return 2; // NOTE: f[3] != 0 since b != 0. Thus alpha is a triple root.
      else
        return (1 + count_quadratic_roots(f)); // alpha is a double root.
    }
  }
  
  // The key trick is to quickly compute the GCD with x^q-x, then check
  // the degree.

  make_monic_nodiv(f, deg);
  
  // q = 2^n.
  // Compute x^q mod f(x). 
  ff2k_t powx[3];
  ff2k_t C = 0;
  
  // Initialize a couple steps into the loop.
  if (n == 1) {
    powx[0] = 0;
    powx[1] = 0;
    powx[2] = 1;
  } else {
    powx[0] = f[0];
    powx[1] = f[1];
    powx[2] = f[2];
  }
  
  // Square-Reduce loop
  for (int i=2; i < n; i++) {

    // Square, in characteristic 2.
    C = ff2k_square(powx[2]);
    powx[0] = ff2k_square(powx[0]);
    powx[2] = ff2k_square(powx[1]);

    // Reduce
    powx[0] ^= ff2k_mult(C, f[0]);
    powx[1]  = ff2k_mult(C, f[1]); // Note: squaring made powx[1] = 0. (Theoretically)
    powx[2] ^= ff2k_mult(C, f[2]);    
  }

  // Subtract x
  powx[1] ^= 1;

  // Compute the GCD degree normally from this point.
  return gcd_degree(f, powx, deg, deg-2);
}


int count_poly_roots(unsigned* poly, int deg) {
  // Probably should put this in the header file, so that the compiler can
  // see that the degree is bounded by 4. Maybe even template this...

  if (deg==1) {
    if ((poly[deg] == 0) and (poly[0] == 0))
      return q; // Identically zero
    else
      return (poly[deg] == 0);
  }
  
  // The key trick is to quickly compute the GCD with x^q-x, then check
  // the degree.

  if (poly[deg] == 0) return count_poly_roots(poly, deg-1);
  if (deg == 4) return count_quartic_roots(poly); // Optimized routine for quartics.

  
  // Declare f. Note this mutates poly.
  const unsigned n = FINITEFIELDBITSIZE;
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


////////////////////////////////////////////////////////////////////////////////
//
// Polynomial root finding.
//
////////////////////////////////////////////////////////////////////////////////

unsigned* quadratic_roots(unsigned* f, unsigned* roots) {
  // Return an array {n, r1, r2}. The first number indicates the number
  // of distinct roots. The remaining values are the actual roots.
  //
  // It is assumed that f is genuinely a quadratic.
  const int n = FINITEFIELDBITSIZE;
  
  // ax^2 + bx + c.
  unsigned c = f[0];
  unsigned b = f[1];
  unsigned a = f[2];

  if (a==0) {
    throw std::invalid_argument("Leading coefficient must be non-zero.");
  }
  
  if (b==0) {
    // The quadratic has a double root.
    roots[0] = 1;
    roots[1] = ff2k_sqrt(ff2k_divi(c, a));
    return roots;
  }

  // First, compute the roots of z^2 + z + ac/b^2.
  // For this, we use linear algebra.
  
  unsigned acbb = ff2k_divi(ff2k_mult(a, c), ff2k_square(b));
  unsigned pivot = (1 << (n-1));
  unsigned rt = 0;
  const unsigned* trace_basis = trace_bases[n];
  const unsigned* pretrace_basis = pretrace_bases[n];
  
  for (int i = 0; i < n-1; i++) {
    if (acbb & trace_basis[i]) {
       rt ^= pretrace_basis[i];
       acbb ^= trace_basis[i];
    }    
    pivot >>= 1;
  }

  // The only thing left at the end of the loop should be the trace.
  if (acbb != 0)
    return roots;
  else {
    // z^2 + z + ac/b^2 =>  ax = bz
    unsigned bdiva = ff2k_divi(b, a);
    rt = ff2k_mult(rt, bdiva);
    roots[0] = 2;
    roots[1] = rt;
    roots[2] = rt ^ bdiva;
    return roots;
  }
}
