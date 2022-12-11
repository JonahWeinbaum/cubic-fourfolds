#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <assert.h>
#include "tableio.h"
#include "Fq_tables.h"

#ifndef COEFFSFILE
#include "coeffs.h"
#else
#include COEFFSFILE
#endif

// we're in F_q. On a modern machine we can assume 32 or 64 bit integers.
unsigned q;
//const unsigned PLACEHOLDER_Fq_elt = 1 << 16;
unsigned NULL_Fq_elt;

// lookup tables
unsigned **mult, **divi,

  ***quadratic_roots,
  // quadratic_roots[i][a][b] is
  // the i'th root of x^2 + ax + b, 
  // or q if we're out of roots

  ***depressed_cubic_roots;
  // depressed_cubic_roots[i][s][t] is
  // the i'th root of x^3 + sx + t,
  // or q if we're out of roots



inline unsigned ff2k_mult(unsigned a, unsigned b) {
  return mult[a][b];
}

inline unsigned ff2k_divi(unsigned a, unsigned b) {
  return divi[a][b];
}


// function prototypes
int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned);
int contribution_at_P3_point(unsigned, unsigned, unsigned, unsigned);


int main(int argc, char **argv) {

  const char argument_error[] = "Argument should be n, for F_{2^n}, where 1 <= n <= 15.\n";
  if (argc != 2) {
    std::cerr << argument_error << std::endl;
    return 1;
  }
  int n = atoi(argv[1]);
  if (n < 1 || n > 15) {
    std::cerr << argument_error << std::endl;
    return 1;
  }
  
  q = 1<<n;
  std::string qq = std::to_string(q);

  // Placeholder value for q. (Global variable)
  NULL_Fq_elt = q;
  
  #ifdef WITHCACHE
  
  // Loading of LOCAL variables
  int* orbit_size = read_table(q, "orbit_size_" + qq, 0);
  unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);
    
  // Loading of GLOBAL variables
  mult = read_table(q, q, "mult_" + qq);
  divi = read_table(q, q, "divi_" + qq);
  
  quadratic_roots = read_table(q, q, 2, "quadratic_roots_" + qq);  
  depressed_cubic_roots = read_table(q, q, 3, "depressed_cubic_roots_" + qq);  

    
  #else

  // NOTE: These variables are declared as globals above.
  // unsigned*** depressed_cubic_roots;
  // unsigned*** quadratic_roots;
  // unsigned** mult;
  // unsigned** divi;
  
  unsigned* orbit_rep;
  int* orbit_size;
  std::tie(depressed_cubic_roots,
           quadratic_roots,
           mult, divi, orbit_rep, orbit_size) = generate_Fq_tables(n);
  
  #endif

  
  //////////////////////////////////////////
  // Main calculation
  // 
  // According to Section~3 of Addington-Auel, it suffices to compute the number of points
  // on the cubic surface X via counting the quantity
  //
  //     diff = #double cover of Delta - #Delta
  //
  // where Delta is the discriminantal hypersurface for the conic fibration.
  // The precise formula is
  //
  //     #X(Fq) = q^4 + q^3 + q * diff + q + 1.
  //
  // In order to compute diff, we need to find all of the rational points on Delta. The
  // equation fed to the code ensures there is a singularity on Delta at (0:0:0:1). Projecting
  // away from this point gives Delta the structure of a generically 3:1 cover of P2.
  //
  // Thus, in order to count points, it suffices to enumerate over the points in P2, compute the
  // set of rational points in the fibre, and then check a criterion on the associated conic
  // (see contribution_at_P3_point) in order to determine the correct contribution.
  
  // The contribution to the point count from the distinguished singularity.
  int diff = contribution_at_P3_point(0,0,0,1);

  // We next look at the open subscheme Delta - (0:0:0:1). 
  // We enumerate over the points in P2. Because the equation for X is defined over F2, we may
  // speed up the computation by using the natural symmetry provided by Galois action.
  // 
  // We use the standard decomposition of P2 = {pt} + A1 + A2.

  // The contribution from the fibre over (0:0:1).
  diff += contribution_of_fibre_over_P2_point(0, 0, 1);

  // The contribution over the hyperplane at infinity.
  for (unsigned y_2 = 0; y_2 < q; y_2++)
    if (orbit_rep[y_2] == y_2)
      diff += contribution_of_fibre_over_P2_point(0, 1, y_2) * orbit_size[y_2];
  
  // The contribution from the A2 part. 
  for (unsigned y_1 = 0; y_1 < q; y_1++)
    if (orbit_rep[y_1] == y_1)
      for (unsigned y_2 = 0; y_2 < q; y_2++)
        diff += contribution_of_fibre_over_P2_point(1, y_1, y_2) * orbit_size[y_1];

  
  //////////////////////////////////////
  // Final calculation of the formula.
  long Q = q; // to avoid overflows
  printf("%ld\n", Q*Q*Q*Q + Q*Q*Q + Q*diff + Q + 1); // Send output to magma pipe.
  return 0;
}

// contribution_of_fibre_over_P2_point
//
// iterate over the three sheets of the quintic,
// i.e. the three roots of a y_3^3 + b y_3^2 + c y_3 + d
//
// Use Cardano's formula to determine the roots of the cubic.
// Then return the number of points on the singular conic associated to that root.
//
int contribution_of_fibre_over_P2_point(unsigned y_0, unsigned y_1, unsigned y_2) {
  unsigned abcde; // coefficients of 3:1 cover, defined in coeffs.h

  assert(a==0); // Ensure the cover is actually degree 3, for this version of the code.

  // For ease of reading.
  unsigned poly[] = {e,d,c,b,a};
  
  // iterate over the three sheets of the quintic,
  // i.e. the three roots of a y_3^3 + b y_3^2 + c y_3 + d

  int ret = 0;

  if (poly[3] != 0) {
    // The cubic is in fact a cubic.

    // Make the cubic monic.
    poly[2] = ff2k_divi(poly[2], poly[3]);
    poly[1] = ff2k_divi(poly[1], poly[3]);
    poly[0] = ff2k_divi(poly[0], poly[3]);

    // Make the cubic depressed.
    unsigned s = ff2k_square(poly[2]) ^ poly[1];
    unsigned t = ff2k_mult(poly[2], poly[1]) ^ poly[0];
    
    for (int i = 0; i < 3; i++) {
      unsigned r = depressed_cubic_roots[s][t][i];
      if (r == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, r ^ poly[2]);
    }
    return ret;
    
  } else if (poly[2] != 0) {
    // The leading coefficient is zero, so there are only two roots.
    // This corresponds to the case the projecting line is tangent at the
    // distinguished point (0:0:0:1).

    // Make it monic.
    poly[1] = ff2k_divi(poly[1], poly[2]); poly[0] = ff2k_divi(poly[0], poly[2]);
    
    for (int i = 0; i < 2; i++) {
      unsigned y_3 = quadratic_roots[poly[1]][poly[0]][i];
      if (y_3 == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, y_3);
    }
    return ret;
    
  } else if (poly[1] != 0) {
    // The cubic degenerates to a linear polynomial.
    return contribution_at_P3_point(y_0, y_1, y_2, ff2k_divi(poly[0], poly[1]));    
    
  } else if (poly[0] == 0) {
    // The cubic degenerates to the identically zero polynomial.
    // That is, the fibre over the input point is an entire line.

    // We loop over points on said line. (Except the distinguished points (0:0:0:1).)
    for (unsigned y_3 = 0; y_3 < q; y_3++)
      ret += contribution_at_P3_point(y_0, y_1, y_2, y_3);

    return ret;
    
  } else {
    // If we avoid all the clauses above, then 
    // a = b = c = 0 and d != 0, so there are no roots and thus no contribution.
    return 0;
  }
}

bool Arf_invariant(unsigned X, unsigned Y, unsigned Z) {
  unsigned XY = ff2k_mult(X, Y);
  unsigned ZZ = ff2k_mult(Z, Z);
  unsigned XYdivZZ = ff2k_divi(XY, ZZ);
  return quadratic_roots[1][XYdivZZ][0] == NULL_Fq_elt;
}

int Arf_invariant_mu2(unsigned X, unsigned Y, unsigned Z) {
  return Arf_invariant(X, Y, Z) == 1 ? -1 : 1;
}

// contribution_at_P3_point (Should be named "points on parametrized conic" or something.)
//
// Determines how many points are on the double cover
// associated to the singular conic
//
//     Ax^2 + Bxy + Cy^2 + Dxz + Eyz + Fz^2 = 0,
//
// minus 1 point for the quintic downstairs
//
int contribution_at_P3_point(unsigned y_0, unsigned y_1, unsigned y_2, unsigned y_3) {
  unsigned ABCDEF; // coefficients of the conic, defined in coeffs.h

  // The conic has rank 0. It is an entire copy of P2.
  // Note that (q^2 + q + 1) - 1 = q(q+1)
  // TODO: Explain why q is the right number.
  if (A == 0 && B == 0 && C == 0 && D == 0 && E == 0 && F == 0)
    return q;
  
  // The conic is rank 1. It is a double line in characteristic 2.
  if (B == 0 && D == 0 && E == 0)
    return 0;

  // The conic is rank 2. We look at Arf invariants.
  if (B != 0)
    return Arf_invariant_mu2(A, C, B);
  if (D != 0)
    return Arf_invariant_mu2(A, F, D);
  if (E != 0)
    return Arf_invariant_mu2(C, F, E);

  
  // We should never run into this case.
  throw std::runtime_error("Conic is not singular.");
}
