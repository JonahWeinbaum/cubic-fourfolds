#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <tuple>
#include "tableio.h"
#include "Fq_tables.h"

#ifndef COEFFSFILE
#include "coeffs.h"
#else
#include COEFFSFILE
#endif

// NOTE: We expect to see a -N flag at compile time.
#ifndef N
println("Error: No N provided at compile time.");
return 1;
#endif

// TODO: Update with Jonah's multiplication function, if it does OK.

// function prototypes
int contribution_at_P3_point(unsigned, unsigned, unsigned, unsigned);
int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned);
int Arf_invariant(unsigned, unsigned, unsigned);
int Arf_invariant_mu2(unsigned, unsigned, unsigned);


// TODO: Perhaps move to Fq.h
unsigned* quadratic_roots(unsigned*);

int main(int argc, char **argv) {

  const char argument_error[] =
    "No argument should be provided for n. Instead, compile with the -N flag.";

  if (argc != 1) {
    std::cerr << argument_error << std::endl;
    return 1;
  }
  
  std::string qq = std::to_string(q);

  // Start timer
  std::clock_t cputime = std::clock();

  // NOTE: I still need to generate these tables. They'll fit in memory.
  // Loading of LOCAL variables
  unsigned* orbit_rep;
  unsigned* orbit_size;
  std::tie(orbit_rep, orbit_size) = generate_orbit_tables(N);
    
  std::cerr << "Orbit table generation: " << (std::clock() - cputime) * 1./CLOCKS_PER_SEC
            << std::endl;
  cputime = std::clock();
  
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
  // In order to compute diff, we project away from a particular point on Delta (assumed to
  // be at (0:0:0:1)). This projection is further assumed to be given by (A : B : C).
  // In other words, that the Arf in variant is *constant along fibres*!
  // 
  // This gives Delta the structure of a generically 4:1 cover of P2.
  //
  // Thus, in order to count points, it suffices to enumerate over the points in P2, compute the
  // Arf invariant (as a function of the base), and count the number of rational points in the
  // fibre. 
  // See contribution_at_P3_point in order to determine the correct contribution.
  
  // The contribution to the point count from the distinguished singularity.
  int diff = contribution_at_P3_point(0,0,0,1);

  
  // We next look at the open subscheme Delta - (0:0:0:1). 
  // We enumerate over the points in P2. Because the equation for X is defined over F2, we may
  // speed up the computation by using the natural symmetry provided by Galois action.
  // 
  // We use the standard decomposition of P2 = {pt} + A1 + A2.

  
  // The contribution from the fibre over (1:0:0).
  diff += contribution_of_fibre_over_P2_point(1, 0, 0);

  // The contribution over the hyperplane at infinity.
  for (unsigned y_2 = 0; y_2 < q; y_2++)
   if (orbit_rep[y_2] == y_2)
     diff += contribution_of_fibre_over_P2_point(y_2, 0, 1) * orbit_size[y_2];
  
  // The contribution from the locus {A != 0, B != 0}. 
  // TODO: Artin-Schrier trickery should speed this up a lot...
  //       Main issue is how to make this compatible with conjugates...
  for (unsigned y_1 = 0; y_1 < q; y_1++)
    if (orbit_rep[y_1] == y_1)
      for (unsigned y_2 = 0; y_2 < q; y_2++) {
        diff += contribution_of_fibre_over_P2_point(y_1, 1, y_2) * orbit_size[y_1];
      }

  // The contribution from the locus {A = 0, B != 0}.
  for (unsigned y_2 = 0; y_2 < q; y_2++)
   if (orbit_rep[y_2] == y_2)
     diff += contribution_of_fibre_over_P2_point(y_2, 0, 1) * orbit_size[y_2];

  
  std::cerr << "Iterate time: " << (std::clock() - cputime) * 1./CLOCKS_PER_SEC
            << std::endl;

  
  //////////////////////////////////////
  // Final calculation of the formula.
  long long Q = q; // to avoid overflows
  std::cout <<  Q*Q*Q*Q + Q*Q*Q + Q*diff + Q + 1 << std::endl; // Send output to magma pipe.
  return 0;
}

// TODO: This only works over the A2 patch. 
int contribution_of_fibre_over_P2_point(unsigned y_0, unsigned y_1, unsigned y_2) {

  // Magma will give us a fibration where A,B,C = y0,y1,y2. In particular, though
  // the degree is higher, the Arf invariant is *constant* along the fibres. Thus,
  // we need only determine the number of roots lying over a given point.

  // TODO: We should only ever call this function when the arf invariant is 1.
  int arf=0;
  if (y_1 == 1) {
    arf = Arf_invariant_mu2(y_0, y_1, y_2);
  } else {
    throw std::runtime_error("Not implemented.");
  }
  
  unsigned abcde; // coefficients of 4:1 cover, defined in coeffs.h
  unsigned f[] = {e,d,c,b,a};

  return arf * count_poly_roots(f, 4);
}

// contribution_of_fibre_over_point_at_infinity() {
//   1. Determine A, C, D, E.
//   2. Compute the roots of the right quadratic equation. (Delta=0, B=0 causes factorization).
//   3. Compute the Arf invariant of the resulting expression.
// }

/*
// contribution_of_fibre_over_P2_point
//
// iterate over the three sheets of the quintic,
// i.e. the three roots of a y_3^3 + b y_3^2 + c y_3 + d
//
// Use Cardano's formula to determine the roots of the cubic.
// Then return the number of points on the singular conic associated to that root.
//
int contribution_of_fibre_over_P2_point(unsigned y_0, unsigned y_1, unsigned y_2) {
  unsigned abcd; // coefficients of 3:1 cover, defined in coeffs.h

  // iterate over the three sheets of the quintic,
  // i.e. the three roots of a y_3^3 + b y_3^2 + c y_3 + d

  int ret = 0;

  if (a != 0) {
    // The cubic is in fact a cubic.

    // Make the cubic monic.
    b = ff2k_divi(b, a); c = ff2k_divi(c, a); d = ff2k_divi(d, a);

    // Make the cubic depressed.
    unsigned s = ff2k_mult(b, b) ^ c, t = ff2k_mult(b, c) ^ d;
    
    for (int i = 0; i < 3; i++) {
      unsigned r = depressed_cubic_roots[s][t][i];
      if (r == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, r ^ b);
    }
    return ret;
    
  } else if (b != 0) {
    // The leading coefficient is zero, so there are only two roots.
    // This corresponds to the case the projecting line is tangent at the
    // distinguished point (0:0:0:1).

    // Make it monic.
    c = ff2k_divi(c, b); d = ff2k_divi(d, b);
    
    for (int i = 0; i < 2; i++) {
      unsigned y_3 = quadratic_roots[c][d][i];
      if (y_3 == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, y_3);
    }
    return ret;
    
  } else if (c != 0) {
    // The cubic degenerates to a linear polynomial.
    return contribution_at_P3_point(y_0, y_1, y_2, ff2k_divi(d, c));    
    
  } else if (d == 0) {
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

*/

// Arf_invariant(unsigned b, unsigned c)
// 
// Given a quadratic polynomial ax^2 + bx + c with b != 0, determine the
// Arf invariant (Equal to 0 or 1).
int Arf_invariant(unsigned a, unsigned b, unsigned c) {

  const int n = FINITEFIELDBITSIZE;
  
  // Use the well-known fact that Arf(a,b,c) = tr(ac/b^2).
  unsigned binv2 = ff2k_square(ff2k_inv(b));
  unsigned acbb  = ff2k_mult(ff2k_mult(a, c), binv2);
  unsigned pow   = acbb;
  unsigned tr    = 0;
  
  for (int i=0; i < n; i++) {
    tr ^= pow;
    pow = ff2k_square(pow);
  }
  return (int)tr;  
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

  // The conic has rank 0, so the correct difference factor is q.
  if (A == 0 && B == 0 && C == 0 && D == 0 && E == 0 && F == 0)
    return q;

  // The conic is rank 1. It is a double line in characteristic 2.
  if (B == 0 && D == 0 && E == 0)
    return 0;

  // The conic is rank 2. We look at Arf invariants.
  if (A == 0 && B == 0 && C == 0) {
    return 1; // The conic is indeed split.
  }

  // We should never run into this case.
  throw std::runtime_error("Not implemented past this point.");

  /*
  if (B != 0)
    return Arf_invariant_mu2(A, C, B);
  if (D != 0)
    return Arf_invariant_mu2(A, F, D);
  if (E != 0)
    return Arf_invariant_mu2(C, F, E);
  */
  
  // We should never run into this case.
  throw std::runtime_error("Conic is not singular.");
}


unsigned* quadratic_roots(unsigned* f) {
  // Return an array {n, r1, r2}. The first number indicates the number
  // of distinct roots. The remaining values are the actual roots.
  //
  // It is assumed that f is genuinely a quadratic.
  
  // ax^2 + bx + c.
  unsigned c = f[0];
  unsigned b = f[1];
  unsigned a = f[2];

  unsigned* roots = (unsigned*) malloc(3*sizeof(unsigned));

  for (int i=0; i<3; i++) {
    roots[i] = 0;
  }

  if (a==0) {
    throw std::runtime_error("Function called with bad arguments.");
  }
  
  if (b==0) {
    // The quadratic has a double root.
    roots[0] = 1;
    roots[1] = ff2k_divi(c, a);
    return roots;
  }

  // First, compute the roots of z^2 + z + ac/b^2.
  // For this, we use linear algebra.
  
  unsigned acbb = ff2k_divi(ff2k_mult(a, c), ff2k_square(b));


  
  return roots;
}
