#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <tuple>
#include <assert.h>
//#include <omp.h>
#include "tableio.h"
#include "Fq_tables.h"

#ifndef COEFFSFILE
#include "coeffs.h"
#else
#include COEFFSFILE
#endif

// NOTE: We expect to see a -N flag at compile time.
#ifndef N
throw std::invalid_argument("Error: No N provided at compile time.");
#endif


// function prototypes
int contribution_at_P3_point(unsigned, unsigned, unsigned, unsigned);
int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned);


int main(int argc, char **argv) {

  const char argument_error[] =
    "No argument should be provided for n. Instead, compile with the -N flag.";

  if (argc != 1) {
    throw std::invalid_argument(argument_error);
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
  // In other words, that the Arf invariant is *constant along fibres*!
  // 
  // This gives Delta the structure of a generically 4:1 cover of P2.
  //
  // Thus, in order to count points, it suffices to enumerate over the points in P2, compute the
  // Arf invariant (as a function of the base), and count the number of rational points in the
  // fibre. 
  // See contribution_at_P3_point in order to determine the correct contribution.
  
  // The contribution to the point count from the distinguished singularity.
  long diff = contribution_at_P3_point(0,0,0,1);

  
  // We next look at the open subscheme Delta - (0:0:0:1). 
  // We enumerate over the points in P2. Because the equation for X is defined over F2, we may
  // speed up the computation by using the natural symmetry provided by Galois action.
  // 
  // We use the standard decomposition of P2 = {pt} + A1 + A2.

  
  // The contribution from the fibre over (1:0:0).
  diff += (long)contribution_of_fibre_over_P2_point(1, 0, 0);

  // The contribution over the hyperplane at infinity.
  for (unsigned y_2 = 0; y_2 < q; y_2++)
   if (orbit_rep[y_2] == y_2)
     diff += (long)contribution_of_fibre_over_P2_point(y_2, 0, 1) * (long)orbit_size[y_2];
  
  // The contribution from the locus {B != 0}. 
  for (unsigned y_1 = 0; y_1 < q; y_1++)
    if (orbit_rep[y_1] == y_1)
      for (unsigned y_2 = 0; y_2 < q; y_2++) {
        diff += (long)contribution_of_fibre_over_P2_point(y_1, 1, y_2) * (long)orbit_size[y_1];
      }
  
  std::cerr << "Iterate time: " << (std::clock() - cputime) * 1./CLOCKS_PER_SEC
            << std::endl;

  
  //////////////////////////////////////
  // Final calculation of the formula.
  __int128 Q = q; // to avoid overflows
  std::string output = int128_to_string(Q*Q*Q*Q + Q*Q*Q + Q*diff + Q + 1);
  std::cout << output  << std::endl; // Send output to magma pipe.
  return 0;
}


int contribution_of_fibre_over_P2_point(unsigned y_0, unsigned y_1, unsigned y_2) {

  // Magma will give us a fibration where A,B,C = y0,y1,y2. In particular, though
  // the degree is higher, the Arf invariant is *constant* along the fibres. Thus,
  // we need only determine the number of roots lying over a given point.

  unsigned abcde; // coefficients of 4:1 cover, defined in coeffs.h
  unsigned f[] = {e,d,c,b,a};

  // TODO: If B is identically zero, then we have an efficiency issue.
  // Otherwise, we are OK.
  
  if (y_1 == 1) {
    int arf = Arf_invariant_mu2_b_equals_1(y_0, y_2);
    return arf * count_poly_roots(f, 4);
  }

  if (y_1 != 0) throw std::invalid_argument("Function called with y_1 not either 0 or 1.");
    
  // Once y_1 = 0, theory guarantees that b = d = 0. Sanity check.
  assert((b==0) && (d==0));

  // Our quadratic is now a(z^2)^2 + c(z^2) + e. We can solve for the roots.
  if (a != 0) {
    unsigned g[] = {e, c, a};
    unsigned rts[3] = {0,0,0};
    quadratic_roots(g, rts);

    int ret = 0;
    for (int i=1; i <= rts[0]; i++) {
      ret += contribution_at_P3_point(y_0, y_1, y_2, ff2k_sqrt(rts[i]));
    }
    
    // free(rts); // Cleanup.
    return ret;

  } else if (c != 0) {
    // Exactly one root.
    return contribution_at_P3_point(y_0, y_1, y_2, ff2k_sqrt(ff2k_divi(e, c)));
    
  } else if (e != 0) {
    // Nothing is a root.
    return 0;
    
  } else {
    // Everything is a root.
    int ret = 0;
    for (int i=0; i<q; i++) {
      ret += contribution_at_P3_point(y_0, y_1, y_2, i);
    }
    return ret;
  }

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
  if (A == 0 && B == 0 && C == 0)
    return 1; // The conic is indeed split.
  if (B != 0)
    return Arf_invariant_mu2(A, B, C);
  if (D != 0)
    return Arf_invariant_mu2(A, D, F);
  if (E != 0)
    return Arf_invariant_mu2(C, E, F);

  
  // We should never run into this case.
  throw std::runtime_error("Conic is not singular.");
}



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
*/
