#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "tableio.h"
#include "coeffs.h"

// we're in F_q
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
  (1<<15) + (1<<1) + 1, };
/* with n > 15 you'll just run out of memory
  (1<<16) + (1<<5) + (1<<3) + (1<<1) + 1,
  (1<<17) + (1<<3) + 1,
  (1<<18) + (1<<3) + 1,
  (1<<19) + (1<<5) + (1<<2) + (1<<1) + 1,
  (1<<20) + (1<<3) + 1,
  (1<<21) + (1<<2) + 1,
  (1<<22) + (1<<1) + 1 }; */

// function prototypes
int contribution_of_fibre_over_P2_point(unsigned, unsigned, unsigned);
int contribution_at_P3_point(unsigned, unsigned, unsigned, unsigned);
void print_it(unsigned);

int main(int argc, char **argv) {
  std::clock_t cputime = std::clock();
  const char syntax_error[] = "Argument should be n, for F_{2^n}, where 1 <= n <= 15.\n";
  if (argc != 2) { printf(syntax_error); return 1; }
  int n = atoi(argv[1]);
  if (n < 1 || n > 15) { printf(syntax_error); return 1; }

  q = 1<<n;
  std::string qq = std::to_string(q);

  // TODO: We are temporarily distinguishing the placeholder value from q. In order to not
  // introduce new bugs, we are keeping the values the same for the time being.
  NULL_Fq_elt = q;
  
  // Polynomial set.
  unsigned p = polynomials[n];
  

  #ifdef WITHCACHE 
  // Loading of LOCAL variables
  int* orbit_size = read_table(q, "orbit_size_" + qq, 0);
  unsigned* orbit_rep = read_table(q, "orbit_rep_" + qq);
    
  // Loading of GLOBAL variables
  mult = read_table(q, q, "mult_" + qq);
  divi = read_table(q, q, "divi_" + qq);
  
  quadratic_roots = read_table(2, q, q, "quadratic_roots_" + qq);  
  depressed_cubic_roots = read_table(3, q, q, "depressed_cubic_roots_" + qq);  

    
  #else
  // lookup tables for multiplication and division
  // allocate
  mult = new unsigned*[q];
  divi = new unsigned*[q];
  for (unsigned i = 0; i < q; i++) {
    mult[i] = new unsigned[q];
    divi[i] = new unsigned[q];
  }

  // fill the tables
  for (unsigned i = 0; i < q; i++)
    for (unsigned j = 0; j <= i; j++) {
      // main multiplication algorithm
      unsigned a = i, b = j, ij = 0;
      while (b != 0) {
	// b = const + higher-deg part
	// if const = 1 then ret += a
	if (b & 1) ij ^= a;
	// kill const, divide higher-deg part by x
	b >>= 1;
	// multiply a by x
	a <<= 1;
	if (a & q) a ^= p;
      }

      mult[i][j] = mult[j][i] = ij;
      divi[ij][i] = j;
      divi[ij][j] = i;
    }
  
  // lookup table for roots of quadratics
  // allocate and initialize
  quadratic_roots = new unsigned**[2];
  for (int i = 0; i < 2; i++) {
    quadratic_roots[i] = new unsigned*[q];
    for (int j = 0; j < q; j++) {
      quadratic_roots[i][j] = new unsigned[q];
      for (int k = 0; k < q; k++)
        quadratic_roots[i][j][k] = NULL_Fq_elt;
    }
  }
  
  // fill the table
  for (int i = 0; i < q; i++)
    for (int j = i; j < q; j++) {
      // x^2 + ax + b = (x+i)(x+j)
      unsigned a = i ^ j, b = mult[i][j];
      quadratic_roots[0][a][b] = i;
      if (j > i)
        quadratic_roots[1][a][b] = j;
    }

  
  // lookup table for roots of depressed cubics
  // allocate and initialize
  depressed_cubic_roots = new unsigned**[3];
  for (int i = 0; i < 3; i++) {
    depressed_cubic_roots[i] = new unsigned*[q];
    for (int j = 0; j < q; j++) {
      depressed_cubic_roots[i][j] = new unsigned[q];
      for (int k = 0; k < q; k++)
        depressed_cubic_roots[i][j][k] = NULL_Fq_elt;
    }
  }

  // Hark! This is at least one bug. Roots of x^3 not counted properly.
  // fill the table
  for (int a = 0; a < q; a++)
    for (int b = 0; b < q; b++) {
      // (x^2 + ax + b)(x+a) = x^3 + sx + t
      unsigned s = mult[a][a] ^ b, t = mult[a][b];
      if (b == 0) { // if a is already a root of x^2 + ax + b
        depressed_cubic_roots[0][s][t] = quadratic_roots[0][a][b];
        depressed_cubic_roots[1][s][t] = quadratic_roots[1][a][b];
      } else {
        depressed_cubic_roots[0][s][t] = a;
        depressed_cubic_roots[1][s][t] = quadratic_roots[0][a][b];
        depressed_cubic_roots[2][s][t] = quadratic_roots[1][a][b];
      }
    }
  
  
  // Galois orbits
  // allocate and initialize
  unsigned *orbit_rep = new unsigned[q];
  for (unsigned i = 0; i < q; i++)
    orbit_rep[i] = NULL_Fq_elt;
  
  int *orbit_size = new int[q]; // only valid if i == orbit_rep[i]

  
  // fill the tables
  // for (unsigned i = 0; i < q; i++) {
  //   if (orbit_rep[i] != NULL_Fq_elt) continue;
    
  //   unsigned j = i;
  //   int size = 1;
  //   while (true) {
  //     orbit_rep[j] = i;
  //     j = mult[j][j];   // Apply Frobenius x->x^2.
  //     if (j == i) {
  //       orbit_size[i] = size;
  //       break;
  //     }
  //     size++;
  //   } 
  // }

  for (unsigned i = 0; i < q; i++) {
    if (orbit_rep[i] != NULL_Fq_elt) continue;

    // Compute the orbit under Gal(F2bar/F2) Frobenius.
    // Once we see the same element twice, we've completed the orbit.

    unsigned start = i;
    orbit_rep[start] = i;

    unsigned mark = mult[start][start];
    int size = 1;

    while (mark != start) {
      orbit_rep[mark] = start;
      mark = mult[mark][mark];
      size++;
    }
    // Set the size of the orbit.
    orbit_size[i] = size;
  }
  #endif

  #ifdef COMPARE
  unsigned*** X = read_table(2, q, q, "quadratic_roots_" + qq);  
  unsigned*** Y = read_table(3, q, q, "depressed_cubic_roots_" + qq);  

  compare_3_table(X, quadratic_roots, 2, q, q);
  compare_3_table(Y, depressed_cubic_roots, 3, q, q);
  #endif
  
  cputime = std::clock() - cputime;
  std::cerr << "precomputation took " << cputime * 1./ CLOCKS_PER_SEC << std::endl;

  //////////////////////////////////////////
  // Main calculation
  // 
  // According to Section~3 of Addington-Auel, it suffices to compute the number of points
  // on the cubic surface X via counting the quantity
  //
  //     diff = #double cover of Delta - #Delta
  //
  // where Delta is the discriminantal hypersurface for the (flat) conic bundle.
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
  unsigned abcd; // coefficients of 3:1 cover, defined in coeffs.h

  // iterate over the three sheets of the quintic,
  // i.e. the three roots of a y_3^3 + b y_3^2 + c y_3 + d

  int ret = 0;

  if (a != 0) {
    // The cubic is in fact a cubic.

    // Make the cubic monic.
    b = divi[b][a]; c = divi[c][a]; d = divi[d][a];

    // Make the cubic depressed.
    unsigned s = mult[b][b] ^ c, t = mult[b][c] ^ d;
    
    for (int i = 0; i < 3; i++) {
      unsigned r = depressed_cubic_roots[i][s][t];
      if (r == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, r ^ b);
    }
    return ret;
    
  } else if (b != 0) {
    // The leading coefficient is zero, so there are only two roots.
    // This corresponds to the case the projecting line is tangent at the
    // distinguished point (0:0:0:1).

    // Make it monic.
    c = divi[c][b]; d = divi[d][b];
    
    for (int i = 0; i < 2; i++) {
      unsigned y_3 = quadratic_roots[i][c][d];
      if (y_3 == NULL_Fq_elt) break;
      ret += contribution_at_P3_point(y_0, y_1, y_2, y_3);
    }
    return ret;
    
  } else if (c != 0) {
    // The cubic degenerates to a linear polynomial.
    return contribution_at_P3_point(y_0, y_1, y_2, divi[d][c]);    
    
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
  unsigned XY = mult[X][Y];
  unsigned ZZ = mult[Z][Z];
  unsigned XYdivZZ = divi[XY][ZZ];
  return quadratic_roots[0][1][XYdivZZ] == NULL_Fq_elt;
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

  // TODO: We need to account for the case where the conic has rank 0.
  if (A == 0 && B == 0 && C == 0 && D == 0 && E == 0 & F == 0)
    throw std::runtime_error("Not implemented if the conic bundle is not flat.");
  
  // The conic is rank 1. It is a double line in characteristic 2.
  if (B == 0 && D == 0 && E == 0)
    return 0;

  // The conic is rank 2. We look at Arf invariants.
  // if (B != 0)
  //   return quadratic_roots[0][1][divi[mult[A][C]][mult[B][B]]] == NULL_Fq_elt ? -1 : 1;
  // if (D != 0)
  //   return quadratic_roots[0][1][divi[mult[A][F]][mult[D][D]]] == NULL_Fq_elt ? -1 : 1;
  // if (E != 0)
  //   return quadratic_roots[0][1][divi[mult[C][F]][mult[E][E]]] == NULL_Fq_elt ? -1 : 1;

  if (B != 0)
    return Arf_invariant_mu2(A, C, B);
  if (D != 0)
    return Arf_invariant_mu2(A, F, D);
  if (E != 0)
    return Arf_invariant_mu2(C, F, E);

  
  // We should never run into this case.
  throw std::runtime_error("Conic is not singular.");
  return -2;
}


// just for debugging
void print_it(unsigned a) {
  if (a == 0) {
    printf("0\n");
    return;
  }

  int exp = 0;
  int first_term = 1;
  while (a) {
    if (a & 1) {
      if (first_term)
        first_term = 0;
      else
        printf(" + ");
      if (exp == 0)
        printf("1");
      else if (exp == 1)
        printf("x");
      else
        printf("x^%d", exp);
    }
    a >>= 1;
    exp++;
  }
  printf("\n");
}
