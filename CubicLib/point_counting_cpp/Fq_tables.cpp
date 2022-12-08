#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <tuple>

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

// Use a compile flag to unlock functionality.
#ifdef N
#include "Fq.h"
#else
// Placeholder function
static unsigned mult_ff2k(unsigned a, unsigned b) {
  return 0;
}
#endif

std::tuple<unsigned*, unsigned*> generate_orbit_tables(int n) {

  unsigned q = (1<<n);
  
  // Galois orbits
  // allocate and initialize
  unsigned *orbit_rep = new unsigned[q];
  for (unsigned i = 0; i < q; i++)
    orbit_rep[i] = q; // means null
  
  unsigned *orbit_size = new unsigned[q]; // only valid if i == orbit_rep[i]

  // fill the tables
  for (unsigned i = 0; i < q; i++) {
    if (orbit_rep[i] != q) continue; // already filled
    unsigned j = i;
    unsigned size = 1;
    while (1) {
      orbit_rep[j] = i;
      j = ff2k_square(j);
      if (j == i) {
        orbit_size[i] = size;
        break;
      }
      size++;
    } 
  }
  
  return std::make_tuple(orbit_rep, orbit_size);
}

using all_table_t = std::tuple<unsigned***, unsigned***,
                               unsigned**, unsigned**,
                               unsigned*, int*>;

all_table_t generate_Fq_tables(int n) {

  // Define q.
  unsigned q = 1<<n;

  // Table variables
  // lookup tables
  unsigned **mult, **divi;

  unsigned ***quadratic_roots;
  // quadratic_roots[i][a][b] is
  // the i'th root of x^2 + ax + b, 
  // or q if we're out of roots

  unsigned ***depressed_cubic_roots;
  // depressed_cubic_roots[i][s][t] is
  // the i'th root of x^3 + sx + t,
  // or q if we're out of roots

  
  // lookup tables for multiplication and division
  // allocate
  mult = new unsigned*[q];
  divi = new unsigned*[q];
  for (unsigned i = 0; i < q; i++) {
    mult[i] = new unsigned[q];
    divi[i] = new unsigned[q];
  }
  
  // fill the tables
  unsigned p = polynomials[n];
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

      mult[i][j] = ij;
      mult[j][i] = ij;
      divi[ij][i] = j;
      divi[ij][j] = i;
    }
  
  // lookup table for roots of quadratics
  // allocate and initialize
  quadratic_roots = new unsigned**[q];
  for (int i = 0; i < q; i++) {
    quadratic_roots[i] = new unsigned*[q];
    for (int j = 0; j < q; j++) {
      quadratic_roots[i][j] = new unsigned[2];
      for (int k = 0; k < 2; k++)
        quadratic_roots[i][j][k] = q;
    }
  }
  // fill the table
  for (int i = 0; i < q; i++)
    for (int j = i; j < q; j++) {
      // x^2 + ax + b = (x+i)(x+j)
      unsigned a = i ^ j, b = mult[i][j];
      quadratic_roots[a][b][0] = i;
      if (j > i)
        quadratic_roots[a][b][1] = j;
    }

  // lookup table for roots of depressed cubics
  // allocate and initialize
  depressed_cubic_roots = new unsigned**[q];
  for (int i = 0; i < q; i++) {
    depressed_cubic_roots[i] = new unsigned*[q];
    for (int j = 0; j < q; j++) {
      depressed_cubic_roots[i][j] = new unsigned[3];
      for (int k = 0; k < 3; k++)
        depressed_cubic_roots[i][j][k] = q;
    }
  }
  // fill the table
  for (int a = 0; a < q; a++)
    for (int b = 0; b < q; b++) {
      // (x^2 + ax + b)(x+a) = x^3 + sx + t
      unsigned s = mult[a][a] ^ b, t = mult[a][b];
      if (b == 0) { // if a is already a root of x^2 + ax + b
        depressed_cubic_roots[s][t][0] = quadratic_roots[a][b][0];
        depressed_cubic_roots[s][t][1] = quadratic_roots[a][b][1];
      } else {
        depressed_cubic_roots[s][t][0] = a;
        depressed_cubic_roots[s][t][1] = quadratic_roots[a][b][0];
        depressed_cubic_roots[s][t][2] = quadratic_roots[a][b][1];
      }
    }

  // Galois orbits
  // allocate and initialize
  unsigned *orbit_rep = new unsigned[q];
  for (unsigned i = 0; i < q; i++)
    orbit_rep[i] = q; // means null
  
  int *orbit_size = new int[q]; // only valid if i == orbit_rep[i]

  // fill the tables
  for (unsigned i = 0; i < q; i++) {
    if (orbit_rep[i] != q) continue; // already filled
    unsigned j = i;
    int size = 1;
    while (1) {
      orbit_rep[j] = i;
      j = mult[j][j];
      if (j == i) {
        orbit_size[i] = size;
        break;
      }
      size++;
    } 
  }

  return std::make_tuple(depressed_cubic_roots, quadratic_roots,
                         mult, divi, orbit_rep, orbit_size);
}


// REMARK: For Conway polynomials, x always reduces to a primitive element. 
int* log_table(int n) {

  // Once we get the fast Fq multiplication working, we'll use that.
  throw std::runtime_error("Error: log_table is not implemented.");
  
  const unsigned q = 1 << n;
  const unsigned halfq = 1 << n;
  int* table = new int[q]; // The zero entry is a placeholder.

  int x = (1<<1);
  int xpow = 1;
  for (int pow=0; pow < halfq; pow++) {
    table[xpow] = pow;
    xpow *= x; // TODO: Wrong multiplication.
  }

  // Flip the sign in the second half for fast inversion.
  for (int pow=halfq; pow < q-1; pow++) {
    table[xpow] = (q - pow);
    xpow *= x; // TODO: Wrong multiplication.
  }
    
  return table;
}


int* cubic_root_count_table(int n) {
  // In order to generate the cubic lookup table. We need to solve the system
  //    (x^2 + ax + b)(x+a) = x^3 + sx + t
  //    a^2 + b = 1.
  //
  // This is trivial.
  //
  //    b = a^2 + 1
  //    s = 1
  //    t = a (a^2 + 1) = a^3 + a.
  //
  // Furthermore, note that the quadratic
  //
  //    x^2 + ax + a^2
  //
  // has a root if and only if either a=0 or Fq has a 3rd root of unity.
    
  unsigned q = 1 << n;
  int* cubic_table = new int[q];
  bool has_cube_root = ( (n % 2) == 0 ? true : false);

  int num_roots_if_solvable = (has_cube_root ? 3 : 1);
  
  // Initialize. Default is zero, since if we don't encounter it, the cubic
  // will have no roots.
  for (unsigned i=0; i<q; i++) {
    cubic_table[i] = 0;
  }

  // Populate
  for (unsigned a=0; a<q; a++) {
    unsigned aaa = ff2k_mult(ff2k_square(a), a);
    cubic_table[aaa ^ a] = num_roots_if_solvable;
  }
  
  return cubic_table;
}
