#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include "tableio.h"

/////////////////////////////////
// Global variables

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

/////////////////////////////////


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

int generate_tables(int n) {
  std::clock_t cputime = std::clock();

  // Define q.
  unsigned q = 1<<n;

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

  cputime = std::clock() - cputime;
  std::cerr << "precomputation took " << cputime * 1./ CLOCKS_PER_SEC << std::endl;

  // Write data to file.
  std::string qq = std::to_string(q);
    
  write_table(depressed_cubic_roots, q, q, 3, "depressed_cubic_roots_" + qq);
  write_table(quadratic_roots, q, q, 2, "quadratic_roots_" + qq);

  write_table(mult, q, q, "mult_" + qq);
  write_table(divi, q, q, "divi_" + qq);

  write_table(orbit_rep, q, "orbit_rep_" + qq);
  write_table(orbit_size, q, "orbit_size_" + qq);

  if (false) {
    compare_2_table(mult, read_table(q, q, "mult_4"), q, q);
    compare_3_table(quadratic_roots, read_table(q, q, 2, "quadratic_roots_4"),
                    q, q, 2);
  }
  
  return 0;
}

int main() {
  for (int n = 1; n < 15; n++)
    generate_tables(n);

  return 0;
}
