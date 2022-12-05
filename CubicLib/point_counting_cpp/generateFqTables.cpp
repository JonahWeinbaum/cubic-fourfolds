#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <tuple>
#include "tableio.h"
#include "Fq_tables.h"


int main() {
  
  for (int n = 1; n < 12; n++) {
    std::clock_t cputime = std::clock();

    unsigned*** depressed_cubic_roots;
    unsigned*** quadratic_roots;
    unsigned** mult;
    unsigned** divi;
    unsigned* orbit_rep;
    int* orbit_size;
    std::tie(depressed_cubic_roots,
             quadratic_roots,
             mult, divi, orbit_rep, orbit_size) = generate_Fq_tables(n);
    
    cputime = std::clock() - cputime;
    std::cerr << "precomputation took " << cputime * 1./ CLOCKS_PER_SEC << std::endl;
  
    // Write data to file.
    unsigned q = 1<<n;    
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
  
  }
  
  return 0;
}
