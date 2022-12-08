using all_table_t = std::tuple<unsigned***, unsigned***,
                               unsigned**, unsigned**,
                               unsigned*, int*>;

all_table_t generate_Fq_tables(int n);

/* std::tuple<unsigned***, unsigned***, */
/*   unsigned**, unsigned**, */
/*   unsigned*, int*> generate_Fq_tables(int n); */

std::tuple<unsigned*, unsigned*> generate_orbit_tables(int);

// Defining polynomials for Fq arithmetic. 
extern const unsigned Fq_polynomials[];
