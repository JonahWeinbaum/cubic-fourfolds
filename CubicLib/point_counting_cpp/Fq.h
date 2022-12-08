#ifndef N
static_assert(false, "No value of N was specified during compilation.");
#endif


using ff2k_t = unsigned;

int initialized_ff_bitsize();

ff2k_t ff2k_mult(ff2k_t a, ff2k_t b);
ff2k_t ff2k_square(ff2k_t a);
ff2k_t ff2k_inv(ff2k_t a);
ff2k_t ff2k_divi(ff2k_t a, ff2k_t b);

int gcd_degree(unsigned* A, unsigned* B, int dA, int dB);

/*
int count_linear_roots(unsigned a, unsigned b);
int count_quadratic_roots(unsigned a, unsigned b, unsigned c);
int count_cubic_roots(unsigned a, unsigned b, unsigned c, unsigned d);
int count_quartic_roots(unsigned a, unsigned b, unsigned c, unsigned d, unsigned e);

int count_quadratic_roots(unsigned* f);
int count_cubic_roots(unsigned* f);
int count_quartic_roots(unsigned* f);
*/

int count_poly_roots(unsigned* poly, int deg);
