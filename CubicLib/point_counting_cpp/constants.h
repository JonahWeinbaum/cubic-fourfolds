
using ff2k_t = unsigned;

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

// Setup aliasing for the value of q.
#if defined N
const unsigned FINITEFIELDBITSIZE = N;
const unsigned q = 1 << FINITEFIELDBITSIZE;
const unsigned p = polynomials[FINITEFIELDBITSIZE];
#endif
