#include <wmmintrin.h>
#include <ctype.h>

const n = 14;

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
  (1<<15) + (1<<1) + 1,
  (1<<16) + (1<<5) + (1<<3) + (1<<1) + 1,
  (1<<17) + (1<<3) + 1,
  (1<<18) + (1<<3) + 1,
  (1<<19) + (1<<5) + (1<<2) + (1<<1) + 1,
  (1<<20) + (1<<3) + 1,
  (1<<21) + (1<<2) + 1,
  (1<<22) + (1<<1) + 1 };


#ifdef WITHCACHE

unsigned mult(unsigned a, unsigned b) {
  return mult_table[a][b];
}

# elif defined FINITEFIELDBITSIZE

const p = polynomials[FINITEFIELDBITSIZE];


static inline ff2k_t _ff2k_pclmul (ff2k_t a, ff2k_t b)
{
     register __m128i A, B, C;
     A[0]=a; B[0] = b;
     C = _mm_clmulepi64_si128 (A,B,0);
     return (ff2k_t)C[0];
}

// NOTE: Need to enable compile flags to actually compile this code.
/*
unsigned mult_new(unsigned a, unsigned b) {
      // main multiplication algorithm
      unsigned a = i, b = j, ij = 0;

      // Convert into 64 bit types
      __m64i a64 = _mm_cvtsi32_si64(a);
      __m128i b128 = _mm_cvtsi32_si128(b);
      
      // Do carryless multiply.
      __m128i polymult = _mm_clmulepi64_si128(a64,   // Input 64-bit "a"
                                              b128,  // Input 128-bit "b'(full of 1 bits)
                                              0      // Multiply the lower 64-bits inputs
                                              );
      
      // Extract lower 64 bit word.
      __m64i ij = _mm_cvtsi128_si64(polymult);

      // Reduce.
      pxi  = p << n;
      for (i= 2*n; i > n; i--) {
        // If the i-th bit is 1, reduce by the polynomial.
        if (ij & (1<<i)) ij ^= pxi;
        pxi = pxi >> 1;
      }

      // Convert down to 32 bits
      unsigned res = _mm_cvtsi64_si32(ij);
      return res;
}
*/

/*
// NOTE: Code below requires intel's compiler.
unsigned mult_new(unsigned a, unsigned b) {
      // main multiplication algorithm
      unsigned a = i, b = j, ij = 0;

      // Convert into 64 bit types
      __m64i a64 = _mm_cvtsi32_si64(a);
      __m128i b128 = _mm_cvtsi32_si128(b);
      
      // Do carryless multiply.
      __m128i polymult = _mm_clmulepi64_si128(a64,   // Input 64-bit "a"
                                              b128,  // Input 128-bit "b'(full of 1 bits)
                                              0      // Multiply the lower 64-bits inputs
                                              );
      
      // Extract lower 64 bit word.
      __m64i ij = _mm_cvtsi128_si64(polymult);

      // Reduce.
      pxi  = p << n;
      for (i= 2*n; i > n; i--) {
        // If the i-th bit is 1, reduce by the polynomial.
        if (ij & (1<<i)) ij ^= pxi;
        pxi = pxi >> 1;
      }

      // Convert down to 32 bits
      unsigned res = _mm_cvtsi64_si32(ij);
      return res;
}
*/

#endif

#ifdef TEST
// TODO: Add a simple test.

#endif

int main() {
  return 0;
}

  /*
    https://learn.microsoft.com/en-us/cpp/intrinsics/x86-intrinsics-list?view=msvc-170

    _mm_cvtsi32_si64() // Convert to 64 bit.
    _mm_cvtsi64_si32() // Convert to 32 bit.

    __m128i Result = _mm_clmulepi64_si128(
		_mm_set_epi64x(0, Bits), // Input 64-bit "a"
		_mm_set1_epi8(0xFF),     // Input 128-bit "b'(full of 1 bits)
		0                        // Multiply the lower 64-bits of each of the inputs
	);
	const std::uint64_t Lo = _mm_cvtsi128_si64(Result);
        int _mm_cvtsi128_si32(__m128i);
  */

