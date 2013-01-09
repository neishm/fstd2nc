#include "compress.h"
#include "io.h"
#include <assert.h>

// A bitwise file interface
typedef struct {
  FILE *file;
  unsigned int current_word;
  unsigned char usable_bits;
} BITFILE;

BITFILE as_bitfile (FILE *f) {
  BITFILE b;
  b.file = f;
  b.current_word = 0;
  b.usable_bits = 0;
  return b;
}

// Read up to 32 bits from a stream
unsigned int bitread (BITFILE *b, unsigned char nbits) {
  unsigned int data = 0;
/*
  if (b->usable_bits == 0) assert (b->current_word == 0);
  if (b->usable_bits == 1) assert (b->current_word < 2);
  if (b->usable_bits == 2) assert (b->current_word < 4);
  if (b->usable_bits == 3) assert (b->current_word < 8);
  if (b->usable_bits == 4) assert (b->current_word < 16);
  if (b->usable_bits == 5) assert (b->current_word < 32);
  if (b->usable_bits == 6) assert (b->current_word < 64);
  if (b->usable_bits == 7) assert (b->current_word < 128);
  if (b->usable_bits == 8) assert (b->current_word < 256);
*/

  // Can we use up the rest of this current word?
  if (nbits > b->usable_bits) {
    data |= b->current_word;
    nbits -= b->usable_bits;
    data <<= nbits;
    b->current_word = fread32(b->file);
//    printf ("?? new word: %08X\n", b->current_word);
    b->usable_bits = 32;
  }

  // Read part of a word
  unsigned char leftover_bits = b->usable_bits - nbits;
  data |= (b->current_word >> leftover_bits);
  b->current_word &= (0xFFFFFFFF >> (32-leftover_bits));  // fails for >>32
  if (leftover_bits == 0) b->current_word = 0; // account for shift failure
  b->usable_bits = leftover_bits;
//  if (b->usable_bits == 0) printf ("0 == %d?\n", b->current_word);
//  if (b->usable_bits == 0) assert (b->current_word == 0);

  return data;
}


// The parameters for a compressed field
typedef struct {

  unsigned char predictor_type;
  unsigned char degree;
  unsigned char step;
  unsigned char mantissa_nbits;
  unsigned char levels;
  unsigned char version;
  unsigned char reserved3;

  unsigned short exp_min;
  unsigned char exp_nbits;

  unsigned char sign_code;
  unsigned char exp_code;
  unsigned char mantissa_code;

} ZParams;


// Read the parameters for a compressed field
int read_compress32_params (FILE *f, ZParams *p) {

  unsigned int unknown, zfstzip, zieee_info;

  unknown = fread32 (f);  // compressed size?
  zfstzip = fread32 (f);
  zieee_info = fread32 (f);

//  printf ("raw params: %08X %08X\n", zfstzip, zieee_info);

  p->reserved3      = (zfstzip & 0xFF000000) >> 24;
  p->version        = (zfstzip & 0x00FC0000) >> 18;
  p->levels         = (zfstzip & 0x00038000) >> 15;
  p->mantissa_nbits = (zfstzip & 0x00007C00) >> 10;
  p->step           = (zfstzip & 0x00000380) >>  7;
  p->degree         = (zfstzip & 0x00000070) >>  4;
  p->predictor_type = (zfstzip & 0x0000000F) >>  0;

  p->exp_min       = (zieee_info & 0xFFFF0000) >> 16;
  p->exp_nbits     = (zieee_info & 0x0000FF00) >>  8;

  p->sign_code     = (zieee_info & 0x00000030) >>  4;
  p->exp_code      = (zieee_info & 0x0000000C) >>  2;
  p->mantissa_code = (zieee_info & 0x00000003) >>  0;

  return 0;
}

int print_compress32_params (ZParams *p) {
  printf ("predictor_type = %2i   degree = %2i   step = %2i\n",
        p->predictor_type,     p->degree,     p->step);

  printf ("mantissa_nbits = %2i   levels = %2i   version = %2i\n",
        p->mantissa_nbits,     p->levels,     p->version);

  printf ("exp_min = %3i   exp_nbits = %2i   sign_code = %2i\n",
        p->exp_min,     p->exp_nbits,     p->sign_code);

  printf ("exp_code = %2i   mantissa_code = %2i\n",
        p->exp_code,     p->mantissa_code);
}

// Read a compressed field
int read_compress32 (FILE *f, RecordHeader *h, int recsize, float *out) {

  ZParams p;

  assert (h->ni * h->nj == recsize);

  // Fill the output
  // TODO - remove this, or set to NaN?
  for (int i = 0; i < recsize; i++) out[i] = 0.0;

  assert (h->datyp == 133);

  read_compress32_params (f, &p);
//  print_compress32_params (&p);

  // Extract the sign bits
  assert (p.sign_code == 0 || p.sign_code == 1 || p.sign_code == 2);
  unsigned char sign[recsize];
  // All positive?
  if (p.sign_code == 0) {
    for (int i = 0; i < recsize; i++) sign[i] = 0;
  }
  // All negative?
  else if (p.sign_code == 1) {
    for (int i = 0; i < recsize; i++) sign[i] = 1;
  }
  // Packed?
  else if (p.sign_code == 2 || p.sign_code == 3) {
    unsigned int len_sign = fread32(f);  // not used?
    BITFILE b = as_bitfile(f);
    unsigned char lastval = 0;
    int i = 0;
    while (i < recsize) {
      unsigned char seq_type = bitread(&b,1);
//      assert (seq_type == 0 || seq_type == 1);
      // sequence?
      if (seq_type == 0) {
//        printf ("?? %4s %d <bit sequence>\n", h->nomvar, h->ip1);
        unsigned char runlength = 7;
        // Truncate if near the end of the record?
        if (i+runlength > recsize) runlength = recsize-i;
        for (int j = 0; j < runlength; j++) sign[i+j] = bitread(&b,1);
        i += runlength; continue;
      }
      // count?
      if (seq_type == 1) {
        unsigned char val = bitread(&b,1);
        unsigned char count = bitread(&b,6);
//        printf ("?? %4s %d bit = %d, i=%d, count=%d, recsize=%d\n", h->nomvar, h->ip1, val, i, count, recsize);
        // Special case: really big run
        if (count == 63) {
//          assert (val == lastval);  // Hope this is consistent!
          for (int j = 0; j < 255; j++) sign[i+j] = lastval;
          i += 255; continue;
        }
        // Normal case
        else {
          assert (i+count <= recsize);
          for (int j = 0; j < count; j++) sign[i+j] = val;
          lastval = val;
          i += count; continue;
        }
      }  // end of "COUNT" encoding

    }  // end of loop over field elements
  }  // end of packed sign condition


  // Extract the exponents
  assert (p.exp_code == 0 || p.exp_code == 1 || p.exp_code == 2);
  unsigned short exp[recsize];
  // Identical exponents?
  if (p.exp_code == 0) {
    for (int i = 0; i < recsize; i++) exp[i] = p.exp_min;
  }
  // Otherwise, have packed exponents
  if (p.exp_code == 1 || p.exp_code == 2) {
    // Temporary variable to hold the exponent differences
    short exp_diff[recsize];
//    for (int k = 0; k < recsize; k++) exp_diff[k] = 0;

    unsigned int len_exp = fread32(f);  // Not used?
//    printf ("len_exp = %d\n", len_exp);
    BITFILE b = as_bitfile(f);
    unsigned char nbits_req_container = bitread(&b,3); // 3 b/c of p.step??
    int ni = h->ni;
    int nj = h->nj;
    assert (ni * nj == recsize);
    // Fill first row
    for (int i = 0; i < ni; i++) {
      exp[i] = bitread(&b,p.exp_nbits);
      exp_diff[i] = 0;
    }
    // Fill first column
    for (int k = ni; k < recsize; k+=ni) {
      exp[k] = bitread(&b,p.exp_nbits);
      exp_diff[k] = 0;
    }

    // Loop over each tile
    for (int j0 = 1; j0 < nj; j0 += p.step) {
      for (int i0 = 1; i0 < ni; i0 += p.step) {
        unsigned char nbits_needed = bitread(&b,nbits_req_container);
//        printf ("nbits_needed: %d\n", nbits_needed);
        for (int dj = 0; dj < p.step && j0+dj < nj; dj++) {
          for (int di = 0; di < p.step && i0+di < ni; di++) {
            int k = ni*(j0+dj) + (i0+di);
            // Special case: no bits needed for the encoding?
            if (nbits_needed == 0) exp_diff[k] = 0;
            // Otherwise, extract the diff
            else {
              int nbits2 = nbits_needed + 1; // minimum of 2 bits for encoding
              int tmp = bitread(&b,nbits2);
              // Apply a sign
              exp_diff[k] = (tmp<<(32-nbits2)) >> (32-nbits2);
//              assert (exp_diff[k] > -1000 && exp_diff[k] < 1000);
            }
          }
        }
      }
    } // end of loop over tiles

    // Integrate the differences to restore the exponent field
    for (int j = 1; j < nj; j++) {
      for (int i = 1; i < ni; i++) {
        int k22 = ni*j     + i;
        int k12 = ni*j     + (i-1);
        int k21 = ni*(j-1) + i;
        int k11 = ni*(j-1) + (i-1);
        exp[k22] = exp_diff[k22] + exp[k21] + exp[k12] - exp[k11];
      }
    }
    // Add the minimum value
    for (int k = 0; k < recsize; k++) exp[k] += p.exp_min;
    // Test the diff
//    for (int k = 0; k < recsize; k++) {
//      assert (exp_diff[k] >= 0);
//      exp[k] = exp_diff[k];
//    }

  } // End of packed exponent case


  // Reconstruct the field from the information we have
//  for (int i = 0; i < recsize; i++) out[i] = 1 - 2*sign[i];
  for (int k = 0; k < recsize; k++) out[k] = exp[k];

  return 0;
}
