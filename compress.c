#include "compress.h"
#include "io.h"
#include <assert.h>

// A bitwise file interface
typedef struct {
  unsigned char *buf_start;
  unsigned char **buf;
  int length;
  unsigned int current_word;
  unsigned char usable_bits;
} BITFILE;

BITFILE as_bitfile (unsigned char **buf) {
  BITFILE b;
  b.length = read32(*buf);
  (*buf) += 4;
  b.buf = buf;
  b.buf_start = *buf;
  b.current_word = 0;
  b.usable_bits = 0;
  return b;
}

// Read up to 32 bits from a stream
unsigned int bitread (BITFILE *b, unsigned char nbits) {
  if (nbits == 0) return 0;

  unsigned int data = 0;

  // Can we use up the rest of this current word?
  if (nbits > b->usable_bits) {
    data |= b->current_word;
    nbits -= b->usable_bits;
    data <<= nbits;
    b->current_word = read32(*(b->buf));
    (*b->buf) += 4;
//    printf ("?? new word: %08X\n", b->current_word);
    b->usable_bits = 32;
  }

  // Read part of a word
  unsigned char leftover_bits = b->usable_bits - nbits;
  data |= (b->current_word >> leftover_bits);
  b->current_word &= (0xFFFFFFFF >> (32-leftover_bits));  // fails for >>32
  if (leftover_bits == 0) b->current_word = 0; // account for shift failure
  b->usable_bits = leftover_bits;

  return data;
}

// Finish with a bit stream
int close_bitstream (BITFILE *b) {
  int used = *(b->buf) - b->buf_start;
//  printf ("used: %d length %d\n", used, b->length);
  if (b->length > 0) {  // skip for the bad length encoded in the mantissa
    *(b->buf) += (b->length - used);
  }
  return 0;
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
int read_compress32_params (unsigned char **buf, ZParams *p) {

  unsigned int unknown, zfstzip, zieee_info;

  unknown = read32 (*buf);  // compressed size?
  (*buf) += 4;
  zfstzip = read32 (*buf);
  (*buf) += 4;
  zieee_info = read32 (*buf);
  (*buf) += 4;

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

  return 0;
}

// Lorenzo predictor decoder
int lorenzo_decoder (BITFILE *b, int ni, int nj, int step, int nbits, int *out) {
  int recsize = ni * nj;
  // Temporary variable to hold the differences
  int diff[recsize];

  // Number of bits required to read the tiles' number of bits
  unsigned char nbits_req_container = bitread(b,3);
  assert (nbits_req_container > 0);

  // Fill first row
  for (int i = 0; i < ni; i++) {
    out[i] = bitread(b,nbits);
    diff[i] = 0;
  }
  // Fill first column
  for (int k = ni; k < recsize; k+=ni) {
    out[k] = bitread(b,nbits);
    diff[k] = 0;
  }

  // Loop over each tile
  for (int j0 = 1; j0 < nj; j0 += step) {
    for (int i0 = 1; i0 < ni; i0 += step) {
      unsigned char nbits_needed = bitread(b,nbits_req_container);
//      printf ("nbits_needed: %d\n", nbits_needed);
      for (int dj = 0; dj < step && j0+dj < nj; dj++) {
        for (int di = 0; di < step && i0+di < ni; di++) {
          int k = ni*(j0+dj) + (i0+di);
          // Special case: no bits needed for the encoding?
          if (nbits_needed == 0) diff[k] = 0;
          // Otherwise, extract the diff
          else {
            int nbits2 = nbits_needed + 1; // minimum of 2 bits for encoding
            int tmp = bitread(b,nbits2);
            // Apply a sign
            diff[k] = (tmp<<(32-nbits2)) >> (32-nbits2);
          }
        }
      }
    }
  } // end of loop over tiles

  // Integrate the differences to restore the field
  for (int j = 1; j < nj; j++) {
    for (int i = 1; i < ni; i++) {
      int k22 = ni*j     + i;
      int k12 = ni*j     + (i-1);
      int k21 = ni*(j-1) + i;
      int k11 = ni*(j-1) + (i-1);
      out[k22] = diff[k22] + out[k21] + out[k12] - out[k11];
    }
  }

  return 0;
}


// Read a compressed field
int read_compress32 (unsigned char *buf, RecordHeader *h, int recsize, float *out) {

  ZParams p;

  assert (h->ni * h->nj == recsize);
  assert (h->datyp == 133);

  read_compress32_params (&buf, &p);
//  print_compress32_params (&p);

//  printf ("?? %4s %d sign_code = %d exp_nbits = %d\n", h->nomvar, h->ip1, p.sign_code, p.exp_nbits);

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
    BITFILE b = as_bitfile(&buf);
//    printf ("len_sign = %d\n", b.length);
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
      else if (seq_type == 1) {
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

    close_bitstream(&b);

  }  // end of packed sign condition


  // Extract the exponents
  assert (p.exp_code == 0 || p.exp_code == 1 || p.exp_code == 2);
  int exp[recsize];
  // Identical exponents?
  if (p.exp_code == 0) {
    for (int i = 0; i < recsize; i++) exp[i] = p.exp_min;
  }
  // Otherwise, have packed exponents
  else if (p.exp_code == 1 || p.exp_code == 2) {
    // Call the Lorenzo predictor decoder
    BITFILE b = as_bitfile(&buf);
//    printf ("len_exp = %d\n", b.length);
//    assert (b.length > 0);
    lorenzo_decoder (&b, h->ni, h->nj, p.step, p.exp_nbits, exp);
    close_bitstream (&b);

    // Add the minimum value
    for (int k = 0; k < recsize; k++) exp[k] += p.exp_min;

  } // End of packed exponent case


  // Extract the mantissas
  unsigned int mantissa[recsize];
  // First case - nothing fancy applied
  assert (p.mantissa_code == 0 || p.mantissa_code == 1);
//  printf ("mantissa code: %d\n", p.mantissa_code);
  if (p.mantissa_code == 1) {
    BITFILE b = as_bitfile(&buf);
    for (int k = 0; k < recsize; k++) mantissa[k] = bitread(&b,p.mantissa_nbits);
    close_bitstream(&b);
  }
  // More fancy case - Lorenzo predictor
  else if (p.mantissa_code == 0) {
    BITFILE b = as_bitfile(&buf);
    lorenzo_decoder (&b, h->ni, h->nj, p.step, p.mantissa_nbits, (int*)mantissa);
    close_bitstream(&b);
  }

  // Account for truncated mantissa
  if (p.mantissa_nbits < 23) {
    for (int k = 0; k < recsize; k++) mantissa[k] <<= (23-p.mantissa_nbits);
  }

  // Reconstruct the field from the information we have
//  for (int i = 0; i < recsize; i++) out[i] = 1 - 2*sign[i];
//  for (int k = 0; k < recsize; k++) out[k] = exp[k];
//  for (int k = 0; k < recsize; k++) out[k] = mantissa[k];
  for (int k = 0; k < recsize; k++) {
    ((unsigned int *)out)[k] = (sign[k]<<31) | (exp[k]<<23) | mantissa[k];
  }

  return 0;
}
