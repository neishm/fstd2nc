#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

// pseudo-RPN interface

typedef unsigned char byte;

unsigned int read32 (byte *b) {
  return (((unsigned int)(b[0]))<<24) | (((int)(b[1]))<<16) | (((int)(b[2]))<<8) | (((int)(b[3]))<<0);
}

int read24 (byte *b) {
  return (((int)(b[0]))<<16) | (((int)(b[1]))<<8) | (((int)(b[2]))<<0);
}

int read16 (byte *b) {
  return (((int)(b[0]))<<8) | (((int)(b[1]))<<0);
}

// Read packed characters
void readchar (byte *dest, byte *src, int n) {
  for (int i = 0; i < n; i++) {
    dest[i] = 0;
    // first byte needed
    byte j = (i*6)/8;
    // shift from the beginning
    byte sh = (i*6)%8;
    dest[i] |= ((byte)(src[j] << sh) >> 2);
    // do we need a second byte?
    if (sh > 2) dest[i] |= (src[j+1]>>(10-sh));
    dest[i] += 32;
  }
}



typedef struct {
  long long int file_size;
  int num_overwrites;
  int nrecs_all;
  int nchunks;
  unsigned long long int last_chunk;
  int max_data_length;
  int num_edits;
  int nrecs;
} FileHeader;

void read_file_header (FILE *f, FileHeader *h) {
  int nbytes;
  byte buf[208];
  char fixed1[16] = "\0\0\0\x1a\0\0\0\0XDF0STDR";
  char fixed2[8] = "\0\x10\0\x09\0\x02\0\x01";
  char fixed3[8] = "\0\0\0\0\0\0\0\0";
  nbytes = fread (buf, 1, 208, f);
  assert (nbytes == 208);
  assert (memcmp(buf, fixed1, 16) == 0);
  assert (memcmp(buf+40, fixed2, 8) == 0);
  assert (memcmp(buf+56, fixed3, 8) == 0);
//  for (int i = 0; i < 16; i++) assert (buf[64+i*8] == 'S');
//  for (int i = 0; i < 16; i++) assert (buf[64+i*8+1] == 'F');
  //TODO: verify the 144-byte fixed thing??
  h->file_size = read32(buf+16) * 8L;
  h->num_overwrites = read32(buf+20);
  h->nrecs_all = read32(buf+24);
  h->nchunks = read32(buf+28);
  h->last_chunk = read32(buf+32) * 8L;
  h->max_data_length = read32(buf+36) * 8;
  h->num_edits = read32(buf+48);
  h->nrecs = read32(buf+52);
}

void print_file_header (FileHeader *h) {
  printf ("file_size: %lld, num_overwrites: %d, nrecs (including erased): %d, nchunks: %d, last_chunk: %llx, max_data_length: %d, num_edits: %d, nrecs: %d\n", h->file_size, h->num_overwrites, h->nrecs_all, h->nchunks, h->last_chunk, h->max_data_length, h->num_edits, h->nrecs);
}

typedef struct {
  unsigned int this_chunk_words;
  unsigned long long int this_chunk;
  unsigned int next_chunk_words;
  unsigned long long int  next_chunk;
  int nrecs;
  unsigned int checksum;
} ChunkHeader;

void read_chunk_header (FILE *f, ChunkHeader *h) {
  int nbytes;
  byte buf[32];
  char fixed1[4] = "\0\0\x09\x04";
  char fixed2[8] = "\0\0\0\0\0\0\0\0";
  char fixed3[4] = "\0\0\0\0";
  nbytes = fread (buf, 1, 32, f);
  assert (nbytes == 32);
  assert (memcmp(buf, fixed1, 4) == 0);
  assert (memcmp(buf+8, fixed2, 8) == 0);
  assert (memcmp(buf+28, fixed3, 4) == 0);
  h->this_chunk_words = read32(buf+4);
  h->this_chunk = h->this_chunk_words * 8L - 8;
  h->next_chunk_words = read32(buf+16);
  h->next_chunk = h->next_chunk_words * 8L;
  if (h->next_chunk > 0) h->next_chunk -= 8; // Rewind a bit
  h->nrecs = read32(buf+20);
  h->checksum = read32(buf+24);
}

void print_chunk_header (ChunkHeader *h) {
  printf ("this_chunk: %llx, next_chunk: %llx, nrecs: %d, checksum: %x\n", h->this_chunk, h->next_chunk, h->nrecs, h->checksum);
}

typedef char Nomvar[5];
typedef char Etiket[13];
typedef char Typvar[3];

typedef struct {
  int status;
  int size;
  unsigned long long int data;
  int deet;
  int npak;
  int ni;
  char grtyp;
  int nj;
  int datyp;
  int nk;
  int npas;
  int ig4;
  int ig2;
  int ig1;
  int ig3;
  Etiket etiket;
  Typvar typvar;
  Nomvar nomvar;
  int ip1;
  int ip2;
  int ip3;
  long long int dateo;
  unsigned int checksum;
} RecordHeader;

void read_record_header (FILE *f, RecordHeader *h) {
  int nbytes;
  byte buf[72];
  nbytes = fread (buf, 1, 72, f);
  assert (nbytes == 72);
  h->status = buf[0];
  h->size = read24(buf+1);
  h->data = read32(buf+4) * 8L;
  assert (h->data > 8);
  h->data -= 8;  // rewind a bit to get the proper start of the data
  h->deet = read16(buf+9);
  h->npak = buf[11];
  h->ni = read24(buf+12);
  h->grtyp = buf[15];
  h->nj = read24(buf+16);
  h->datyp = buf[19];
  h->nk = read24(buf+20)>>4;
  h->npas = (read32(buf+24)>>4)/4;
  h->ig4 = read24(buf+28);
//  h->ig2 = buf[27]*65536 + buf[35]*256 + buf[39];  // this isn't part of ig2???
  h->ig2 = buf[35]*256 + buf[39];
  h->ig1 = read24(buf+32);
  h->ig3 = read24(buf+36);
  readchar (h->etiket, buf+40, 5);
  readchar (h->etiket+5, buf+44, 5);
  readchar (h->etiket+10, buf+48, 2);
  h->etiket[12] = 0;
  h->typvar[0] = ((buf[49]&0x0F)<<2) + (buf[50]>>6) + 32;
  h->typvar[1] = (buf[50]&0x3F) + 32;
  h->typvar[2] = 0;
  readchar (h->nomvar, buf+52, 4);
  h->nomvar[4] = 0;
  h->ip1 = read32(buf+56) >> 4;
  h->ip2 = read32(buf+60) >> 4;
  h->ip3 = read32(buf+64) >> 4;
  h->dateo = read32(buf+68);

  // Checksum
  h->checksum = 0;
  for (int i = 0; i < 72; i+= 4) {
    h->checksum ^= ((buf[i]<<24) | (buf[i+1]<<16) | (buf[i+2]<<8) | buf[i+3]);
  }
}

void print_record_header (RecordHeader *h) {
  printf ("status: %d, size: %d, data pointer: %llx, deet: %d, npak: %d, ni: %d, grtyp: '%c', nj: %d, datyp: %d nk: %d, npas: %d, ig4: %d, ig2: %d, ig1: %d, ig3: %d, etiket: '%s', typvar: '%s', nomvar: '%s', ip1: %d, ip2: %d, ip3: %d, \ndateo: %lld\n", h->status, h->size, h->data, h->deet, h->npak, h->ni, h->grtyp, h->nj, h->datyp, h->nk, h->npas, h->ig4, h->ig2, h->ig1, h->ig3, h->etiket, h->typvar, h->nomvar, h->ip1, h->ip2, h->ip2, h->dateo);
}


// Get the number of records in a file
int get_num_records (char *filename) {
  FileHeader h;
  FILE *file = fopen (filename, "rb");
  read_file_header (file, &h);
  fclose (file);
  return h.nrecs;
}

// Read the records into a pre-allocated array
int get_record_headers (char *filename, RecordHeader *h) {
  FILE *f = fopen(filename, "rb");
  FileHeader fileheader;
  fseek (f, 0, SEEK_SET);
  read_file_header (f, &fileheader);
  int rec = 0;
  for (int c = 0; c < fileheader.nchunks; c++) {
    ChunkHeader chunkheader;
    read_chunk_header (f, &chunkheader);
    unsigned int checksum = 0;
    for (int r = 0; r < chunkheader.nrecs; r++) {
      read_record_header (f, h+rec);
      // Skip erased records
      if (h[rec].status == 255) continue;
      checksum ^= h[rec].checksum;
      rec++;
    }
    checksum ^= chunkheader.nrecs;
    checksum ^= (chunkheader.next_chunk_words);
//    printf ("chunk %d checksum: 0x%08X   computed checksum: 0x%08X, xor: %08X\n", c, chunkheader.checksum, checksum, chunkheader.checksum ^ checksum);
    assert (checksum == chunkheader.checksum);
    // Go to next chunk
    fseek (f, chunkheader.next_chunk, SEEK_SET);
  }
  fclose (f);
  return 0;
}


// Compare RecordHeaders for equality
// NOTE: this used to compare a RecordHeader with a Varinfo_entry,
// hence the 'h' and 'v' names.
int receq (RecordHeader *h, RecordHeader *v) {
  if (strncmp(h->nomvar, v->nomvar, sizeof(Nomvar)) != 0) {
    printf ("'%s' != '%s'\n", h->nomvar, v->nomvar);
    return 0;
  }
  if (strncmp(h->etiket, v->etiket, sizeof(Etiket)) != 0) return 0;
  if (strncmp(h->typvar, v->typvar, sizeof(Typvar)) != 0) return 0;
  if (h->grtyp != v->grtyp) return 0;
  if (h->ip1 != v->ip1) return 0;
  if (h->ip2 != v->ip2) return 0;
  if (h->ip3 != v->ip3) return 0;
  if (h->ig1 != v->ig1) return 0;
  if (h->ig2 != v->ig2) return 0;
  if (h->ig3 != v->ig3) return 0;
  if (h->ig4 != v->ig4) return 0;
  if (h->ni != v->ni) return 0;
  if (h->nj != v->nj) return 0;
  if (h->nk != v->nk) return 0;
  return 1;

}


// Read a chunk of data
int read_data (char *filename, int nrecs, RecordHeader *headers, int recsize, float *out) {

  FILE *file = fopen (filename, "rb");

  // Loop over all records
  for (int rec = 0; rec < nrecs; rec++) {

      unsigned long long offset = headers[rec].data;
//      printf ("offset: %llx\n", offset);
      RecordHeader h;
      fseek (file, offset, SEEK_SET);
      read_record_header (file, &h);
//      print_record_header (&h);
      // Make sure the header matches
      assert (receq(headers+rec, &h) == 1);
      assert (h.datyp == 1 || h.datyp == 5); //currently only support packed floating-point/IEEE floating point

      byte b[4];
      fread (b, 1, 4, file);
      assert (read32(b) == 0);
      fread (b, 1, 4, file);
      assert (read32(b) == 0);

      // Easiest case: 32-bit IEEE float
      if (h.datyp == 5) {
        assert (h.npak == 32);  // must be 32-bit?!
        assert (sizeof(float) == sizeof(uint32_t));
        byte *raw = malloc(4*recsize);
        fread (raw, 4, recsize, file);
        for (int i = 0; i < recsize; i++) {
          ((uint32_t*)(out))[i] = read32(raw+4*i);
        }
        free(raw);
//        for (int i = 0; i < recsize; i++) {
//          printf ("%g  ", out[i]);
//        }
//        printf ("\n");
      }
      // Other supported case: packed float
      else if (h.datyp == 1) {

        // Get and validate record size
        fread (b, 1, 4, file);
        int recsize_ = read32(b) - 0x7ff00000;
        assert (recsize_ == recsize);

        // Get encoding info
        fread (b, 1, 4, file);
        int32_t info = read32(b);
//        printf ("info: 0x%08x\n", info);
        int diff_shift = (info >> 16) - 0x0ff0;
//        printf ("diff_shift: %d\n", diff_shift);

        // Get min value
        float min;
        fread (b, 1, 4, file);
        assert (b[3] == 0);
        int mincode = read24(b);
        if ((info & 0x0000ff00) == 0x00001100) {
          min = 0;
//          printf ("min is zero\n");
        }
        else {
          int min_shift = ((info & 0x0000fff0) >> 4) - 0x3d0 + 1;
//          printf ("min_shift: %d\n", min_shift);
//          printf ("mincode: %d\n", mincode);
          min = mincode / 16777216.;
          while (min_shift > 0) { min *= 2; min_shift--; }
          while (min_shift < 0) { min /= 2; min_shift++; }
          if ((info & 0x0000000f) == 1) min *= -1;
//          printf ("min: %g\n", min);
        }

        // Get and verify npak
        fread (b, 1, 3, file);
        int npak = read24(b);
        assert (npak == h.npak);
//        printf ("npak: %d\n", npak);

        // Fast case: pack=16
        if (npak == 16) {
          byte *raw = malloc(2*recsize);
          fread (raw, 2, recsize, file);
          for (int i = 0; i < recsize; i++) {
            out[i] = 1. * read16(raw+2*i) / (1<<16);
          }
          free (raw);
        }

        // Fast case: pack=24
        else if (npak == 24) {
          byte *raw = malloc(3*recsize);
          fread (raw, 3, recsize, file);
          for (int i = 0; i < recsize; i++) {
            out[i] = 1. * read24(raw+3*i) / (1<<24);
          }
          free (raw);
        }

        // Fast case: pack=32
        else if (npak == 32) {
          byte *raw = malloc(4*recsize);
          fread (raw, 4, recsize, file);
          for (int i = 0; i < recsize; i++) {
            out[i] = 1. * read32(raw+4*i) / (1LL<<32);
          }
          free (raw);
        }

        // Slow case: other packing density
        else {
          printf ("warning: slow unpacking!\n");
          // Read the data into a buffer
          byte *raw = malloc(recsize * npak / 8);
          byte *bits = malloc(recsize * npak);
          unsigned int *codes = malloc(sizeof(unsigned int) * recsize);
          // First, read in bytes
          fread (raw, 1, recsize*npak/8, file);
          // Next, expand bytes into bit array
          for (int i = 0; i < recsize * npak / 8; i++) {
            byte x = raw[i];
            for (int j = 7; j >= 0; j--) {
              bits[i*8+j] = x & 0x01;
              x >>= 1;
            }
          }
          // Now, collapse this back into an integer code of size <npak>
          for (int i = 0; i < recsize; i++) {
            unsigned int x = 0;
            for (int j = 0; j < npak; j++) x = (x<<1) + bits[i*npak + j];
            codes[i] = x;
//            printf ("%d ", x);
          }
//          printf ("\n");
          // Decode this into a float
          for (int i = 0; i < recsize; i++) {
            out[i] = codes[i];
            // convert to a value from 0 to 1
            for (int j = 0; j < npak; j++) out[i] /= 2;

          }
//          printf ("\n");
          free (codes);
          free (bits);
          free (raw);
        }

        // Finish decoding the values
        for (int i = 0; i < recsize; i++) {
          //TODO: more robust diff shift
          if (diff_shift >= 0) out[i] *= (1<<diff_shift);
          else out[i] /= (1<<(-diff_shift));

          // apply min
          out[i] += min;
        }

      }


      out += recsize;
  }

  fclose (file);

  return 0;
}

