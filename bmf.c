#include <stdio.h>
#include <assert.h>

typedef unsigned char byte;

unsigned int read32 (byte *b) {
  return (((unsigned int)(b[0]))<<24) | (((int)(b[1]))<<16) | (((int)(b[2]))<<8) | (((int)(b[3]))<<0);
}

unsigned int fread32 (FILE *f) {
  byte b[4];
  int n = fread (b, 1, 4, f);
  assert (n == 4);
  return read32(b);
}

float freadfloat (FILE *f) {
  int x = fread32(f);
  return *((float*)(&x));
}

void ftn_section (FILE *f, int n) {
  int m = fread32(f);
  assert (m == n);
}

#define ftn_start ftn_section
#define ftn_end ftn_section

typedef struct {
  char nom[5];
  int ni, istart, iend;
  int nj, jstart, jend;
  int nk, kstart, kend;
  int time1, time2;
  int hgrid, vgrid;
  int dtyp, scat;
  int ndata;
  long long data_start;
} BMF_Header;

int read_header (FILE *f, BMF_Header *h) {
//  ftn_start(f,8);
  byte test[4];
  fread (test, 1, 4, f);
  if (feof(f)) return 0;
  assert (read32(test) == 8);
    int head_size = fread32(f);
    assert (head_size == 56);
    fread (h->nom, 1, 4, f);
    h->nom[4] = 0;
  ftn_end(f,8);
  ftn_start(f,12);
    h->ni = fread32(f);
    h->istart = fread32(f);
    h->iend = fread32(f);
  ftn_end(f,12);
  ftn_start(f,12);
    h->nj = fread32(f);
    h->jstart = fread32(f);
    h->jend = fread32(f);
  ftn_end(f,12);
  ftn_start(f,12);
    h->nk = fread32(f);
    h->kstart = fread32(f);
    h->kend = fread32(f);
  ftn_end(f,12);
  ftn_start(f,8);
    h->time1 = fread32(f);
    h->time2 = fread32(f);
  ftn_end(f,8);
  ftn_start(f,8);
    h->hgrid = fread32(f);
    h->vgrid = fread32(f);
  ftn_end(f,8);
  ftn_start(f,8);
    h->dtyp = fread32(f);
    h->scat = fread32(f);
  ftn_end(f,8);
  ftn_start(f,4);
    h->ndata = fread32(f);
  ftn_end(f,4);
  ftn_start(f,4);
    head_size = fread32(f);
    assert (head_size == 56);
  ftn_end(f,4);
  h->data_start = ftello64(f);
  return 1;
}

void print_header (BMF_Header *h) {
  printf ("%s (%d,%d,%d)  [%d..%d, %d..%d, %d..%d]\n", h->nom, h->ni, h->nj, h->nk, h->istart, h->iend, h->jstart, h->jend, h->kstart, h->kend);
  printf ("  %08d %08d\n", h->time1, h->time2);
  printf ("  %d %d %d %d\n", h->hgrid, h->vgrid, h->dtyp, h->scat);
}

void skip_data (FILE *f) {
  ftn_start(f,4);
    int data_size = fread32(f);
  ftn_end(f,4);
  ftn_start(f,data_size*4);
    fseek(f, data_size*4, SEEK_CUR);
  ftn_end(f,data_size*4);
  ftn_start(f,4);
    int data_size2 = fread32(f);
    assert (data_size2 == data_size);
  ftn_end(f,4);
}

int num_records (FILE *f) {
  int nrec = 0;
  BMF_Header h;
  fseek(f,0,SEEK_SET);
  while (read_header (f, &h)) {
    nrec ++;
    skip_data(f);
  }
  return nrec;
}

int read_all_headers (FILE *f, BMF_Header *h) {
  fseek(f,0,SEEK_SET);
  while (read_header(f, h++)) skip_data(f);
  return 1;
}

int read_data (FILE *f, BMF_Header *h, float *out) {
  fseeko64(f, h->data_start, SEEK_SET);
  ftn_start(f,4);
    int data_size = fread32(f);
    assert (data_size == h->ndata);
  ftn_end(f,4);
  ftn_start(f,data_size*4);
    for (int i = 0; i < h->ndata; i++) out[i] = freadfloat(f);
//    fread(out, 4, h->ndata, f);
  ftn_end(f,data_size*4);
  ftn_start(f,4);
    int data_size2 = fread32(f);
    assert (data_size2 == data_size);
  ftn_end(f,4);
  return 1;
}


int main (int argc, char *argv[]) {
  char filename[] = "bm20090102000000-00-00";
  FILE *f = fopen(filename, "rb");
  printf ("nrec: %d\n", num_records(f));
  fclose(f);
  return 0;
}
