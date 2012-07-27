#include <stdio.h>
#include <assert.h>
#include "io.h"

int get_file_info (FILE *f, int *p_bgeo_top, int *p_bgeo_siz) {
  fseek (f, 0, SEEK_SET);
  ftn_start (f, 8);
  *p_bgeo_top = fread32(f);
  *p_bgeo_siz = fread32(f);
  ftn_end (f, 8);
  return 1;
}

int read_metadata (FILE *f, int p_bgeo_top, char *geonm1, char *geonm2, char *geonm5, int *geopar1, int *geopar2, int *geopar3) {
  fseek (f, 16, SEEK_SET);  // Skip p_bgeo_top and p_bgeo_siz
  // geonm
  ftn_start (f, p_bgeo_top*16*3);
  fread (geonm1, 16, p_bgeo_top, f);
  fread (geonm2, 16, p_bgeo_top, f);
  fread (geonm5, 16, p_bgeo_top, f);
  ftn_end (f, p_bgeo_top*16*3);
  // geopar
  ftn_start (f, p_bgeo_top*4*3);
  for (int i = 0; i < p_bgeo_top; i++) geopar1[i] = fread32(f);
  for (int i = 0; i < p_bgeo_top; i++) geopar2[i] = fread32(f);
  for (int i = 0; i < p_bgeo_top; i++) geopar3[i] = fread32(f);
  ftn_end (f, p_bgeo_top*4*3);
  return 1;
}

int test_base_offset (FILE *f, int offset, int p_bgeo_siz) {
  fseek (f, offset, SEEK_SET);
  ftn_start (f, p_bgeo_siz*4);
  fseek (f, p_bgeo_siz*4, SEEK_CUR);
  ftn_end (f, p_bgeo_siz*4);
  return 1;
}

int read_data(FILE *f, int offset, int size, float *out) {
  fseek(f, offset, SEEK_SET);
  for (int i = 0; i < size; i++) out[i] = freadfloat(f);
  return 1;
}
