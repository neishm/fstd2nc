#include <stdio.h>

typedef unsigned char byte;

// Read a 32-bit integer from a buffer
unsigned int read32 (byte *b);

// Read a 24-bit integer from a buffer
int read24 (byte *b);

// Read a 16-bit integer from a buffer
int read16 (byte *b);

// Read a 32-bit integer from a stream
unsigned int fread32 (FILE *f);

// Read a 32-bit float from a stream
float freadfloat (FILE *f);

// Validate the byte count at the beginning/end of an unformatted Fortran stream
void ftn_section (FILE *f, int n);

#define ftn_start ftn_section
#define ftn_end ftn_section


