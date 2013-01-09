#include <stdio.h>
#include "rpn.h"

// Read a compressed field
int read_compress32 (FILE *f, RecordHeader *header, int recsize, float *out);
