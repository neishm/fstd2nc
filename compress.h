#include <stdio.h>
#include "rpn.h"

// Read a compressed field
int read_compress32 (unsigned char *buf, RecordHeader *header, int recsize, float *out);
