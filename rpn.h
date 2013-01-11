typedef struct {
  long long int file_size;
  int num_overwrites;
  int num_extensions;
  int nchunks;
  unsigned long long int last_chunk;
  int max_data_length;
  int num_erasures;
  int nrecs;
} FileHeader;


typedef struct {
  unsigned int this_chunk_words;
  unsigned long long int this_chunk;
  unsigned int next_chunk_words;
  unsigned long long int  next_chunk;
  int nrecs;
  unsigned int checksum;
} ChunkHeader;


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

