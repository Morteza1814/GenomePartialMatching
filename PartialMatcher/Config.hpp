#include<math.h>

#define ALPHABET_SIZE 4
#define KMER_SIZE 4
#define NUMBER_OF_BITS_FOR_EACH_KMER 4
#define NUMBER_OF_REGIONS 3

int NUMBER_OF_BIT_VECTOR_BYTES = ((4*pow(4, KMER_SIZE))/8);
