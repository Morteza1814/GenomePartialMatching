#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Config.hpp"

using namespace std;

unsigned long getIndexOfKmerInDictionary(const char* kmer)
{
    if((unsigned)strlen(kmer) != KMER_SIZE)
    {
        cout << "Error : kmer length is: " << (unsigned)strlen(kmer) << endl;
        return -1;
    }

    unsigned long alphabetInTwoBits = 0;
    unsigned long kmerInBits = 0 ;
    unsigned int j=0;
    for(int i=KMER_SIZE-1 ; i>=0 ; i--)
    {
        switch (kmer[i])
        {
        case 'A':
            alphabetInTwoBits = 0;
            break;
        case 'C':
            alphabetInTwoBits = 1;
            break;
        case 'G':
            alphabetInTwoBits = 2;
            break;
        case 'T':
            alphabetInTwoBits = 3;
            break;
        default:
            cout << "Error : Alphabet :" << kmer[i] << " deos not exist!" << endl;
            return -1;
            break;
        }
        kmerInBits |= alphabetInTwoBits << (j*2);
        j++;
    }
    return kmerInBits;
}

void setBitVector(const char* kmer, unsigned char* bitVector, int kmerPos)
{
    // printf("kmer = %s \t", kmer);
    int indexOfKmer = getIndexOfKmerInDictionary(kmer);
    // cout << indexOfKmer << endl;
    int indexInBitVector = indexOfKmer / 2;
    int firstOrSecondHalf = indexOfKmer % 2;
    unsigned int kmerRegionBitPattern = (unsigned int)pow(2, (kmerPos/(CONTIG_SIZE/NUMBER_OF_REGIONS)));
    bitVector[indexInBitVector] |= (unsigned char)(((kmerRegionBitPattern & 0x07) << (firstOrSecondHalf*NUMBER_OF_BITS_FOR_EACH_KMER)));
    // cout << kmerRegionBitPattern << endl;
    // cout << indexInBitVector << endl;
    // cout << (unsigned int)bitVector[indexInBitVector] << endl;
}

int getBitVector(const char* kmer, unsigned char* bitVector, int kmerPos)
{
    printf("kmer = %s \t", kmer);
    int indexOfKmer = getIndexOfKmerInDictionary(kmer);
    cout << indexOfKmer << endl;
    int indexInBitVector = indexOfKmer / 2;
    int firstOrSecondHalf = indexOfKmer % 2;
    return (unsigned char)((bitVector[indexInBitVector] >> (firstOrSecondHalf*NUMBER_OF_BITS_FOR_EACH_KMER)) & 0x07);
}

void kmerizeAndMakeBitVector(const char* contig, unsigned char* bitVector)
{
    char kmer[KMER_SIZE];
    if(strlen(contig) != CONTIG_SIZE)
    {
        cout << "Error : Contig size :" << strlen(contig) << " is not vlaid!" << endl;
    }

    for(int i=0 ; i <= CONTIG_SIZE-KMER_SIZE ; i++)
    {
        strncpy(kmer, contig+i, KMER_SIZE);
        kmer[KMER_SIZE]='\0';
        setBitVector(kmer, bitVector, i+(KMER_SIZE/2));
    }
}

void kmerizeAndReadBitVector(const char* contig, unsigned char* bitVector)
{
    char kmer[KMER_SIZE];
    if(strlen(contig) != CONTIG_SIZE)
    {
        cout << "Error : Contig size :" << strlen(contig) << " is not vlaid!" << endl;
    }

    for(int i=0 ; i <= CONTIG_SIZE-KMER_SIZE ; i++)
    {
        strncpy(kmer, contig+i, KMER_SIZE);
        kmer[KMER_SIZE]='\0';
        cout << "i = " << i <<"\t : " << getBitVector(kmer, bitVector, i+(KMER_SIZE/2)) << endl;
    }
}

int main()
{
    char contig[150];
    unsigned char bitVector[NUMBER_OF_BIT_VECTOR_BYTES];

    for(int i=0 ; i<NUMBER_OF_BIT_VECTOR_BYTES ; i++)
    {
        bitVector[i] = 0;
    }
    cout << "kmer size is : " << KMER_SIZE << endl;
    cout << "number of bit vector bytes is : " << NUMBER_OF_BIT_VECTOR_BYTES << endl;
    
    strcpy(contig, "GCCCTGATGAATTACCTCGTCTTTTCTCATATAACATGTCCTGGGAAGCCACAACATTGTGGTAAAGCTGTTCAACTACCACCGATACCATAGCAAGATGCTCATCTAGACTGTGACGACAATTACATCGTGAGAGATTGTGCTCTAGGC");
    kmerizeAndMakeBitVector(contig, bitVector);
    // kmerizeAndReadBitVector(contig, bitVector);
    // for(int i=0 ; i<NUMBER_OF_BIT_VECTOR_BYTES ; i++)
    // {
    //     cout << i << " = " << (unsigned int)bitVector[i] << endl;
    // }
    
    return 0;
}
