#include <iostream>
#include <fstream>
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
    char contig[151];
    string line;
    // ifstream readFile ("/bigtemp/rgq5aw/naiveBitVectorData/sampledDataFromSRR.fa");
    // ifstream queryFile ("/bigtemp/rgq5aw/naiveBitVectorData/sampledQueryDataFromSRR.fa");
    ifstream readFile("in.txt");
    ifstream queryFile("q.txt");
    int readCount = 0;
    int queryCount = 0;
    unsigned char **readsBitVector = new unsigned char*[NUMBER_OF_READS_CONTIGS];
    unsigned char **queriesBitVector = new unsigned char*[1000];
    
    for (int i = 0; i < NUMBER_OF_READS_CONTIGS; i++)
    {
        readsBitVector[i] = new unsigned char[NUMBER_OF_BIT_VECTOR_BYTES];
    }
        
    for (int i = 0; i < NUMBER_OF_QUERIES_CONTIGS; i++)
    {
        queriesBitVector[i] = new unsigned char[NUMBER_OF_BIT_VECTOR_BYTES];
    }

    if (readFile.is_open())
    {
        while (getline (readFile, line))
        {
            // cout << line << '\n';
            if (line.find("chr") == std::string::npos && !line.empty())
            {
                strncpy(contig, line.c_str(), CONTIG_SIZE);
                contig[CONTIG_SIZE] = '\0';
                kmerizeAndMakeBitVector(contig, readsBitVector[readCount]);
                readCount++;
            }
        }
        readFile.close();
    } else cout << "Unable to open file"; 

    if (queryFile.is_open())
    {
        while (getline (queryFile, line))
        {
            // cout << line << '\n';
            if (line.find("chr") == std::string::npos && !line.empty())
            {
                strncpy(contig, line.c_str(), CONTIG_SIZE);
                contig[CONTIG_SIZE] = '\0';
                kmerizeAndMakeBitVector(contig, queriesBitVector[queryCount]);
                queryCount++;
            }
        }
        queryFile.close();
    } else cout << "Unable to open file"; 

    cout << "kmer size is : " << KMER_SIZE << endl;
    cout << "number of reads : " << readCount << " number of queries : " << queryCount << endl;
    unsigned char query[NUMBER_OF_BIT_VECTOR_BYTES];
    int theSameKmers=0;
    for (int i = 0; i < queryCount; i++)
    {
        // for (int j = 0; j < NUMBER_OF_BIT_VECTOR_BYTES; j++)
        // {
        //     query[i] = queriesBitVector[i][j];
        // }
        for (int k = 0; k < readCount; k++)
        {
            theSameKmers=0;
            for (int l = 0; l < NUMBER_OF_BIT_VECTOR_BYTES; l++)
            {
                if(l==0)  
                    cout << "read : " << (int)readsBitVector[k][l] << "\t query : " << (int)query[l] << endl;    
                if (readsBitVector[k][l] & queriesBitVector[i][l] != 0)
                {
                    theSameKmers++;
                }
            }
            cout << "the same kmers:" <<theSameKmers << endl;
            if (theSameKmers == (CONTIG_SIZE - KMER_SIZE + 1))
            {
                cout << "read : " << k << "\t query : " << i << endl;
            }
                      
        }
        cout << " query : " << i << endl;
        
    }
    
    // cout << "number of bit vector bytes is : " << NUMBER_OF_BIT_VECTOR_BYTES << endl;
    
    // strcpy(contig, "GCCCTGAT?GAATTACCTCGTCTTTTCTCATATAACATGTCCTGGGAAGCCACAACATTGTGGTAAAGCTGTTCAACTACCACCGATACCATAGCAAGATGCTCATCTAGACTGTGACGACAATTACATCGTGAGAGATTGTGCTCTAGGC");
    // kmerizeAndMakeBitVector(contig, bitVector);
    // kmerizeAndReadBitVector(contig, bitVector);
    // for(int i=0 ; i<NUMBER_OF_BIT_VECTOR_BYTES ; i++)
    // {
    //     cout << i << " = " << (unsigned int)bitVector[i] << endl;
    // }
    
    return 0;
}
