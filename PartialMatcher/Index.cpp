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
    int nonzeroPatterns = 0;
    if(strlen(contig) != CONTIG_SIZE)
    {
        cout << "Error : Contig size :" << strlen(contig) << " is not vlaid!" << endl;
    }

    for(int i=0 ; i <= CONTIG_SIZE-KMER_SIZE ; i++)
    {
        strncpy(kmer, contig+i, KMER_SIZE);
        kmer[KMER_SIZE]='\0';
        int pattern = getBitVector(kmer, bitVector, i+(KMER_SIZE/2));
        // if (pattern != 0)
        // {
        //     cout << "i = " << i <<"\t pattern : " << pattern << endl;
        // }
        // nonzeroPatterns++;   
    }
    // cout << "nonzeroPatterns : " << nonzeroPatterns << endl;
}

unsigned char ** createBitVector(int numberOfContigs, int initialValue, int numberOfBytesInBitVector)
{
    unsigned char **bitVector = new unsigned char*[numberOfContigs];
    for (int i = 0; i < numberOfContigs; i++)
    {
        bitVector[i] = new unsigned char[numberOfBytesInBitVector];
        for (int j = 0; j < numberOfBytesInBitVector; ++j) {
            bitVector[i][j] = initialValue;
        }
    }
    return bitVector;
}


int main()
{
    char contig[151];
    string line;
    ifstream readFile ("/bigtemp/rgq5aw/naiveBitVectorData/sampledDataFromSRR.fa");
    ifstream queryFile ("/bigtemp/rgq5aw/naiveBitVectorData/sampledQueryDataFromSRR.fa");
    // ifstream readFile("in.txt");
    // ifstream queryFile("q.txt");
    int readCount = 0;
    int queryCount = 0;
    unsigned char **readsBitVector = createBitVector(NUMBER_OF_READS_CONTIGS, 0, NUMBER_OF_BIT_VECTOR_BYTES);
    unsigned char **queriesBitVector = createBitVector(NUMBER_OF_QUERIES_CONTIGS, 0, NUMBER_OF_BIT_VECTOR_BYTES);
    string readsArr[NUMBER_OF_READS_CONTIGS];
    string queriesArr[NUMBER_OF_QUERIES_CONTIGS];
    ofstream out("output.txt");
    
 
    if (readFile.is_open())
    {
        while (getline (readFile, line))
        {
            readsArr[readCount] = line;
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
            queriesArr[queryCount] = line;
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

    cout << "number of reads : " << readCount << " number of queries : " << queryCount << endl;
    out << "number of reads : " << readCount << " number of queries : " << queryCount << endl;
    
    unsigned char query[NUMBER_OF_BIT_VECTOR_BYTES];
    int similarKmers = 0;
    int similarPairs = 0;
    int similarReadsForQuery = 0;
    int trueSimilarPairs = 0;
    int trueSimilarReadsForQuery = 0;
    int firstHalfPatternCheckResult = 0;
    int secondHalfPatternCheckResult = 0;
    
    for (int i = 0; i < queryCount; i++)
    {
        out << "query : " << i << "\t["<< queriesArr[i] << "]" << endl;
        firstHalfPatternCheckResult = 0;
        secondHalfPatternCheckResult = 0;
        similarReadsForQuery = 0;
        trueSimilarReadsForQuery = 0;
        for (int k = 0; k < readCount; k++)
        {
            similarKmers=0;
            for (int l = 0; l < NUMBER_OF_BIT_VECTOR_BYTES; l++)
            {
                firstHalfPatternCheckResult = (readsBitVector[k][l] & 0x07) & (queriesBitVector[i][l] & 0x07);
                secondHalfPatternCheckResult = ((readsBitVector[k][l] >> NUMBER_OF_BITS_FOR_EACH_KMER) & 0x07) & ((queriesBitVector[i][l] >> NUMBER_OF_BITS_FOR_EACH_KMER)& 0x07);
                if (firstHalfPatternCheckResult != 0)
                {
                    similarKmers++;
                    // cout << ", " << 2*l;

                }
                if (secondHalfPatternCheckResult != 0)
                {
                    similarKmers++;
                    // cout << ", " << 2*l+1;
                }
                
            }
            // cout << "the same kmers:" <<theSameKmers << endl;

            if (similarKmers == (CONTIG_SIZE - KMER_SIZE + 1))
            {
                out << "read : " << k << "\t [" << readsArr[k] << "]" <<endl;
                similarReadsForQuery++;
                if (queriesArr[i].compare(readsArr[k]) == 0)
                {
                    trueSimilarReadsForQuery++;
                    // cout << "true!" << endl;
                }
                
            }         
        } 
        out << "query : " << i << "\tsimilarReadsForQuery : " << similarReadsForQuery << "\ttrueSimilarReadsForQuery : " << trueSimilarReadsForQuery <<endl;
        cout << "query : " << i << "\tsimilarReadsForQuery : " << similarReadsForQuery << "\ttrueSimilarReadsForQuery : " << trueSimilarReadsForQuery <<endl;
        similarPairs += similarReadsForQuery;
        trueSimilarPairs += trueSimilarReadsForQuery;
    }
    out << "the same pairs : " << similarPairs << "\t true similar pairs : "<< trueSimilarPairs << endl;
    cout << "the same pairs : " << similarPairs << "\t true similar pairs : "<< trueSimilarPairs << endl;

    //Free each sub-array
    for(int i = 0; i < NUMBER_OF_READS_CONTIGS; i++) {
        delete[] readsBitVector[i];   
    }
    for(int i = 0; i < NUMBER_OF_QUERIES_CONTIGS; i++) {
        delete[] queriesBitVector[i];   
    }
    //Free the array of pointers
    delete[] readsBitVector;
    delete[] queriesBitVector;
    out.close();
    return 0;
}
