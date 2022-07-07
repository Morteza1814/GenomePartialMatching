#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <map>

#define KMER_SIZE 10
#define CONTIG_SIZE 150


using namespace std;

int NUMBER_OF_POSSIBLE_KMERS = ((pow(4, KMER_SIZE)));

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

void kmerizeAndCountKmers(const char* contig, map<int,int> &kmerCounts)
{
    char kmer[KMER_SIZE];
    int indexOfKmer = 0;
    map<int,int>::iterator it;
    if(strlen(contig) != CONTIG_SIZE)
    {
        cout << "Error : Contig size :" << strlen(contig) << " is not vlaid!" << endl;
    }
    for(int i=0 ; i <= CONTIG_SIZE-KMER_SIZE ; i++)
    {
        strncpy(kmer, contig+i, KMER_SIZE);
        kmer[KMER_SIZE]='\0';
        indexOfKmer = getIndexOfKmerInDictionary(kmer);
        it = kmerCounts.find(indexOfKmer);
        if (it != kmerCounts.end())
        {
            kmerCounts[indexOfKmer]++;
        } else
        {
            kmerCounts[indexOfKmer] = 1;
        }        
    }
}

int main()
{
    ifstream readFile ("/bigtemp/rgq5aw/batFiles/SRR1246727.fasta");
    string line;
    map<int,int> kmerCounts;
    char contig[151];
    int contigCount = 0;
    int nonAlphCount = 0;
    int totalNumberOfContigs = 0;
    ofstream out("output.txt");

    if (readFile.is_open())
    {
        while (getline (readFile, line))
        {   
            if (!line.empty() && line.find("SRR") == std::string::npos)
            {
                totalNumberOfContigs++;
                
                // if (totalNumberOfContigs>10000)
                // {
                //     cout << "number of all contigs : " << totalNumberOfContigs << endl;
                //     break;
                // }

                if (line.find_first_not_of("ACGTacgt") != std::string::npos)
                {
                    nonAlphCount++;
                    continue;                
                }
                
                strncpy(contig, line.c_str(), CONTIG_SIZE);
                contig[CONTIG_SIZE] = '\0';
                kmerizeAndCountKmers(contig, kmerCounts);
                contigCount++;
            } else
            {
                continue;
            }
            
            if (contigCount % 1000000 == 0)
            {
                cout << contigCount/1000000 << endl;
            }
            
        }
        readFile.close();
    } else cout << "Unable to open file"; 
    cout << "contig count : " << contigCount << "\tnon Alph Count : " << nonAlphCount << endl;
    cout << "number of different kmers : " << kmerCounts.size() << "\tnumber of possible kmers : " << NUMBER_OF_POSSIBLE_KMERS << "\tpercent of kmers : " << (double)(kmerCounts.size()/NUMBER_OF_POSSIBLE_KMERS) << endl;
    out << "contig count : " << contigCount << "\tnon Alph Count : " << nonAlphCount << endl;
    out << "number of different kmers : " << kmerCounts.size() << "\tnumber of possible kmers : " << NUMBER_OF_POSSIBLE_KMERS << "\tpercent of kmers : " << (double)(kmerCounts.size()/NUMBER_OF_POSSIBLE_KMERS) << endl;
    long long totalNumberOfKmers = 0;
    
    for (map<int, int>::iterator it = kmerCounts.begin(); it != kmerCounts.end(); it++)
    {
        totalNumberOfKmers += it->second;
    }
    cout << "total Number Of Kmers : " << totalNumberOfKmers << endl;
    out << "total Number Of Kmers : " << totalNumberOfKmers << endl;

    int maxCount = 1;
    for (map<int, int>::iterator it = kmerCounts.begin(); it != kmerCounts.end(); it++)
    {
        if (it->second > maxCount)
        {
            maxCount = it->second;
        }   
    }
    cout << "maxCount : " << maxCount << endl;
    int *histogramArray = new int[maxCount+1]();
    for (int i = 0; i < maxCount+1; i++)
    {
        int count = 0;
        for (map<int, int>::iterator it = kmerCounts.begin(); it != kmerCounts.end(); it++)
        {
            if (it->second == i)
            {
                count++;
            }
        }
        histogramArray[i]+=count;
    }
    
    cout << "---------------histo array finished------------------" << endl;
    out << "---------------histo array finished------------------" << endl;
    long long total = 0;
    out << "count" <<  "\t, frequency" << endl;
    for (int i = 0; i < maxCount+1; i++)
    {
        out << i << "\t,\t" << histogramArray[i] << endl;
        total += histogramArray[i];
    }
    cout << "total : " << total << endl;
    kmerCounts.clear();
    delete[] histogramArray;
    out.close();
}
