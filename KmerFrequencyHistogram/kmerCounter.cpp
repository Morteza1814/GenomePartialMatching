#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <map>
#include <vector>
#include <algorithm>

#define KMER_SIZE 16
#define CONTIG_SIZE 150


using namespace std;

unsigned long NUMBER_OF_POSSIBLE_KMERS = ((pow(4, KMER_SIZE)));

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

bool cmp(int a, int b)
{
    return a < b;
}

vector<int> sort(map<int, int>& M, ofstream& out)
{  
    vector<int> A;
    long long totalNumberOfKmers = 0;

    // Copy key-value pair from Map
    // to vector of pairs
    for (map<int, int>::iterator it = M.begin(); it != M.end(); it++) {
        A.push_back(it->second);
        totalNumberOfKmers += it->second;
    }
    cout << "total Number Of Kmers : " << totalNumberOfKmers << endl;
    out << "total Number Of Kmers : " << totalNumberOfKmers << endl;
  
    // Sort using comparator function
    sort(A.begin(), A.end(), cmp);
  
    // Print the sorted value
    // for (auto& it : A) {
  
    //     cout << it.first << ' '
    //          << it.second << endl;
    // }
    return A;
}

int main()
{
    ifstream readFile ("/bigtemp/rgq5aw/batFiles/SRR1246727.fasta");
    string line;
    map<int,int> kmerCounts;
    vector<int> kmerCountsVec;
    char contig[151];
    unsigned int contigCount = 0;
    unsigned int nonAlphCount = 0;
    unsigned int totalNumberOfContigs = 0;
    ofstream out("output.txt");
    cout << "kmer size : " << KMER_SIZE << endl;
    out << "kmer size : " << KMER_SIZE << endl;
    double tStart = (double)(clock());
    if (readFile.is_open())
    {
        while (getline (readFile, line))
        {   
            if (!line.empty() && line.find("SRR") == std::string::npos)
            {
                totalNumberOfContigs++;
                
                // if (totalNumberOfContigs>1000)
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
    cout << "timing for counting kmers : " << (double)(clock() - tStart) << endl;
    cout << "contig count : " << contigCount << "\tnon Alph Count : " << nonAlphCount << endl;
    cout << "number of different kmers : " << kmerCounts.size() << "\tnumber of possible kmers : " << NUMBER_OF_POSSIBLE_KMERS << "\tpercent of kmers : " << (double)(kmerCounts.size()/NUMBER_OF_POSSIBLE_KMERS) << endl;
    //dump to file
    out << "contig count : " << contigCount << "\tnon Alph Count : " << nonAlphCount << endl;
    out << "number of different kmers : " << kmerCounts.size() << "\tnumber of possible kmers : " << NUMBER_OF_POSSIBLE_KMERS << "\tpercent of kmers : " << (double)(kmerCounts.size()/NUMBER_OF_POSSIBLE_KMERS) << endl;
    out << "timing for counting kmers : " << (double)(clock() - tStart) << endl;

    tStart = (double)(clock());
    kmerCountsVec = sort(kmerCounts, out);
    kmerCounts.clear();
    int numberOfKmers = kmerCountsVec.size();
    cout << "numberOfKmers : " << numberOfKmers << endl;
    cout << "maximum vector count : " << kmerCountsVec[numberOfKmers-1] << endl;
    cout << "minimum vector count : " << kmerCountsVec[0] << endl;
    //dump to file
    out << "numberOfKmers : " << numberOfKmers << endl;
    out << "maximum vector count : " << kmerCountsVec[numberOfKmers-1] << endl;
    out << "minimum vector count : " << kmerCountsVec[0] << endl;

    map<int, int> histogramMap;
    int currentCount = kmerCountsVec[0];
    int count = 1; 
    histogramMap[kmerCountsVec[0]] = 1;
    for (int i = 1; i < numberOfKmers; i++)
    {
        if (kmerCountsVec[i] == kmerCountsVec[i-1])
        {
            histogramMap[kmerCountsVec[i]]++;
        }
        else
        {
            histogramMap[kmerCountsVec[i]] = 1;
        }     
    } 
    out << "-------------------------------------" << endl;
    cout << "timing for counting frequeincies : " << (double)(clock() - tStart) << endl;
    out << "timing for counting frequeincies : " << (double)(clock() - tStart) << endl;

    long long total = 0;
    out << "count" <<  "\t, frequency" << endl;
    for (map<int, int>::iterator it = histogramMap.begin(); it != histogramMap.end(); it++)
    {
        out << it->first << "\t,\t" << it->second << endl;
        total += (it->first * it->second);
    }
    cout << "total histogram kmers: " << total << endl;
    cout << "different kmer counts : " << histogramMap.size() << endl;
    histogramMap.clear();
    out.close();
}
