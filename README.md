# GenomePartialMatching
This repository includes all the source code files that are used in DNA partial matching project. The pogramming languages used are C++ and python.

## KmerFrequencyHistogram
 This project was about analyzing the kmer frequency in a real genome. It consists of 3 phases:
   - Kmerizing each read in the genome file and counted the (global) number of each specific kmer in all the reads. 
   - Counting similar counts to attain the frequency of each count. 
   - Creating a diagram based on the range of the count frequencies, i.e., 0-2k, 2k-5k, 5k-25k, 25k-50k, and more than 50k. 

## NaivePartialMatcher
This project created to evaluate the idea of Bit Vector Matchig. It has 3 steps:
   - Sampling 10k reads from the Metagenomic database and sampling 1000 of them as queries. 
   - Kmerizing each query and searched for kmers in each read of the reference database. 
   - If the number of matched kmers equals the total number of kmers in a query (length_of_query - kmer_length + 1), the query is accepted as a "Bit Vector Method" match
   - If found a match, compare the query and the read bp by bp to see if it is a True match. 
   - Finally, I counted the number of "Bit Vector Method" matches and True matches. By subtracting these two numbers, we will have the number of false positives.

## BitVectorPartialMatcher
This project is the exact implementation of Bit Vector Idea.
   - Checking the number of the matched kmers between the query and the read. 
   - If the [number of matched kmers] equals [the number of kmers in a read - (edit_distance*kmer_length)], accept it as a partial match.
   - To calculate the number of True partial matches, I compared the read and the query bp to bp. 

## PartialMatchWithPosition
This project conducts the partial matching between read database and queries by finding the same adjacent kmers in a read and queries using the position of kmers.
- Making the Index Table (hash map) whose keys are kmers and the values are the positions of the kmer in all the reads in the read database.
- If the number of adjacent kmers in read and query are more than a threshold, it is partial match
