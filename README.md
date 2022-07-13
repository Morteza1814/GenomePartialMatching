# GenomePartialMatching
This repository includes all the source code files that are used in DNA partial matching project. The pogramming languages used are C++ and python.

1- KmerFrequencyHistogram : This project was about analyzing the kmer frequency in a real genome. It consists of 3 phases:
  - Kmerizing each read in the genome file and counted the (global) number of each specific kmer in all the reads. 
  - Counting similar counts to attain the frequency of each count. 
  - Creating a diagram based on the range of the count frequencies, i.e., 0-2k, 2k-5k, 5k-25k, 25k-50k, and more than 50k. 
