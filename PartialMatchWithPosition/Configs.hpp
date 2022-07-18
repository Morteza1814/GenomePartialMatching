#ifndef CONFIGS_H
#define CONFIGS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#define KMER_SZ 32
#define READ_SZ 150
#define ACCEPTABLE_PARTIAL_MATCH_REGION 75
#define NUMBER_OF_KMERS (READ_SZ - KMER_SZ + 1)
#define ADJACENT_KMERS_WINDOW_SZ (ACCEPTABLE_PARTIAL_MATCH_REGION - KMER_SZ + 1)
#define NUMBER_OF_WINDOWS_IN_READ (NUMBER_OF_KMERS - ADJACENT_KMERS_WINDOW_SZ + 1)
#define NUMBER_OF_ALLOWED_EDIT_DISTANCES 2
#define UNKNOWN_REF_FILE_ADDRESS "/bigtemp/rgq5aw/SRR.csv"
// #define UNKNOWN_REF_FILE_ADDRESS "/bigtemp/rgq5aw/testCases/tmpsrc.txt"
#define NUMBER_OF_READS_FROM_UNKNOWN_REF_FILE 999999
#define QUERY_FILE_ADDRESS "/bigtemp/rgq5aw/RMY.csv"
// #define QUERY_FILE_ADDRESS "/bigtemp/rgq5aw/testCases/1_query_30th60th100th_BP_different"
#define NUMBER_OF_QUERIES_FROM_QUERY_FILE 180

using namespace std;

#endif