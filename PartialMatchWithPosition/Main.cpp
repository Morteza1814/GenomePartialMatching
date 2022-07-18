#include "PartialMatcher.hpp"
#include "MakeIndexTable.hpp"
#include <time.h>

map<string, map<int, vector<int> > > indexTable; //a map of kmers, positions, indexes

void printConfigs ()
{
	cout << "**********************CONFIGURATIONS*****************************" << endl;
	cout << "READ_SZ : " << READ_SZ << endl;
	cout << "KMER_SZ : " << KMER_SZ << endl;
	cout << "NUMBER_OF_KMERS : " << NUMBER_OF_KMERS << endl;
	cout << "NUMBER_OF_ALLOWED_EDIT_DISTANCES : " << NUMBER_OF_ALLOWED_EDIT_DISTANCES << endl;
	cout << "ACCEPTABLE_PARTIAL_MATCH_REGION : " << ACCEPTABLE_PARTIAL_MATCH_REGION << endl;
	cout << "ADJACENT_KMERS_WINDOW_SZ : " << ADJACENT_KMERS_WINDOW_SZ << endl;
	cout << "NUMBER_OF_WINDOWS_IN_READ : " << NUMBER_OF_WINDOWS_IN_READ << endl;
	cout << "UNKNOWN_REF_FILE_ADDRESS : " << UNKNOWN_REF_FILE_ADDRESS << endl;
	cout << "QUERY_FILE_ADDRESS : " << QUERY_FILE_ADDRESS << endl;
	cout << "NUMBER_OF_READS_FROM_UNKNOWN_REF_FILE : " << NUMBER_OF_READS_FROM_UNKNOWN_REF_FILE << endl;
	cout << "NUMBER_OF_QUERIES_FROM_QUERY_FILE : " << NUMBER_OF_QUERIES_FROM_QUERY_FILE << endl;
	cout << "*****************************************************************" << endl;
}

int main(){
	printConfigs();
	clock_t timeOfIndexTable = clock();
	makeIndexTable(NUMBER_OF_READS_FROM_UNKNOWN_REF_FILE);
	cout << "--------------------------------------" << endl;
	printf("\ntimeOfIndexTable: %.2fs\n", (double)(clock() - timeOfIndexTable)/CLOCKS_PER_SEC);
	cout << "--------------------------------------" << endl;
	clock_t timeOfPartialMatch = clock();
	partialMatcher(NUMBER_OF_QUERIES_FROM_QUERY_FILE);
	cout << "--------------------------------------" << endl;
	printf("timeOfPartialMatch: %.2fs\n", (double)(clock() - timeOfPartialMatch)/CLOCKS_PER_SEC);
	cout << "--------------------------------------" << endl;
	return 0;
}