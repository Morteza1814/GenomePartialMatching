#include "PartialMatcher.hpp"

//log parameters
int numOfKmersFromQueryFoundInIndexTable = 0;

void printindexesForEachKmerInQuery(int start, int end, map<int, map <string, vector<int> > > *indexesForEachKmerInQuery)
{
	cout << "\n\n*****index array*****" << endl;
	for (int i = start; i < end; i++)
	{
		cout << "position = " << i << " : ";
		for (map<int, map <string, vector<int> > >::iterator it = indexesForEachKmerInQuery[i].begin(); it != indexesForEachKmerInQuery[i].end(); it++)
		{
			cout << "index [" << it->first << "] = ";
			for (map <string, vector<int> >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			{
				cout << "kmer (" << it2->first << ") -> positions = ";
				for (vector<int>::iterator it3 = it2->second.begin(); it3 !=it2->second.end(); it3++)
				{
					cout << *it3 << ", ";
				}
			}
			cout << endl;
		}	
	}
}

void fillRelevantIndexsAndPositionsArrayForQuery(string query, map<int, map <string, vector<int> > > *indexesAndPositionsPairForTheKmersInQuery)
{
	int startPos = 0;
	string kmer;
	int numberOfItemsInindexForQuery = 0, numberOfItemsForQuery = 0; 
	while (startPos < NUMBER_OF_KMERS)
	{
		kmer = query.substr(startPos, KMER_SZ);
		if (indexTable.find(kmer) != indexTable.end())
		{		
			numOfKmersFromQueryFoundInIndexTable++;
			for (map<int, vector<int> >::iterator itr = indexTable[kmer].begin(); itr != indexTable[kmer].end(); itr++) 
			{
				map <string, vector<int> > tmp;
				tmp.insert(make_pair(kmer, itr->second));
  				indexesAndPositionsPairForTheKmersInQuery[startPos].insert(make_pair(itr->first, tmp));
				numberOfItemsInindexForQuery++;
    		}
			numberOfItemsForQuery += indexesAndPositionsPairForTheKmersInQuery[startPos].size();
		}
		startPos++;
	}
	// printindexesForEachKmerInQuery(0, 4, indexesForEachKmerInQuery);
	// cout << "\n\n\n should be equal : " << numberOfItemsInindexForQuery << " = " << numberOfItemsForQuery << "\n\n\n";
}


vector<int> getAllInvolvedIndexes(map<int, map <string, vector<int> > > *indexesForEachKmerInQuery)
{
	vector<int> indexes;
	for (int i = 0; i < NUMBER_OF_KMERS; i++)
	{
		if (indexesForEachKmerInQuery[i].empty() == true)
		{
			continue;
		}
		for(map<int, map <string, vector<int> > >::iterator it = indexesForEachKmerInQuery[i].begin(); it != indexesForEachKmerInQuery[i].end(); it++) 
		{
			vector<int>::iterator itt = find(indexes.begin(), indexes.end(), it->first);
			// cout << "\n the index is : " << it->first;
			if (itt == indexes.end())
			{
				// cout << " = added" << endl;
				indexes.push_back(it->first);
			}
 		}
	}
	return indexes;
}

// void fillPositionsFromIndexTable(map<int, map <string, vector<int> > > *indexesAndPositionsPairForTheKmer, int index, int* positions)
// {
// 	map <string, int > repeatedIndexes; 
// 	for (int i = 0; i < NUMBER_OF_KMERS; i++)
// 	{
// 		for(map<int, map <string, vector<int> > >::iterator it = indexesAndPositionsPairForTheKmer[i].begin(); it != indexesAndPositionsPairForTheKmer[i].end(); it++) 
// 		{
// 			if(it->first == index)
// 			{
// 				map <string, vector<int> > kmerPositionsInRead = it->second;
// 				string kmer;
// 				if(kmerPositionsInRead.size() != 1)
// 				{
// 					cout << "\n\nwow.... errrrrrrrr" << endl;
// 					exit(0);
// 				}
// 				kmer = kmerPositionsInRead.begin()->first;
// 				int numberOfTheSameKmers = kmerPositionsInRead[kmer].size();
// 				// cout << "number of the same kmers : " << numberOfTheSameKmers;
// 				if(kmerPositionsInRead[kmer].size() > 1)
// 				{
// 					int numberOfPops = 0;
// 					while(numberOfPops < kmerPositionsInRead[kmer].size())
// 					{
// 						if(repeatedIndexes.find(kmer) != repeatedIndexes.end())
// 						{
// 							numberOfPops = repeatedIndexes[kmer];
// 							repeatedIndexes[kmer] = ++numberOfPops;
// 						}
// 						else
// 						{
// 							repeatedIndexes.insert(make_pair(kmer, 0));
// 						}
// 						positions[i] = kmerPositionsInRead[kmer][numberOfPops];
// 						if (i==0)
// 						{
// 							break;
// 						}
// 						else
// 						{
// 							if(positions[i] >= positions[i-1])
// 							{
// 								break;
// 							}
// 							else
// 							{
// 								positions[i] = -1;
// 								continue;
// 							}
// 						}
// 					}
// 					// cout << "\nposition = " << kmerPositionsInRead[kmer][numberOfPops];
// 				}
// 				else
// 				{
// 					// cout << "\nposition = " << kmerPositionsInRead[kmer].back();
// 					positions[i] = kmerPositionsInRead[kmer].back();
// 				}
// 			}
// 		}
// 	}
// }

// bool checkWindow(int* window)
// {
// 	int dissimilarKmers = 0;
// 	for (int i = 0; i < NUMBER_OF_ADJACENT_MATCHES-1; i++)
// 	{
// 		if((window[i] == -1) || (window[i] != window[i+1]-1))
// 		{
// 			dissimilarKmers++;
// 		}
// 	}
// 	if (dissimilarKmers <= NUMBER_OF_EDIT_DISTANCES)
// 	{
// 		cout << "\n\n*******a match found*********";
// 		cout << "\ndissimilarKmers = " << dissimilarKmers << endl;
// 		return true;
// 	}
// 	return false;
// }

void fillPositionsArrayForEachIndexFromIndexesAndPositionsArray(map<int, map <string, vector<int> > > *indexesAndPositionsPairForTheKmer, int index, vector<int> *positionsForAnIndexInQuery)
{
	//For each kmer position in query :
	for (int i = 0; i < NUMBER_OF_KMERS; i++)
	{
		for(map<int, map <string, vector<int> > >::iterator it = indexesAndPositionsPairForTheKmer[i].begin(); it != indexesAndPositionsPairForTheKmer[i].end(); it++) 
		{
			if(it->first == index)
			{
				map <string, vector<int> > positionsOfKmerFromQueryInRead = it->second;
				string kmer;
				if(positionsOfKmerFromQueryInRead.size() != 1)
				{
					cout << "\n\nwow.... errrrrrrrr" << endl;
					exit(0);
				}
				for(vector<int>::iterator itt = (positionsOfKmerFromQueryInRead.begin()->second).begin(); itt != (positionsOfKmerFromQueryInRead.begin()->second).end(); itt++) 
				{
        			positionsForAnIndexInQuery[i].push_back(*itt);
				}		
			}
		}
	}
}

bool checkPartialMatchForTheIndex(vector<int> *positionsForAnIndexInQuery)
{
	int editDistance = 0;
	int matchedKmers = 0;
	int startPositionOfMatchInQuery = -1;
	int startPositionOfMatchInRead = -1;
	int currentPosition;
	int previousPosition;
	//for each kmer position in query :

	for (int i = 0; i < NUMBER_OF_WINDOWS_IN_READ + NUMBER_OF_ALLOWED_EDIT_DISTANCES; i++)
	{
		//reset
		// repeatedPositions.clear();
		editDistance = 0;
		matchedKmers = 0;
		startPositionOfMatchInQuery = i;

		if (positionsForAnIndexInQuery[i].empty())
		{
			continue;
		}
		//check a window of NUMBER_OF_ADJACENT_MATCHES adjacent kmers starting with position i in query
		for(vector<int>::iterator it = positionsForAnIndexInQuery[i].begin(); it != positionsForAnIndexInQuery[i].end(); it++)
		{
			// repeatedPositions.clear();
			editDistance = 0;
			previousPosition = -1;
			matchedKmers = 1;
			startPositionOfMatchInRead = *it;
			previousPosition = *it;
			// repeatedPositions.push_back(*it);
			//check for adjacency of next kmers in query
			for (int j = i+1; (j < i + ADJACENT_KMERS_WINDOW_SZ) && (j < NUMBER_OF_KMERS); j++)
			{
				currentPosition = -1;
				if(editDistance > NUMBER_OF_ALLOWED_EDIT_DISTANCES)
				{
					break;
				}
				//for each position in query, check all positions of the unknown DB for a possible match to previous one
				for(vector<int>::iterator itt = positionsForAnIndexInQuery[j].begin(); itt != positionsForAnIndexInQuery[j].end(); itt++)
				{
					if(editDistance > NUMBER_OF_ALLOWED_EDIT_DISTANCES)
					{
						break;
					}
					currentPosition = (*itt);
					if ((currentPosition <= previousPosition) || (currentPosition > previousPosition + NUMBER_OF_ALLOWED_EDIT_DISTANCES))
					{
						continue;
					}
					else
					{
						editDistance += currentPosition - previousPosition - 1;
						previousPosition = currentPosition;
						break;
					}
				}
				if (previousPosition != currentPosition)
				{
					editDistance++;
				}
				else
				{
					matchedKmers++;
				}
			}
			if (matchedKmers >= (ADJACENT_KMERS_WINDOW_SZ - NUMBER_OF_ALLOWED_EDIT_DISTANCES))
			{
				if (matchedKmers < ADJACENT_KMERS_WINDOW_SZ && editDistance == 0)
				{
					editDistance = ADJACENT_KMERS_WINDOW_SZ - matchedKmers;
				}
				
				cout << "\nedit distance : " << editDistance << ", number of matched kmers : " << matchedKmers << ", start Position Of Match in Query : " << startPositionOfMatchInQuery << ", start Position Of Match in Read : " << startPositionOfMatchInRead;
				return true;
			}
		}
	}
	return false;
}

int partialMatcher(int numberOfQueries)
{
	ifstream file(QUERY_FILE_ADDRESS);
	string kmer, query;
	int queryIndex = 0, numOfPartialMatches = 0, numberOfQueriesMatched = 0;
	map<int, map <string, vector<int> > > indexesAndPositionsPairForTheKmersInQuery[NUMBER_OF_KMERS];
	vector<int> positionsForAnIndexInQuery[NUMBER_OF_KMERS];
	bool partialMatch = false;
	cout << "\n*********************PARTIAL MATCH***********************" << endl;
	//For each query : 
	while (queryIndex < numberOfQueries)
	{
		getline(file, query);
		if (file.eof())
		{
			cout << "queries file reached the end!!" << endl;
			break;
		}
		//reset
		for (int i = 0; i < NUMBER_OF_KMERS; i++)
		{
			indexesAndPositionsPairForTheKmersInQuery[i].clear();
		}

		// 	cout<<"main str: " << str << "\n";
		if (query.length() < READ_SZ)
		{
			cout << "\n\nquery error";
			file.close();
			return 0;
		}
		//for each kmer in the query: search in the index table for indexes and positions and 
		//make an array containing indexes and positions from index table to the positions of query
		fillRelevantIndexsAndPositionsArrayForQuery(query, indexesAndPositionsPairForTheKmersInQuery);
		//Make an array containing all the different indexes of unknown read DB that involved in this query
		vector<int> indexes = getAllInvolvedIndexes(indexesAndPositionsPairForTheKmersInQuery);
		//for each index:
		for(vector<int>::iterator it = indexes.begin(); it != indexes.end(); it++) 
		{
			//reset
			for (int i = 0; i < NUMBER_OF_KMERS; i++)
			{
				positionsForAnIndexInQuery[i].clear();
			}
			//For each index fill the corresponding positions of query (from index table)
			fillPositionsArrayForEachIndexFromIndexesAndPositionsArray(indexesAndPositionsPairForTheKmersInQuery, *it, positionsForAnIndexInQuery);

			if(checkPartialMatchForTheIndex(positionsForAnIndexInQuery))
			{
				numOfPartialMatches++;
				partialMatch = true;
				cout << " -> index = " << (*it);
			}
		}
		if (partialMatch)
		{
			cout << " => query : " << query << endl;
			numberOfQueriesMatched++;
			partialMatch = false;
		}
		queryIndex++;
	}
	cout << "**************************************************************" << endl;
	cout << "number of partial matches : " << numOfPartialMatches << endl;
	cout << "number of queries matched : " << numberOfQueriesMatched << endl;
	cout << "num Of Kmers In Query Found In Index Table : " << numOfKmersFromQueryFoundInIndexTable << endl;
	cout << "index table size : " << indexTable.size() << endl;
	cout << "**************************************************************" << endl;
	file.close();
	return 1;
}
