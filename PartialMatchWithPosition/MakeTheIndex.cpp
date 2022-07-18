#include "MakeIndexTable.hpp"

int makeIndexTable(int numberOfReads)
{
	ifstream file(UNKNOWN_REF_FILE_ADDRESS);
	//  ofstream MyExcelFile;
	string read;
	int startPos = 0, readIndex = 0;
	string kmer;
	int repeatedKmers = 0, numberOfKmersInIndexTable = 0, similarKmersInOneRead = 0;
	int count=0;
	cout << "\n***********************MAKE INDEX********************" << endl;

	while (readIndex < numberOfReads)
	{
		getline(file, read);
		if (file.eof())
		{
			cout << "reads file reached the end!!" << endl;
			break;
		}
		// 	cout<<"main str: " << str << "\n";
		if (read.length() < READ_SZ)
		{
			cout << "\n\nread error";
			file.close();
			return 0;
		}
		startPos = 0;
		while (startPos < NUMBER_OF_KMERS)
		{
			count++;
			kmer = read.substr(startPos, KMER_SZ);
			//	cout << "\n read : " << read;
			if (indexTable.find(kmer) != indexTable.end())
			{
				repeatedKmers++;
				if(indexTable[kmer].find(readIndex) != indexTable[kmer].end())
				{
					similarKmersInOneRead++;
				}
			}
			else
			{
				numberOfKmersInIndexTable++;
			}
			indexTable[kmer][readIndex].push_back(startPos); //add the position to the related kmer and the index
			startPos++;
		}
		
		readIndex++;
	}
	// cout << "\nreadindex : " << readIndex;
	// cout << "\n count : " << count;
	// cout << "\n" << NUMBER_OF_KMERS;
	// cout << "\nNumber of reads : " << readIndex << ", total number of kmers : " << (readIndex * NUMBER_OF_KMERS);
	cout << "similar kmers in unidentified ref DB : " << repeatedKmers;
	cout << "\nnumber Of Kmers In Index Table : " << numberOfKmersInIndexTable;
	cout << "\nsimilar Kmers in One Read : " << similarKmersInOneRead; 
	cout << "\ntotal number of items in index table : " << repeatedKmers+numberOfKmersInIndexTable;
	cout << "\nsize of index table : " << indexTable.size() << endl;
	int elemCount = 0, kmerCount = 0;
	for (map < string, map <int , vector<int> > >::const_iterator it = indexTable.begin(); it != indexTable.end(); it++)
	{
		//	cout << "\nread = " << it->first << ": ";
		kmerCount++; //number of different reads in the index table
		for (map<int, vector<int> >::const_iterator b = it->second.begin(); b != it->second.end(); b++)
		{
			//		cout << "index = " << b->first << "[";
			for (vector<int>::const_iterator c = b->second.begin(); c != b->second.end(); ++c)
			{
				//			cout << *c << ", ";
				elemCount++;
			}
			//		cout << " ], ";
		}
		//	cout << endl;
	}
	file.close();
	cout << "elem count in index table = " << elemCount << endl;
	cout << "kmer count in index table = " << kmerCount << endl;
	cout << "**************************************************************" << endl;

}
