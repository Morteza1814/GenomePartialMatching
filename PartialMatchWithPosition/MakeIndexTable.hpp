#ifndef MAKEINDEXTABLE_H
#define MAKEINDEXTABLE_H

#include "Configs.hpp"

extern map<string, map<int, vector <int> > > indexTable; //a map of kmers, indexes, and positions

int makeIndexTable(int numberOfReads);

#endif