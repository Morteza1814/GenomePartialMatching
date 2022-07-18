#ifndef PARTIALMATCHER_H
#define PARTIALMATCHER_H

#include "Configs.hpp"

extern map<string, map<int, vector<int> > > indexTable; //a map of kmers, positions, indexes

int partialMatcher(int numberOfQueries);


#endif