#ifndef _Individual_H_
#define _Individual_H_

#include <vector>
#include "Chromosome.h"

using namespace std;

template<typename T> class Individual : public std::vector < Chromosome<T> * > {

public:
	// designed to support co-evolution, self-adaptive which requires more than one type of chromosome
	double fitness;

	Individual(const size_t size) : std::vector< Chromosome<T>* >(size, NULL)
	{
	}

	Individual<T>* clone();
	// reserved for constraint problems
	// vector<double> equalityConstraints;
	// vector<double> inequalityConstraints;	
};

#endif
