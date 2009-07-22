/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */
#ifndef _CellularGA_H_
#define _CellularGA_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../Elements/Chromosome_Real.h"
#include "../Elements/Population.h"
#include "../Operators/Mutations/Mutation.h"
#include "../Operators/Crossovers/Crossover.h"

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "GlobalSearch.h"

class CellularGA : public GlobalSearch<double>
{
public:
	CellularGA(void);	
	CellularGA(int populationWidth, int populationHeight, ObjectiveFunction* objectiveFunction, Mutation<double>* mutationOperator,
						Crossover<double>* crossoverOperator);
	~CellularGA() { }
	virtual void evolve(unsigned int nGeneration = 1);
public:
	Mutation<double>* mutation;
	Crossover<double>* crossover;
	int popW, popH;
};

#endif

