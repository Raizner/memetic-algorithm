/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Selection_RouletteWheel_H_
#define _Selection_RouletteWheel_H_

#include <vector>
#include "../../Elements/Chromosome.h"
#include "Selection_RouletteWheel.h"
#include "Selection.h"

using namespace std;

/*!
 * \brief
 * Implementation of the biased roulette wheel selection.
 *  
 * Chance that an individual is selected is proportioned to its scaled fitness value.
 *
 * \see
 * Selection
 */

template<typename T> class Selection_RouletteWheel : public Selection<T> {

public:
	Selection_RouletteWheel(Scaling* scaling);
	virtual vector< Chromosome<T>* > select(vector< Chromosome<T>* >& individuals, unsigned int nSelect);

private:
};


template <typename T>
Selection_RouletteWheel<T>::Selection_RouletteWheel(Scaling* scaling) : Selection<T>(scaling)
{
}

template <typename T>
vector< Chromosome<T>* > Selection_RouletteWheel<T>::select(vector< Chromosome<T>* >& individuals, unsigned int nSelect)
{
	vector< Chromosome<T>* > selectedPool(nSelect);

	unsigned int nIndv = individuals.size();
	vector <double> fitness(nIndv);
	unsigned int i;

	// extract all the fitness values
	for(i=0; i<nIndv; i++)
	{
		fitness[i] = individuals[i]->fitness;
	}
	
	// scale the fitness in order to avoid negative fitness values
	(*this->scaling)(fitness);
	
	// calculate the total sum of all fitness values, this will be
	// the roulette wheel's size
	double fitnessSum = 0.0;
	for(i=0; i<nIndv; i++)
	{
		fitnessSum+=fitness[i];
	}

	// if converged
	if (fitnessSum < 1e-10)
	{
		for(i=0; i<nSelect; i++)
		{
			selectedPool[i] = individuals[0]->clone();
		}
		return selectedPool;
	}

	// vector of random numbers, representing chosen individuals
	// we are generating nSelect random points from 0 to fitnessSum
	vector<double> rndPoints(nSelect);
	for(i=0; i<nSelect; i++)
	{
		rndPoints[i] = Rng::uni(0.0, fitnessSum);
	}

	// append to the pool
	// we have nSelect points and nIndv segments on the roulette wheel
	// we need to find the segment that each selected point lies in
	// below is an O(nSelect*log(nSelect)) algorithm
	// we first sort the values and scan the individuals array only once

	sort(rndPoints.begin(), rndPoints.end());

	unsigned int selPos = 0; // position in the rndPoints vector
	unsigned int indPos = 0; // position in the fitness/individuals vectors
	double selectedSum = 0.0;
	
	while(selPos < nSelect)
	{
		// if the considering point rndPoints[selPos] lies within
		// the current segment of individuals[indPos]
		if (selectedSum + fitness[indPos] > rndPoints[selPos])
		{
			selectedPool[selPos] = individuals[indPos]->clone();
			selPos++;
		}
		else
		{
			selectedSum += fitness[indPos];
			indPos++;
		}
	}

	return selectedPool;
}

#endif
