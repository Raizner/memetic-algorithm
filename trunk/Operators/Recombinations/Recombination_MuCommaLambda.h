/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Recombination_MuCommaLambda_H_
#define _Recombination_MuCommaLambda_H_

#include "../../Global.h"
#include "../../Elements/Chromosome.h"
#include "Recombination.h"


/*!
 * \brief
 * Implementation of mu comma lambda recombination
 *  
 * The current population size is mu, each round lambda offspring will be generated
 * the best mu from these lambda offsprings are selected for the next generations.
 * the entire current population is discarded.
 *
 * \see
 * Recombination | Recombination_KeepBest | Recombination_MuPlusLambda
 */
template<typename T> class Recombination_MuCommaLambda: public Recombination<T> {

public:
	Recombination_MuCommaLambda()
	{
	}

	virtual bool recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings);
};


template<typename T>
bool Recombination_MuCommaLambda<T>::recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings)
{
	if (offsprings.size() < currentPopulation.size()) return false;
	unsigned int i, j = 0;	

	for(i=0; i<offsprings.size(); i++)
	{
		for(j=i+1; j<offsprings.size(); j++)
		{
			if (offsprings[i]->fitness < offsprings[j]->fitness)
			{
				swap(offsprings[i], offsprings[j]);
			}
		}
	}
		
	for(i=0; i<currentPopulation.size(); i++)
	{
		currentPopulation[i] = offsprings[i];
	}	

	return true;
}

#endif
