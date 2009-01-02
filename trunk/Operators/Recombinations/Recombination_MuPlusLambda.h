/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Recombination_MuPlusLambda_H_
#define _Recombination_MuPlusLambda_H_

#include "../../Global.h"
#include "../../Elements/Chromosome.h"
#include "Recombination.h"

/*!
 * \brief
 * Implementation of mu plus lambda recombination
 *  
 * The current population size is mu, each round lambda offspring will be generated
 * the best \mu from (\mu+\lambda) individuals are selected for the next generations.
 *
 * \see
 * Recombination | Recombination_KeepBest | Recombination_MuCommaLambda
 */
template<typename T> class Recombination_MuPlusLambda: public Recombination<T> {
public:
	Recombination_MuPlusLambda()
	{
	}

	virtual bool recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings);
};


template<typename T>
bool Recombination_MuPlusLambda<T>::recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings)
{
	unsigned int i, j = 0;
	vector< Chromosome<T>* > tmp;

	
	for(i=0; i<currentPopulation.size(); i++)
	{
		tmp.push_back(currentPopulation[i]);
	}

	for(i=0; i<offsprings.size(); i++)
	{
		tmp.push_back(offsprings[i]);
	}

	for(i=0; i<tmp.size(); i++)
	{
		for(j=i+1; j<tmp.size(); j++)
		{
			if (tmp[i]->fitness < tmp[j]->fitness)
			{
				swap(tmp[i], tmp[j]);
			}
		}
	}
		
	for(i=0; i<currentPopulation.size(); i++)
	{
		currentPopulation[i] = tmp[i];
	}	

	return true;
}

#endif
