/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Recombination_KeepBest_H_
#define _Recombination_KeepBest_H_

#include "../../Global.h"
#include "../../Elements/Chromosome.h"
#include "Recombination.h"

/*!
 * \brief
 * Implementation of elitism recombination
 *  
 * The best individuals in the current population are retained and 
 * the remaining part of population is replaced by the top individuals in the offpring.
 *
 * \see
 * Recombination | Recombination_MuCommaLambda | Recombination_MuPlusLambda
 */

template<typename T> class Recombination_KeepBest: public Recombination<T> {

public:
	/*!
	 * \brief
	 * Constructor
	 * 
	 * \param nElitists
	 * Number of elitists
	 * 	 
	 */
	Recombination_KeepBest(unsigned int nElitists) : nElitists(nElitists)
	{
	}

	virtual bool recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings);
	
protected:	
	unsigned int nElitists;
};


template<typename T>
bool Recombination_KeepBest<T>::recombine(vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings)
{
	if (nElitists >= currentPopulation.size()) return false;
	if (currentPopulation.size() - nElitists > offsprings.size()) return false;

	unsigned int i, j = 0;

	for(i=0; i<currentPopulation.size(); i++)
	{
		for(j=i+1; j<currentPopulation.size(); j++)
		{
			if (currentPopulation[i]->fitness < currentPopulation[j]->fitness)
			{
				swap(currentPopulation[i], currentPopulation[j]);
			}
		}
	}

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
	

	j = 0;
	for(i=nElitists; i<currentPopulation.size(); i++)
	{
		swap(offsprings[j++], currentPopulation[i]);
	}

	for(i=0; i<offsprings.size(); i++) delete offsprings[i];

	return true;
}

#endif
