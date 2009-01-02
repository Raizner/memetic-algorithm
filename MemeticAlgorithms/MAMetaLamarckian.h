/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _MAMetaLamarckian_H_
#define _MAMetaLamarckian_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../Elements/Population.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "../LocalSearches/LocalSearch.h"
#include "../EvolutionaryAlgorithms/GlobalSearch.h"
#include "MemeticAlgorithm.h"

/*!
 * \brief
 * Implemetation of Meta-Lamarckian MA.
 * 
 * \see
 * MemeticAlgorithm<T>
 */
template<typename T>
class MAMetaLamarckian: public MemeticAlgorithm<T> {
public:
	MAMetaLamarckian();
	MAMetaLamarckian(GlobalSearch<T>* globalSearch);
	void evolve(unsigned int nGenerations = 1);
	
	/*!
	 * \brief
	 * Local search pool.
	 * 	 
	 * \remarks
	 * To replace the local search parameter in the superclass, i.e., MemeticAlgorithm::ls.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	vector<LocalSearch*> lsPool;
protected:	
	double rewards[1000];
	unsigned counter;
	virtual void initialize();
};

/*!
 * \brief
 * Default constructor.
 * 
 */
template<typename T>
MAMetaLamarckian<T>::MAMetaLamarckian()
{	
	initialize();
}

/*!
 * \brief
 * Constructor.
 * 
 * \param globalSearch
 * Global search method of type GlobalSearch.
 *  
 * 
 * \remarks
 * No local search as we need to maintain a pool of local searchers.
 * 
 * \see
 * GlobalSearch
 */
template<typename T>
MAMetaLamarckian<T>::MAMetaLamarckian(GlobalSearch<T>* globalSearch)
{
	MemeticAlgorithm<T>(globalSearch, NULL);
	this->gs = globalSearch;
	initialize();
}

template<typename T>
void MAMetaLamarckian<T>::initialize()
{
	MemeticAlgorithm<T>::initialize();
	counter = 0;
	for (int i=0; i<1000; i++)
	{
		rewards[i] = 0;
	}
}

/*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGenerations
 * Number of generation to evolve.
 *
 * We have 2 phases: local searcher selection to find the most appropriate local serach to apply
 * and individual selection which is similar to one in canonical MA.
 */
template<typename T>
void MAMetaLamarckian<T>::evolve(unsigned int nGenerations)
{

	for(unsigned int i=0; i<nGenerations; i++)
	{
		if (this->gs != NULL) this->gs->evolve();
		

		// Local search phase
		if (lsPool.size() == 0) continue;
		vector<int> v = this->getLSIndividuals();

		// for every individual
		for(unsigned i=0;i<v.size();i++)
		{
			// if not all the LS has been run at least 5 times
			// then they will take turn to run
			int curLS = 0;
			if (counter < 5*lsPool.size())
			{
				// Give a single opportunity to each method			
				curLS = counter % lsPool.size();
			}
			else
			{			
				// if all the methods have been chosen at least 4 times,
				// use biased roulette wheel to choose based on the previous performance
					
				double sum = 0;
				for (unsigned j=0; j<lsPool.size(); j++) 
				{
					sum += rewards[j];
				}			


				double tmp;

				//tmp  = ((double)(rand() % 11)) / 10 * sum;
				tmp  = sum * Rng::uni(0, 1);

				double sum2=0;
				for (unsigned j=0; j<lsPool.size(); j++)
				{
					sum2 += rewards[j];
					if (sum2 > tmp) 
					{
						curLS = j;
						break;
					}
				}
			}

			counter++;
			// run the LS
			this->ls = lsPool[curLS];
			this->ls->evaluationLimit = this->tLS;
			//cout << "Running " << ls->className() << endl;

			Chromosome<T>* ch = this->gs->pop[v[i]];	
			double oldFitness = ch->fitness;
			unsigned oldEval = this->gs->fObj->nEvaluations;
			this->ls->search(*ch);
			double newFitness = ch->fitness = this->gs->evaluate(*ch);
			unsigned newEval = this->gs->fObj->nEvaluations;
			
			// calculate reward
			int ret = 0;		
			double nheta = 0, beta;
			
			beta = newFitness / this->gs->bestEvaluation();			
			nheta = beta * (newFitness - oldFitness);
			nheta /= (newEval - oldEval);

			rewards[curLS] += nheta;
		}
	}
}

#endif
