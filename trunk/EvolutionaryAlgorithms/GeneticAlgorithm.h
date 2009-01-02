/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _GeneticAlgorithm_H_
#define _GeneticAlgorithm_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../Elements/Population.h"
#include "../Operators/Mutations/Mutation.h"
#include "../Operators/Crossovers/Crossover.h"
#include "../Operators/Selections/Selection.h"
#include "../Operators/Recombinations/Recombination.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "GlobalSearch.h"

/*!
 * \brief
 * This class implement standard Genetic Algorithm
 * 
 * \see
 * GlobalSearch
 */
template<typename T>
class GeneticAlgorithm: public GlobalSearch<T>
{
public:
	GeneticAlgorithm() { }
	GeneticAlgorithm(Population<T> population, ObjectiveFunction* objectiveFunction, Mutation<T>* mutationOperator,
						Crossover<T>* crossoverOperator, Selection<T>* selectionScheme, Recombination<T>* recombinationScheme);
	~GeneticAlgorithm() { }
	virtual void evolve(unsigned int nGenerations = 1);

public:
	Mutation<T>* mutation;
	Crossover<T>* crossover;
	Selection<T>* selection;
	Recombination<T>* recombination;
};


/*!
 * \brief
 * Constructor
 * 
 * \param population
 * Input population of type Population
 * 
 * \param objectiveFunction
 * Objective function of type ObjectiveFunction
 * 
 * \param mutationOperator
 * Mutation operator
 * 
 * \param crossoverOperator
 * Crossover operator
 * 
 * \param selectionScheme
 * Selection operator
 * 
 * \param recombinationScheme
 * Recombination operator
 *  
 * \remarks
 * The mutation, crossover, selection and recombination must be implemented for the same coding type with population
 * 
 * \see
 * Population | ObjectiveFunction | Crossover | Mutation | Selection | Recombination
 */
template<typename T>
GeneticAlgorithm<T>::GeneticAlgorithm(Population<T> population, ObjectiveFunction* objectiveFunction, Mutation<T>* mutationOperator,
									  Crossover<T>* crossoverOperator, Selection<T>* selectionScheme, Recombination<T>* recombinationScheme):
mutation(mutationOperator), crossover(crossoverOperator), selection(selectionScheme), recombination(recombinationScheme)
{	
	this->pop = population;
	this->fObj = objectiveFunction;
	unsigned int i;
	// evaluation the population
	for(i=0; i<this->pop.size(); i++)
	{
		//Chromosome<T>*& id = *(pop[j]);
		vector<double> tmp = this->pop[i]->toDoubleVector();

		this->pop[i]->fitness = this->evaluate(tmp);
	}
}

/*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGenerations
 * Number of generation to evolve.
 */
template<typename T>
void GeneticAlgorithm<T>::evolve(unsigned int nGenerations)
{
	unsigned int i, j;

	for(i=0; i<nGenerations; i++)
	{
		this->nGen++;

		// selection
		vector< Chromosome<T>* >  matingPool = selection->select(this->pop, this->pop.size());

		/*for(j=0; j<matingPool.size(); j++)
			cout << matingPool[j]->fitness << " ";

		cout << endl;*/
		
		// scramble the mating pool
		vector<int>	order(matingPool.size());

		for(j=0; j<order.size(); j++) order[j] = j;
		for(j=0; j<order.size(); j++)
		{
			int k = (int)Rng::uni(0, (double)matingPool.size() - 1 - 1e-9);
			swap(order[j], order[k]);
		}
		
		// crossover
		for(j=0; j<matingPool.size()-1; j+=2)
		{				
			crossover->cross(*matingPool[order[j]], *matingPool[order[j+1]]);		
		}				

		// mutation
		for(j=0; j<matingPool.size(); j++)
		{
			mutation->mutate(*matingPool[j]);
		}

		// evaluation the offspring
		for(j=0; j<matingPool.size(); j++)
		{			
			vector<double> tmp = matingPool[j]->toDoubleVector();

			matingPool[j]->fitness = this->evaluate(tmp);
		}

		recombination->recombine(this->pop, matingPool);
	}
}

#endif
