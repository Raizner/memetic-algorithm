/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _MemeticAlgorithm_H_
#define _MemeticAlgorithm_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../Elements/Population.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "../LocalSearches/LocalSearch.h"
#include "../EvolutionaryAlgorithms/GlobalSearch.h"

/*!
 * \brief
 * The super class for Memetic Algorithms.
 * 
 * \param T
 * Type T can be bool (binary problems), int (combinatorial problems) or real (continuous problems).
 *  
 * \see
 * GeneticAlgorithm | LocalSearch
 */

template<typename T>
class MemeticAlgorithm {
public:	
	MemeticAlgorithm();
	MemeticAlgorithm(GlobalSearch<T>* globalSearch, LocalSearch* localSearch = NULL);
	void evolve(unsigned int nGenerations = 1);
	
	/*!
	 * \brief
	 * Different strategy to select individuals undergoing local search.
	 *
	 * Best: individuals with best fitnesses are selected.
	 * Stratified: individuals are selected so as to maintain highest diversity.
	 * Random: individuals are selected randomly.
	 */
	typedef enum {maLSBest, maLSStratified, maLSRandom} MA_SELECTION_STRATEGY;
	typedef enum {maLSLamarckian, maLSBaldwinian} MA_LEARNING_STRATEGY;
	int nGenerations();
	bool done();
public:
	
	/*!
	 * \brief
	 * Global search method	 
	 *
	 * \see
	 * GlobalSearch
	 */
	GlobalSearch<T>* gs;

	/*!
	 * \brief
	 * Local search method	 
	 *
	 * \see
	 * LocalSearch
	 */
	LocalSearch* ls;

	/*!
	 * \brief
	 * Maximum number of evaluations allocated for each individual learning.
	 */
	int tLS;

	/*!
	 * \brief
	 * Percentage of the population should undergo individual learning.
	 */
	double pLS;
	
	/*!
	 * \brief
	 * Strategy to select individuals undergoing individual learning.
	 * 
	 * \see
	 * MA_SELECTION_STRATEGY
	 */
	MA_SELECTION_STRATEGY maSelectionStrategy;

	/*!
	 * \brief
	 * Strategy for individual learning.
	 * 
	 * \see
	 * MA_LEARNING_STRATEGY
	 */
	MA_LEARNING_STRATEGY maLearningStrategy;

	/*!
	 * \brief
	 * (Stopping criteria) Maximum number of evaluations allowed.	 
	 */
	unsigned int maxEvaluations;

	/*!
	 * \brief
	 * (Stopping criteria) Maximum number of generations allowed.	 
	 */
	unsigned int maxGenerations;
protected:
	vector<int> getLSIndividuals();	
	virtual void initialize();
};


/*!
 * \brief
 * Default constructor
 */ 
template<typename T>
MemeticAlgorithm<T>::MemeticAlgorithm()
{
	gs = NULL;
	ls = NULL;
	initialize();
}

/*!
 * \brief
 * Constructor.
 * 
 * \param globalSearch
 * Global search method of type GlobalSearch.
 * 
 * \param localSearch
 * Local search method of type LocalSearch.
 *    
 * \remarks
 * Global search must be implemented for the same coding type.
 * 
 * \see
 * GlobalSearch | LocalSearch
 */
template<typename T>
MemeticAlgorithm<T>::MemeticAlgorithm(GlobalSearch<T>* globalSearch, LocalSearch* localSearch):gs(globalSearch),ls(localSearch)
{
	initialize();
}

/*!
 * \brief
 * Initialization.
 *  
 */
template<typename T>
void MemeticAlgorithm<T>::initialize()
{
	pLS = 1;	
	maSelectionStrategy = maLSBest;
	maLearningStrategy = maLSLamarckian;
	if (gs != NULL)
	{
		tLS = 10*gs->fObj->nDimensions();
		maxEvaluations = 10000*gs->fObj->nDimensions();
		maxGenerations = 10000*gs->fObj->nDimensions() / gs->pop.size();
	}
	if(ls != NULL) ls->evaluationLimit = tLS;
}

/*!
 * \brief
 * Number of generations elapsed.
 */
template<typename T>
int MemeticAlgorithm<T>::nGenerations()
{
	return gs->nGenerations();
}

/*!
 * \brief
 * Check if stopping criteria are satisfied.
 * 
 * \returns
 * True if one of the stopping conditions is true, false otherwise.
 *  
 */
template<typename T>
bool MemeticAlgorithm<T>::done()
{
	return (gs->nGenerations() >= maxGenerations || gs->fObj->nEvaluations >= maxEvaluations);
}

// MA - best
/*!
 * \brief
 * Generate the set of individuals undergoing local search, based on the setting of
 * maStrategy.
 * 
 * \returns
 * Vector of positions of the individuals.
 *  
 *  
 * \see
 * MA_STRATEGY
 */
template<typename T>
vector<int> MemeticAlgorithm<T>::getLSIndividuals()
{
	vector<int> LSIndividualSet;
	vector<int> order(gs->pop.size());
	switch (maSelectionStrategy)
	{
	case maLSBest:
		for(unsigned i=0; i < (gs->pop.size()*pLS); i++) LSIndividualSet.push_back(i);
		break;
	case maLSRandom:		
		for(unsigned i=0; i<gs->pop.size(); i++) order[i] = i;
		for(unsigned i=0; i<gs->pop.size() * pLS; i++)
		{
			int pos = (int)Rng::uni(i+1, gs->pop.size() - 1e-5);
			swap(order[i], order[pos]);
		}

		for(unsigned i=0; i < (gs->pop.size()*pLS); i++)
		{			
			LSIndividualSet.push_back(order[i]);
		}
		break;
	case maLSStratified:
		int step = (int)(1.0/pLS);
		for(unsigned i=0; i < gs->pop.size(); i+=step)
		{			
			LSIndividualSet.push_back(i);
		}
		break;
	};
	return LSIndividualSet;
}

 /*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGenerations
 * Number of generation to evolve.
 *
 * Typically, one of more generations of global search is followed by an individual learning phase.
 *  
 */
template<typename T>
void MemeticAlgorithm<T>::evolve(unsigned int nGenerations)
{
	if (ls != NULL) ls->evaluationLimit = tLS;

	for(unsigned int i=0; i<nGenerations; i++)
	{
		/*for(unsigned j=0; j<gs->pop.size(); j++)
		{
			cout << gs->pop[j]->fitness << " ";
		}
		cout << endl;*/

		if (gs != NULL) gs->evolve();

		/*for(unsigned j=0; j<gs->pop.size(); j++)
		{
			cout << gs->pop[j]->fitness << " ";
		}
		cout << endl;*/


		if (ls != NULL)
		{
			
			vector<int> v = getLSIndividuals();
			for(unsigned j=0; j<v.size(); j++)
			{
				//cout << -gs->pop[v[j]]->fitness << " --> ";
				vector<double> ch = gs->pop[v[j]]->toDoubleVector();
				//cout << ls->search(ch) << " == ";
				ls->search(ch);

				// if lamarckian learning then copy back
				if (maLearningStrategy == maLSLamarckian) gs->pop[v[j]]->fromDoubleVector(ch);	
				
				gs->pop[v[j]]->fitness = gs->evaluate(ch);

				//cout << -gs->pop[v[j]]->fitness << endl;
			}
		}
	}
}

#endif
