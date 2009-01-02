/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _MALamarcBaldwin_H_
#define _MALamarcBaldwin_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "../LocalSearches/LocalSearch.h"
#include "../EvolutionaryAlgorithms/GeneticAlgorithm.h"

/*!
 * \brief
 * The super class for Memetic Algorithms.
 *  
 * \see
 * GeneticAlgorithm | LocalSearch
 */

class MALamarcBaldwin {
public:	
	MALamarcBaldwin(ObjectiveFunction* f);	
	void evolve(unsigned int nGenerations = 1);
	
	/*!
	 * \brief
	 * Different strategy to select individuals undergoing local search.
	 *
	 * Best: individuals with best fitnesses are selected.
	 * Stratified: individuals are selected so as to maintain highest diversity.
	 * Random: individuals are selected randomly.
	 */
	typedef enum {maLSBest, maLSStratified, maLSRandom} MA_STRATEGY;
	int nGenerations();
	bool done();
public:
	
	ObjectiveFunction* f;

	/*!
	 * \brief
	 * Global search methods
	 *
	 * \see
	 * GlobalSearch
	 */
	GeneticAlgorithm<double>* gs1;
	GeneticAlgorithm<bool>* gs2;

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
	 * MA_STRATEGY
	 */
	MA_STRATEGY maStrategy;

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
#endif
