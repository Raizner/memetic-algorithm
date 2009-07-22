/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */
#ifndef _CellularMA_H_
#define _CellularMA_H_

#include "../Global.h"
#include "../Rng/GlobalRng.h"

#include "../Elements/Population.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "../LocalSearches/LocalSearch.h"
#include "../EvolutionaryAlgorithms/CellularGA.h"

#include "MemeticAlgorithm.h"

/*!
 * \brief
 * Implementation of Cellular Memetic Algorithms.  
 *  
 * \see
 * MemeticAlgorithm | CellularGA
 */

class CellularMA : public MemeticAlgorithm<double>{
public:	
	CellularMA();
	CellularMA(CellularGA* globalSearch, LocalSearch* localSearch = NULL);
	void evolve(unsigned int nGeneration = 1);
protected:
	vector<int> getLSIndividuals();
	void initialize();
};
#endif

