/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Recombination_H_
#define _Recombination_H_

#include "../../Global.h"
#include "../../Elements/Chromosome.h"

/*!
 * \brief
 * The super class for the recombination operators.
 *  
 * \see
 * Recombination_KeepBest | Recombination_MuCommaLambda | Recombination_MuPlusLambda
 */
template<typename T> class Recombination {

public:	
	
	/*!
	 * \brief
	 * Wrapper operator for the recombine function.
	 *
	 * Offsprings will be combined into the current population.
	 */
	virtual bool operator() (vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings)
	{
		return recombine(currentPopulation, offsprings);
	}

	/*!
	 * \brief
	 * Implementation of the recombination strategy.
	 * 	
	 */
	virtual bool recombine (vector< Chromosome<T>* >& currentPopulation, vector< Chromosome<T>* >& offsprings)
	{
		return false;
	}
	
protected:	
};

#endif
