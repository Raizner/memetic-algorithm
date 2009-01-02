/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Crossover_H_
#define _Crossover_H_

#include "../../Global.h"
#include "../../Elements/Chromosome.h"


/*!
 * \brief
 * The super class for crossover operators.
 *  
 * \see
 * Crossover_nPoint | Crossover_Uniform
 */

template<typename T> class Crossover {
public:
	
	/*!
	 * \brief
	 * Constructor.
	 * 
	 * \param pCrossover
	 * Crossover probability.	 
	 */
	Crossover(double pCrossover);
	
	/*!
	 * \brief
	 * Wrapper operator for the cross function.
	 */
	virtual bool operator()(Chromosome<T>& dad, Chromosome<T>& mom) { return cross(dad, mom); }
	
	/*!
	 * \brief
	 * Implementation of the crossover strategy.
	 * 	
	 */
	virtual bool cross(Chromosome<T>& dad, Chromosome<T>& mom) { return false; }

protected:
	double pCrossover;
};

/********************************************************************
Implementation
********************************************************************/
template<typename T> Crossover<T>::Crossover(double pCrossover) : pCrossover(pCrossover)
{	
}

#endif
