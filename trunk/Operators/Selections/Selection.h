/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Selection_H_
#define _Selection_H_

#include <vector>
#include "../../Elements/Chromosome.h"
#include "../Scalings/Scaling.h"

using namespace std;

/*!
 * \brief
 * The super class for the selection operators.
 *  
 * \see
 * Selection_RouletteWheel | Selection_UniversalStochastic 
 */

template<typename T> class Selection {

public:	

	/*!
	 * \brief
	 * Constructor.
	 * 
	 * \param scalingOperator
	 * The scaling operator used to scale the population's fitness values into a certain range.
	 *
	 * \see
	 * Scaling
	 */
	Selection(Scaling* scalingOperator):scaling(scalingOperator)
	{
	}

	/*!
	 * \brief
	 * Wrapper operator for the select function.
	 */
	virtual vector< Chromosome<T>* > operator () (vector< Chromosome<T>* >& individuals, unsigned int nSelect)
	{
		return select(individuals, nSelect);
	}

	/*!
	 * \brief
	 * Implementation of the selection strategy.
	 * 	
	 * \param individuals
	 * Set of individuals to be selected from
	 *
	 * \param nSelect
	 * Number of individuals to be selected
	 *
	 * \see
	 * Chromosome
	 */

	virtual vector< Chromosome<T>* > select(vector< Chromosome<T>* >& individuals, unsigned int nSelect)
	{
		return vector< Chromosome<T>* >(nSelect);
	}
	
protected:
	Scaling* scaling;
};

#endif
