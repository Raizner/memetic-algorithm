/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Mutation_H_
#define _Mutation_H_


#include "../../Elements/Chromosome.h"

/********************************************************************
Headers
********************************************************************/
/*!
 * \brief
 * The super class for mutation operators.
 *  
 * \see
 * Mutation_BitFlip | Mutation_Gaussian
 */
template<typename T> 
class Mutation {
public:
	/*!
	 * \brief
	 * Default constructor	 
	 */
	Mutation() { pMutation = 0; }

	Mutation(double pMutation);

		
	/*!
	 * \brief
	 * Wrapper operator for the mutate function.
	 */
	virtual bool operator ()(Chromosome<T>& chromosome) { return mutate(chromosome); }

	/*!
	 * \brief
	 * Implementation of the mutation strategy.
	 * 	
	 */
	virtual bool mutate(Chromosome<T>& chromosome) { return false; }

protected:
	double pMutation;
};

/********************************************************************
Implementation
********************************************************************/
/*!
 * \brief
 * Constructor.
 * 
 * \param pMutation
 * Mutation probability
 */
template<typename T>
Mutation<T>::Mutation(double pMutation) : pMutation(pMutation) 
{
} 

#endif
