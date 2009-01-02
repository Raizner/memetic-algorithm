/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Mutation_Gaussian_H_
#define _Mutation_Gaussian_H_

#include "../../Elements/Chromosome.h"
#include "Mutation.h"

#include "../../Rng/GlobalRng.h"

using namespace std;

/********************************************************************
Headers
********************************************************************/
/*!
 * \brief
 * Implementation of Gaussian mutation (applicable for real-coded chromosomes).
 *  
 * \see
 * Mutation | Mutation_BitFlip
 */
class Mutation_Gaussian : public Mutation<double> {

public:
	Mutation_Gaussian(double pMutation);

	virtual bool mutate(Chromosome<double>& chromosome);
	
	/*!
	 * \brief
	 * The sigma value for the Gaussian distribution.
	 */
	double sigma;
private:
	
};

#endif
