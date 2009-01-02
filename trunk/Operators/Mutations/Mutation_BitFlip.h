/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Mutation_BitFlip_H_
#define _Mutation_BitFlip_H_

#include "../../Elements/Chromosome.h"
#include "Mutation.h"

#include "../../Rng/GlobalRng.h"

using namespace std;

/********************************************************************
Headers
********************************************************************/
/*!
 * \brief
 * Implementation of bit-flip mutation (applicable binary-coded chromosomes).
 *  
 * \see
 * Mutation | Mutation_Gaussian
 */
class Mutation_BitFlip : public Mutation<bool> {
public:
	Mutation_BitFlip(double pMutation);

	virtual bool mutate(Chromosome<bool>& chromosome);

private:	
};

#endif
