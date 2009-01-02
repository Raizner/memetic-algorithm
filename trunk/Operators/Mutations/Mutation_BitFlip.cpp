/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "Mutation_BitFlip.h"

/********************************************************************
Implementation
********************************************************************/
Mutation_BitFlip::Mutation_BitFlip(double pMutation) : Mutation<bool>(pMutation)
{
}


bool Mutation_BitFlip::mutate(Chromosome<bool>& chromosome)
{
	//double randomNumber = Rng::uni();
	//if (randomNumber > this->pMutation) return false;

	for (unsigned int i=0; i<chromosome.size(); i++)
	{
		double randomNumber = Rng::uni();

		if (randomNumber <= pMutation) chromosome[i] = ~chromosome[i];
	}

	return true;
}
