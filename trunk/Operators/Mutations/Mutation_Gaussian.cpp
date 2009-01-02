/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "Mutation_Gaussian.h"

/********************************************************************
Implementation
********************************************************************/
Mutation_Gaussian::Mutation_Gaussian(double pMutation) : Mutation<double>(pMutation)
{
	sigma = 1.0;
}


bool Mutation_Gaussian::mutate(Chromosome<double>& chromosome)
{
	//double randomNumber = Rng::uni();
	//if (randomNumber > this->pMutation) return false;

	for (unsigned int i=0; i<chromosome.size(); i++)
	{
		double randomNumber = Rng::uni();

		if (randomNumber <= pMutation) 
		{
			double tmpSigma = sigma;
            double tmp;
            do {
                    tmp = Rng::gauss(0, tmpSigma);
					// we need to reduce sigma
					// to avoid the case of very small range
					// causing the loop to run infinitely
                    tmpSigma /= 2;
			} while (chromosome[i] + tmp < chromosome.lowerBounds[i] || chromosome[i] + tmp > chromosome.upperBounds[i]);

			chromosome[i] += tmp;
		}
	}

	return true;
}

