/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Crossover_Uniform_H_
#define _Crossover_Uniform_H_

#include "../../Elements/Chromosome.h"
#include "Crossover.h"

#include "../../Rng/GlobalRng.h"

/********************************************************************
Headers
********************************************************************/
/*!
 * \brief
 * Implementation of uniform crossover (applicable for both binary/real-coded)
 *  
 */
template <typename T> class Crossover_Uniform : public Crossover<T> {
public:
	Crossover_Uniform(double pCrossover);
	virtual bool cross(Chromosome<T>& dad, Chromosome<T>& mom);

private:	
};


/********************************************************************
Implementation
********************************************************************/
template <typename T> 
Crossover_Uniform<T>::Crossover_Uniform(double pCrossover): Crossover<T>(pCrossover)
{	
}


template <typename T>
bool Crossover_Uniform<T>::cross(Chromosome<T>& dad, Chromosome<T>& mom)
{
	double randomNumber = Rng::uni();

	// should we do crossover?
	if (randomNumber > this->pCrossover) return false;
		
	// the following code is slow, block copy is faster
	for(unsigned i=0; i<dad.size(); i++)
	{
		randomNumber = Rng::uni();
		if (randomNumber >= 0.5)
		{
			T tmp = dad[i];
			dad[i] =  mom[i];
			mom[i] = tmp;
		}
	}
	return true;
}

#endif
