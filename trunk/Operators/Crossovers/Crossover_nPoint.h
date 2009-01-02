/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Crossover_nPoint_H_
#define _Crossover_nPoint_H_

#include "../../Elements/Chromosome.h"
#include "Crossover.h"

#include "../../Rng/GlobalRng.h"

/********************************************************************
Headers
********************************************************************/
/*!
 * \brief
 * Implementation of n-point crossover (applicable for both binary/real-coded).
 *  
 */
template <typename T> class Crossover_nPoint : public Crossover<T> {
public:
	Crossover_nPoint(int nPoints, double pCrossover);
	virtual bool cross(Chromosome<T>& dad, Chromosome<T>& mom);

private:
	int n;
};


/********************************************************************
Implementation
********************************************************************/
/*!
 * \brief
 * Write brief comment for Crossover_nPoint here.
 * 
 * \param nPoints
 * Number of crossing points.
 * 
 * \param pCrossover
 * Crossover probability.
 *  
 * 
 * \see
 * CrossOver | Crossover_Uniform
 */
template <typename T> 
Crossover_nPoint<T>::Crossover_nPoint(int nPoints, double pCrossover): n(nPoints), Crossover<T>(pCrossover)
{	
}

template <typename T>
bool Crossover_nPoint<T>::cross(Chromosome<T>& dad, Chromosome<T>& mom)
{
	double randomNumber = Rng::uni();

	// should we do crossover?
	if (randomNumber > this->pCrossover) return false;
	
	// generate n randoms crossover points
	vector<int> crossSites(n);
	int chroLen = dad.size();

	int i, j;

	for(i=0; i<n; i++)
	{
		crossSites[i] = (int)Rng::uni(0.0,  (double)chroLen - 1 - 1E-9);
	}
	
	sort(crossSites.begin(), crossSites.end());

	// we add-in one safe-guard element
	crossSites.push_back(chroLen);	

	// the following code is slow, block copy is faster
	for(i=0; i<n; i++)
	{
		if (i % 2 == 0)
		{
			for(j=crossSites[i]; j<crossSites[i+1]; j++)
			{
				T tmp = dad[j];
				dad[j] =  mom[j];
				mom[j] = tmp;
			}
		}
	}
	return true;
}

#endif
