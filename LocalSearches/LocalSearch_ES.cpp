/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */
#include "LocalSearch_ES.h"
#include "../Rng/GlobalRng.h"

LocalSearch_ES::LocalSearch_ES( ObjectiveFunction *f) : LocalSearch( f )
{
	lambda = f->nDimensions();
	sigma = 1.0;
}

LocalSearch_ES::~LocalSearch_ES()
{
}

/*!
 * \brief
 * Search
 * 
 * \param x
 * Initial point.
 * 
 * \returns
 * Final finess value.
 *  
 */
double LocalSearch_ES::search(vector<double>& x) 
{
	unsigned evalCount = 0;

	vector<double> xtmp, xbest = x;
	double fx = evaluate(x), fbest = fx, ftmp;
	int nGood = 0, nBad = 0;
	while(evalCount < evaluationLimit)
	{
		// each iteration we generate lambda offspring
		// every n offspring (n = number of dimension), we check the adaptation rule
		for(unsigned i=0; i<fObj->nDimensions(); i++)
		{
			xtmp[i] = xbest[i] + Rng::gauss() * sigma;
		}		
		ftmp = evaluate(xtmp);
		evalCount++;

		// checking if the new point is better than the current point
		if (ftmp <= fx)
		{
			nGood++;			
		}
		else
		{
			nBad++;
		}

		// store the best point
		if (ftmp < fbest) 
		{ 
			xbest = xtmp;
			fbest = ftmp;
		}

		// 1/5 adaptive rule
		if (nGood + nBad == fObj->nDimensions())
		{
			if (nGood >= 0.2 * fObj->nDimensions()) sigma*=2; else sigma/=2;
			nGood = nBad = 0;
		}

		if (evalCount % lambda == 0 || evalCount == evaluationLimit)
		{
			if (fbest < fx)
			{
				x = xbest;
				fx = fbest;
			}
		}

	}
	
	return fx;
}
