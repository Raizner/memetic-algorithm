/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FSchwefel102.h"

FSchwefel102::FSchwefel102( int nDimensions ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
}

FSchwefel102::FSchwefel102( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FSchwefel102::FSchwefel102( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{

}

FSchwefel102::~FSchwefel102()
{
}


double FSchwefel102::evaluate_( vector<double>& x )
{
	int i, j;


	ObjectiveFunction::evaluate_( x );

		
	double res = 0.0;

	for( i=0; i<nDim; i++)
	{
		double tmp = 0;
		for( j=0; j<=i; j++)
		{
			tmp += x[j];
		}
		res += tmp*tmp;
	}

	return res;
}

vector<double> FSchwefel102::gradient_( vector<double>& x )
{	
	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
