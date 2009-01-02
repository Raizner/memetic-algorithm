/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FStep.h"

FStep::FStep( int nDimensions ):ObjectiveFunction( -1, nDimensions, -5.12, 5.12)
{
}

FStep::FStep( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FStep::FStep( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{

}

FStep::~FStep()
{
}


double FStep::evaluate_( vector<double>& x )
{
	int i;
	double sum1;

	ObjectiveFunction::evaluate_( x );
	
	sum1 = 0;	

	for( i = 0; i < nDim; i++ ) {
		sum1 += floor(x[i]);
	}

	return sum1 + 6*nDim;
}

vector<double> FStep::gradient_( vector<double>& x )
{		
	return vector<double>(x.size(), 0.0);
}
