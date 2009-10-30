/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FBump.h"

FBump::FBump( int nDimensions ):ObjectiveFunction(1, nDimensions, 0, 10)
{
}

FBump::FBump( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction(1, nDimensions, lowerBound, upperBound)
{
}

FBump::FBump( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction(1, nDimensions, lowerBounds, upperBounds)
{
}

double FBump::evaluate_( vector<double>& x )
{
	int i;

	ObjectiveFunction::evaluate_( x );

	// check constraints
	double sumx = 0, prodx = 1;
	for (i=0; i<nDim; i++)
	{
		sumx += x[i];
		prodx *= x[i];
	}

	if ((sumx >= 7.5 * nDim) || (prodx <= 0.75)) return 0; // penalty


	double sumc4=0, prodc2=1, sumixi2=0;
	double tmp;

	for( i = 0; i < nDim; i++ ) {
		tmp = cos(x[i]);
		sumc4 += tmp*tmp*tmp*tmp;;
		prodc2 *= tmp*tmp;
		sumixi2 += (i+1)*x[i]*x[i];
	}
	
	return fabs(sumc4 - 2*prodc2) / sqrt(sumixi2);
}

vector<double> FBump::gradient_( vector<double>& x )
{	
	return ObjectiveFunction::finiteDifference(x);
}
