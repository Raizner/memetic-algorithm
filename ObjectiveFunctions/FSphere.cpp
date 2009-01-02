/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FSphere.h"

FSphere::FSphere( int nDimensions ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
}

FSphere::FSphere( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FSphere::FSphere( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FSphere::~FSphere()
{
}

double FSphere::evaluate_( vector<double>& x )
{
	int i;
	double res;

	ObjectiveFunction::evaluate_( x );

	res = 0;

	for( i = 0; i < nDim; i++ ) {
		res += (x[i] * x[i]);
	}
		
	return res;
}

vector<double> FSphere::gradient_( vector<double>& x )
{
	vector<double> grad(x.size());
	//<Manhtung>
	//Rewrite later//
	int i;
	for(i=0;i<nDim;i++){
		grad [ i ] = 2 * x [ i ];
	}
	//</Manhtung>

	return grad;
}
