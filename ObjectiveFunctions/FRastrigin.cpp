/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FRastrigin.h"

FRastrigin::FRastrigin( int nDimensions ):ObjectiveFunction( -1, nDimensions, -5.12f, 5.12f)
{
}

FRastrigin::FRastrigin( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FRastrigin::FRastrigin( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FRastrigin::~FRastrigin()
{
}

double FRastrigin::evaluate_( vector<double>& x )
{
	int i;
	double res;

	ObjectiveFunction::evaluate_( x );

	res = 0;

	for( i = 0; i < nDim; i++ ) {
		res += ( x[ i ] * x[ i ] - 10 * cos ( 2 * PI * x[ i ] ) + 10 );
	}


	return res;
}

vector<double> FRastrigin::gradient_( vector<double>& x )
{
	vector<double> grad(x.size());
	//<Manhtung>
	//Rewrite later//
	int i;	

	for(i=0;i<nDim;i++){
		grad [ i ] = 2 * x[ i ] - 10 * -sin ( 2 * PI * x[ i ] ) * 2 * PI;
	}
	//</Manhtung>

	return grad;
}
