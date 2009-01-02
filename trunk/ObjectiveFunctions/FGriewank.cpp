/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FGriewank.h"

FGriewank::FGriewank( int nDimensions ):ObjectiveFunction( -1, nDimensions, -600, 600)
{
}

FGriewank::FGriewank( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FGriewank::FGriewank( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FGriewank::~FGriewank()
{
}

double FGriewank::evaluate_( vector<double>& x )
{
	int i;
	double sum, prod;

	ObjectiveFunction::evaluate_( x );

	sum = 0;
	prod = 1;

	for( i = 0; i < nDim; i++ ) {
		sum += ( x[ i ] * x[ i ] );
		prod *= cos( x[ i ] / sqrt( 1.0 + i ) );
	}	

	return sum / 4000 - prod + 1;
}

vector<double> FGriewank::gradient_( vector<double>& x )
{
	vector <double> grad(x.size());

	//<Manhtung>
	//Rewrite later//
	double prod;
	int i,j;

	for(i=0;i<nDim;i++){
		grad[i] = x[i]/2000;
		prod = 1;
		for(j=0;j<i;j++){
			prod *= cos ( x[ j ] / sqrt ( 1.0+j ) );
		}
		for(j=i+1;j<nDim;j++){
			prod *= cos ( x[j] / sqrt( 1.0+j ) );
		}
		prod *= ( 1 / sqrt(1.0+i) ) * -sin ( x[ i ] / sqrt( 1.0+i ) );
		grad[i] = grad[i]-prod;
	}
	//</Manhtung>

	return grad;
}
