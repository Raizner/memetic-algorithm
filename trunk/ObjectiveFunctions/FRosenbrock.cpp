/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FRosenbrock.h"

FRosenbrock::FRosenbrock( int nDimensions ):ObjectiveFunction( -1, nDimensions, -3, 3)
{
}

FRosenbrock::FRosenbrock( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FRosenbrock::FRosenbrock( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FRosenbrock::~FRosenbrock()
{
}

double FRosenbrock::evaluate_( vector<double>& x )
{
	int i;
	double res;
	double _x[1000];

	for( i=0; i<nDim; i++) _x[i] = x[i] + 1;	

	res = 0;
	
	for( i = 0; i < nDim - 1; i++ ) {
		res += ( 100 * pow( _x[ i] * _x[ i] - _x[ i + 1], 2 ) + pow( _x[ i] - 1, 2 ) );
	}

	return res;
}

vector<double> FRosenbrock::gradient_( vector<double>& x )
{
	
	vector<double> grad(x.size());
	vector<double> _x(x.size());
	
	//<Manhtung>
	//Rewrite later//
	int i;

	for(i=0; i<nDim; i++) _x[i] = x[i] + 1;

	grad[0] = 100.0 * 2.0 * ( _x[0] * _x[0] - _x[1] ) * 2 * _x[0] + 2 * (_x[0] - 1);
	grad[nDim-1] = 100.0 * 2.0 * ( _x [nDim - 2] * _x[nDim - 2] - _x[nDim - 1] ) * -1.0;

	for(i=1; i<=nDim-2; i++){
		
		grad[i] = 100.0 * 2.0 * ( _x[i - 1] * _x[i - 1] - _x[i]) * -1
			+ 100.0 * 2.0 * ( _x[i] * _x[i] - _x[i + 1] ) * 2 * _x[i] + 2 * (_x[i] - 1);
	}
	//</Manhtung>

	return grad;
}
