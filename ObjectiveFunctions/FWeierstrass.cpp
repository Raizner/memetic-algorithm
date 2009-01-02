/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FWeierstrass.h"

FWeierstrass::FWeierstrass( int nDimensions ):ObjectiveFunction( -1, nDimensions, -0.5, 0.5)
{
	init();
}

FWeierstrass::FWeierstrass( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
	init();
}

FWeierstrass::FWeierstrass( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
	init();
}

FWeierstrass::~FWeierstrass()
{
}

void FWeierstrass::init()
{	
	a = 0.5;
	b = 3;
	k = 20;
	sum2 = 0;
	pi = PI;
	
	for( int j = 0; j <= k; j++ ) {
		pa[j] =  pow( a, j );
		pb[j] = pow( b, j );
		sum2 += pa[j] * cos(pi*pb[j]);
    }

	sum2 *= nDim;
}

double FWeierstrass::evaluate_( vector<double>& x )
{
	int i, j;
	double sum1;

	ObjectiveFunction::evaluate_( x );	

	sum1 = 0;

	for( i = 0; i < nDim; i++ ) {
		for( j = 0; j <= k; j++ ) {
			sum1 += (  pa[j] * cos( 2 * pi * pb[j] * (x[ i ] + 0.5) ) );
		}
	}

	return sum1 - sum2;
}

vector<double> FWeierstrass::gradient_( vector<double>& x )
{
	vector<double> grad(x.size());
	//<Manhtung>
	//Rewrite later//
	int i, j;
	double a, b , k;
	a = 0.5;
	b = 3;
	k = 20;

	for ( i = 0 ; i < nDim ; i++ ){
		
		double sum = 0;

		for ( j = 0 ; j<= k; j++ ){

			sum += pow ( a , j ) * -sin ( 2 * PI * pow ( b, j) * ( x [ i ] + 0.5 ) ) 
				* 2 * pi * pow ( b , j );
		}

		grad [ i ] = sum;
	}
	//</Manhtung>

	return grad;
}
