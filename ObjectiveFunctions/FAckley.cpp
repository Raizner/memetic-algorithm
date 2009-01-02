/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FAckley.h"

FAckley::FAckley( int nDimensions ):ObjectiveFunction( -1, nDimensions, -32, 32)
{
}

FAckley::FAckley( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FAckley::FAckley( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{

}

FAckley::~FAckley()
{
}


double FAckley::evaluate_( vector<double>& x )
{
	int i;
	double sum1, sum2;

	ObjectiveFunction::evaluate_( x );
	
	sum1 = 0.0;
	sum2 = 0.0;

	for(i=0; i<nDim; i++) {
		sum1 += ( x[ i ] * x[ i ] );
		sum2 += cos( 2 * PI * x[ i ] );
	}

	return -20 * exp( -0.2 * sqrt( sum1/nDim ) ) - exp( sum2/nDim ) + 20 + exp(1.0);
}

vector<double> FAckley::gradient_( vector<double>& x )
{
	vector<double> grad(x.size(), 0.0);
	//<Manhtung>
	//Rewrite later//
	int i;
	double sum1 = 0, sum2 = 0;

	for(i=0; i<nDim; i++){
		sum1 += ( x[ i ] * x[ i ] );
		sum2 += cos( 2 * PI * x[ i ] );
	}

	for(i=0;i<nDim;i++){		
		if(sum1!=0)
        {
                grad[i] = -20 * exp( -0.2 * sqrt( sum1/nDim ) ) * ( -0.2 / sqrt( 1.0*nDim ) ) * ( 1/(2 * sqrt( sum1 )) ) * 2 * x[i]
                - exp( sum2 / nDim ) * ( 2 * acos( -1.0 ) / nDim ) * ( -sin ( 2 * acos( -1.0 ) * x[i] ) );
        }
        else
        {
                grad[i] = - exp( sum2 / nDim ) * ( 2 * acos( -1.0 ) / nDim ) * ( -sin ( 2 * acos( -1.0 ) * x[i] ) );
        }
	}
	//</Manhtung>

	return grad;
}
