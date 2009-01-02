/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FElliptic.h"

FElliptic::FElliptic( int nDimensions ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
}

FElliptic::FElliptic( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FElliptic::FElliptic( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{

}

FElliptic::~FElliptic()
{
}


double FElliptic::evaluate_( vector<double>& x )
{
	int i;

	ObjectiveFunction::evaluate_( x );	

	double res = 0.0;
	for( i=0; i<nDim; i++)
	{
		res += pow(1.0e6, (double)i/(nDim-1)) * x[i] * x[i];
	}

	return res;
}

vector<double> FElliptic::gradient_( vector<double>& x )
{

	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
