/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FGriewankRosenbrock.h"

FGriewankRosenbrock::FGriewankRosenbrock( int nDimensions ):ObjectiveFunction( -1, nDimensions, -5, 5)
{
}

FGriewankRosenbrock::FGriewankRosenbrock( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FGriewankRosenbrock::FGriewankRosenbrock( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FGriewankRosenbrock::~FGriewankRosenbrock()
{
}

double FGriewankRosenbrock::evaluate_( vector<double>& x )
{
	int i;	
	ObjectiveFunction::evaluate_( x );
	
	double res = 0.0;	
	double* _x = new double[nDim];

	for (i=0; i<nDim; i++) _x[i] = x[i] + 1;

    for (i=0; i<nDim-1; i++)
    {
        double temp = 100.0*pow((_x[i]*_x[i]-_x[i+1]),2.0) + pow((_x[i]-1.0),2.0);
	        res += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }

	double temp = 100.0*pow((_x[nDim-1]*_x[nDim-1]-_x[0]),2.0) + pow((_x[i]-1.0),2.0);
    res += (temp*temp)/4000.0 - cos(temp) + 1.0;

    return res;	
}

vector<double> FGriewankRosenbrock::gradient_( vector<double>& x )
{	
	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
