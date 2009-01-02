/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FScaffer.h"

FScaffer::FScaffer( int nDimensions ):ObjectiveFunction( -1, nDimensions, -5, 5)
{
}

FScaffer::FScaffer( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FScaffer::FScaffer( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
}

FScaffer::~FScaffer()
{
}

double FScaffer::evaluate_( vector<double>& x )
{
	int i;

	ObjectiveFunction::evaluate_( x );

	double res = 0.0;

    for (i=0; i<nDim-1; i++)
    {
        double temp1 = pow((sin(sqrt(pow(x[i],2.0)+pow(x[i+1],2.0)))),2.0);
        double temp2 = 1.0 + 0.001*(pow(x[i],2.0)+pow(x[i+1],2.0));
        res += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    double temp1 = pow((sin(sqrt(pow(x[nDim-1],2.0)+pow(x[0],2.0)))),2.0);
    double temp2 = 1.0 + 0.001*(pow(x[nDim-1],2.0)+pow(x[0],2.0));

    res += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    return res;	
}

vector<double> FScaffer::gradient_( vector<double>& x )
{
	
	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
