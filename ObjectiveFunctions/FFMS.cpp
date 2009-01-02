/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FFMS.h"

FFMS::FFMS():ObjectiveFunction( -1, 6, -6.4, 6.35)
{
	for(int t=0;t<=100;t++)
	{
		FMS_val[t] = FMS_helper(t,1,5,-1.5,4.8,2,4.9);
	}
}

FFMS::FFMS(double lowerBound, double upperBound ):ObjectiveFunction( -1, 6, lowerBound, upperBound)
{
	for(int t=0;t<=100;t++)
	{
		FMS_val[t] = FMS_helper(t,1,5,-1.5,4.8,2,4.9);
	}
}

FFMS::FFMS(vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, 6, lowerBounds, upperBounds)
{
	for(int t=0;t<=100;t++)
	{
		FMS_val[t] = FMS_helper(t,1,5,-1.5,4.8,2,4.9);
	}
}

FFMS::~FFMS()
{
}

double FFMS::FMS_helper(int t,double a1,double w1,double a2,double w2,double a3,double w3)
{	
	double theta = 2*PI/100;
	double res = a1*sin(w1*t*theta + a2*sin(w2*t*theta+a3*sin(w3*t*theta)));
	return res;
}

double FFMS::evaluate_( vector<double>& x )
{
	ObjectiveFunction::evaluate_( x );
	
	double fitness = 0;

	double theta = 2*PI/100;
	double y, y0;
	for(int t=0;t<=100;t++)
	{		
		y = x[0]*sin(x[1]*t*theta+x[2]*sin(x[3]*t*theta+x[4]*sin(x[5]*t*theta)));
		y0 = FMS_val[t];		
		fitness += (y-y0)*(y-y0);
	}	
	
	return fitness;
}

vector<double> FFMS::gradient_( vector<double>& x )
{		
	return finiteDifference(x);
}
