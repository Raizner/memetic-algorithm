/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #8: Shifted rotated Ackley function
******************************************************************************/

#ifndef _FACKLEY_H
#define _FACKLEY_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Ackley function.
 * f(x) = -20 * exp(-0.2 * sqrt(sum{i = 1..n}(x(i) ^ 2) / n)) - exp(sum{i = 1..n}(cos(2 * PI * x(i))) / n) + 20 + exp(1)
 */
class FAckley:public ObjectiveFunction
{
public:
	FAckley( int nDimensions );
	FAckley( int nDimensions, double lowerBound, double upperBound );
	FAckley( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FAckley();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = -20 * exp(-0.2 * sqrt(sum{i = 1..n}(x(i) ^ 2) / n)) - exp(sum{i = 1..n}(cos(2 * PI * x(i))) / n) + 20 + exp(1)";
	}
	virtual const string classname()
	{
		return "Ackley";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	
};

#endif
