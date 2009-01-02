/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #6: Shifted Rosenbrock function
******************************************************************************/

#ifndef _FROSENBROCK_H
#define _FROSENBROCK_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Rosenbrock function.
 * f(x) = sum{i = 1..(n - 1)}(100 * (x(i + 1) - x(i) ^ 2) ^ 2 + (x(i) - 1) ^ 2)
 */
class FRosenbrock:public ObjectiveFunction
{
public:
	FRosenbrock( int nDimensions );
	FRosenbrock( int nDimensions, double lowerBound, double upperBound );
	FRosenbrock( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FRosenbrock();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = sum{i = 1..(n - 1)}(100 * (x(i + 1) - x(i) ^ 2) ^ 2 + (x(i) - 1) ^ 2)";
	}
	virtual const string classname()
	{
		return "Rosenbrock";
	}
// </identifier>
	
private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
};

#endif
