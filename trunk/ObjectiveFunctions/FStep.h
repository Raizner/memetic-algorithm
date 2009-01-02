/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #8: Shifted rotated Ackley function
******************************************************************************/

#ifndef _FStep_H
#define _FStep_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Step function.
 * f(x) = 6n + \\sum_{i=1}^n \\floor{x_i}
 */
class FStep:public ObjectiveFunction
{
public:
	FStep( int nDimensions );
	FStep( int nDimensions, double lowerBound, double upperBound );
	FStep( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FStep();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = 6n + \\sum_{i=1}^n \\floor{x_i}";
	}
	virtual const string classname()
	{
		return "Step";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	
};

#endif
