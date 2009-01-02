/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #9, #10: Shifted rotated Rastrigin function
******************************************************************************/

#ifndef _FRASTRIGIN_H
#define _FRASTRIGIN_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Rastrigin function.
 * f(x) = sum{i = 1..n}(x(i) ^ 2 - 10 * cos(2 * PI * x(i)) + 10)
 */
class FRastrigin:public ObjectiveFunction
{
public:
	FRastrigin( int nDimensions );
	FRastrigin( int nDimensions, double lowerBound, double upperBound );
	FRastrigin( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FRastrigin();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = sum{i = 1..n}(x(i) ^ 2 - 10 * cos(2 * PI * x(i)) + 10)";
	}
	virtual const string classname()
	{
		return "Rastrigin";
	}
// </identifier>
	
private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
};

#endif
