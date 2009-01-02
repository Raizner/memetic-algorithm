/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #1
******************************************************************************/

#ifndef _FSPHERE_H
#define _FSPHERE_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Sphere function.
 * f(x) = sum{i = 1..n}(x(i) ^ 2)
 */
class FSphere:public ObjectiveFunction
{
public:
	FSphere( int nDimensions );
	FSphere( int nDimensions, double lowerBound, double upperBound );
	FSphere( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FSphere();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = sum{i = 1..n}(x(i) ^ 2)";
	}
	virtual const string classname()
	{
		return "Sphere";
	}
// </identifier>
	
private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
};

#endif
