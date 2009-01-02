/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #14: Shifted rotated expanded Scaffer's F6 function
******************************************************************************/

#ifndef _FScaffer_H
#define _FScaffer_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Scaffer function.
 */
class FScaffer:public ObjectiveFunction
{
public:
	FScaffer( int nDimensions );
	FScaffer( int nDimensions, double lowerBound, double upperBound );
	FScaffer( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FScaffer();

// <identifier>
	virtual const string equation()
	{
		return "";
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
