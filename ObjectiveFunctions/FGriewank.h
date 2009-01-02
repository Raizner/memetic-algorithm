/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #7: Shifted rotated Griewank function
******************************************************************************/

#ifndef _FGRIEWANK_H
#define _FGRIEWANK_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Griewank function.
 * f(x) = sum{i = 1..n}(x(i) ^ 2) / 4000 - prod{i = 1..n}(cos(x(i) / sqrt(i))) + 1
 */
class FGriewank:public ObjectiveFunction
{
public:
	FGriewank( int nDimensions );
	FGriewank( int nDimensions, double lowerBound, double upperBound );
	FGriewank( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FGriewank();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = sum{i = 1..n}(x(i) ^ 2) / 4000 - prod{i = 1..n}(cos(x(i) / sqrt(i))) + 1";
	}
	virtual const string classname()
	{
		return "Griewank";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
};

#endif
