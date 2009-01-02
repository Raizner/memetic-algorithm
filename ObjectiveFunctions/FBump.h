/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _FBUMP_H
#define _FBUMP_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Bump function.
 * f(x) = frac{abs(sum(cos(x_i)^4) - 2prod(cos(x_i)^2))}{sqrt(sum{i=1}^n i*(x_i)^2)}; prod x_i > 0.75; sum x_i < 15*n/2
 */
class FBump:public ObjectiveFunction
{
public:
	FBump( int nDimensions );
	FBump( int nDimensions, double lowerBound, double upperBound );
	FBump( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FBump();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = frac{abs(sum(cos(x_i)^4) - 2prod(cos(x_i)^2))}{sqrt(sum{i=1}^n i*(x_i)^2)}; prod x_i > 0.75; sum x_i < 15*n/2";
	}
	virtual const string classname()
	{
		return "Bump";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
};

#endif
