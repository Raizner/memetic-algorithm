/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #13: Shifted expanded Griewank plus Rosenbrock function
******************************************************************************/

#ifndef _FGriewankRosenbrock_H
#define _FGriewankRosenbrock_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the hybrid Griewank-Rosenbrock function.
 */
class FGriewankRosenbrock:public ObjectiveFunction
{
public:
	FGriewankRosenbrock( int nDimensions );
	FGriewankRosenbrock( int nDimensions, double lowerBound, double upperBound );
	FGriewankRosenbrock( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FGriewankRosenbrock();

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
