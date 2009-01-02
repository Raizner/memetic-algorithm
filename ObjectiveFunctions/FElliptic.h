/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #3: Shifted rotated high conditioned elliptic function
******************************************************************************/

#ifndef _FELLIPTIC_H
#define _FELLIPTIC_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Elliptic function
 */
class FElliptic:public ObjectiveFunction
{
public:
	FElliptic( int nDimensions );
	FElliptic( int nDimensions, double lowerBound, double upperBound );
	FElliptic( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FElliptic();

// <identifier>
	virtual const string equation()
	{
		return "";
	}
	virtual const string classname()
	{
		return "Elliptic";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	
};

#endif
