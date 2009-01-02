/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #2: Schwefel's problem 1.2
******************************************************************************/

#ifndef _FSchwefel102_H
#define _FSchwefel102_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Schwefel function 1.02.
 */
class FSchwefel102:public ObjectiveFunction
{
public:
	FSchwefel102( int nDimensions );
	FSchwefel102( int nDimensions, double lowerBound, double upperBound );
	FSchwefel102( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FSchwefel102();

// <identifier>
	virtual const string equation()
	{
		return "";
	}
	virtual const string classname()
	{
		return "Schwefel";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	
};

#endif
