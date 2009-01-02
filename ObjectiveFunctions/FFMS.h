/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function: FMS
******************************************************************************/

#ifndef _FFMS_H
#define _FFMS_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Frequency Modulation of Sound function.
 * f(x) = -20 * exp(-0.2 * sqrt(sum{i = 1..n}(x(i) ^ 2) / n)) - exp(sum{i = 1..n}(cos(2 * PI * x(i))) / n) + 20 + exp(1)
 */
class FFMS:public ObjectiveFunction
{
public:
	FFMS();
	FFMS(double lowerBound, double upperBound );
	FFMS(vector<double>& lowerBounds, vector<double>& upperBounds );
	~FFMS();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = -20 * exp(-0.2 * sqrt(sum{i = 1..n}(x(i) ^ 2) / n)) - exp(sum{i = 1..n}(cos(2 * PI * x(i))) / n) + 20 + exp(1)";
	}
	virtual const string classname()
	{
		return "FMS";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	

	double FMS_helper(int t,double a1,double w1,double a2,double w2,double a3,double w3);

	double FMS_val[101];
};

#endif
