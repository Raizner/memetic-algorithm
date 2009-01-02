/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #11: shifted rotated Weierstrass
******************************************************************************/

#ifndef _FWEIERSTRASS_H
#define _FWEIERSTRASS_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Weierstrass function.
 * f(x) = sum{i = 1..n}(sum{j = 0..k}(a ^ j * cos(2 * PI * b ^ j * (x(i) + 0.5)))) - n * sum{j = 0..k}(a ^ j * cos(PI * b ^ j)), a = 0.5, b = 3, k = 20
 */
class FWeierstrass:public ObjectiveFunction
{
public:
	FWeierstrass( int nDimensions );
	FWeierstrass( int nDimensions, double lowerBound, double upperBound );
	FWeierstrass( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FWeierstrass();

// <identifier>
	virtual const string equation()
	{
		return "f(x) = sum{i = 1..n}(sum{j = 0..k}(a ^ j * cos(2 * PI * b ^ j * (x(i) + 0.5)))) - n * sum{j = 0..k}(a ^ j * cos(PI * b ^ j)), a = 0.5, b = 3, k = 20";
	}
	virtual const string classname()
	{
		return "Weierstrass";
	}
// </identifier>

private:
	void init();
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
private:
	double pa[100];
	double pb[100];
	double sum2;
	double a;
	double b;
	double pi;
	double k;
};

#endif
