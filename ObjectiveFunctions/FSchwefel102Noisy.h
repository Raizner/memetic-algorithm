/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */


/******************************************************************************
Function #4: shifted schwefel's problem 1.2 with noise
******************************************************************************/

#ifndef _FSchwefel102Noisy_H
#define _FSchwefel102Noisy_H

#include "ObjectiveFunction.h"

/*!
 * \brief
 * Implementation of the Schwefel function 1.02 with Gaussian noise.
 */
class FSchwefel102Noisy:public ObjectiveFunction
{
public:
	FSchwefel102Noisy( int nDimensions );
	FSchwefel102Noisy( int nDimensions, double lowerBound, double upperBound );
	FSchwefel102Noisy( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );
	~FSchwefel102Noisy();

// <identifier>
	virtual const string equation()
	{
		return "";
	}
	virtual const string classname()
	{
		return "Schwefel with noise";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );	

private:
	/* Variable declarations for the random number generator */
	double seed;
	double oldrand[55];
	int jrand;
	int rndcalcflag;
	double rndx1, rndx2;

	/* Function declarations for the random number generator */
	void randomize();
	void warmup_random (double seed);
	void advance_random ();
	double randomperc();
	//int rnd (int low, int high);
	//double rndreal (double low, double high);
	void initrandomnormaldeviate();
	//double noise (double mu, double sigma);
	double randomnormaldeviate();

};

#endif
