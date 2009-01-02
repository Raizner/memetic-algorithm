/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "FSchwefel102Noisy.h"
#include "../Rng/GlobalRng.h"

FSchwefel102Noisy::FSchwefel102Noisy( int nDimensions ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
	seed = 0.0;
	randomize();
	initrandomnormaldeviate();
}

FSchwefel102Noisy::FSchwefel102Noisy( int nDimensions, double lowerBound, double upperBound ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
}

FSchwefel102Noisy::FSchwefel102Noisy( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{

}

FSchwefel102Noisy::~FSchwefel102Noisy()
{
}

// random generator
/* Compute the noise */

/* Get seed number for random and start it up */
void FSchwefel102Noisy::randomize()
{
      int j1;
      for(j1=0; j1<=54; j1++)
      {
            oldrand[j1] = 0.0;
      }
      jrand=0;
      warmup_random (seed);
      return;
}

/* Get randomize off and running */
void FSchwefel102Noisy::warmup_random (double seed)
{
      int j1, ii;
      double new_random, prev_random;
      oldrand[54] = seed;
      new_random = 0.000000001;
      prev_random = seed;
      for(j1=1; j1<=54; j1++)
      {
            ii = (21*j1)%54;
            oldrand[ii] = new_random;
            new_random = prev_random-new_random;
            if(new_random<0.0)
            {
                  new_random += 1.0;
            }
            prev_random = oldrand[ii];
      }
      advance_random ();
      advance_random ();
      advance_random ();
      jrand = 0;
      return;
}

/* Initialize the randome generator for normal distribution */
void FSchwefel102Noisy::initrandomnormaldeviate()
{
    rndcalcflag = 1;
    return;
}

/* Create next batch of 55 random numbers */
void FSchwefel102Noisy::advance_random ()
{
      int j1;
      double new_random;
      for(j1=0; j1<24; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1+31];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
      for(j1=24; j1<55; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1-24];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
}

/* Fetch a single random number between 0.0 and 1.0 */
double FSchwefel102Noisy::randomperc()
{
      jrand++;
      if(jrand>=55)
      {
            jrand = 1;
            advance_random();
      }
      return((double)oldrand[jrand]);
}

/* Compute the noise */
double FSchwefel102Noisy::randomnormaldeviate()
{
    double t;
    if(rndcalcflag)
    {
        rndx1 = sqrt(- 2.0*log(randomperc()));
        t = 6.2831853072*randomperc();
        rndx2 = sin(t);
        rndcalcflag = 0;
        return(rndx1*cos(t));
    }
    else
    {
        rndcalcflag = 1;
        return(rndx1*rndx2);
    }
}
// end of random generator

double FSchwefel102Noisy::evaluate_( vector<double>& x )
{
	int i, j;

	ObjectiveFunction::evaluate_( x );

	double res = 0.0;	

	for( i=0; i<nDim; i++)
	{
		double tmp = 0;
		for( j=0; j<=i; j++)
		{
			tmp += x[j];
		}
		res += tmp*tmp;
	}

	double r = randomnormaldeviate();
	return res*(1.0+0.4*fabs(r));
	//return res+bias;
}

vector<double> FSchwefel102Noisy::gradient_( vector<double>& x )
{
	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
