/*
The Lennard-Jones potential (also referred to as the L-J potential, 6-12 potential or, 
less commonly, 12-6 potential) is a mathematically simple model that describes the 
interaction between a pair of neutral atoms or molecules. A form of the potential was 
first proposed in 1924 by John Lennard-Jones.

For more information, please refer to:
http://en.wikipedia.org/wiki/Lennard-Jones_potential
 
This code implement the Lennard-Jones potential function

Copyright (c) 2010 by <Quang Huy / NTU>
*/

#include "FLennardJones.h"


FLennardJones::FLennardJones(int nDimensions, double epsilon, double sigma) : ObjectiveFunction(-1, nDimensions, INT_MIN, INT_MAX)
{
	this->epsilon = epsilon;
	this->sigma = sigma;
	eps6 = pow(epsilon, 6.0);
	eps12 = pow(epsilon, 12.0);
}

double FLennardJones::evaluate_( vector<double>& x )
{
	int i, j, k;	

	double res = 0;
	int nAtom = x.size() / 3;

	RP(i, nAtom)
	{
		FR(j, i+1, nAtom - 1)
		{
			double r2 = 0;
			RP(k, 3) r2+=SQR(x[i*3+k]-x[j*3+k]);
			// cut-off at rc = 2.5 * \sigma
			if (r2 > 6.25 * sigma * sigma) continue;
			res += 4 * epsilon * (eps12/pow(r2, 6.0) - eps6/pow(r2, 3.0));
		}
	}	

	return res;
}

vector<double> FLennardJones::gradient_( vector<double>& x )
{	
	int i, j, k;
	vector<double> grad(x.size());

	// doing type casting here - I hate this
	int nAtom = x.size() / 3;
	
	RP(i, nAtom) RP(j, 3) grad[i*3+j] = 0;

	RP(i, nAtom)
	{
		FR(j, i+1, nAtom - 1)
		{
			double r2 = 0;
			RP(k, 3) r2+=SQR(x[i*3+k]-x[j*3+k]);
			// cut-off at rc = 2.5 * \sigma
			if (r2 > 6.25 * sigma * sigma) continue;
			
			RP(k, 3)
			{
				double tmp = 4*epsilon*(-eps12/pow(r2, 6.5)+eps6/pow(r2, 3.5)) * (2*(x[i*3+k]-x[j*3+k]));			
				grad[i*3+k] += tmp;
				grad[j*3+k] -= tmp;
			}
		}
	}	

	return grad;
}