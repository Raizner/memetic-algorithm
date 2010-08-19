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

#ifndef _FLennardJones_H
#define _FLennardJones_H

#define NMAX 40

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class FLennardJones:public ObjectiveFunction
{
public:
	FLennardJones(int nDimensions, double epsilon, double sigma);	
	~FLennardJones() {}	

	// <identifier>
	virtual const string equation()
	{
		return "V(r) = 4 * epsilon * [(sigma / r)^12 - (sigma / r)^6]";
	}
	virtual const string classname()
	{
		return "Lennard Jones potential";
	}
	// </identifier>

private:
	// parameters
	double epsilon, sigma;

	// epsilon ^ 6, epsilon ^ 12
	double eps6, eps12;
	
private:
	virtual double evaluate_( vector<double>& x ); 
	virtual vector<double> gradient_( vector<double>& x );	
};

#endif
