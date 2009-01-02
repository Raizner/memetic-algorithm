/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "Scaling_Linear.h"

/*
Given non-negative fitness f(x), the scaled
function: g=af+b
Use heuristics to obtain values for a and b
– maintaining average fitness value before and after scaling
g_{avg} = f_{avg}.
– g_{max}=C * f_{avg} to restrict the fittest candidate to have atmost 2
samples in the mating pool by setting C=2.
- a = f_{avg}(C-1) / (f_{max}-f_{avg}) b = f_{avg}(1-a)
*/

bool Scaling_Linear::operator() (vector<double> & v)
{
	double sum = 0;
	double vmax = v[0];
	double vmin = v[0];
	double C = 2;

	unsigned int i;

	for(i=0; i<v.size(); i++)
	{
		vmax = max(vmax, v[i]);
		vmin = min(vmin, v[i]);
		sum += v[i];
	}

	if (vmax-vmin < 1e-12) return false;

	// now we are to make sure that all the value in v is > 0
	if (vmax <= 0)
	{
		for(i=0; i<v.size(); i++) v[i] = v[i] - vmax - vmin;
		double tmp = vmax;
		vmax = -vmin;
		vmin = -tmp;
		sum = sum + v.size() *(vmax + vmin);
	}	

	
	double vavg = sum / v.size();
	
	double a = vavg*(C-1) / (vmax - vavg);
	double b = vavg*(1-a);
	
	if (a*vmin + b < 0) 
	{
		a = vavg/(vavg-vmin);
		b = vavg*(1-a);
	}

	for(i=0; i<v.size(); i++)
	{
		v[i] = v[i]*a + b;
	}

	return true;
}