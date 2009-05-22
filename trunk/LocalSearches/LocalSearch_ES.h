/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */
#ifndef _LocalSearch_ES_H
#define _LocalSearch_ES_H

#include "LocalSearch.h"

/*!
 * \brief
 * Implement the (1+\lambda) Evolutionary Strategy with 1/5-rule adaptation.
 *  
 * \remarks
 * A stochastic local search method
 * 
 * \see
 * LocalSearch
 */
class LocalSearch_ES:public LocalSearch
{
public:
	LocalSearch_ES( ObjectiveFunction *f );	
	~LocalSearch_ES();

	virtual double search(vector<double>& x);
	virtual const char* className() { return "LS_ES"; }

	int lambda; /* number of offspring */
	double sigma; /* mutation strength */
	double multiplier; /* used in 1/5 rule */
private:	
};

#endif
