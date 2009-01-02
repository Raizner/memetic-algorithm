/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _LocalSearch_DFP_H
#define _LocalSearch_DFP_H

#include "LocalSearch.h"

/*!
 * \brief
 * Implement the Davidon-Fletcher-Powell Method
 *  
 * 
 * \remarks
 * A quasi-newton method
 * 
 * \see
 * LocalSearch
 */
class LocalSearch_DFP:public LocalSearch
{
public:
	LocalSearch_DFP( ObjectiveFunction *f );	
	~LocalSearch_DFP();

	virtual double search(vector<double>& x);
	virtual const char* className() { return "LS_DFP"; }
	
private:
	void lnsrch(vector<double> &xold, const double fold, vector<double> &g, vector<double> &p, 
		vector<double> &x, double &f, const double stpmax, bool &check);
	double dfpmin(vector<double>& p);
};

#endif
