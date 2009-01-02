/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _LocalSearch_DSCG_H
#define _LocalSearch_DSCG_H

#include "LocalSearch.h"

/*!
 * \brief
 * Implement the strategy of Davis, Swan and Campey, with Gram-Schmidt orthogonalization.
 *  
 * \remarks
 * A direct search method
 * 
 * \see
 * LocalSearch
 */
class LocalSearch_DSCG:public LocalSearch
{
public:
	LocalSearch_DSCG( ObjectiveFunction *f );	
	~LocalSearch_DSCG();

	virtual double search(vector<double>& x);
	virtual const char* className() { return "LS_DSCG"; }

private:
	void dscg(vector<double>& x, double& fx);
	void linesearch(vector<double>& x, double& fx, vector<double>& dir, double step, vector<double>& xend, double& fxend);
};

#endif
