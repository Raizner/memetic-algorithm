/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2010 by <Quang Huy / NTU>
 */

#ifndef _FOss2_H
#define _FOss2_H

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class FOrthodontist:public ObjectiveFunction
{
public:
	FOrthodontist(int nDimensions);	
	~FOrthodontist() {}	

	// <identifier>
	virtual const string equation()
	{
		return "<See the paper>";
	}
	virtual const string classname()
	{
		return "Oss2";
	}
	// </identifier>

	virtual int isInfeasible( vector<double>& x);
private:
	virtual double evaluate_( vector<double>& x ); 	
};

#endif
