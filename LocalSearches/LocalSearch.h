/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _LocalSearch_H_
#define _LocalSearch_H_

#include "../Global.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"

/*!
 * \brief
 * The super class for local search methods.
 *  
 * 
 * \remarks
 * Note that the nature of global search is to minimize
 * due to the search strategy generally based on gradients / hessian matrix / quadratic approximation.
 * Therefore, if the objective function is to maximize,
 * the evaluation should return the negative value of it.
 * 
 * \see
 * ObjectiveFunction
 */
class LocalSearch
{
public:
	// stopping criteria
	double accuracy;     /* required accuracy */
	unsigned int evaluationLimit; /* computational budget allocated */
	unsigned int timeLimit; /* time budget allocated */
	unsigned int iterationLimit; /* number of iterations allocated */

	// objective function
	ObjectiveFunction* fObj;
	
	// initial step length
	vector<double> stepLength;
	vector<double> lowerBounds;
	vector<double> upperBounds;

protected:	
	unsigned int nEvalStart; /* function evaluation count at the beginning of the search */
	time_t nTimeStart; /* clock value at the beginning of the search */

	bool done(); /* check weather termination condition is satisfied */
	double evaluate(vector<double>& x);	
public:
	LocalSearch(ObjectiveFunction* objectiveFunction);	
	~LocalSearch();

	double operator() (vector<double>& x) { return search(x); }
	virtual double search(vector<double>& x);
	
	virtual const char* className() { return "LS base class"; }
private:
	void initialize();
	
};

#endif
