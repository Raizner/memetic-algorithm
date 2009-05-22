/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _ObjectiveFunction_H_
#define _ObjectiveFunction_H_

#include "../Global.h"
#include "../Utilities/Statistics.h"

/*!
 * \brief
 * The super class for objective functions. 
 */
class ObjectiveFunction {

	//variables
public:
	enum { MINIMIZE = -1, MAXIMIZE = 1 };
protected:
// function parameters
	
	int nDim;					/* number of dimensions */
	int nFunc;					/* number of basic functions */

// algorithm variables	
	double bestEval;			/* best evaluation value */
	vector<double> bestSol;    /* best solution so far */

//methods
public:
	int minimaxi;				/* nature of function: minimizing = -1 maximizing = +1 */
	bool isMaximizing();
			
	ObjectiveFunction( int minmax, unsigned int nDimensions, double lowerBound, double upperBound );
	ObjectiveFunction( int minmax, unsigned int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds );	
	~ObjectiveFunction();

	// <identifier>
	virtual const string equation()
	{
		return "Objective Function base class";
	}
	virtual const string classname()
	{
		return "";
	}
	// </identifier>

	unsigned int nDimensions();
	unsigned int nBasicFunctions();	

	double bestEvaluation();
	vector<double> bestSolution();

	bool isInBound( vector<double> x);

public:
	unsigned int nEvaluations; /* number of evaluations */	
	vector<double> translationVector;

	vector<double> lowerBounds;		/* lower bound vector */
	vector<double> upperBounds;		/* upper bound vector */
	
	double operator()( vector<double>& x );
	double evaluate(vector<double>& x);
	vector<double> gradient(vector<double>& x);

	Statistics* statModule;

protected:
	virtual double evaluate_( vector<double>& x );
	virtual vector<double> gradient_( vector<double>& x );
	vector<double> finiteDifference( vector<double>& x );
	void initialize();
};

#endif
