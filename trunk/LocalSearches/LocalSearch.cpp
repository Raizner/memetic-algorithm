/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "LocalSearch.h"

/*!
 * \brief
 * Constructor.
 * 
 * \param objectiveFunction
 * Objective function to be optimized.
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
LocalSearch::LocalSearch( ObjectiveFunction* objectiveFunction )
{
	this->fObj = objectiveFunction;	
	upperBounds = objectiveFunction->upperBounds;
	lowerBounds = objectiveFunction->lowerBounds;	
	initialize();
}

/*!
 * \brief
 * Destructor.
 *  
 */
LocalSearch::~LocalSearch()
{	
}

/*!
 * \brief
 * Initializaion of step lengths, accuracy, etc.
 */

void LocalSearch::initialize()
{
	unsigned int i;
	unsigned int nDim = fObj->nDimensions();
	
	// default values for the stopping criteria
	accuracy = (double)1e-5;	
	evaluationLimit = INT_MAX;
	timeLimit = INT_MAX;
	iterationLimit = INT_MAX;

	stepLength.resize(nDim);

	for(i=0; i<nDim; i++) {
		stepLength[i] = (upperBounds[i] - lowerBounds[i]) / 100;
	}

	// function evaluation count at the beginning of the search
	nEvalStart = fObj->nEvaluations;
	nTimeStart = time(NULL);
}

/*!
 * \brief
 * Check if stopping criteria are satisfied.
 * 
 * \returns
 * True if one of the stopping conditions is true, false otherwise.
 *  
 */
bool LocalSearch::done()
{
	if (fObj->nEvaluations - nEvalStart >= this->evaluationLimit) return true;
	if (time(NULL) - this->nTimeStart >= this->timeLimit ) return true;

	return false;
}

/*!
 * \brief
 * Evaluate the objective function.
 * 
 * \param x
 * Input phenotype.
 * 
 * \returns
 * Fitness value.
 * 
 * 
 * A wrapper for the objective function evaluation
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
double LocalSearch::evaluate(vector<double>& x)
{
	double res = (*fObj)(x);

	//for(unsigned i=0; i<x.size(); i++) cout << x[i] << " ";
	//cout << res << endl;

	if (fObj->isMaximizing()) return -res;

	return res;
}

/*!
 * \brief
 * Actual search strategy
 * 
 * \param x
 * Initial point.
 * 
 * \returns
 * Final finess value.
 *  
 *  
 * 
 * \remarks
 * To be inherited by different search strategies.
 * 
 * \see
 * LocalSearch_DFP | LocalSearch_DSCG
 */
double LocalSearch::search(vector<double>& x)
{ 
	nEvalStart = fObj->nEvaluations;
	nTimeStart = time(NULL);
	return 0; 
}
