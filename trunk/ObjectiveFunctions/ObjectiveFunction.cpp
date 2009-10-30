/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#include "ObjectiveFunction.h"


/*********************************************************/
/* Constructors and destructor                           */
/*********************************************************/
//constructor
/*!
 * \brief
 * Constructor
 * 
 * \param minmax
 * Minimize or maximize(0 = minimize, 1 = maximize). 
 * 
 * \param nDimensions
 * Number of optimization variables.
 * 
 * \param lowerBound
 * Lower bound value, applied for all variables.
 * 
 * \param upperBound
 * Upper bound value, applied for all variables.
 *  
 */
ObjectiveFunction::ObjectiveFunction( int minmax, unsigned int nDimensions, double lowerBound, double upperBound )
{
	minimaxi = minmax;
	nDim = nDimensions;

	lowerBounds.assign(nDimensions, lowerBound);
	upperBounds.assign(nDimensions, upperBound);
	
	initialize();
}

//constructor
/*!
 * \brief
 * Write brief comment for ObjectiveFunction here.
 * 
 * \param minmax
 * Minimize or maximize(0 = minimize, 1 = maximize). 
 * 
 * \param nDimensions
 * Number of optimization variables.
 * 
 * \param lowerBounds
 * Vector of lower bound values.
 * 
 * \param upperBounds
 * Vector of upper bound values.
 *  
 */
ObjectiveFunction::ObjectiveFunction( int minmax, unsigned int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds )
{
	minimaxi = minmax;
	nDim = nDimensions;

	lowerBounds.assign(lowerBounds.begin(), lowerBounds.end());
	upperBounds.assign(upperBounds.begin(), upperBounds.end());
		
	initialize();
}

//destructor
/*!
 * \brief
 * Destructor.
 */
ObjectiveFunction::~ObjectiveFunction()
{	
}

void ObjectiveFunction::initialize()
{
	statModule = NULL;
	nEvaluations = 0;
	bestEval = 0;
	bestSol.resize(nDim);

	translationVector.assign(nDim, 0.0);
	noiseLevel = 0;
}
/*********************************************************/
/* Evaluating					                         */
/*********************************************************/
/*!
 * \brief
 * Bound checking.
 * 
 * \returns
 * True if all the variables is inbound, false otherwise.
 * 	 
 * 
 * Check if the chromsome is within the given bounds.
 * 
 * \remarks
 * Used in many local search / constrained optimization problem
 * 	 
 */
bool ObjectiveFunction::isInBound( vector<double> x )
{
	for( int i=0; i<nDim; i++)
	{
		if (x[i] < lowerBounds[i] || x[i] > upperBounds[i]) return false;
	}

	return true;
}

/*!
 * \brief
 * Number of optimization variables.
 * 
 */
unsigned int ObjectiveFunction::nDimensions()
{
	return nDim;
}

unsigned int ObjectiveFunction::nBasicFunctions()
{
	return nFunc;
}

/*!
 * \brief
 * Best result so far.
 * 
 * \returns
 * Best objective value obtained.
 *  
 * \remarks
 * Dependent on the minmax value.
 * 
 * \see
 * ObjectiveFunction::bestSolution()
 */
double ObjectiveFunction::bestEvaluation()
{
	return bestEval;
}

/*!
 * \brief
 * Best solution so far.
 * 
 * \returns
 * The soultion with the best fitness value.
 *  
 * \remarks
 * Dependent on the minmax value.
 * 
 * \see
 * ObjectiveFunction::bestEvaluation()
 */
vector<double> ObjectiveFunction::bestSolution()
{
	return bestSol;
}

double ObjectiveFunction::evaluate_(vector<double>& x)
{
	return 0.0;
}

vector<double> ObjectiveFunction::gradient_(vector<double>& x)
{
	return vector<double>(x.size(), 0.0);
}

/*!
 * \brief
 * Implementation of the finite differencing method  to calculate the 1st order derivative.
 * 
 * \param x
 * Input values.
 * 
 * \returns
 * Gradient vector.
 * 
 * \remarks
 * Machine dependent. Should only be used when analytical gradient is not available.
 * 
 * \see
 * ObjectiveFunction::gradient
 */


vector<double> ObjectiveFunction::finiteDifference(vector<double>& x)
{
	return vector<double>(x.size(), 0.0);
}

/*!
 * \brief
 * Wrapper of the evaluate function.
 */
 
double ObjectiveFunction::operator() ( vector<double>& x )
{
	return evaluate(x);
}

/*!
 * \brief
 * Evaluate the fitness value.
 * 
 * \param x
 * Input values.
 * 
 * \returns
 * Fitness value.
 * 
 * \see
 * ObjectiveFunction::gradient()
 */
double ObjectiveFunction::evaluate( vector<double>& x )
{
	// increase the number of evaluations
	nEvaluations++;
	// bound checking
	if (!isInBound(x)) return INT_MIN * minimaxi;

	// everything's ok, let's evaluate the fitness of x
	// translate it first
	vector<double> xtmp(x);
	for(unsigned i=0; i<x.size(); i++)
	{
		xtmp[i] -= translationVector[i];
	}

	double res = evaluate_(xtmp);
	// enable noise
	res *= (1.0 + noiseLevel * fabs(Rng::gauss()));

	// update the best
	if (nEvaluations == 1) 
	{
		bestEval = res;
		bestSol = x;
	}
	else 
	{
		if (minimaxi == MINIMIZE) 
		{
			if (bestEval > res)
			{
				bestEval = res;
				bestSol = x;
			}
		}
		else 
		{
			if (bestEval < res)
			{
				bestEval = res;
				bestSol = x;
			}
		}
	}

	// recording
	if (statModule != NULL) statModule->addEntry(nEvaluations, x, res);
	
	return res;
}

bool ObjectiveFunction::isMaximizing()
{
	return (this->minimaxi == MAXIMIZE);
}

/*!
 * \brief
 * Calculate gradient (first order derivative).
 * 
 * \param x
 * Input values.
 * 
 * \returns
 * Gradient vector.  
 * 
 * \see
 * ObjectiveFunction::evaluate()
 */
vector<double> ObjectiveFunction::gradient(vector<double>& x)
{
	// increase the number of evaluations by the number of dimensions
	nEvaluations += nDim;

	// translate it first
	vector<double> xtmp(x);
	for(int i=0; i<nDim; i++)
	{
			xtmp[i] -= translationVector[i];
	}
	
	vector<double> res = this->gradient_(xtmp);

	// enable noise
	for(int i=0; i<nDim; i++) { res[i] *= (1.0 +  noiseLevel * fabs(Rng::gauss()));  }
	
	return res;
}
