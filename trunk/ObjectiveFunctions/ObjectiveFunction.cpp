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
 * Best fitness value so far.
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

/*!
 * \brief
 * Evaluate the fitness value
 * 
 * \returns
 * Fitness value
 *  
 * \remarks
 * Virtual function, to be implemented by subclass
 * 
 * \see
 * ObjectiveFunction::evaluate()
 */
double ObjectiveFunction::evaluate_(vector<double>& x)
{
	return 0.0;
}

/*!
 * \brief
 * Calculate the gradient vector
 * 
 * \returns
 * Gradient vector
 *  
 * \remarks
 * Virtual function, to be implemented by subclass, finite differencing is used by default
 * 
 * \see
 * ObjectiveFunction::gradient()
 */
vector<double> ObjectiveFunction::gradient_(vector<double>& x)
{
	return finiteDifference(x);
}

/*!
 * \brief
 * Check if the given solution is feasible, used for constrained objective function
 * 
 * \returns
 * 0 if the solution is <b>feasible</b>, otherwise returns the number of violated constraints 
 *  
 * \remarks
 * Virtual function, to be implemented by subclass, perform bound checking by default;
 * 
 * \see
 * ObjectiveFunction::gradient()
 */
int ObjectiveFunction::isInfeasible(vector<double>& x)
{
	if (isInBound(x)) return 0; else return 1;
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


vector<double> ObjectiveFunction::finiteDifference(vector<double> &x)
{
	double eps = 1.0;	

	// first we find the smallest eps we can handle
	double one = 1.0;
	double two = 2.0;
	double test = 0.0;
	do
	{
		eps = eps/two;
		test = eps * one + one;
	} while (test > 1.0);
	eps = eps*two*two;	

	// now let's go for the finite different	
	vector<double> result(nDim, 0);
	double fx = evaluate(x);

	double delta, xi, udelta = 0, rsteps = sqrt(eps);

	for(int i=0; i<nDim; i++) 
	{
		xi = x[i];
		double tmp = 1.0;
		if (fabs(xi) >= 1.0) tmp = fabs(xi);

		if (udelta > rsteps * tmp)
		{
			delta = udelta;
		}
		else
		{
			delta = rsteps * tmp;
		}

		if (xi < 0) delta = -delta;

		tmp = x[i];
		x[i] += delta;
		double ft = evaluate(x);
		result[i] = (ft - fx) / delta;
		x[i] = tmp;
	}

	return result;
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
	// feasibility checking
	if (isInfeasible(x)) return 1.0 * INT_MIN * minimaxi;

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
