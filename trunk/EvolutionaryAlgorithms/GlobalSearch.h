/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _GlobalSearch_H_
#define _GlobalSearch_H_

#include "../Global.h"
#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include "../Elements/Population.h"

/*!
 * \brief
 * The super class for population-based search methods
 * 
 * \param T
 * Type T can be bool (binary problems), int (combinatorial problems) or real (continuous problems).
 * 
 * This class proposes the general interface for a population-based search algorithms,
 *  
 * \remark
 * Note that the nature of global search is to maximize
 * due to the natural selection process (one with larger fitness value has 
 * better chance to survive). Therefore, if the objective function is to minimize,
 * the evaluation should return the negative value of it.
 *
 * \see
 * GeneticAlgorithm | DifferentialEvolution | EvolutionaryStrategy
 */
template<typename T> 
class GlobalSearch
{
public:
	GlobalSearch();
	~GlobalSearch() { }
	virtual void evolve(unsigned int nGeneration = 1) { nGen += nGeneration; }
	unsigned nGenerations();
	double evaluate(vector<double>& x);
	double bestEvaluation();
public:
	/*!
	 * \brief
	 * Population to be evolved	 	 
	 * 
	 * \remarks
	 * Write remarks for pop here.
	 * 	 
	 */
	Population<T> pop;
	
	
	/*!
	 * \brief
	 * Objective function to optimize	 	 
	 * 	 	 
	 */
	ObjectiveFunction* fObj;
protected:
	
	/*!
	 * \brief
	 * Number of generation evolved so far in the search process
	 */
	unsigned nGen;
};

/*!
 * \brief
 * Default constructor
 */
template<typename T> 
GlobalSearch<T>::GlobalSearch()
{
	pop = NULL;
	fObj = NULL;
	nGen = 0;
}

/*!
 * \brief
 * Number of generation evolved so far in the search process. 
 *
 */
template<typename T> 
unsigned GlobalSearch<T>::nGenerations()
{
	return nGen;
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
 * Wrapper operator for the objective function evaluation
 * 
 * \remarks
 * Note that the nature of global search is to maximize
 * due to the natural selection process (one with larger fitness value has 
 * better chance to survive). Therefore, if the objective function is to minimize,
 * the evaluation should return the negative value of it.
 * 
 * \see
 * ObjectiveFunction
 */
template<typename T> 
double GlobalSearch<T>::evaluate(vector<double>& x)
{
	if (fObj->isMaximizing())
	{
		return (*fObj)(x);
	}

	return -(*fObj)(x);
}

/*!
 * \brief
 * Best(highest) result so far.
 * 
 * \returns
 * Best(highest) objective value obtained.
 * 
 * \throws <exception class>
 * Description of criteria for throwing this exception.
 * 
 * Return the best evaluation obtained so far during the search process
 * 
 * \remarks
 * If the objective function is to minimize, this will return the negative value of
 * the fitness.
 * 
 * \see
 * ObjectiveFunction
 */
template<typename T> 
double GlobalSearch<T>::bestEvaluation()
{
	if (fObj->isMaximizing())
	{
		return fObj->bestEvaluation();
	}

	return -fObj->bestEvaluation();
}

#endif

