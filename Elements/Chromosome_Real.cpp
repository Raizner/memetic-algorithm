/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "Chromosome_Real.h"


/*!
 * \brief
 * Implement a real encoding chromosome
 * 
 * \param nDimension
 * Number of optimization(search) variables
 * 
 * \param lowerBounds
 * Vector of lower bound values
 * 
 * \param upperBounds
 * Vector of upper bound values
 *    
 * \remarks
 * Widely used
 * 
 * \see
 * Chromosome<T> | Chromosome_Binary
 */
Chromosome_Real::Chromosome_Real(const unsigned int nDimension, vector<double> lowerBounds, vector<double> upperBounds):Chromosome<double>(nDimension, lowerBounds, upperBounds)
{
}

/*!
 * \brief
 * Not used. The chromosome is a vector of doubles itself 
 * 
 */
void Chromosome_Real::fromDoubleVector(vector<double>& phenotype)
{
	this->assign(phenotype.begin(), phenotype.end());
}

/*!
 * \brief
 * Not used. The chromosome is a vector of doubles itself
 *  
 *  
 */
vector<double> Chromosome_Real::toDoubleVector()
{	
	return *this;
}
