/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "Chromosome_Binary.h"

/*!
 * \brief
 * Constructor.
 * 
 * \param nBitPerDimension
 * Number of bit used for each optimization(search) variable.
 * 
 * \param nDimension
 * Number of optimization(search) variables.
 * 
 * \param lowerBounds
 * Vector of lower bound values.
 * 
 * \param upperBounds
 * Vector of upper bound values.
 *  
 * Implement a binary encoding chromosome.
 * 
 * \remarks
 * Note that in binary chromosome, we actually descretize the search space.
 * So be careful with precision loss.
 * 
 * \see
 * Chromosome<T> | Chromosome_Real
 */
Chromosome_Binary::Chromosome_Binary(const unsigned int nBitPerDimension, const unsigned int nDimension, vector<double> lowerBounds, vector<double> upperBounds):nBitPerDim(nBitPerDimension), nDim(nDimension), Chromosome<bool>(nBitPerDimension*nDimension, lowerBounds, upperBounds)
{
}


/****************************************************
Encoding and Decoding function
***************************************************/

/*!
 * \brief
 * Convert a vector of double to a binary chromosome.
 * 
 * \param phenotype
 * The input phenotype.
 *  
 * 
 * For binary-encoded chromosome, we discretize the range
 * [upperBound, lowerBound] to 2^n steps where n is the 
 * number of bits allocated for one dimension.
 * 
 * \remarks
 * Be careful of precision loss.
 * 
 * \see
 * Chromosome_Binary::toDoubleVector
 */
void Chromosome_Binary::fromDoubleVector(vector<double>& phenotype)
{
	// note that we won't do range checking here
	// make sure the phenotype values are within the range

	// number of steps
	int64 nRanges = (1<<nBitPerDim);
	int i, j;

	int pos = 0;

	for(i=0; i<nDim; i++)
	{
		// step size of each step
		double stepSize = (upperBounds[i] - lowerBounds[i]) / nRanges;

		// discretize the phenotype[i]
		int64 nSteps = (int64)((phenotype[i] - lowerBounds[i]) / stepSize);

		// convert nSteps to binary
		for(j=0; j<nBitPerDim; j++)
		{
			(*this)[pos] = ((nSteps & (1LL << j)) > 0);
			pos++;
		}
	}
}

/*!
 * \brief
 * Decode a binary string into a vector of double values.
 * 
 * \returns
 * The decoded vector of doubles.
 *  
 * Convert the binary string to decimal, then scale it back to the range.
 * 
 * \remarks
 * Optimized for fast calculation (using bit-shift operators).
 * 
 * \see
 * Chromosome_Binary::fromDoubleVector()
 */
vector<double> Chromosome_Binary::toDoubleVector()
{
	vector<double> res(nDim);
	
	int i, j;
	int pos = 0;
	
	for(i=0; i<nDim; i++)
	{	
		long long tmp = 0;
		for(j=0; j<nBitPerDim; j++)
		{
			if(this->at(pos)) tmp |= (1LL << j);
			pos++;
		}

		res[i] = lowerBounds[i] + ((upperBounds[i] - lowerBounds[i]) / (1<<nBitPerDim))*tmp;
	}
	return res;
}
