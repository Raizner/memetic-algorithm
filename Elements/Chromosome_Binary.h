/*!
 * <Binary chromosome implemetation>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Chromosome_Binary_H_
#define _Chromosome_Binary_H_


#include "Chromosome.h"
#include "../Global.h"

/*!
 * \brief
 * Implement a binary-coded chromosome 
 * 
 * \see
 * Chromosome | Chromosome_Real
 */
class Chromosome_Binary: public Chromosome<bool>{

public:
	Chromosome_Binary(const unsigned int nBitPerDimension, const unsigned int nDimension, vector<double> lowerBounds, vector<double> upperBounds);
	
		
	virtual void fromDoubleVector(vector<double>& phenotype);
	virtual vector<double> toDoubleVector();

	/*!
	 * \brief
	 * Clone function
	 * 
	 * \returns
	 * A "twin" of the current chromosome	 
	 * 	 
	 * 	 	
	 */
	virtual Chromosome_Binary* clone()
	{
		Chromosome_Binary* ch = new Chromosome_Binary(*this);
		return ch;
	}
private:	
	int nBitPerDim;
	int nDim;
};

#endif
