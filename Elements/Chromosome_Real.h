/*!
 * <Real chromosome implemetation>
 * 
 * Copyright (c) 2008 by <Quang Huy/ NTU>
 */
#ifndef _Chromosome_Real_H_
#define _Chromosome_Real_H_


#include "Chromosome.h"
#include "../Global.h"


/*!
 * \brief
 * Implement a real-coded chromosome
 *  
 * \see
 * Chromosome | Chromosome_Binary
 */
class Chromosome_Real: public Chromosome<double>{

public:
	Chromosome_Real(const unsigned int nDimension, vector<double> lowerBounds, vector<double> upperBounds);
	
		
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
	virtual Chromosome_Real* clone()
	{
		Chromosome_Real* ch = new Chromosome_Real(*this);
		return ch;
	}
private:	
	
};

#endif
