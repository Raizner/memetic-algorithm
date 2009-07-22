/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Population_H_
#define _Population_H_

#include <vector>
#include "Chromosome.h"

using namespace std;

/*!
 * \brief
 * Population class.
 * 
 * \param T
 * Type T can be bool (binary problems), int (combinatorial problems) or real (continuous problems).
 * 
 * Implementation of population class.
 * 
 * \remarks
 * Used in population-based search method
 * 
 * \see
 * GlobalSearch< T >
 */
template<typename T> class Population : public std::vector< Chromosome<T>* > {

public:
	/*!
	 * \brief
	 * Constructor.
	 * 	 
	 * 
	 * Create a population from a vector of Chromosome.
	 * 
	 * \remarks
	 * Actually a copy constructor
	 * 
	 * \see
	 * Chromosome<T>
	 */
	Population():std::vector< Chromosome<T>* >()
	{
	}

	/*!
	 * \brief
	 * Constructor.
	 * 
	 * \param size
	 * Population size.
	 * 	 
	 * 
	 * Create a population with given size.
	 * 	 
	 */
	Population(const size_t size):std::vector< Chromosome<T>* >(size) 
	{ 
	}

	virtual Population* clone()
	{
		Population* p = new Population();
		for(unsigned i = 0; i < this->size(); i++)
		{
			p->push_back(this->at(i)->clone());
		}
		return p;
	}
	
	/*!
	 * \brief
	 * Constructor.
	 * 
	 * \param size
	 * Population size.
	 * 
	 * \param chromosomeLength
	 * Length of chromosome (number of optimization variables).
	 * 	 
	 * 
	 * Create a population with given size, and fill it with chromosomes of given length.
	 * 	 
	 * 
	 * \see
	 * Chromosome<T>
	 */
	Population(const size_t size, const size_t chromosomeLength):std::vector< Chromosome<T>* >(size, new Chromosome<T>*(chromosomeLength)) 
	{ 
	}
	
	void sort()
	{
		unsigned i, j;

		for(i=0; i<this->size(); i++)
		{
			for(j=i+1; j<this->size(); j++)
			{
				if (this->at(i)->fitness < this->at(j)->fitness)
				{
					std::swap(this->at(i), this->at(j));
				}
			}
		}
	}
		
private:
	// To do list
	// Statistical functions
	// Assignment and comparison operators
	// Data structure-related methods: sort, insert, etc.
};

#endif
