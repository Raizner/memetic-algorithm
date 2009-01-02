/*!
 * <Chromosome class>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Chromosome_H_
#define _Chromosome_H_

#include <vector>
#include <sstream>

using namespace std;

/*!
 * \brief
 * The super class for different type of chromosomes.
 * 
 * \param T
 * Type T can be bool (binary problems), int (combinatorial problems) or real (continuous problems).
 *  
 *  
 * \see
 * Chromosome_Binary | Chromosome_Real
 */
template<typename T> class Chromosome: public vector<T>{

public:	
	/*!
	 * \brief
	 * Default constructor.
	 * 	 
	 * 
	 * Create an instance of chromosome of type T.
	 * 	 	 
	 * 	 
	 */
	Chromosome() : std::vector<T>()
	{
	}

	/*!
	 * \brief
	 * Constructor	 
	 * 
	 * \param v
	 * Input vector
	 * 	 	 
	 * Create the chromosome from a vector of type T.
	 * 
	 * \remarks
	 * Actually a copy constructor.
	 * 	 
	 */
	Chromosome(const vector<T>& v):std::vector<T>(v)
	{
	}


	/*!
	 * \brief
	 * Constructor
	 * 
	 * \param size
	 * Chromosome size, normally equal to the number of optimization variables.
	 * 
	 * \param lowerBounds
	 * Vector of lower bound values.
	 * 
	 * \param upperBounds
	 * Vector of upper bound values.
	 * 	 
	 * 
	 * Create the chromosome of given size and bounds.
	 */
	Chromosome(const size_t size, vector<double> lowerBounds, vector<double> upperBounds):std::vector<T>(size), lowerBounds(lowerBounds), upperBounds(upperBounds)
	{
	}

	/*!
	 * \brief
	 * Constructor
	 * 
	 * \param phenotype
	 * Input phenotype.
	 * 	 
	 * \param lowerBounds
	 * Vector of lower bound values.
	 * 
	 * \param upperBounds
	 * Vector of upper bound values.
	 * 
	 * Generate the chromosome from a vector of doubles.
	 * 
	 * \remarks
	 * Same as Chromosome<T>::fromDoubleVector(v).
	 * 
	 * \see
	 * Chromosome<T>::fromDoubleVector
	 */
	Chromosome(vector<double> phenotype, vector<double> lowerBounds, vector<double> upperBounds):lowerBounds(lowerBounds), upperBounds(upperBounds)
	{
		this->fromDoubleVector(phenotype);
	}

	/*!
	 * \brief
	 * Fitness value.
	 * 	 
	 * 	 	 
	 * \see
	 * ObjectiveFunction
	 */
	double fitness;
	
	/*!
	 * \brief
	 * Vector of upper bound values.
	 * 	 
	 * 	 
	 * \see
	 * lowerBounds
	 */
	vector<double> upperBounds;
	
	/*!
	 * \brief
	 * Vector of lower bound values.
	 * 	 
	 * 	 
	 * 
	 * \see
	 * upperBounds
	 */
	vector<double> lowerBounds;

	
	/*!
	 * \brief
	 * Encode function
	 * 
	 * \param phenotype
	 * Phenotype (optimization variables).
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Encode a phennotype to its genotype representation
	 * 
	 * \remarks
	 * Write remarks for fromDoubleVector here.
	 * 
	 * \see
	 * toDoubleVector()
	 */
	virtual void fromDoubleVector(vector<double>& phenotype) {}
	
	
	
	/*!
	 * \brief
	 * Decode function.
	 * 
	 * \returns
	 * The decoded vector of doubles.
	 * 		 
	 * Decode a genotype to a double vector	 	 
	 * 
	 * \see
	 * fromDoubleVector()
	 */
	virtual vector<double> toDoubleVector() { return vector<double>() ;}

	/*!
	 * \brief
	 * Clone function
	 * 
	 * \returns
	 * A "twin" of the current chromosome	 
	 * 	 
	 * 	 	
	 */
	virtual Chromosome* clone() { Chromosome* ch = new Chromosome(*this); return ch; }
	
	bool operator < (Chromosome<T>& chromosome);
	
	virtual string toString();

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
	bool isInBound()
	{
		vector<double> vd = this->toDoubleVector();
		unsigned i;
		for(i=0; i<vd.size(); i++)
		{
			if (vd[i] < lowerBounds[i] || vd[i] > upperBounds[i]) return false;
		}

		return true;
	}

protected:
	
};

/*!
 * \brief
 * Compare 2 chromosome.
 * 
 * \param chromosome
 * The chromosome to compare.
 * 
 * \returns
 * True if the current chromosome has lower fitness value, false otherwise.
 *  
 * 
 * Compare the fitness value with another chromosome.
 *  
 * 
 */
template<typename T>
bool Chromosome<T>::operator < (Chromosome<T>& chromosome)
{
	return this->fitness < chromosome.fitness;
}

/*!
 * \brief
 * Convert the chromosome to string.
 * 
 * \returns
 * The string representing the chromosome.
 *  
 * 
 * Convert the chromosome to string.
 * 
 * \remarks
 * To be inherited by different chromosome classes.
 * E.g.: we want a concatenated binary strings, but a space-delimited or comma-delimited double string,
 * i.e., 0100110010010 vs. {1.2, 0.3, 5.9}
 *  
 */
template<typename T>
string Chromosome<T>::toString()
{
	ostringstream oss;
	for(unsigned int i=0; i<this->size(); i++)
	{
		if (i > 0) oss << " ";
		oss << this->at(i);
	}

	return oss.str();
}

#endif
