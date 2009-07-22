/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */

#include "CellularGA.h"
/*!
 * \brief
 * Constructor
 * 
 * \param populationWidth, populationHeight
 * Define the 2 dimensions of the population (width x height)
 * 
 * \param objectiveFunction
 * Objective function of type ObjectiveFunction
 * 
 * \param mutationOperator
 * Mutation operator
 * 
 * \param crossoverOperator
 * Crossover operator  
 *  
 * \remarks
 * The mutation, crossover must be implemented for the same coding type with population
 * 
 * \see
 * Population | ObjectiveFunction | Crossover | Mutation
 */
CellularGA::CellularGA(int populationWidth, int populationHeight, ObjectiveFunction* objectiveFunction, Mutation<double>* mutationOperator,
								Crossover<double>* crossoverOperator): popW(populationWidth), popH(populationHeight), mutation(mutationOperator), crossover(crossoverOperator)
{	
	this->fObj = objectiveFunction;	
	
	for(int i=0; i<popW*popH; i++)
	{
		vector<double> indv;
		for(unsigned j=0; j<fObj->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(fObj->lowerBounds[j], fObj->upperBounds[j]));
		}

		Chromosome_Real* ch = new Chromosome_Real(fObj->nDimensions(), fObj->lowerBounds, fObj->upperBounds);
		ch->fromDoubleVector(indv);
		ch->fitness = this->evaluate(indv);
		pop.push_back(ch);		
	}	
}

/*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGeneration
 * Number of generation to evolve.
 */
void CellularGA::evolve(unsigned int nGeneration )
{
	GlobalSearch<double>::evolve(nGeneration);
	int dx[] = {1, 0, -1, 0};
	int dy[] = {0, 1, 0, -1};

	int ndim = fObj->nDimensions();	
	Population<double> nextPop(pop.size());

	for(int i=0; i<popH; i++)
	{
		for(int j=0; j<popW; j++)
		{
			// parent 1 at (i, j)
			// parent 2 to be selected
			int xpar1 = j, ypar1 = i;
			int xpar2 = j, ypar2 = i;
			
			for(int d=0; d<4; d++)
			{
				// wrap around here 
				// neighbor of the last row/column is the first row/column
				int y = (i+dy[d] + popH) % popH;
				int x = (j+dx[d] + popW) % popW;

				// note that we are maximizing the fitness
				if (d == 0 || (pop[y*popH + x]->fitness > pop[ypar2*popH + xpar2]->fitness))
				{
					xpar2 = x;
					ypar2 = y;
				}
			}

			Chromosome_Real* parent1 = (Chromosome_Real*)pop[ypar1*popH + xpar1]->clone();
			Chromosome_Real* parent2 = (Chromosome_Real*)pop[ypar2*popH + xpar2]->clone();

			// crossover and mutation
			// note that we only use parent1 as it's the current position
			crossover->cross(*parent1, *parent2);
			mutation->mutate(*parent1);
		
			vector<double> tmp = parent1->toDoubleVector();
			parent1->fitness = GlobalSearch<double>::evaluate(tmp);

			// again, not that we are maximizing the fitness
			if (parent1->fitness > pop[ypar1*popH + xpar1]->fitness)
			{
				nextPop[ypar1*popH + xpar1] = parent1->clone();
			}
			else
			{
				nextPop[ypar1*popH + xpar1] = pop[ypar1*popH + xpar1]->clone();
			}

			delete parent1;
			delete parent2;
		}
	}
	
	// discard the current population
	for(unsigned i=0; i<pop.size(); i++)
	{
		delete pop[i];
	}
	
	// update with the next one
	pop = nextPop;
			
}
