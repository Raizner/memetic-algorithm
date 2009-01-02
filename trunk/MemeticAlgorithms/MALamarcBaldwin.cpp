/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#include "MALamarcBaldwin.h"

#include "../Elements/Population.h"
#include "../Elements/Chromosome_Real.h"
#include "../Elements/Chromosome_Binary.h"

// Mutations
#include "../Operators/Mutations/Mutation_BitFlip.h"
#include "../Operators/Mutations/Mutation_Gaussian.h"

// Crossovers
#include "../Operators/Crossovers/Crossover_nPoint.h"
#include "../Operators/Crossovers/Crossover_Uniform.h"

// Scaling & Selections 
#include "../Operators/Scalings/Scaling_Linear.h"
#include "../Operators/Selections/Selection_RouletteWheel.h"

// Recombinations
#include "../Operators/Recombinations/Recombination_KeepBest.h"
#include "../Operators/Recombinations/Recombination_MuPlusLambda.h"
#include "../Operators/Recombinations/Recombination_MuCommaLambda.h"


/*!
 * \brief
 * Default constructor
 */ 

MALamarcBaldwin::MALamarcBaldwin(ObjectiveFunction* f)
{
	this->f = f;
	gs1 = NULL;
	gs2 = NULL;
	ls = NULL;
	initialize();
}


/*!
 * \brief
 * Initialization.
 *  
 */

void MALamarcBaldwin::initialize()
{
	unsigned popsz1 = 10;
	unsigned popsz2 = 10;
	unsigned i, j;
	Population<double> pop1;
	Population<bool> pop2;

	// we have 2 populations here
	// one binary, one real
	// The real-coded population will have 
	// half of the population be the same as the binary one
	for(i=0; i<popsz1; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->lowerBounds[j], f->upperBounds[j]));
		}

		Chromosome_Real* ch = new Chromosome_Real(f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch->fromDoubleVector(indv);
				
		pop1.push_back(ch);		
	}

	for(i=0; i<popsz2; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->lowerBounds[j], f->upperBounds[j]));
		}

		Chromosome_Real* ch1 = new Chromosome_Real(f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch1->fromDoubleVector(indv);
		Chromosome_Binary* ch2 = new Chromosome_Binary(20, f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch2->fromDoubleVector(indv);
		
		pop1.push_back(ch1);
		pop2.push_back(ch2);		
	}
	
	Mutation<double>* mut1 = new Mutation_Gaussian(0.01);		
	Crossover<double>* crs1 = new Crossover_Uniform<double>(0.8);	
	Selection<double>* sel1 = new Selection_RouletteWheel<double>(new Scaling_Linear());
	Recombination<double>* rec1 = new Recombination_KeepBest<double>(1);
	
	gs1 = new GeneticAlgorithm<double>(pop1, f, mut1, crs1, sel1, rec1);
	
	
	Mutation<bool>* mut2 = new Mutation_BitFlip(0.002);
	Crossover<bool>* crs2 = new Crossover_Uniform<bool>(0.8);
	Selection<bool>* sel2 = new Selection_RouletteWheel<bool>(new Scaling_Linear());
	Recombination<bool>* rec2 = new Recombination_KeepBest<bool>(0);

	gs2 = new GeneticAlgorithm<bool>(pop2, f, mut2, crs2, sel2, rec2);

	pLS = 1;	
	maStrategy = maLSBest;

	tLS = 10*f->nDimensions();
	maxEvaluations = 10000*f->nDimensions();
	maxGenerations = 100000;
	if(ls != NULL) ls->evaluationLimit = tLS;
}

/*!
 * \brief
 * Number of generations elapsed.
 */

int MALamarcBaldwin::nGenerations()
{
	return gs1->nGenerations();
}

/*!
 * \brief
 * Check if stopping criteria are satisfied.
 * 
 * \returns
 * True if one of the stopping conditions is true, false otherwise.
 *  
 */

bool MALamarcBaldwin::done()
{
	return (gs1->nGenerations() >= maxGenerations || f->nEvaluations >= maxEvaluations);
}


 /*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGenerations
 * Number of generation to evolve.
 *
 * Typically, one of more generations of global search is followed by an individual learning phase.
 *  
 */

void MALamarcBaldwin::evolve(unsigned int nGenerations)
{
	if (ls != NULL) ls->evaluationLimit = tLS;

	for(unsigned int i=0; i<nGenerations; i++)
	{

		cout << "gs1: " << endl;
		gs1->evolve();
		for(unsigned j=0; j<gs1->pop.size(); j++)
			cout << gs1->pop[j]->fitness << " ";
		cout << endl;

		cout << "gs2: " << endl;
		gs2->evolve();
		for(unsigned j=0; j<gs2->pop.size(); j++)
			cout << gs2->pop[j]->fitness << " ";
		cout << endl;

		gs1->pop.sort();

		if (ls != NULL)
		{		
			// for the 1st population - the real coded
			// we apply lamarckian learning
			for(unsigned j=0; j<(gs1->pop.size()-gs2->pop.size()); j++)
			{	
				cout << -gs1->pop[j]->fitness << " --> ";
				vector<double> ch = gs1->pop[j]->toDoubleVector();
				cout << ls->search(ch) << " == ";
				gs1->pop[j]->fromDoubleVector(ch);	
				gs1->pop[j]->fitness = gs1->evaluate(ch);
				cout << -gs1->pop[j]->fitness << endl;
			}


			// for the 2nd population - the binary coded
			// we apply baldwinian learning
			for(unsigned j=0; j<gs2->pop.size(); j++)
			{
				cout << -gs2->pop[j]->fitness << " --> ";
				vector<double> ch = gs2->pop[j]->toDoubleVector();
				cout << ls->search(ch) << " == ";				
				gs2->pop[j]->fitness = gs2->evaluate(ch);
				cout << -gs2->pop[j]->fitness << endl;
			}			
		}

		// finally, we replace the second part of the 1st population
		// with the chromosomes in the 2nd one to prepare for the next
		// generation

		cout << "merge: " << endl;
		for(unsigned j=0; j<gs2->pop.size(); j++)
		{
			vector<double> tmp = gs2->pop[j]->toDoubleVector();
			gs1->pop[gs1->pop.size()-gs2->pop.size()+j]->fromDoubleVector(tmp);
			gs1->pop[gs1->pop.size()-gs2->pop.size()+j]->fitness = gs2->pop[j]->fitness;
		}

		for(unsigned j=0; j<gs1->pop.size(); j++)
			cout << gs1->pop[j]->fitness << " ";
		cout << endl;
	}
}
