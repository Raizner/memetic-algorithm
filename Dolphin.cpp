/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

// Dolphin.cpp : Defines the entry point for the console application.
//

#include <iostream>

#include "Includes/Functions.h"
#include "Includes/GlobalSearchers.h"
#include "Includes/LocalSearchers.h"
#include "Includes/Operators.h"

#include "Elements/Chromosome_Binary.h"
#include "Elements/Chromosome_Real.h"

#include "MemeticAlgorithms/MemeticAlgorithm.h"
#include "MemeticAlgorithms/MAMetaLamarckian.h"
#include "MemeticAlgorithms/MALamarcBaldwin.h"

using namespace std;

void testLS()
{
	ObjectiveFunction* f = new FSphere(2, -100, 100);
	vector<double> v(2, 4.0);

	/*
	for(v[0] = -4; v[0] <= 4; v[0]+=2)
		for(v[1] = -4; v[1] <= 4; v[1]+=2)
		{
			vector<double> g = f->gradient(v);
			cout << v[0] << " " << v[1] << " " << g[0] << " " << g[1] << endl;
		}
	*/
	/*
	for(int i=0; i<30; i++)
	{
		v[i] = Rng::uni(-5, 5);
	}
	*/


	//LocalSearch_DFP ls(f);
	//LocalSearch_DSCG ls(f);
	LocalSearch_ES ls(f);
	ls.evaluationLimit = 100000;
	//ls.stepLength = vector<double>(2, 0.9);	
	ls.accuracy = 1e-5;

	int tmp = 0;
	for(v[0] = -4; v[0] <= 4; v[0]+=2)
	{
		for(v[1] = -4; v[1] <= 4; v[1]+=2)			
		{
			tmp++;
			//v[0] = -4;
			//v[1] = 2;
			//v[2] = 4;
			vector<double> v2 = v;
			//cout << tmp << ": (" << v[0] << "," << v[1] << " " << v[2] << ") - " << (*f)(v) << " --> ";
			//cout << tmp << ": (" << v[0] << "," << v[1] << ") - " << (*f)(v) << " --> ";
			printf("%.12lf\n", (*f)(v));
			//f->nEvaluations = 0;
			double res = ls(v2);
			//cout << "(" << v2[0] << " " << v2[1] << " " << v2[2] << ") - " << res << " " << f->nEvaluations << endl;	
			printf("%.12lf\n", res);
		}
	}
}

void testGA()
{
	
	unsigned int i, j;

	//ObjectiveFunction* f = new FSchwefel102(30);
	ObjectiveFunction* f = new FRastrigin(30);
	//ObjectiveFunction* f = new FSchwefel102Noisy(30);

	for(i=0; i<f->nDimensions(); i++)
	{
//		f->upperBounds[i] = 100.0 - 2*i;
//		f->lowerBounds[i] = -100.0 + 2*i;
	}

	cout << "Translation vector:" << endl;
	for(i=0; i<f->nDimensions(); i++)
	{
//		f->translationVector[i] = Rng::uni();
		cout << f->translationVector[i] << " ";
	}
	cout << endl;

	cout << "Test evaluation: " << (*f)(f->translationVector) << endl;

	f->nEvaluations = 0;

	Population<double> pop;


	for(i=0; i<20; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->lowerBounds[j], f->upperBounds[j]));
		}

		Chromosome_Real* ch = new Chromosome_Real(f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch->fromDoubleVector(indv);
				
		pop.push_back(ch);		
	}
	
	Mutation<double>* mut = new Mutation_Gaussian(0.01);		
	Crossover<double>* crs = new Crossover_Uniform<double>(0.8);	
	Selection<double>* sel = new Selection_RouletteWheel<double>(new Scaling_Linear());
	Recombination<double>* rec = new Recombination_KeepBest<double>(1);


	GeneticAlgorithm<double> ga(pop, f, mut, crs, sel, rec);

	while(f->nEvaluations < f->nDimensions() * 10000)
	{					
		ga.evolve();
		//	if (ma.nGenerations() % 50 == 0)
		{
			cout << f->nEvaluations << " evals: " << f->bestEvaluation() << endl;
		}
	}

	cout << "Best solution so far: " << endl;
	for(unsigned int j=0; j<f->bestSolution().size(); j++) cout << f->bestSolution()[j] << " ";
	cout << endl << "Fitness: " << f->bestEvaluation() << endl;
}

void testMA()
{
	unsigned int i, j;

	//ObjectiveFunction* f = new FSchwefel102(30);
	ObjectiveFunction* f = new FAckley(30);
	//ObjectiveFunction* f = new FSchwefel102Noisy(30);

	for(i=0; i<f->nDimensions(); i++)
	{
//		f->upperBounds[i] = 100.0 - 2*i;
//		f->lowerBounds[i] = -100.0 + 2*i;
	}

	cout << "Translation vector:" << endl;
	for(i=0; i<f->nDimensions(); i++)
	{
//		f->translationVector[i] = Rng::uni();
		cout << f->translationVector[i] << " ";
	}
	cout << endl;

	cout << "Test evaluation: " << (*f)(f->translationVector) << endl;


	f->nEvaluations = 0;

	/* binary-coded GA here
	Population<bool>* pop = new Population<bool>();

	unsigned int i, j;
	for(i=0; i<50; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->getLowerBounds(j), f->getUpperBounds(j)));
		}

		Chromosome<bool>* ch = new Chromosome_Binary(20, f->nDimensions(), f->getLowerBounds(), f->getUpperBounds());
		ch->fromDoubleVector(indv);
//		vector<double> vd = ch->toDoubleVector();
//		ch->fitness = (*f)(indv);
		pop->push_back(ch);
	}

	Scaling* scl = new Scaling_Linear();
	Mutation<bool>* mut = new Mutation_BitFlip(0.002);
	Crossover<bool>* crs = new Crossover_Uniform<bool>(0.8);
	Selection<bool>* sel = new Selection_RouletteWheel<bool>(scl);
	Recombination<bool>* rec = new Recombination_KeepBest<bool>(0);
	*/

	Population<double> pop;


	for(i=0; i<50; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->lowerBounds[j], f->upperBounds[j]));
		}

		Chromosome_Real* ch = new Chromosome_Real(f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch->fromDoubleVector(indv);
				
		pop.push_back(ch);		
	}
	
	Mutation<double>* mut = new Mutation_Gaussian(0.01);		
	Crossover<double>* crs = new Crossover_Uniform<double>(0.8);	
	Selection<double>* sel = new Selection_RouletteWheel<double>(new Scaling_Linear());
	Recombination<double>* rec = new Recombination_KeepBest<double>(1);


	GeneticAlgorithm<double>* ga = new GeneticAlgorithm<double>(pop, f, mut, crs, sel, rec);
	
	//LocalSearch* ls = NULL;
	
	
	MemeticAlgorithm<double> ma(ga);
	//MAMetaLamarckian<double> ma(ga);

	LocalSearch* ls1 = new LocalSearch_DSCG(f);
	LocalSearch* ls2 = new LocalSearch_DFP(f);
	LocalSearch* ls3 = new LocalSearch_ES(f);
	ls1->stepLength = vector<double>(f->nDimensions(), 1.0);
	ls2->stepLength = vector<double>(f->nDimensions(), 0.7);


	//ma.lsPool.push_back(ls1);
	//ma.lsPool.push_back(ls2);
	ma.ls = ls3;
	ls3->evaluationLimit = 300;
	
	ma.pLS = 0.1;
	ma.maSelectionStrategy = ma.maLSBest;
	//ma.maLearningStrategy = ma.maLSBaldwinian;
	ma.maLearningStrategy = ma.maLSLamarckian;

	//ma.maxEvaluations = 100000;
	while(!ma.done())
	{					
		ma.evolve();
		//	if (ma.nGenerations() % 50 == 0)
		{
			cout << f->nEvaluations << " evals: " << f->bestEvaluation() << endl;
		}
	}

	cout << "Best solution so far: " << endl;
	for(unsigned int j=0; j<f->bestSolution().size(); j++) cout << f->bestSolution()[j] << " ";
	cout << endl << "Fitness: " << f->bestEvaluation() << endl;
}

void testMALamaBald()
{
	unsigned int i, j;

	ObjectiveFunction* f = new FRastrigin(30);

	for(i=0; i<f->nDimensions(); i++)
	{
//		f->upperBounds[i] = 100.0 - 2*i;
//		f->lowerBounds[i] = -100.0 + 2*i;
	}

	cout << "Translation vector:" << endl;
	for(i=0; i<f->nDimensions(); i++)
	{
//		f->translationVector[i] = Rng::uni();
		cout << f->translationVector[i] << " ";
	}
	cout << endl;

	cout << "Test evaluation: " << (*f)(f->translationVector) << endl;


	f->nEvaluations = 0;

	
	Population<double> pop;

	for(i=0; i<20; i++)
	{
		vector<double> indv;
		for(j=0; j<f->nDimensions(); j++)
		{
			indv.push_back(Rng::uni(f->lowerBounds[j], f->upperBounds[j]));
		}

		Chromosome_Real* ch = new Chromosome_Real(f->nDimensions(), f->lowerBounds, f->upperBounds);
		ch->fromDoubleVector(indv);
		ch->fitness = f->evaluate(indv);

		
		pop.push_back(ch);
	}
	
	
	for(i=0; i<20; i++) cout << pop[i]->fitness << " ";
	cout << endl;

	pop.sort();

	for(i=0; i<20; i++) cout << pop[i]->fitness << " ";
	cout << endl;

	MALamarcBaldwin ma(f);
	LocalSearch* ls1 = new LocalSearch_DSCG(f);
	ma.ls = ls1;

	while(!ma.done())
	{					
		ma.evolve();
		//	if (ma.nGenerations() % 50 == 0)
		{
			cout << f->nEvaluations << " evals: " << f->bestEvaluation() << endl;
		}
	}

}

int main(int argc, char* argv[])
{
	Rng::seed((long)time(NULL));

	int option = 0;
	if (argc == 1)
	{
		cout << argv[0] << " GA|MA|MA2|MS" << endl;
	}
	else
	{
		if (strcmp(argv[1], "MA") == 0) option = 0;
		else if (strcmp(argv[1], "MA2") == 0) option = 1;
		else if (strcmp(argv[1], "MS") == 0) option = 2;
		else if (strcmp(argv[1], "GA") == 0) option = 3;
	}

	switch(option)
	{
		case 0: testMA(); break;		
		case 1: testMALamaBald(); break;
		case 2: testLS(); break;
		case 3: testGA(); break;
	}
	return 0;
}

