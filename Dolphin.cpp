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

#include "EvolutionaryAlgorithms/CMAES.h"

#include "Utilities/Statistics.h"

using namespace std;

/* the objective (fitness) function to be minized */
double fitfun(double const *x, int N) { /* function "cigtab" */
  int i; 
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)  
    sum += x[i]*x[i]; 
  return sum;
}

void testES()
{
	CMAES::cmaes_t evo; /* an CMA-ES type struct or "object" */
  double *arFunvals, *const*pop, *xfinal;
  int i; 

  CMAES cma;
  /* Initialize everything into the struct evo, 0 means default */
  arFunvals = cma.cmaes_init(&evo, 0, NULL, NULL, 0, 0, "initials.par"); 
  printf("%s\n", cma.cmaes_SayHello(&evo));
  cma.cmaes_ReadSignals(&evo, "signals.par");  /* write header and initial values */

  /* Iterate until stop criterion holds */
  while(!cma.cmaes_TestForTermination(&evo))
    { 
      /* generate lambda new search points, sample population */
      pop = cma.cmaes_SamplePopulation(&evo); /* do not change content of pop */

      /* Here you may resample each solution point pop[i] until it
	 becomes feasible, e.g. for box constraints (variable
	 boundaries). function is_feasible(...) needs to be
	 user-defined.  
	 Assumptions: the feasible domain is convex, the optimum is
	 not on (or very close to) the domain boundary, initialX is
	 feasible and initialStandardDeviations are sufficiently small
	 to prevent quasi-infinite looping.
      */
      /* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) 
	   while (!is_feasible(pop[i])) 
	     cmaes_ReSampleSingle(&evo, i); 
      */

      /* evaluate the new search points using fitfun from above */ 
      for (i = 0; i < cma.cmaes_Get(&evo, "lambda"); ++i) {
			arFunvals[i] = fitfun(pop[i], (int) cma.cmaes_Get(&evo, "dim"));
      }

      /* update the search distribution used for cmaes_SampleDistribution() */
      cma.cmaes_UpdateDistribution(&evo, arFunvals);  

      /* read instructions for printing output or changing termination conditions */ 
      cma.cmaes_ReadSignals(&evo, "signals.par");   
      fflush(stdout); /* useful in MinGW */
    }
  printf("Stop:\n%s\n",  cma.cmaes_TestForTermination(&evo)); /* print termination reason */
  cma.cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

  /* get best estimator for the optimum, xmean */
  xfinal = cma.cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
  cma.cmaes_exit(&evo); /* release memory */ 

  /* do something with final solution and finally release memory */
  free(xfinal);   
}

void testLS()
{	
	string filename;
	Statistics s;
	for(int mul = 1; mul < 1000; mul*=10)
	{
		for(int ndim = 2; ndim<=10; ndim++)
		{
			cout << "Runing on " << ndim*mul << " dimensions...";
			ObjectiveFunction* f = new FSphere(ndim*mul, -100, 100);
			f->statModule = &s;
			vector<double> v(f->nDimensions());

			LocalSearch_ES ls(f);
			ls.evaluationLimit = 100000;
			ls.accuracy = 1e-5;
			
			for(int samples = 0; samples < 100; samples++)
			{
				ostringstream oss;
				oss << "D:/HUYNQ/Works/My papers/Gecco 2009/Data/Sphere/ES/output_" << ndim*mul << "_" << samples;
				s.startRecording(oss.str().c_str());
				for(int i=0; i<f->nDimensions()-1; i++) v[i] = Rng::uni(-100, 100);
				//printf("%.12lf --> ", (*f)(v));
				double res = ls(v);
				//printf("%.12lf\n", res);
				s.stopRecording();
			}
			
			delete f;
			cout << endl;
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
	ObjectiveFunction* f = new FScaffer(30, -100, 100);
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
	ls1->stepLength = vector<double>(f->nDimensions(), 0.7);
	ls2->stepLength = vector<double>(f->nDimensions(), 0.7);


	//ma.lsPool.push_back(ls1);
	//ma.lsPool.push_back(ls2);
	ma.ls = ls1;
	ls1->evaluationLimit = 1000;
	//ls3->evaluationLimit = 300;
	
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

	//testES();
	testLS();
	return 0;
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

