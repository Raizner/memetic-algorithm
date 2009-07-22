/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2009 by <Quang Huy / NTU>
 */
#include "CellularMA.h"

/*!
 * \brief
 * Default constructor
 */ 

CellularMA::CellularMA()
{
	gs = NULL;
	ls = NULL;
	initialize();
}

/*!
 * \brief
 * Constructor.
 * 
 * \param globalSearch
 * Global search method of type GlobalSearch.
 * 
 * \param localSearch
 * Local search method of type LocalSearch.
 *    
 * \remarks
 * Global search must be implemented for the same coding type.
 * 
 * \see
 * GlobalSearch | LocalSearch
 */

CellularMA::CellularMA(CellularGA* globalSearch, LocalSearch* localSearch)
{
	gs = globalSearch;
	ls = localSearch;
	initialize();
}

/*!
 * \brief
 * Initialization.
 *  
 */
void CellularMA::initialize()
{
	pLS = 1;		
	maLearningStrategy = maLSLamarckian;
	if (gs != NULL)
	{
		tLS = 10*gs->fObj->nDimensions();
		maxEvaluations = 10000*gs->fObj->nDimensions();
		maxGenerations = 10000*gs->fObj->nDimensions() / gs->pop.size();
	}
	if(ls != NULL) ls->evaluationLimit = tLS;
}

// MA - best
/*!
 * \brief
 * Generate the set of individuals undergoing local search, based on the setting of
 * maStrategy.
 * 
 * \returns
 * Vector of positions of the individuals.
 *  
 *  
 * \see
 * MA_STRATEGY
 */

vector<int> CellularMA::getLSIndividuals()
{
	int popSize = gs->pop.size();
	int nLS = (int)(pLS * popSize);

	if (nLS == 0) return vector<int>();

	vector<int> LSIndividualSet;
	vector< pair<double, int> > individuals;	

	for(int i=0; i<popSize; i++)
	{
		individuals.push_back(make_pair(gs->pop[i]->fitness, i));
	}

	// sort all the fitness from max to min 
	// (note that the global search is to maximize)
	sort(individuals.rbegin(), individuals.rend());

	double range = (individuals[0].first - individuals[popSize-1].first) / nLS;

	// if the range is too small, we return an empty set
	if (range < 1e-5) return LSIndividualSet;

	// just simply choose the first individual in each range
	int pos = 0;
	for(int i=0; i<nLS; i++)
	{
		while (true)
		{
			if (pos == popSize || individuals[pos].first <= individuals[0].first - (i+1)*range) break;
			if (individuals[pos].first <= individuals[0].first - i*range)
			{
				LSIndividualSet.push_back(individuals[pos].second);
				break;
			}
			pos++;
		}
	}	

	return LSIndividualSet;
}

 /*!
 * \brief
 * Evolving algorithm
 * 
 * \param nGeneration
 * Number of generation to evolve.
 *
 * Typically, one of more generations of global search is followed by an individual learning phase.
 *  
 */

void CellularMA::evolve(unsigned int nGeneration)
{
	if (ls != NULL) ls->evaluationLimit = tLS;

	for(unsigned int i=0; i<nGeneration; i++)
	{		
		if (gs != NULL) gs->evolve();

		if (ls != NULL)
		{
			int nevaltmp = gs->fObj->nEvaluations;
			vector<int> v = getLSIndividuals();
			for(unsigned j=0; j<v.size(); j++)
			{
				//cout << -gs->pop[v[j]]->fitness << " --> ";
				vector<double> ch = gs->pop[v[j]]->toDoubleVector();
				//cout << ls->search(ch) << " == ";
				ls->search(ch);

				// if lamarckian learning then copy back
				if (maLearningStrategy == maLSLamarckian) gs->pop[v[j]]->fromDoubleVector(ch);	
				
				gs->pop[v[j]]->fitness = gs->evaluate(ch);

				//cout << -gs->pop[v[j]]->fitness << endl;
			}
			cout << "LS: " << v.size() << " times -> " << gs->fObj->nEvaluations - nevaltmp << " evals. " << endl;
		}
	}
}

