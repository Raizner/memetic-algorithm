/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#include "LocalSearch_DSCG.h"

LocalSearch_DSCG::LocalSearch_DSCG( ObjectiveFunction *f):LocalSearch( f )
{
}

LocalSearch_DSCG::~LocalSearch_DSCG()
{
}

/*!
 * \brief
 * Search
 * 
 * \param x
 * Initial point.
 * 
 * \returns
 * Final finess value.
 *  
 */
double LocalSearch_DSCG::search(vector<double>& x) 
{
	double fx = this->evaluate(x);
	this->dscg(x, fx);	
	return fx;
}

void LocalSearch_DSCG::dscg(vector<double>& x, double& fx)
{
	// number of dimensions
	unsigned n = fObj->nDimensions();	

	// initialize direction vectors v with unit vector
	vector< vector<double> > v(n+1, vector<double>(n, 0.0));

	// copy the step length
	vector<double> stepLen = stepLength;

	unsigned i, j, k;

	for(i=0; i<n; ++i) {
		v[i][i] = 1.0;
	}		

	// direction counter
	unsigned int firstdir = 0;

	// iteration counter
	int iteration = 0;

	vector< vector<double> > xsearch(n+2, vector<double>(n));
	vector<double> fxsearch(n+2);
	xsearch[0] = x;
	fxsearch[0] = fx;

	while(true)
	{
		// iteration count
		iteration++;

		// we are going to perform n line search
		// along the direction v[i] (i=1..n)
		// initial points x0, final points xn	

		for(i=firstdir; i<n; i++)
		{			
			linesearch(xsearch[i], fxsearch[i], v[i], stepLen[i], xsearch[i+1], fxsearch[i+1]);			
		}

		// after the n line search, we check the improvement
		double impv = 0;
		for(i=0; i<n; i++)
		{
			v[n][i] = xsearch[n][i] - xsearch[0][i]; 
			impv += v[n][i]*v[n][i];
		}
		impv = sqrt(impv);
		
		// significant improvement?
		// I'm not quite sure if we should replace 0.0 with some epsilon value
		if (impv > 0.0)
		{
			// one more line search
			// normalize the direction
			for(i=0; i<n; i++) v[n][i] /= impv;

			// we search with step length of the first direction 
			// since xsearch[n] will be come xsearch[0] in the next iteration
			linesearch(xsearch[n], fxsearch[n], v[n], stepLen[0], xsearch[n+1], fxsearch[n+1]);

			// and we have to calculate the improvement again
			impv = 0;
			for(i=0; i<n; i++)
			{			
				impv += (xsearch[n+1][i] - xsearch[0][i])*(xsearch[n+1][i] - xsearch[0][i]);
			}
			impv = sqrt(impv);
		}
		else 
		{
			xsearch[n+1] = xsearch[n];
			fxsearch[n+1] = fxsearch[n];
		}
		
		// now we need to check if the improvement is bigger than step length
		// note that here all the directions have the same value of step length 
		// (I don't know why, should they be different?)

		if (impv > stepLen[0])
		{
			// that's great, let's perform orthogonalization
			// but we should check first if we have more than one dimension to search (sounds weird huh?)
			if (n > 1)
			{
				// I'm going to write the orthogonalization here, but I'm too hungry
				
				// before the orthogonalization step, we need to order the direction so that
				// the improvement di > epsilon
								
				vector<double> d(n, 0.0);
				for(i=0; i<n; i++)
				{
					for (j=0; j<n; j++) 
					{
						if (fabs(v[i][j]) > 1e-16) 
						{ 
							d[i] = (xsearch[i+1][j] - xsearch[i][j]) / v[i][j]; 
							break;
						}
					}				
				}

				// we will sort the directions such that those with significant improvements are shifted
				// to the lower indices
				vector<int> indices;
				for(i=0; i<n; i++)
				{
					if (fabs(d[i]) >= this->accuracy) indices.push_back(i);
				}

				if (indices.size() <= 1)
				{
					// only one of the directions get improvement, we should reduce the step length
					// I reduce the whole stepLength array for 2 purposes:
					// 1. To be consistent
					// 2. For future use, if we want to use different steplengths for different directions.

					bool quitting = false;
					for(i=0; i<n; i++) 
					{
						stepLen[i] /= 10;
						// if steplength becomes too small, then we say bye bye
						if (stepLen[i] < this->accuracy)
						{
							quitting = true;
							break;
						}
					}
					

					if (quitting) 
					{ 
						// copy the last best values
						x = xsearch[n+1];
						fx = fxsearch[n+1];
						break;
					}

					xsearch[0] = xsearch[n+1];			
					fxsearch[0] = fxsearch[n+1];
					firstdir = 0;
				}
				else
				{
					// shift all the good director upfront
					unsigned pos = 0;

					unsigned nDirections = indices.size();
					for(i=0; i<indices.size(); i++)
					{						
						swap(d[indices[i]], d[i]);
						swap(v[indices[i]], v[i]);
					}
				
					// here we perform the Gram-Schmidt orthogonalization
					// no change in step length				
					vector< vector<double> > a ( n, vector<double>(n, 0.0) );
					vector< vector<double> > w ( n, vector<double>(n, 0.0) );

					// first of all, we calculate
					// a_i = sum_{j=i}^n d_j * v_j

					i = nDirections;
					do
					{
						i--;
						for(j=0; j<n; j++)
						{
							if (i == n-1)
							{
								a[i][j] = d[i]*v[i][j];
							}
							else
							{
								a[i][j] = a[i+1][j] + d[i]*v[i][j];
							}
						}
					} while( i > 0 );				

					// now the tough part
					// we calculate
					// w_0 = a_0
					// w_i = a_i - sum_{j=0}^{i-1}(a_i^T * v_j) * v_j (i = 0 .. n-1)

					for(i=0; i<nDirections; i++)
					{					
						for(j=0; j<i; j++)
						{						
							double tmp = 0;
							for(k=0; k<n; k++)
								tmp += a[i][k]*v[j][k];

							for(k=0; k<n; k++)
								w[i][k] += tmp*v[j][k];
						}
						
						double lenw = 0;
						for(j=0; j<n; j++)
						{
							w[i][j] = a[i][j] - w[i][j];
							lenw += w[i][j]*w[i][j];
						}

						// normalize
						lenw = sqrt(lenw);
						for(j=0; j<n; j++)
						{
							v[i][j] = w[i][j] / lenw;
						}
					}

					xsearch[0] = xsearch[n];
					fxsearch[0] = fxsearch[n];
					xsearch[1] = xsearch[n+1];					
					fxsearch[1] = fxsearch[n+1];
					firstdir = 1;
					//lastdir = indices.size() - 1;	
				}
			}
		}
		else
		{
			// no signification improvement achieved, we should reduce the step length
			// I reduce the whole stepLength array for 2 purposes:
			// 1. To be consistent
			// 2. For future use, if we want to use different steplengths for different directions.

			bool quitting = false;
			for(i=0; i<n; i++) 
			{
				stepLen[i] /= 10;
				// if steplength becomes too small, then we say bye bye
				if (stepLen[i] < this->accuracy)
				{
					quitting = true;
					break;
				}
			}
			
			if (quitting) 
			{ 
				// copy the last best values
				x = xsearch[n+1];
				fx = fxsearch[n+1];
				break;
			}

			
			xsearch[0] = xsearch[n+1];
			fxsearch[0] = fxsearch[n+1];			
			firstdir = 0;
		}
		if (done())
		{
			x = xsearch[n+1];
			fx = fxsearch[n+1];			
			break;
		}
	}
}



void LocalSearch_DSCG::linesearch(vector<double>& x, double& fx, vector<double>& dir, double step, vector<double> &xend, double& fxend)
{
	bool bothTrialsFailed = false;

	unsigned i;
	unsigned n = fObj->nDimensions();
	
	// previous, current and next search points.	
	vector<double> prev(n), current(n), next(n);
	double fprev, fcurrent, fnext;

	current = x;
	fcurrent = fx;

	// first trial
	for(i=0; i<n; i++)
	{
		next[i] = current[i] + dir[i]*step;
	}	
	fnext = this->evaluate(next);
	
	// if next is not better than current
	if (fcurrent < fnext)
	{
		// second trial, in the opposite direction
		step = -step;
		prev = next;
		fprev = fnext;

		for(i=0; i<n; i++)
		{
			next[i] = current[i] + dir[i]*step;
		}	
		fnext = this->evaluate(next);

		if (fcurrent < fnext) bothTrialsFailed = true;
	}

	if (!bothTrialsFailed)
	{
		// start moving in the direction
		// after each successful move, we double the step length
		while (fnext <= fcurrent)
		{
			// double the step length, set the current point to new point
			step *= 2;
			prev = current; fprev = fcurrent;
			current = next; fcurrent = fnext;

			for(i=0; i<n; i++)
			{
				next[i] = current[i] + dir[i]*step;
			}	
			fnext = evaluate(next);
		}
		
		// here we have fnext > fcurrent, we are going doing interpolation
		// but it requires a bit of preparation
		// since next - current = 2(current - prev)
		// we are considering another point in the middle of (current, next)
		step /= 2;

		vector<double> tmpPoint(n);		
		double ftmpPoint;

		for(i=0; i<n; i++)
		{		
			tmpPoint[i] = current[i] + dir[i]*step;
		}	
		ftmpPoint = this->evaluate(tmpPoint);
		
		if (ftmpPoint <= fcurrent)
		{
			prev = current; fprev = fcurrent;
			current = tmpPoint; fcurrent = ftmpPoint;
		}
		else
		{
			next = tmpPoint; fnext = ftmpPoint;
		}
	}

	// Here, we either have both trials failed or we have finished searching in the line
	// Whatever the case is, we have 3 points prev, current, next satisfying
	// 1. fcurrent <= fprev and fcurrent <= fnext
	// 2. current - prev = next - current i.e. current is the mid-point
	// We are going to do a quadratic interpolation
	// result = current + dir*step*(fprev - fnext) / (2*(fprev + fnext - 2*fcurrent))	
	if (fprev == fcurrent && fnext == fcurrent)
	{
		xend = next;
		fxend = fnext;
	}
	else
	{
		double ratio = step*(fprev - fnext) / (2*(fprev + fnext - 2*fcurrent));
		for(i=0; i<n; i++)
		{
			xend[i] = current[i] + dir[i]*ratio;
		}
		fxend = evaluate(xend);

		if (fxend > fcurrent) { xend = current; fxend = fcurrent; }
	}	
}