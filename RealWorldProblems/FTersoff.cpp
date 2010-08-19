/*! Potential model for Ge(n)Si(m) molecule.
Reference

Tersoff, New empirical model for the structural properties of silicon, Phys. Rev. Lett. 56, 632 (1986).

Copyright (c) 2010 by <Quang Huy / NTU>
*/

#include "FTersoff.h"

FTersoff::FTersoff(int nAtoms, int nSi, double parameters[11][2]) : ObjectiveFunction(-1, nAtoms * 3, INT_MIN, INT_MAX)
{
	this->nAtoms = nAtoms;
	this->nSi = nSi;

	for(int i=0; i<11; i++)
	{		
		this->params[i][0] = parameters[i][0];
		this->params[i][1] = parameters[i][1];
	}
}

double FTersoff::evaluate_( vector<double>& x )
{
	int i, j;
	double _x[100][3];

	int cnt = 0;
	// number of atoms
	int nAtoms = x.size() / 3;
	int nO = (nAtoms + 1) / 3;

	RP(i, nAtoms)
	{
		RP(j, 3) _x[i][j] = x[cnt++];
	}

	double res = energy(_x);

	return res;
}

void FTersoff::init(double _r[][3])
{
	int tmp;
	for(int i=0; i<nAtoms; i++)
	{
		if (i<nSi) tmp = 0; else tmp = 1;
		
		_A[i] = params[A_][tmp];
		_B[i] = params[B_][tmp];
		_S[i] = params[S_][tmp];
                _R[i] = params[R_][tmp];
		_lambda[i] = params[lambda_][tmp];
                _mu[i] = params[mu_][tmp];
		c[i] = params[c_][tmp];
                d[i] = params[d_][tmp];
		n[i] = params[n_][tmp];
                h[i] = params[h_][tmp];
		beta[i] = params[A_][tmp];
		for(int j=0; j<nAtoms; j++) 
		{
			chi[i][j] = 1.00061;
			r[i][j] = sqrt(pow(_r[i][0]-_r[j][0],2)+pow(_r[i][1]-_r[j][1],2)+pow(_r[i][2]-_r[j][2],2));
		}
	}

	cout << _A[0] << " " << _A[1] << endl;

	
}

double FTersoff::energy(double _r[][3])
{
	double res = 0;
	init(_r);

	for(int i=0; i<nAtoms; i++)	
	{
		for(int j=0; j<nAtoms; j++)
		{
			if (j == i) continue;
			cout << i << " " << j << endl;

			lambda[i][j] = (_lambda[i] + _lambda[j]) / 2;
			mu[i][j] = (_mu[i] + _mu[j]) / 2;
			A[i][j] = sqrt(_A[i]*_A[j]);
			B[i][j] = sqrt(_B[i]*_B[j]);
			R[i][j] = sqrt(_R[i]*_R[j]);
            S[i][j] = sqrt(_S[i]*_S[j]);			
			
			fr[i][j] = A[i][j]*exp(-lambda[i][j]*r[i][j]);
			fa[i][j] = -B[i][j]*exp(-mu[i][j]*r[i][j]);

			if (r[i][j] < R[i][j])
			{
				fc[i][j] = 1;
			}
			else if (r[i][j] < S[i][j])
			{
				fc[i][j] = 0.5 + 0.5*cos(PI*(r[i][j] - R[i][j])/(S[i][j] - R[i][j]));
			}
			else
			{
				fc[i][j] = 0;
			}

			zeta[i][j] = 0.0;
			for(int k=0; k<nAtoms; k++)
			{
			
				if (k!=i && k!=j) 
				{
					theta[i][j][k] = r[i][j]*r[i][j] + r[j][k]*r[j][k] - r[i][k]*r[i][k] / 2*r[i][j]*r[j][k];
					g[i][j][k] = 1+c[i]*c[i]*( 1 / d[i]*d[i] - 1 / (d[i]*d[i] + pow(h[i]-cos(theta[i][j][k]), 2)) );
					zeta[i][j] += fc[i][k]*omega[i][k]*g[i][j][k];
				}
			}
			b[i][j] = chi[i][j]*(1+pow(beta[i], n[i]) * pow(zeta[i][j], n[i]));
			V[i][j] = fc[i][j]*(fr[i][j] + b[i][j]*fa[i][j]);
			res += V[i][j];
		}
	}

	return res / 2;
}
