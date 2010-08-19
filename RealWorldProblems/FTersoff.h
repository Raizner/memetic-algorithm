/*! Potential model for Ge(n)Si(m) molecule.
Reference

Tersoff, New empirical model for the structural properties of silicon, Phys. Rev. Lett. 56, 632 (1986).

Copyright (c) 2010 by <Quang Huy / NTU>
*/


#ifndef _FTersoff_H
#define _FTersoff_H


#define A_ 0
#define B_ 1
#define lambda_ 2
#define mu_ 3
#define beta_ 4
#define n_ 5
#define c_ 6
#define d_ 7
#define h_ 8
#define R_ 9
#define S_ 10

#define NMAX 40

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class FTersoff:public ObjectiveFunction
{
public:
	FTersoff(int nAtoms, int nSi, double parameters[11][2]);
	~FTersoff() {}	

	// <identifier>
	virtual const string equation()
	{
		return "<See the paper>";
	}
	virtual const string classname()
	{
		return "Tersoff";
	}
	// </identifier>


private:
	// parameters
	double params[11][2];

	// temporary variables
	double _A[NMAX], _B[NMAX], _R[NMAX], _S[NMAX], _lambda[NMAX], _mu[NMAX];
	double c[NMAX], d[NMAX], n[NMAX], h[NMAX], beta[NMAX];
	double A[NMAX][NMAX], B[NMAX][NMAX],  S[NMAX][NMAX], R[NMAX][NMAX], lambda[NMAX][NMAX], mu[NMAX][NMAX], omega[NMAX][NMAX], chi[NMAX][NMAX], zeta[NMAX][NMAX];
	double b[NMAX][NMAX], r[NMAX][NMAX];
	double fr[NMAX][NMAX], fa[NMAX][NMAX], fc[NMAX][NMAX];
	double g[NMAX][NMAX][NMAX], theta[NMAX][NMAX][NMAX];
	double V[NMAX][NMAX];

	// number of atoms and silic
	int nAtoms, nSi;

private:
	virtual double evaluate_( vector<double>& x ); 		
	double energy(double x[][3]);
	void init(double x[][3]);
};

#endif
