/*
Ojamäe, Shavitt and Singer designed a family of potentials describing both 
intramolecular and intermolecular interactions, and allowing dissociation 
and formation of solvated proton (H2O)nH+ complexes in water called the OSSn 
potentials where n = 1 to 3 are fitted to ab initio Moller Plesset (MP2) 
calculations. The potentials are parametrized in the rm of interactions 
between H+ and O2- ions, with additionalree-body (H-O-H) interaction terms and 
self-consistent treatment of the polarizability of the oxygen ions.

For more information, please refer to 

Ojamme, L.; Shavitt, I. & Singer, S. J. (1998), 'Potential models for simulations of the solvated proton in water', The Journal of Chemical Physics 109, 5547.
 
This code implement the OSS2 function for water molecule (protonated, deprotonated and neutral)

Copyright (c) 2010 by <Quang Huy / NTU>
*/

#ifndef _FOss2_H
#define _FOss2_H

//wikipedia scale
//#define RESCALE 0.5291772108
//fortran scale
#define RESCALE 0.529177249

//#define E_H -0.500272784180
//#define E_O -74.9663616401

//definition of parameters
#define p_r0 alpha[0]
#define p_theta0 alpha[1]
#define p_a1 alpha[2]
#define p_a2 alpha[3]
#define p_b1 alpha[4]
#define p_b2 alpha[5]
#define p_c1 alpha[6]
#define p_c2 alpha[7]

#define p_h1 alpha[8]
#define p_h2 alpha[9]
#define p_h3 alpha[10]
#define p_h4 alpha[11]
#define p_h5 alpha[12]

#define p_o1 alpha[13]
#define p_o2 alpha[14]
#define p_o3 alpha[15]
#define p_o4 alpha[16]
#define p_o5 alpha[17]
#define p_o6 alpha[18]
#define p_o7 alpha[19]

#define p_k1 alpha[20]
#define p_k2 alpha[21]
#define p_k3 alpha[22]
#define p_k4 alpha[23]
#define p_k5 alpha[24]
#define p_k6 alpha[25]
#define p_k7 alpha[26]
#define p_k8 alpha[27]
#define p_k9 alpha[28]
#define p_k10 alpha[29]
#define p_k11 alpha[30]
#define p_k12 alpha[31]
#define p_k13 alpha[32]
#define p_k14 alpha[33]
#define p_k15 alpha[34]
#define p_k16 alpha[35]

#define p_m1 alpha[36]
#define p_m2 alpha[37]
#define p_m3 alpha[38]

#define p_alpha alpha[39]

#define p_EH alpha[40]
#define p_EO alpha[41]

#define NMAX 40

#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class FOss2:public ObjectiveFunction
{
public:
	FOss2(int nAtoms, int nO, double alpha[42]);	
	~FOss2() {}	

	// <identifier>
	virtual const string equation()
	{
		return "<See the paper>";
	}
	virtual const string classname()
	{
		return "Oss2";
	}
	// </identifier>

private:
	// parameters
	double alpha[42];

	// charges
	double q[NMAX];

	// distance matrices
	double rv[NMAX][NMAX][3];
	double r[NMAX][NMAX];
	double r2[NMAX][NMAX];

	// distance matrices's derivation
	double dr[NMAX][NMAX][3];

	// field cutoff functions
	double Scd[NMAX][NMAX];
	double dScd[NMAX][NMAX];
	double Sdd[NMAX][NMAX];
	double dSdd[NMAX][NMAX];

	// T matrix
	double T[NMAX][NMAX];

	// D matrix
	double D[NMAX][NMAX];
	double Dtmp[NMAX][NMAX];
	double Dinv[NMAX][NMAX];

	// E vector
	double E[NMAX];

	// mu vector
	double mu[NMAX];

	// number of atoms and number of oxi
	int nAtoms, nO;

private:
	virtual double evaluate_( vector<double>& x ); 
	virtual vector<double> gradient_( vector<double>& x );	

	void inv(double A[NMAX][NMAX], double Ainv[NMAX][NMAX], int N);	
	double energy(double x[][3]);	
	void energy_grad(double x[][3], double& energy, double grad[100][3]);
};

#endif
