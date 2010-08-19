#ifndef _FENERGY_H
#define _FENERGY_H

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
//addtional parameters
#define p_EH alpha[40]
#define p_EO alpha[41]

#define NMAX 40
//#define NMAX 30

#define EPS 1.0e-12
#define TOLX (4*EPS) 
//Convergence criterion on x values.
#define STPMX 3.0 
//Scaled maximum step length allowed in line searches.
#define FMAX(A,B) ((A)>(B)?(A):(B))
#define SQR(x) (x)*(x)
#define ALF 1.0e-5 


#include "../ObjectiveFunctions/ObjectiveFunction.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class FEnergy:public OFObjectiveFunction
{
public:
	FEnergy::FEnergy(const char* boundfile, const char* neureffile, const char* posreffile, const char* negreffile);
	
	~FEnergy() {}	

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

public:
	struct input{
		double energy;
		int nO;
		int nAtoms;
		double **r;
		double **g;
	};


private:
	double oss1best[42];
	double oss2best[42];
public:
	double* getBest();
	void printBestAlpha(double* x, const char* fname);
	void evaluateAndPrint(double *x);
	void setWeight(int pos, double w);
	int getNParam();
	void loadInputs(const char* filename);
	
private:
	int nparam, nfixed;
	double ub[100];
	double lb[100];
	int fixed[100];

	vector<double> weights;
	vector< vector<input> > inputs;
	
	input ref_inp[3];

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

private:
	double rmsEnergy(double* alpha, vector<input>& inputs, bool print);
	void loadReference(const char* filename, input& ref_input);

	double localopt(double **x0,int nAtom,int nO,double *fileAlpha,double acc,int maxStep,double &energy,double &rmsGrad);

	virtual double evaluate_( double const *x );

	void inv(double A[NMAX][NMAX], double Ainv[NMAX][NMAX], int N);
	//double energy(double x[][3], int nAtom, int nO, double alpha[]);
	double energy(double** x, int nAtom, int nO, double alpha[]);
	//void energy_grad(double x[][3], int nAtom, int nO, double alpha[], double& energy, double grad[][3]);
	void energy_grad(double** x, int nAtom, int nO, double alpha[], double& energy, double grad[][3]);	

	double func(double* x, int nAtom, int nO, double alpha[]);

	double dfunc(double* x, int nAtom, int nO, double alpha[], double* grad);

	void lnsrch(double xold[], double fold, double g[], double p[], double x[],double *fret, double stpmax, int *check, int nAtom, int nO, double alpha[]);
	int dfpmin(double p[], double gtol, int iterMax, double *fret,double *rmsGrad, int nAtom, int nO, double alpha[]);

};

#endif
