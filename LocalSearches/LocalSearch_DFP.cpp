/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#include "LocalSearch_DFP.h"

LocalSearch_DFP::LocalSearch_DFP(ObjectiveFunction* objectiveFunction):LocalSearch( objectiveFunction )
{
}

LocalSearch_DFP::~LocalSearch_DFP()
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
double LocalSearch_DFP::search( vector<double>& x )
{
	LocalSearch::search(x);	
	double res = this->dfpmin(x);	
	return res;
}

void LocalSearch_DFP::lnsrch(vector<double> &xold, const double fold, vector<double> &g, vector<double> &p,
	vector<double> &x, double &f, const double stpmax, bool &check)
{
	const double ALF=1.0e-4, TOLX = 4e-12;
	int i;
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,slope,sum,temp,test,tmplam;

	int n=xold.size();
	check=false;
	sum=0.0;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++) p[i] *= stpmax/sum;
	slope=0.0;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) { 
		cout << "Roundoff problem in lnsrch." << endl; 
	}
	test=0.0;
	for (i=0;i<n;i++) {
		temp=fabs(p[i])/max(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		f=evaluate(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=MAX(tmplam,0.1*alam);
	}
}

double LocalSearch_DFP::dfpmin(vector<double>& p)
{
	const int ITMAX = 200;
	const double EPS = 1e-40;
	const double TOLX = 4*EPS, STPMX=100.0;
	bool check;
	int i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test, fret;

	int n = fObj->nDimensions();

	vector<double> dg(n),g(n),hdg(n),pnew(n),xi(n);

	vector< vector<double> > hessin(n, vector<double> (n));

	fp=evaluate(p);
	g = fObj->gradient(p);

	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}

	stpmax=STPMX*MAX(sqrt(sum),double(n));

	for (its=0;its<ITMAX;its++) {
		//iter=its;
		lnsrch(p,fp,g,xi,pnew,fret,stpmax,check);
		
		fp=fret;
		for (i=0;i<n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp= fabs(xi[i]) / max(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return fret;
		for (i=0;i<n;i++) dg[i]=g[i];
		g = fObj->gradient(p);
		test=0.0;
		den=MAX(fret,1.0);
		for (i=0;i<n;i++) {
			temp=fabs(g[i]) * max(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < this->accuracy)
			return fret;
		for (i=0;i< n;i++) dg[i]=g[i]-dg[i];
		for (i=0;i<n;i++) {
			hdg[i]=0.0;
			for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=0;i<n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += (dg[i]*dg[i]);
			sumxi += (xi[i]*xi[i]);
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=0;i<n;i++) {
				for (j=i;j<n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
						-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin[j][i]=hessin[i][j];
				}
			}
		}
		for (i=0;i<n;i++) {
			xi[i]=0.0;
			for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}

	return -1;

}
