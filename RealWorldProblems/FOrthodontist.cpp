/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2010 by <Quang Huy / NTU>
 */

#include "FOrthodontist.h"
#include <limits.h>

FOrthodontist::FOrthodontist(int nDimensions) : ObjectiveFunction(-1, 8, INT_MIN, INT_MAX)
{	
	lowerBounds[0] = lowerBounds[6] = 5e-3; lowerBounds[1] = lowerBounds[2] = lowerBounds[3] = lowerBounds[4] = lowerBounds[5] = 1.5e-3; lowerBounds[7] = 0.5e-3;
	upperBounds[0] = upperBounds[6] = 9e-3; upperBounds[1] = 10e-3; upperBounds[2] = 5e-3; upperBounds[3] = 16e-3; upperBounds[4] = 5e-3; upperBounds[5] = 10e-3; upperBounds[7] = 3e-3; 
}

double FOrthodontist::evaluate_( vector<double>& x )
{
	// Distance of Gap (m)
	double e = 0.5e-3;

	// Total of Length (m)
	double Lt = 20e-3;

	// Total of Height (m)
	double Ht = 15e-3;

	// Offset (m)
	double d = 1e-3;

	// Geometri
	// Radius of T (m)
	double R = x[7]; 

	// Angle of Gable
	double gable = 0.; 	// degree
	double theta = gable*PI/180.; 	// Rad

	// Material : Stainless Steel
	// Modulus Elasticity (N/m^2)
	double E = 2e11;

	// Inersia Momen (m^4)
	double I = (1./12)*(0.0005588)*pow(0.0004064,3);


	// Distance Activation (m) 
	double ux = 1e-3;

	// Force Equation (N)
	double Fx =  
	(12*(x[0]+x[6]+x[1]+x[2]+2*R*PI+x[3]+x[4]+x[5])*E*I/(48*R*R*x[3]*x[1]+12*x[1]*x[1]*x[3]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[4]+4*pow(sin
	(1/180*theta*PI),2)*pow(x[6],3)*x[4]+12*x[5]*x[5]*x[6]*x[3]+24*R*pow(x[0],3)*PI+24*R*R*PI*x[5]*x[5]+48*R*R*x[3]*x[5]+8*pow(x[1],3)*R*PI+12*x[1]*x[1]
	*x[5]*x[6]+12*x[1]*x[1]*x[5]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[3]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[0]-12*x[5]*x[5]*x[1]*x[0]
	+48*R*R*x[3]*x[4]-6*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[1]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[5]+48*x[1]*x[1]*x[2]*x[6]+48*R*R*x[0]*x[0]
	*PI*PI+12*x[1]*x[1]*x[4]*x[0]+48*R*R*x[3]*x[6]+12*x[1]*x[1]*x[2]*x[0]+36*pow(R,3)*PI*x[6]+48*R*R*x[3]*x[0]+48*R*R*x[3]*x[2]+12*x[5]*x[5]*x[6]*x[0]
	+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[3]+48*x[1]*R*x
	[3]*x[0]+96*x[1]*R*x[3]*x[6]+pow(x[1],4)+pow(x[5],4)+24*x[1]*x[1]*R*x[3]-60*x[5]*x[5]*x[1]*x[6]+48*x[1]*x[1]*x[3]*x[6]+48*x[1]*x[1]*x[4]*x[6]+12*x[1]*x[1]*x
	[6]*x[0]+12*x[5]*x[5]*x[6]*x[2]+12*x[5]*x[5]*x[6]*x[4]+4*pow(x[1],3)*x[0]+28*pow(x[1],3)*x[6]+4*pow(x[1],3)*x[2]+4*pow(x[1],3)*x[3]+4*pow(x[1],3)*x[4]+4*pow
	(x[1],3)*x[5]-6*x[1]*x[1]*x[5]*x[5]+4*pow(x[5],3)*x[0]+28*pow(x[5],3)*x[6]+4*pow(x[5],3)*x[1]+4*pow(x[5],3)*x[2]+4*pow(x[5],3)*x[3]+4*pow(x[5],3)*x[4]-24*x
	[5]*x[1]*x[6]*x[0]-48*x[5]*x[1]*x[6]*x[2]-48*x[5]*x[1]*x[6]*x[3]-48*x[5]*x[1]*x[6]*x[4]-48*x[5]*x[6]*R*x[3]+24*R*x[3]*x[5]*x[5]+96*x[1]*x[1]*x[6]
	*R*PI+24*x[5]*x[5]*x[6]*R*PI+8*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*R*PI-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[2]-24*x[0]*x[0]*sin
	(1/180*theta*PI)*R*R*PI-24*R*pow(x[0],3)*PI*pow(cos(1/180*theta*PI),2)-48*R*R*x[0]*x[0]*PI*PI*pow(cos(1/180*theta*PI),2)-12*x[0]*x[0]*sin
	(1/180*theta*PI)*x[1]*x[3]-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[3]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[4]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]
	*x[5]-36*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[6]+12*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*x[6]+36*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[6]-12*x[0]*sin
	(1/180*theta*PI)*x[5]*x[6]*x[6]+72*x[0]*sin(1/180*theta*PI)*x[5]*x[5]*x[6]-16*pow(sin(1/180*theta*PI),2)*pow(x[0],3)*R*PI+48*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[4]*x[6]+24*R*x[0]*x[0]*PI*x[6]+24*R*x[0]*x[0]*PI*x[1]+24*R*x[0]*x[0]*PI*x[2]+24*R*x[0]*x[0]*PI*x[3]+24*R*x[0]*x[0]*PI*x[4]
	+24*R*x[0]*x[0]*PI*x[5]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*x[6]+48*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[5]*x[6]+24*R*x[1]*x[1]*PI*x[0]-72*x[0]*sin(1/180*theta*PI)*x[1]*x[1]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*x[6]
	+48*R*R*x[1]*PI*x[0]+96*R*R*x[1]*PI*x[6]+12*x[6]*x[6]*sin(1/180*theta*PI)
	*x[5]*x[1]+36*pow(R,3)*x[3]*PI+8*pow(x[5],3)*R*PI+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[2]+24*R*R*x[1]*x[1]*PI+28*pow(x[0],3)
	*pow(sin(1/180*theta*PI),2)*x[6]-18*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]*x[6]+6*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[5]-18*x[6]*x[6]*sin
	(1/180*theta*PI)*x[5]*x[5]+18*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[1]+36*pow(R,3)*PI*x[0]+36*pow(R,3)*PI*x[1]+36*pow(R,3)*PI*x[2]+36*pow(R,3)
	*PI*x[4]+36*pow(R,3)*PI*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[2]+24*pow(R,4)*PI*PI+pow
	(sin(1/180*theta*PI),2)*pow(x[6],4)+pow(sin(1/180*theta*PI),2)*pow(x[0],4)-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]*R*PI-24*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[3]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]*R*PI-96*x[5]*x[1]*x
	[6]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*R*PI+72*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]*R*PI-192*x[0]*sin(1/180*theta*PI)*x[1]*x
	[6]*R*PI+96*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*R*PI-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*PI-24*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*R*PI+48*x
	[6]*x[6]*sin(1/180*theta*PI)*x[1]*R*PI-48*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*R*PI-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[2]-12*x[6]*x[6]
	*sin(1/180*theta*PI)*x[5]*x[3]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[4]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[2]+24*x[6]*x[6]*sin
	(1/180*theta*PI)*x[1]*x[3]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[4]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]-24*x[6]*x[6]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[2]-24*pow(x[6],2)*x[0]*pow(sin(1/180*theta*PI),2)*x[3]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]-24*x[6]*x[6]*x[0]*pow
	(sin(1/180*theta*PI),2)*x[5]-48*x[5]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*x[3]-48*R*R*x[0]*x[0]
	*pow(sin(1/180*theta*PI),2)*PI*PI-96*x[0]*sin(1/180*theta*PI)*x[1]*x[2]*x[6]-96*x[0]*sin(1/180*theta*PI)*R*R*PI*x[6]-24*R*x[0]*x[0]*PI*pow(cos
	(1/180*theta*PI),2)*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[1]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[2]-24*R*x[0]*x[0]
	*PI*pow(cos(1/180*theta*PI),2)*x[3]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[4]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[5]-96*x[0]
	*sin(1/180*theta*PI)*x[1]*x[3]*x[6]-96*R*x[0]*sin(1/180*theta*PI)*x[3]*x[6]-96*x[0]*sin(1/180*theta*PI)*x[1]*x[4]*x[6]-48*x[0]*sin(1/180*theta*PI)*x
	[1]*x[5]*x[6]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[2]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[3]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[4])*
	(ux)+12*(-sin(1/180*theta*PI)*x[0]*x[0]+2*x[5]*x[6]-sin(1/180*theta*PI)*x[6]*x[6]+
	2*x[0]*sin(1/180*theta*PI)*x[6]-2*x[1]*x[6]+x[1]*x[1]-2*x[0]*sin(1/180*theta*PI)*x[1]+2*x[1]*x[2]-2*x[0]*sin(1/180*theta*PI)*x[2]
	+4*R*R*PI+4*x[1]*R*PI-4*R*x[0]*sin(1/180*theta*PI)*PI+2*x[1]*x[3]-2*x[0]*sin(1/180*theta*PI)*x[3]+4*R*x[3]+2*x[1]*x[4]-2*x[0]*sin
	(1/180*theta*PI)*x[4]-2*x[0]*sin(1/180*theta*PI)*x[5]-x[5]*x[5]+2*x[5]*x[1])*E*I/(48*R*R*x[3]*x[1]+12*x[1]*x[1]*x[3]*x[0]+4*pow(x[0],3)*pow(sin
	(1/180*theta*PI),2)*x[4]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[4]+12*x[5]*x[5]*x[6]*x[3]+24*R*pow(x[0],3)*PI+24*R*R*PI*x[5]*x[5]+48*R*R*x[3]*x
	[5]+8*pow(x[1],3)*R*PI+12*x[1]*x[1]*x[5]*x[6]+12*x[1]*x[1]*x[5]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[3]+4*pow(sin(1/180*theta*PI),2)*pow(x
	[6],3)*x[0]-12*x[5]*x[5]*x[1]*x[0]+48*R*R*x[3]*x[4]-6*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[1]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[5]+48*x[1]*x[1]
	*x[2]*x[6]+48*R*R*x[0]*x[0]*PI*PI+12*x[1]*x[1]*x[4]*x[0]+48*R*R*x[3]*x[6]+12*x[1]*x[1]*x[2]*x[0]+36*pow(R,3)*PI*x[6]+48*R*R*x[3]*x[0]+48*R*R*x[3]*x[2]
	+12*x[5]*x[5]*x[6]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x
	[6],3)*x[3]+48*x[1]*R*x[3]*x[0]+96*x[1]*R*x[3]*x[6]+pow(x[1],4)+pow(x[5],4)+24*x[1]*x[1]*R*x[3]-60*x[5]*x[5]*x[1]*x[6]+48*x[1]*x[1]*x[3]*x[6]+48*x[1]*x[1]*x
	[4]*x[6]+12*x[1]*x[1]*x[6]*x[0]+12*x[5]*x[5]*x[6]*x[2]+12*x[5]*x[5]*x[6]*x[4]+4*pow(x[1],3)*x[0]+28*pow(x[1],3)*x[6]+4*pow(x[1],3)*x[2]+4*pow(x[1],3)*x[3]
	+4*pow(x[1],3)*x[4]+4*pow(x[1],3)*x[5]-6*x[1]*x[1]*x[5]*x[5]+4*pow(x[5],3)*x[0]+28*pow(x[5],3)*x[6]+4*pow(x[5],3)*x[1]+4*pow(x[5],3)*x[2]+4*pow(x[5],3)*x[3]
	+4*pow(x[5],3)*x[4]-24*x[5]*x[1]*x[6]*x[0]-48*x[5]*x[1]*x[6]*x[2]-48*x[5]*x[1]*x[6]*x[3]-48*x[5]*x[1]*x[6]*x[4]-48*x[5]*x[6]*R*x[3]+24*R*x[3]*x[5]*x[5]+96*x
	[1]*x[1]*x[6]*R*PI+24*x[5]*x[5]*x[6]*R*PI+8*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*R*PI-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[2]-24*x[0]*x[0]
	*sin(1/180*theta*PI)*R*R*PI-24*R*pow(x[0],3)*PI*pow(cos(1/180*theta*PI),2)-48*R*R*x[0]*x[0]*PI*PI*pow(cos(1/180*theta*PI),2)-12*x[0]*x[0]*sin
	(1/180*theta*PI)*x[1]*x[3]-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[3]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[4]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]
	*x[5]-36*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[6]+12*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*x[6]+36*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[6]-12*x[0]*sin
	(1/180*theta*PI)*x[5]*x[6]*x[6]+72*x[0]*sin(1/180*theta*PI)*x[5]*x[5]*x[6]-16*pow(sin(1/180*theta*PI),2)*pow(x[0],3)*R*PI+48*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[4]*x[6]+24*R*x[0]*x[0]*PI*x[6]+24*R*x[0]*x[0]*PI*x[1]+24*R*x[0]*x[0]*PI*x[2]+24*R*x[0]*x[0]*PI*x[3]+24*R*x[0]*x[0]*PI*x[4]
	+24*R*x[0]*x[0]*PI*x[5]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*x[6]+48*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[5]*x[6]+24*R*x[1]*x[1]*PI*x[0]-72*x[0]*sin(1/180*theta*PI)*x[1]*x[1]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*x[6]
	+48*R*R*x[1]*PI*x[0]+96*R*R*x[1]*PI*x[6]+12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[1]+36*pow(R,3)*x[3]*PI+8*pow(x[5],3)*R*PI+4*pow(x[0],3)*pow(sin
	(1/180*theta*PI),2)*x[2]+24*R*R*x[1]*x[1]*PI+28*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[6]-18*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]*x[6]+6*x[0]
	*x[0]*sin(1/180*theta*PI)*x[5]*x[5]-18*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[5]+18*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[1]+36*pow(R,3)*PI*x[0]+36*pow
	(R,3)*PI*x[1]+36*pow(R,3)*PI*x[2]+36*pow(R,3)*PI*x[4]+36*pow(R,3)*PI*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[1]+4*pow(sin
	(1/180*theta*PI),2)*pow(x[6],3)*x[2]+24*pow(R,4)*PI*PI+pow(sin(1/180*theta*PI),2)*pow(x[6],4)+pow(sin(1/180*theta*PI),2)*pow(x[0],4)-24*x[0]*x[0]
	*pow(sin(1/180*theta*PI),2)*x[4]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*R*PI-24*x
	[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]*R*PI-96*x[5]*x[1]*x[6]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*R*PI+72*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[6]*R*PI-192*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*R*PI+96*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*R*PI-24*R*pow(x[0],2)*sin
	(1/180*theta*PI)*x[1]*PI-24*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*R*PI+48*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*R*PI-48*x[6]*x[6]*x[0]*pow(sin
	(1/180*theta*PI),2)*R*PI-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[2]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[3]-12*x[6]*x[6]*sin(1/180*theta*PI)*x
	[5]*x[4]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[2]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[3]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[4]-24*x[6]*x[6]
	*x[0]*pow(sin(1/180*theta*PI),2)*x[1]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]-24*x[6]*x[6]
	*x[0]*pow(sin(1/180*theta*PI),2)*x[4]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]-48*x[5]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]
	*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*x[3]-48*R*R*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*PI*PI-96*x[0]*sin(1/180*theta*PI)*x[1]*x[2]*x[6]-
	96*x[0]*sin(1/180*theta*PI)*R*R*PI*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[1]-
	24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[2]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[3]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)
	*x[4]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[5]-96*x[0]*sin(1/180*theta*PI)*x[1]*x[3]*x[6]-96*R*x[0]*sin(1/180*theta*PI)*x[3]*x[6]-
	96*x[0]*sin(1/180*theta*PI)*x[1]*x[4]*x[6]-48*x[0]*sin(1/180*theta*PI)*x[1]*x[5]*x[6]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[2]+48*x
	[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[3]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[4])*(gable*PI/180));



	// Moment Equation (Nm)
	double Mz = 6*(-sin(1/180*theta*PI)*x[0]*x[0]+2*x[5]*x[6]-sin(1/180*theta*PI)*x[6]*x[6]+2*x[0]*sin(1/180*theta*PI)*x[6]-2*x[1]*x[6]+x[1]*x[1]-2*x[0]*sin(1/180*theta*PI)*x[1]+2*x[1]*x[2]-2*x[0]*sin(1/180*theta*PI)*x[2]+4*R*R*PI+4*x[1]*R*PI-4*R*x[0]*sin(1/180*theta*PI)*PI+2*x[1]*x[3]-2*x[0]*sin(1/180*theta*PI)*x[3]+4*R*x[3]+2*x[1]*x[4]-2*x[0]*sin(1/180*theta*PI)*x[4]-2*x[0]*sin(1/180*theta*PI)*x[5]-x[5]*x[5]+2*x[5]*x[1])*E*I/(48*R*R*x[3]*x[1]+12*x[1]*x[1]*x[3]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[4]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[4]+12*x[5]*x[5]*x[6]*x[3]
	+24*R*pow(x[0],3)*PI+24*R*R*PI*x[5]*x[5]+48*R*R*x[3]*x[5]+8*pow(x[1],3)*R*PI+12*x[1]*x[1]*x[5]*x[6]+12*x[1]*x[1]*x[5]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[3]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[0]-12*x[5]*x[5]*x[1]*x[0]+48*R*R*x[3]*x[4]-6*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[1]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[5]+48*x[1]*x[1]*x[2]*x[6]+48*R*R*x[0]*x[0]*PI*PI+12*x[1]*x[1]*x[4]*x[0]+48*R*R*x[3]*x[6]+12*x[1]*x[1]*x[2]*x[0]+36*R*R*R*PI*x[6]+48*R*R*x[3]*x[0]+48*R*R*x[3]*x[2]+12*x[5]*x[5]*x[6]*x[0]+4*x[0]*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[3]+48*x[1]*R*x[3]*x[0]+96*x[1]*R*x[3]*x[6]+pow(x[1],4)+pow(x[5],4)+24*x[1]*x[1]*R*x[3]-60*x[5]*x[5]*x[1]*x[6]+48*x[1]*x[1]*x[3]*x[6]+48*x[1]*x[1]*x[4]*x[6]+12*x[1]*x[1]*x[6]*x[0]+12*x[5]*x[5]*x[6]*x[2]+12*x[5]*x[5]*x[6]*x[4]+4*pow(x[1],3)*x[0]+28*pow(x[1],3)*x[6]+4*pow(x[1],3)*x[2]+4*pow(x[1],3)*x[3]+4*pow(x[1],3)*x[4]+4*pow(x[1],3)*x[5]-6*x[1]*x[1]*x[5]*x[5]+4*pow(x[5],3)*x[0]+28*pow(x[5],3)*x[6]+4*pow(x[5],3)*x[1]+4*pow(x[5],3)*x[2]+4*pow(x[5],3)*x[3]+4*pow(x[5],3)*x[4]-24*x[5]*x[1]*x[6]*x[0]-48*x[5]*x[1]*x[6]*x[2]-48*x[5]*x[1]
	*x[6]*x[3]-48*x[5]*x[1]*x[6]*x[4]-48*x[5]*x[6]*R*x[3]+24*R*x[3]*pow(x[5],2)+96*x[1]*x[1]*x[6]*R*PI+24*x[5]*x[5]*x[6]*R*PI+8*pow(sin(1/180*theta*PI),2)
	*pow(x[6],3)*R*PI-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[2]-24*x[0]*x[0]*sin(1/180*theta*PI)*R*R*PI-24*R*pow(x[0],3)*PI*pow(cos(1/180*theta*PI),2)-48*R*R*x[0]*x[0]*PI*PI*pow(cos(1/180*theta*PI),2)-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[3]-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[3]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[4]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[5]-36*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[6]+12*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*x[6]+36*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[6]-12*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[6]+72*x[0]*sin(1/180*theta*PI)*x[5]*x[5]*x[6]-16*pow(sin(1/180*theta*PI),2)*pow(x[0],3)*R*PI+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]*x[6]+24*R*x[0]*x[0]*PI*x[6]+24*R*x[0]*x[0]*PI*x[1]+24*R*x[0]*x[0]*PI*x[2]+24*R*x[0]*x[0]*PI*x[3]+24*R*x[0]*x[0]*PI*x[4]+24*R*x[0]*x[0]*PI*x[5]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]*x[6]+24*R*pow(x[1],2)
	*PI*x[0]-72*x[0]*sin(1/180*theta*PI)*pow(x[1],2)*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*x[6]+48*R*R*x[1]*PI*x[0]+96*R*R*x[1]*PI*x[6]
	+12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[1]+36*pow(R,3)*x[3]*PI+8*pow(x[5],3)*R*PI+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[2]+24*R*R*x[1]*x[1]
	*PI+28*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[6]-18*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]*x[6]+6*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[5]-18*x
	[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[5]+18*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[1]+36*pow(R,3)*PI*x[0]+36*pow(R,3)*PI*x[1]+36*pow(R,3)*PI*x[2]
	+36*pow(R,3)*PI*x[4]+36*pow(R,3)*PI*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[2]+24*pow(R,4)*pow
	(PI,2)+pow(sin(1/180*theta*PI),2)*pow(x[6],4)+pow(sin(1/180*theta*PI),2)*pow(x[0],4)-24*pow(x[0],2)*pow(sin(1/180*theta*PI),2)*x[4]*R*PI-24*x[0]*x
	[0]*pow(sin(1/180*theta*PI),2)*x[3]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]*R*PI-
	96*x[5]*x[1]*x[6]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*R*PI+72*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]*R*PI-192*x[0]*sin
	(1/180*theta*PI)*x[1]*x[6]*R*PI+96*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*R*PI-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*PI-24*x[6]*x[6]*sin
	(1/180*theta*PI)*x[5]*R*PI+48*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*R*PI-48*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*R*PI-12*x[6]*x[6]*sin
	(1/180*theta*PI)*x[5]*x[2]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[3]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[4]+24*x[6]*x[6]*sin(1/180*theta*PI)*x
	[1]*x[2]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[3]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[4]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]-24*x
	[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]-24*x
	[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]-48*x[5]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*x[3]-
	48*R*R*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*PI*PI-96*x[0]*sin(1/180*theta*PI)*x[1]*x[2]*x[6]-96*x[0]*sin(1/180*theta*PI)*R*R*PI*x[6]-24*R*x[0]*x
	[0]*PI*pow(cos(1/180*theta*PI),2)*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[1]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[2]-
	24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[3]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[4]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)
	*x[5]-96*x[0]*sin(1/180*theta*PI)*x[1]*x[3]*x[6]-96*R*x[0]*sin(1/180*theta*PI)*x[3]*x[6]-96*x[0]*sin(1/180*theta*PI)*x[1]*x[4]*x[6]-48*x[0]*sin
	(1/180*theta*PI)*x[1]*x[5]*x[6]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[2]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[3]+48*x[0]*sin(1/180*theta*PI)*x
	[5]*x[6]*x[4])*ux+8*(3*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]+pow(sin(1/180*theta*PI),2)*pow(x[0],3)+3*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x
	[3]-3*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)+3*x[6]*x[6]*sin(1/180*theta*PI)*x[1]-3*x[6]*x[6]*sin(1/180*theta*PI)*x[5]+3*x[0]*sin(1/180*theta*PI)
	*x[5]*x[5]+3*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[6]+6*R*x[1]*x[1]*PI+3*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]+12*R*R*x[1]*PI+pow(x[1],3)-3*x
	[0]*sin(1/180*theta*PI)*x[1]*x[1]+3*x[1]*x[1]*x[2]-6*x[0]*sin(1/180*theta*PI)*x[1]*x[2]+3*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]-12*x[0]*sin
	(1/180*theta*PI)*R*R*PI+6*R*x[0]*x[0]*PI-6*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)+9*pow(R,3)*PI-12*R*x[0]*sin(1/180*theta*PI)*x[1]*PI-6*x
	[0]*sin(1/180*theta*PI)*x[1]*x[3]+3*x[1]*x[1]*x[3]+12*x[1]*R*x[3]-12*R*x[0]*sin(1/180*theta*PI)*x[3]+12*R*R*x[3]+3*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[4]+3*x[1]*x[1]*x[4]-6*x[0]*sin(1/180*theta*PI)*x[1]*x[4]-6*x[0]*sin(1/180*theta*PI)*x[1]*x[5]+3*x[1]*x[1]*x[5]+pow(x[5],3)+pow
	(sin(1/180*theta*PI),2)*pow(x[6],3)-6*x[0]*sin(1/180*theta*PI)*x[1]*x[6]+3*x[1]*x[1]*x[6]-6*x[5]*x[1]*x[6]+3*x[5]*x[5]*x[6]+6*x[0]*sin(1/180*theta*PI)
	*x[5]*x[6]-3*x[5]*x[5]*x[1])*E*I/(48*R*R*x[3]*x[1]+12*x[1]*x[1]*x[3]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[4]+4*pow(sin(1/180*theta*PI),2)*pow
	(x[6],3)*x[4]+12*x[5]*x[5]*x[6]*x[3]+24*R*pow(x[0],3)*PI+24*R*R*PI*x[5]*x[5]+48*R*R*x[3]*x[5]+8*pow(x[1],3)*R*PI+12*x[1]*x[1]*x[5]*x[6]+12*x[1]*x[1]*x
	[5]*x[0]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[3]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[0]-12*x[5]*x[5]*x[1]*x[0]+48*R*R*x[3]*x[4]-6*x[0]*x[0]
	*sin(1/180*theta*PI)*x[1]*x[1]+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[5]+48*pow(x[1],2)*x[2]*x[6]+48*R*R*x[0]*x[0]*PI*PI+12*x[1]*x[1]*x[4]*x[0]
	+48*R*R*x[3]*x[6]+12*x[1]*x[1]*x[2]*x[0]+36*pow(R,3)*PI*x[6]+48*R*R*x[3]*x[0]+48*R*R*x[3]*x[2]+12*x[5]*x[5]*x[6]*x[0]+4*pow(x[0],3)*pow(sin
	(1/180*theta*PI),2)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[3]+48*x[1]*R*x[3]*x[0]+96*x[1]*R*x
	[3]*x[6]+pow(x[1],4)+pow(x[5],4)+24*pow(x[1],2)*R*x[3]-60*pow(x[5],2)*x[1]*x[6]+48*pow(x[1],2)*x[3]*x[6]+48*pow(x[1],2)*x[4]*x[6]+12*pow(x[1],2)*x[6]*x[0]
	+12*pow(x[5],2)*x[6]*x[2]+12*pow(x[5],2)*x[6]*x[4]+4*pow(x[1],3)*x[0]+28*pow(x[1],3)*x[6]+4*pow(x[1],3)*x[2]+4*pow(x[1],3)*x[3]+4*pow(x[1],3)*x[4]+4*x[1]*x
	[1]*x[5]-6*x[1]*x[1]*x[5]*x[5]+4*pow(x[5],3)*x[0]+28*pow(x[5],3)*x[6]+4*pow(x[5],3)*x[1]+4*pow(x[5],3)*x[2]+4*pow(x[5],3)*x[3]+4*pow(x[5],3)*x[4]-24*x[5]*x
	[1]*x[6]*x[0]-48*x[5]*x[1]*x[6]*x[2]-48*x[5]*x[1]*x[6]*x[3]-48*x[5]*x[1]*x[6]*x[4]-48*x[5]*x[6]*R*x[3]+24*R*x[3]*pow(x[5],2)+96*x[1]*x[1]*x[6]*R*PI+24*x
	[5]*x[5]*x[6]*R*PI+8*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*R*PI-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[2]-24*x[0]*x[0]*sin(1/180*theta*PI)
	*R*R*PI-24*R*pow(x[0],3)*PI*pow(cos(1/180*theta*PI),2)-48*R*R*x[0]*x[0]*PI*PI*pow(cos(1/180*theta*PI),2)-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]
	*x[3]-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[3]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[4]-12*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*x[5]-36*x[0]*x[0]*sin
	(1/180*theta*PI)*x[1]*x[6]+12*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*x[6]+36*x[0]*x[0]*sin(1/180*theta*PI)*x[5]*x[6]-12*x[0]*sin(1/180*theta*PI)*x[5]*x
	[6]*x[6]+72*x[0]*sin(1/180*theta*PI)*x[5]*x[5]*x[6]-16*pow(sin(1/180*theta*PI),2)*pow(x[0],3)*R*PI+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]*x[6]
	+24*R*x[0]*x[0]*PI*x[6]+24*R*x[0]*x[0]*PI*x[1]+24*R*x[0]*x[0]*PI*x[2]+24*R*x[0]*x[0]*PI*x[3]+24*R*x[0]*x[0]*PI*x[4]+24*R*x[0]*x[0]*PI*x[5]+48*x
	[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*x[6]+48*pow(x[0],2)*pow(sin(1/180*theta*PI),2)*x[5]*x[6]
	+24*R*x[1]*x[1]*PI*x[0]-72*x[0]*sin(1/180*theta*PI)*x[1]*x[1]*x[6]+48*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*x[6]+48*R*R*x[1]*PI*x[0]+96*R*R*x[1]
	*PI*x[6]+12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[1]+36*pow(R,3)*x[3]*PI+8*pow(x[5],3)*R*PI+4*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[2]+24*R*R*x
	[1]*x[1]*PI+28*pow(x[0],3)*pow(sin(1/180*theta*PI),2)*x[6]-18*pow(x[0],2)*pow(sin(1/180*theta*PI),2)*pow(x[6],2)+6*pow(x[0],2)*sin(1/180*theta*PI)
	*pow(x[5],2)-18*pow(x[6],2)*sin(1/180*theta*PI)*x[5]*x[5]+18*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[1]+36*R*R*R*PI*x[0]+36*R*R*R*PI*x[1]
	+36*R*R*R*PI*x[2]+36*R*R*R*PI*x[4]+36*R*R*R*PI*x[5]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[1]+4*pow(sin(1/180*theta*PI),2)*pow(x[6],3)*x[2]
	+24*pow(R,4)*PI*PI+pow(sin(1/180*theta*PI),2)*pow(x[6],4)+pow(sin(1/180*theta*PI),2)*pow(x[0],4)-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[4]
	*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[1]*R*PI-24*x[0]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[5]*R*PI-96*x[5]*x[1]*x[6]*R*PI-24*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]*R*PI+72*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*x
	[6]*R*PI-192*x[0]*sin(1/180*theta*PI)*x[1]*x[6]*R*PI+96*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*R*PI-24*R*x[0]*x[0]*sin(1/180*theta*PI)*x[1]*PI-
	24*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*R*PI+48*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*R*PI-48*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*R*PI-12*x[6]
	*x[6]*sin(1/180*theta*PI)*x[5]*x[2]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[3]-12*x[6]*x[6]*sin(1/180*theta*PI)*x[5]*x[4]+24*x[6]*x[6]*sin
	(1/180*theta*PI)*x[1]*x[2]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[3]+24*x[6]*x[6]*sin(1/180*theta*PI)*x[1]*x[4]-24*x[6]*x[6]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[1]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[2]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[3]-24*x[6]*x[6]*x[0]*pow(sin
	(1/180*theta*PI),2)*x[4]-24*x[6]*x[6]*x[0]*pow(sin(1/180*theta*PI),2)*x[5]-48*x[5]*x[6]*R*R*PI+24*sin(1/180*theta*PI)*x[6]*x[6]*R*R*PI+24*sin
	(1/180*theta*PI)*x[6]*x[6]*R*x[3]-48*R*R*x[0]*x[0]*pow(sin(1/180*theta*PI),2)*PI*PI-96*x[0]*sin(1/180*theta*PI)*x[1]*x[2]*x[6]-96*x[0]*sin
	(1/180*theta*PI)*R*R*PI*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[6]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[1]-24*R*x[0]*x[0]
	*PI*pow(cos(1/180*theta*PI),2)*x[2]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[3]-24*R*x[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[4]-24*R*x
	[0]*x[0]*PI*pow(cos(1/180*theta*PI),2)*x[5]-96*x[0]*sin(1/180*theta*PI)*x[1]*x[3]*x[6]-96*R*x[0]*sin(1/180*theta*PI)*x[3]*x[6]-96*x[0]*sin
	(1/180*theta*PI)*x[1]*x[4]*x[6]-48*x[0]*sin(1/180*theta*PI)*x[1]*x[5]*x[6]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[2]+48*x[0]*sin(1/180*theta*PI)*x
	[5]*x[6]*x[3]+48*x[0]*sin(1/180*theta*PI)*x[5]*x[6]*x[4])*gable*PI/180;	



	// Moment to Force Ratio Equation (mm)
	double Rt = Mz/Fx * 1000;

	// Moment to force Ratio by orthodontist (mm)
	double Ra = 3;

	// Objective Function (dimensionless)
	double res = pow(((Rt/Ra)-1),2);

	return res;
}

int FOrthodontist::isInfeasible( vector<double> &x )
{
	double h[3] = {0};
	double g[2] = {0};

	// Distance of Gap (m)
	double e = 0.5e-3;

	// Total of Length (m)
	double Lt = 20e-3;

	// Total of Height (m)
	double Ht = 15e-3;

	// Offset (m)
	double d = 1e-3;

	// Equality Contraints
	h[0] = x[2]- x[4]; 		//(m)
	h[1] = x[3]- (x[2]+x[4]+e) ; 	//(m)
	h[2] = x[1]- (x[5]-d); 		//(m)
 

	// Inequality Constraints
	g[0] = x[0]+x[6]+e-Lt; //(m)
	g[1] = x[5]+2*x[7]-Ht; //(m)

	int violations = 0;

    for(int i=0; i<2; i++) {
		if(g[i]>EPS) violations++; break;
    }
    
	for(int i=0; i<3; i++) {
		if(fabs(h[i])>EPS) violations++; break;	
    }

	return violations;
}