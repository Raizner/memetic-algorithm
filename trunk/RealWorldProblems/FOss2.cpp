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

#include "FOss2.h"

FOss2::FOss2(int nAtoms, int nO, double alpha[42]) : ObjectiveFunction(-1, nAtoms * 3, INT_MIN, INT_MAX)
{
	this->nAtoms = nAtoms;
	this->nO = nO;

	for(int i=0; i<42; i++)
	{
		this->alpha[i] = alpha[i];
	}
}

double FOss2::evaluate_( vector<double>& x )
{
	int i, j;
	double _x[100][3];

	int cnt = 0;
	// number of atoms

	RP(i, nAtoms)
	{
		RP(j, 3) _x[i][j] = x[cnt++];
	}

	double res = energy(_x);

	return res;
}

vector<double> FOss2::gradient_( vector<double>& x )
{
	double energy;
	int i, j;	
	double _x[100][3];
    int cnt = 0;
	double _grad[100][3];
	vector<double> grad(x.size());

	// doing type casting here - I hate this	
	RP(i, nAtoms)
	{
		RP(j, 3) _x[i][j] = x[cnt++];
	}	

	energy_grad(_x, energy, _grad);

	cnt = 0;
	RP(i, nAtoms)
    {
        RP(j, 3) grad[cnt++] = _grad[i][j];
    }
	
	return grad;
}


void FOss2::inv(double A[NMAX][NMAX], double Ainv[NMAX][NMAX], int N)
{
	int i, j, k;
	// for inversion
	double b[NMAX][NMAX];
	double scale[NMAX];
	int tidx[NMAX];
	
	//* Matrix b is initialized to the identity matrix

	RP(i, N) RP(j, N) b[i][j] = 0.0;
	RP(i, N)
	{		
		b[i][i] = 1.0;
	}

	//* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
	RP(i, N) {
		tidx[i] = i;			  // Initialize row tidx list
		double scalemax = 0.;
		RP(j, N) 
			if (scalemax < ABS(A[i][j])) scalemax = ABS(A[i][j]);
		
		scale[i] = scalemax;
	 }

	//* Loop over rows k = 0, ..., (N-2)
	int signDet = 1;
	RP(k, N-1) {
		//* Select pivot row from max( |a(j,k)/s(j)| )
		double ratiomax = 0.0;
		int jPivot = k;

		FR(i, k, N-1) 
		{
			double ratio = ABS(A[tidx[i]][k])/scale[tidx[i]];
			if( ratio > ratiomax ) {
				jPivot=i;
				ratiomax = ratio;
			}
		}

		//* Perform pivoting using row tidx list
		int tidxJ = tidx[k];

		if( jPivot != k ) {	          // Pivot
			  tidxJ = tidx[jPivot];
			  tidx[jPivot] = tidx[k];   // Swap tidx jPivot and k
			  tidx[k] = tidxJ;
			  signDet *= -1;			  // Flip sign of determinant
		}

		//* Perform forward elimination
		FR(i, k+1, N-1) {
			  double coeff = A[tidx[i]][k]/A[tidxJ][k];
			  FR(j, k+1, N-1)
					A[tidx[i]][j] -= coeff*A[tidxJ][j];

			  A[tidx[i]][k] = coeff;

			  RP(j, N) 
				b[tidx[i]][j] -= A[tidx[i]][k]*b[tidxJ][j];
		}
	}

	//* Compute determinant as product of diagonal elements
	double determ = signDet;	   // Sign of determinant
	RP(i, N) determ *= A[tidx[i]][i];

	//* Perform backsubstitution
	RP(k, N) {
		Ainv[N-1][k] = b[tidx[N-1]][k]/A[tidx[N-1]][N-1];

		for(i=N-2; i>=0; i--) {
			double sum = b[tidx[i]][k];

			FR(j, i+1, N-1)
				sum -= A[tidx[i]][j]*Ainv[j][k];

			Ainv[i][k] = sum/A[tidx[i]][i];
		}
	}
}

double FOss2::energy(double x[][3])
{
    // all kind of distance: dOO, dHH, dOH
    
	int i, j, k, l;

	RP(i, nAtoms) q[i] = ((i<nO)? -2 : 1);

    RP(i, nAtoms)
    {
        FR(j, i+1, nAtoms-1)
        {
            r2[i][j] = 0;
			RP(k,3) r2[i][j] += SQR(x[i][k] - x[j][k]);
            r[i][j] = sqrt(r2[i][j]);
        }
    }

    // <formula 10>
    double VOO = 0;
    RP(i, nO)
    {
        FR(j, i+1, nO-1)
        {
            VOO += p_o1*exp(-p_o2 * r[i][j]);
            VOO += p_o3*exp(-p_o4 * r[i][j]);
            VOO += p_o5*exp(-p_o6 * SQR(r[i][j] - p_o7));

			//VOO += exp(-30.0*(r[i][j]-1.8));
        }
    }
    // </formula 10>


	double VHH = 0;
	for(i=nO;i<nAtoms;i++)
		for(j=i+1;j<nAtoms;j++); 
				//VHH += exp(-50.0*(r[i][j]-1));

    // <formula 9>
    double h5tmp = pow(1.0-p_h5,2) / (pow(1.0-p_h5,2) + p_h5*p_h5);

    double VOH = 0;

    RP(i, nO)
    {
        FR(j, nO, nAtoms-1)
        {
            VOH += p_h1 * pow(1 - h5tmp*exp(-p_h3*(r[i][j]-p_h2))
            - (1-h5tmp)*exp(-p_h4*(r[i][j]-p_h2)), 2) - p_h1;

			//VOH += exp(-70.0*(r[i][j]-0.3));
        }
    }
    // <formula 9>


    // <formula 11>
    double VHOH = 0;
   
	double r11, r21, r12, r22, t, dt, dt2;
	RP(i, nO)
	{
		FR(j, nO, nAtoms-1)
		{
			FR(k, j+1, nAtoms-1)
			{
                r11 = r[i][j] - p_r0;
                r21 = r[i][k] - p_r0;
				r12 = SQR(r11);
				r22 = SQR(r21);
                
				// theta & shifted theta
				t = acos( (r2[i][j] + r2[i][k] - r2[j][k]) / (2.0*r[i][j]*r[i][k]));
                dt = t - p_theta0;
				dt2 = SQR(dt);
                
				double fcutoff = exp(-(p_m1*(r12+r22) +  p_m2*dt2 + p_m3*(r12+r22)*dt2));

				double temp = p_k1+ p_k2*(r11+r21) + p_k3*dt
							+ p_k4*(r12+r22) + p_k5*r11*r21 + p_k6*dt2 + p_k7*(r11+r21)*dt
							+ p_k8*(r12*r11+r22*r21) + p_k9*(r12*r21+r11*r22) + p_k10*dt*dt2 + p_k11*(r12+r22)*dt
							+ p_k12*r11*r21*dt + p_k13*(r11+r21)*dt2
							+ p_k14*(SQR(r12) + SQR(r22)) + p_k15*r12*r22 + p_k16*SQR(dt2);
	
				VHOH+=temp*fcutoff;
			}
		}
    }
    // </formula 11>

// <formula 1>

	// <part 1: V coulomb>
	double Vq = 0;
	
	RP(i, nAtoms)
	{
		FR(j, i+1, nAtoms-1)
		{
			Vq += q[i]*q[j] / r[i][j];
		}
	}
	// </part 1: V coulomb>

	// <part 2: ScdOH, ScdOO, SddOO>
	// formula 6
	RP(i, nO)
	{
		FR(j, nO, nAtoms-1)
		{
			Scd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_a1*exp(-p_a2*r[i][j])));
			
		}
	}

	// formula 7
	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{
			Scd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_b1*exp(-p_b2*r[i][j])));
			
		}
	}

	// formula 8
	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{
			Sdd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_c1*exp(-p_c2*r[i][j])));			
		}
	}
	// </part 2: ScdOH, ScdOO, SddOO>

	// <part 3: matrix T & D>
	double rij[3];

	RP(i, 3*nO)
		RP(j, 3*nO) {  T[i][j] = 0; D[i][j] = 0; }

	RP(i, 3*nO)
	{
		T[i][i] = 0;
		D[i][i] = 0;
	}

	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{			
			RP(k, 3)
			{
				rij[k] = x[j][k] - x[i][k];
			}
			RP(k, 3)
			{
				RP(l, 3)
				{					
					if (k==l)
						T[3*i+k][3*j+l] = 1 - 3*rij[k]*rij[l] / r2[i][j];
					else
						T[3*i+k][3*j+l] = -3*rij[k]*rij[l] / r2[i][j];

					D[3*i+k][3*j+l]  = T[3*i+k][3*j+l] * Sdd[i][j];

					T[3*j+l][3*i+k] = T[3*i+k][3*j+l];
					D[3*j+l][3*i+k] = D[3*i+k][3*j+l];
				}
			}
		}
	}
	// </part 3: matrix T & D>

	// <part 4: matrix E>
	RP(i, nO)
	{
		double tmp[3];
		
		RP(k,3) tmp[k] = 0;
		
		RP(j, i)
		{
			RP(k, 3)
			{
				rij[k] = x[j][k] - x[i][k];
				tmp[k] += -2 * Scd[j][i] * rij[k];
			}
		}

		FR(j, i+1, nO-1)
		{
			RP(k, 3)
			{
				rij[k] = x[j][k] - x[i][k];
				tmp[k] += -2 * Scd[i][j] * rij[k];
			}
		}

		FR(j, nO, nAtoms-1)
		{
			RP(k, 3)
			{
				rij[k] = x[j][k] - x[i][k];
				//cout << Scd[i][j] << endl;
				tmp[k] += 1 * Scd[i][j] * rij[k];
			}
		}

		RP(j, 3)
		{
			E[3*i+j] = -p_alpha*tmp[j];	
		}
	}
	// </part 4: matrix E>


	// <part5: mu>	
	RP(i, 3*nO)
	{
		RP(j, 3*nO)
		{			
			if (i == j)
				Dtmp[i][j] = 1 + p_alpha*D[i][j];
			else
				Dtmp[i][j] = p_alpha*D[i][j];			
		}
	}

	inv(Dtmp, Dinv, 3*nO);

	RP(i, 3*nO)
	{
		mu[i] = 0;
		RP(j, 3*nO)
		{
			mu[i] += Dinv[i][j] * E[j];
		}
		//cout << mu[i] << " ";
	}
	// cout << endl;
	// </part5: mu>

	// <part 6: 3 Vsd, Vcd, Vdd>
	double Vsd = 0;
	RP(i, nO*3) Vsd += mu[i]*mu[i];
	Vsd /= 2*p_alpha;


	double Vdd = 0;
	RP(i, 3*nO)
	{
		FR(j, i+1, 3*nO-1)
		{
			Vdd += mu[i]*D[i][j]*mu[j];
		}
	}

	double Vcd = 0;
	RP(i, 3*nO)
	{
		Vcd -= mu[i] * E[i] / p_alpha;
	}
	// </part 6: the rest>

	// <part 7: Vel & energy>
	double Vel = 0, Ener = 0;

	Vel = (Vdd + Vcd + Vq + Vsd) * RESCALE;

	Ener = Vel + VOH + VHOH + VOO + VHH;

	// </part 7: Vel>

// </formula 1>

	
// results
	/*
	cout << "VOO: " << VOO << " VOH: " << VOH << " VHOH: " << VHOH << endl;
	cout << "Vq: " << Vq << endl;
	cout << "Self-dipole: " << Vsd << endl;
	cout << "Dipole-dipole interaction: " << Vdd << endl;
	cout << "Charge-dipole interaction: " << Vcd << endl;
	cout << "Vel: " << Vel << endl;
	cout << "Energy: " << Ener << endl;
	*/

	return Ener;
}

void FOss2::energy_grad(double x[][3], double& energy, double grad[100][3])
{
	int i, j, k, l, t;
	double pt1, pt2, pt3, dp, dpk;

    // clear the grad matrix
	RP(i, nAtoms) RP(j, 3) grad[i][j] = 0;

	// charges
	RP(i, nAtoms) q[i] = ((i<nO)? -2 : 1);

	// distance vectors, euclidian distances and their derivations
    RP(i, nAtoms)
    {
        FR(j, i+1, nAtoms-1)
        {
			r2[i][j] = 0;

			RP(k, 3)
			{
				rv[i][j][k] = x[i][k] - x[j][k];
				rv[j][i][k] = -rv[i][j][k];
				r2[i][j] += SQR(rv[i][j][k]);
			}
            r[i][j] = sqrt(r2[i][j]);
			r[j][i] = r[i][j];

			//  dij = SQRT(SQR(ri1 - rj1) + SQR(ri2 - rj2) + SQR(ri3 - rj3))'
			// dij'(rik) = 1/(2d) * (2*(rik-rjk)) = (rik - rjk) / d
			RP(k, 3) 
			{ 
				dr[i][j][k] = rv[i][j][k]/r[i][j]; 
				dr[j][i][k] =-dr[i][j][k]; 
			}
        }
    }

	// <formula 1>
	// <part 1: V coulomb>
	double Vq = 0;
	double temp;	
	
	RP(i, nAtoms)
	{
		FR(j, i+1, nAtoms-1)
		{
			Vq += q[i]*q[j] / r[i][j];
			RP(k,3)
			{
				temp=(-q[i]*q[j]/r2[i][j])*dr[i][j][k];
				grad[i][k]+= temp;
				grad[j][k]+=-temp;
			}
		}
	}
	// </part 1: V coulomb>

	// <part 2: ScdOH, ScdOO, SddOO>
	// formula 6
	RP(i, nO)
	{
		FR(j, i+1, nAtoms-1)
		{
			Scd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_a1*exp(-p_a2*r[i][j])));
			Scd[j][i] = Scd[i][j];
			dScd[i][j]=-SQR(Scd[i][j])*(3.0*r2[i][j] + p_a1*exp(-p_a2*r[i][j])*(1.0-p_a2*r[i][j]));
			dScd[j][i]=dScd[i][j];
		}
	}

	// formula 7
	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{
			Scd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_b1*exp(-p_b2*r[i][j])));
			Scd[j][i] = Scd[i][j];
			dScd[i][j]=-SQR(Scd[i][j])*(3.0*r2[i][j] + p_b1*exp(-p_b2*r[i][j])*(1.0-p_b2*r[i][j]));
			dScd[j][i]=dScd[i][j];
		}
	}

	// formula 8
	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{
			Sdd[i][j] = 1 / (r[i][j]*(r2[i][j] + p_c1*exp(-p_c2*r[i][j])));
			Sdd[j][i] = Sdd[i][j];
			dSdd[i][j]=-SQR(Sdd[i][j])*(3.0*r2[i][j] + p_c1*exp(-p_c2*r[i][j])*(1.0-p_c2*r[i][j]));
			dSdd[j][i]=dSdd[i][j];
		}
	}
	// </part 2: ScdOH, ScdOO, SddOO>

	// <part 3: matrix T & D>
	RP(i, 3*nO)
		RP(j, 3*nO) { T[i][j] = 0; D[i][j] = 0; }

	RP(i, 3*nO)
	{
		T[i][i] = 0;
		D[i][i] = 0;
	}

	RP(i, nO)
	{
		FR(j, i+1, nO-1)
		{
			RP(k, 3)
			{
				RP(l, 3)
				{					
					if (k==l)
						T[3*i+k][3*j+l] = 1 - 3*rv[j][i][k]*rv[j][i][l] / r2[i][j];
					else
						T[3*i+k][3*j+l] = -3*rv[j][i][k]*rv[j][i][l] / r2[i][j];

					D[3*i+k][3*j+l]  = T[3*i+k][3*j+l] * Sdd[i][j];

					T[3*j+l][3*i+k] = T[3*i+k][3*j+l];
					D[3*j+l][3*i+k] = D[3*i+k][3*j+l];
				}
			}
		}
	}
	// </part 3: matrix T & D>

	// <part 4: matrix E>
	RP(i, nO)
	{
		double tmp[3];
		
		RP(k,3) tmp[k] = 0;
		
		RP(j, i)
		{
			RP(k, 3)
			{				
				tmp[k] += -2 * Scd[j][i] * rv[j][i][k];
			}
		}

		FR(j, i+1, nO-1)
		{
			RP(k, 3)
			{				
				tmp[k] += -2 * Scd[i][j] * rv[j][i][k];
			}
		}

		FR(j, nO, nAtoms-1)
		{
			RP(k, 3)
			{
				tmp[k] += 1 * Scd[i][j] * rv[j][i][k];
			}
		}

		RP(j, 3)
		{
			E[3*i+j] = -p_alpha*tmp[j];	
		}
	}
	// </part 4: matrix E>

	// <part5: mu>	
	RP(i, 3*nO)
	{
		RP(j, 3*nO)
		{			
			if (i == j)
				Dtmp[i][j] = 1 + p_alpha*D[i][j];
			else
				Dtmp[i][j] = p_alpha*D[i][j];			
		}
	}

	inv(Dtmp, Dinv, 3*nO);

	RP(i, 3*nO)
	{
		mu[i] = 0;
		RP(j, 3*nO)
		{
			mu[i] += Dinv[i][j] * E[j];
		}
		//cout << mu[i] << " ";
	}
	// cout << endl;
	// </part5: mu>

	// <part 6: 3 Vsd, Vcd, Vdd>
	double Vsd = 0;
	RP(i, nO*3) Vsd += mu[i]*mu[i];
	Vsd /= 2*p_alpha;


	double Vdd = 0;
	RP(i, 3*nO)
	{
		FR(j, i+1, 3*nO-1)
		{
			Vdd += mu[i]*D[i][j]*mu[j];
		}
	}
	
	// gradient Vdd'
	// Vdd = mu[i]*T[ij]*mu[j]*Sdd[ij]
	// Vdd'= mu[i]*T'[ij]*mu[j]*Sdd[ij] + mu[i]*T[ij]*mu[j]*S'dr[ij]

	RP(i, nO)
	{
		RP(j, nO)
		{
			if(i!=j){	
				pt1 = 0;
				RP(k,3) RP(l,3) pt1 += mu[3*i+k]*T[3*i+k][3*j+l]*mu[3*j+l];				

				pt2 = 0;
				pt3 = 0;

				RP(k,3)
				{	
					pt2 += mu[i*3+k] * dr[i][j][k];
					pt3 += mu[j*3+k] * dr[i][j][k];
				}
				
				RP(t, 3)
				{
					temp = pt1*dSdd[i][j]*dr[i][j][t];
					temp += 3.0*(2.0*dr[i][j][t]*pt2*pt3
							- mu[j*3+t]*pt2
							- mu[i*3+t]*pt3)/r[i][j]*Sdd[i][j];
			
					grad[i][t]+=temp;
				}
			}
		}			
	}

	double Vcd = 0;
	RP(i, 3*nO)
	{
		Vcd -= mu[i] * E[i] / p_alpha;
	}

	// gradient Vcd'
	// Vcd = q[i]*mu[i]*rji*Scd[ij]
	//Vcd = q[i]*mu[i]*rji*Scd[ij]

	RP(i, nAtoms)
	{
		RP(j, nO)
		{
			if(i!=j){
				
				RP(t,3)				
				{
					double tmp = 0;
					RP(k,3) tmp += mu[3*j+k] * rv[i][j][k];

					tmp =q[i]*tmp*dScd[i][j]*dr[i][j][t]; 
					tmp+=q[i]*mu[3*j+t]*Scd[i][j];

					grad[i][t]+=tmp;
					grad[j][t]-=tmp;
				}
			}
		}
	}
	
	// </part 6: the rest>


	// <part 7: Vel & energy>
	double Vel = 0, Ener = 0;

	Vel = (Vdd + Vcd + Vq + Vsd) * RESCALE;	
	RP(i, nAtoms) RP(j, 3) grad[i][j] *= RESCALE;
	// </part 7: Vel>

// </formula 1>

    // <formula 9: 2-body interaction OH>
    double h5tmp = SQR(1.0-p_h5) / (SQR(1.0-p_h5) + p_h5*p_h5);

    double VOH = 0;
	

    RP(i, nO)
    {
        FR(j, nO, nAtoms-1)
        {
			pt1 = h5tmp*exp(-p_h3*(r[i][j]-p_h2));
			pt2 = (1-h5tmp)*exp(-p_h4*(r[i][j]-p_h2));

            VOH += p_h1 * SQR(1 - pt1 - pt2) - p_h1;
			VOH += exp(-70.0*(r[i][j]-0.3));
			
			// pt3 = pt1' + pt2'
			pt3= -p_h3 * pt1 - p_h4*pt2;			
			dp = -p_h1 * 2 * (1-pt1-pt2) * pt3;
			dp+=-70.0*exp(-70.0*(r[i][j]-0.3)); //prevent unphysically short O-H
			RP(k, 3)
			{				
				dpk = dp * dr[i][j][k];
				grad[i][k] += dpk;
				grad[j][k] -= dpk;
			}
        }
		
    }
    // <formula 9>

	// <formula 10: 2-body interaction OO>
    double VOO = 0;
	
    RP(i, nO)
    {
        FR(j, i+1, nO-1)
        {
            pt1 = p_o1*exp(-p_o2 * r[i][j]);
            pt2 = p_o3*exp(-p_o4 * r[i][j]);
            pt3 = p_o5*exp(-p_o6 * SQR(r[i][j] - p_o7));

			VOO += pt1 + pt2 + pt3;
			VOO += exp(-30.0*(r[i][j]-1.8));
			

			// dp = pt1' + pt2' + pt3'
			dp = (-p_o2*pt1 - p_o4*pt2 - p_o6*2.0*(r[i][j] - p_o7)*pt3);
			dp+=-30.0*exp(-30.0*(r[i][j]-1.8)); //prevent unphysically short O-O
			
			RP(k, 3)
			{
				dpk = dp * dr[i][j][k];
				grad[i][k] += dpk;
				grad[j][k] -= dpk;
			}
        }
    }
    // </formula 10>

	double VHH = 0;
	
	for(i=nO;i<nAtoms;i++)
		for(j=i+1;j<nAtoms;j++){ //H-H
			VHH += exp(-50.0*(r[i][j]-1));
			temp=-50.0*exp(-50.0*(r[i][j]-1.0));
			RP(t, 3)
			{
				temp*=dr[i][j][t];
				grad[i][t]+=  temp;
				grad[j][t]+= -temp;
			}
		}

    // <formula 11>
    double VHOH = 0;
   
	double r11, r21, r12, r22, cost, theta, dt, dt2, fcutoff;

	double dfcutoff_r1=0,dfcutoff_r2=0,dfcutoff_t=0;

	double dtemp_r1=0,dtemp_r2=0,dtemp_t=0;

	double dVHOH_r1=0,dVHOH_r2=0,dVHOH_t=0;	
	
	double dt_ri,dt_rj,dt_rk;

	RP(i, nO)
	{
		FR(j, nO, nAtoms-1)
		{
			FR(k, j+1, nAtoms-1)
			{
				// shifted distance
                r11 = r[i][j] - p_r0;
                r21 = r[i][k] - p_r0;
				r12 = SQR(r11);
				r22 = SQR(r21);
                
				// theta & shifted theta
				cost = (r2[i][j] + r2[i][k] - r2[j][k]) / (2.0*r[i][j]*r[i][k]);
				theta = acos(cost);
                dt = theta - p_theta0;
				dt2 = SQR(dt);
                
				fcutoff = exp(-(p_m1*(r12+r22) +  p_m2*dt2 + p_m3*(r12+r22)*dt2));

				double temp = p_k1+ p_k2*(r11+r21) + p_k3*dt
							+ p_k4*(r12+r22) + p_k5*r11*r21 + p_k6*dt2 + p_k7*(r11+r21)*dt
							+ p_k8*(r12*r11+r22*r21) + p_k9*(r12*r21+r11*r22) + p_k10*dt*dt2 + p_k11*(r12+r22)*dt
							+ p_k12*r11*r21*dt + p_k13*(r11+r21)*dt2
							+ p_k14*(SQR(r12) + SQR(r22)) + p_k15*r12*r22 + p_k16*SQR(dt2);
	

				VHOH+=temp*fcutoff;


				// the damn derivatives follow :((
				// VHOH' = temp'*fcutoff + fcutoff'*temp

				// fcutoff'
				dfcutoff_r1=-2.0*(p_m1*r11 + p_m3*r11*dt2)*fcutoff;
				dfcutoff_r2=-2.0*(p_m1*r21 + p_m3*r21*dt2)*fcutoff;
				dfcutoff_t =-2.0*(p_m2+p_m3*(r12+r22))*dt*fcutoff;

				// temp'
				dtemp_r1 = p_k2 + p_k4*2.0*r11 + p_k5*r21 + p_k7*dt 
							+ p_k8*3.0*r12 + p_k9*(2.0*r11*r21+r22) + p_k11*2.0*r11*dt
							+ p_k12*r21*dt + p_k13*dt2 +p_k14*4.0*r11*r12 + p_k15*2.0*r11*r22;
				dtemp_r2 = p_k2 + p_k4*2.0*r21 + p_k5*r11 + p_k7*dt 
							+ p_k8*3.0*r22 + p_k9*(2.0*r21*r11+r12) + p_k11*2.0*r21*dt
							+ p_k12*r11*dt + p_k13*dt2 +p_k14*4.0*r21*r22 + p_k15*2.0*r21*r12;
				dtemp_t  = p_k3 + p_k6*2.0*dt + p_k7*(r11+r21) + p_k10*3.0*dt2 + p_k11*(r12+r22)
							+ p_k12*r11*r21 + p_k13*(r11+r21)*2.0*dt + p_k16*4.0*dt*dt2;

				// VHOH'
				dVHOH_r1 = fcutoff*dtemp_r1 + dfcutoff_r1*temp;
				dVHOH_r2 = fcutoff*dtemp_r2 + dfcutoff_r2*temp;
				dVHOH_t  = fcutoff*dtemp_t  + dfcutoff_t *temp;

				// r11', r21' and theta'
				// r11' = dr[i][j] && r21' = dr[i][k]
				double sint = sqrt(1-SQR(cost));
				
				RP(t, 3)
				{
					dt_rj = (-1.0/sint) * (cost*dr[i][j][t] - dr[i][k][t])/r[i][j];
					dt_rk = (-1.0/sint) * (cost*dr[i][k][t] - dr[i][j][t])/r[i][k];
					dt_ri = -(dt_rj + dt_rk);
																		
					grad[i][t]+= dVHOH_r1*dr[i][j][t] + dVHOH_r2*dr[i][k][t] + dVHOH_t * dt_ri;
					grad[j][t]+= dVHOH_r1*dr[j][i][t]  						 + dVHOH_t * dt_rj;
					grad[k][t]+=                        dVHOH_r2*dr[k][i][t] + dVHOH_t * dt_rk;	
				}


			}
		}
    }
    // </formula 11>	

	energy = Vel + VOH + VHOH + VOO + VHH;	
}