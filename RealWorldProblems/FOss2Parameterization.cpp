#include "FOss2Parameterization.h"


FOss2Parameterization::FOss2Parameterization(const char* boundfile, const char* neureffile, const char* posreffile, const char* negreffile) : ObjectiveFunction(-1, 42, INT_MIN, INT_MAX)
{

	cout << "==========================================\n";
	cout << "Loading bounds..." << endl;
	ifstream f(boundfile);

	nparam=0;
	int i, j, k;
	while(f >> oss2best[nparam] >> lb[nparam] >> ub[nparam] >> fixed[nparam])
	{
		cout << "  " << oss2best[nparam] << " " << lb[nparam] << " " << ub[nparam] << " " <<  fixed[nparam] << endl;
		nparam++;
	}

	nfixed = 0;
	for(i=0; i<nparam; i++) nfixed+=fixed[i];

	ndim = nparam - nfixed;
	
	int t = 0;

	this->lbound = new double[ ndim ];
	this->ubound = new double[ ndim ];

	for(i=0; i<nparam; i++)
	{
		if (!fixed[i])
		{
			lbound[t] = lb[i] - fabs(ub[i]-lb[i]) * 0.01;
			ubound[t] = ub[i] + fabs(ub[i]-lb[i]) * 0.01;			
			t++;
		}
	}
	
	f.close();

	cout << nparam << " variables, " << nfixed << " fixed values...";
	cout << "done." << endl;

    	cout << "==========================================\n";

	cout << "Initialize...";
	initialize(-1, ndim);	
	cout << "done." << endl;

    	cout << "==========================================\n";

	loadReference(neureffile, ref_inp[0]);
	loadReference(posreffile, ref_inp[1]);
	loadReference(negreffile, ref_inp[2]);

	cout << "==========================================\n";
}

double FOss2Parameterization::rmsEnergy(double* alpha, vector<input>& inputs, bool print)
{
	double rmsEner=0;
    double ener_shift[3];

	
	//	cout << "Alpha:" << endl;
		int j;
	//	RP(j, 40) cout << alpha[j] << endl;

	for(int i=0; i<3; i++) 
	{
	
		/*
		// optimize the reference point
		double tmp[30];
		int ccc = 0;

		for(int ii=0; ii<ref_inp[i].nAtoms; ii++)
		for(int jj=0; jj<3; jj++) tmp[ccc++] = ref_inp[i].r[ii][jj];

		// dfpmin(double p[], double gtol, int iterMax, double *fret,double *rmsGrad, int nAtom, int nO, double alpha[])

		double fret;
		double rmsGrad;

		

		cout << "Optimizing " << ref_inp[i].nAtoms << " atoms" << endl;
		cout << "Initial energy: " <<  energy(ref_inp[i].r, ref_inp[i].nAtoms, ref_inp[i].nO, alpha) << endl;

		int nsteps = dfpmin(tmp, 4e-9, 100, &fret, &rmsGrad, ref_inp[i].nAtoms, ref_inp[i].nO, alpha);

		cout << "Energy: " << ref_inp[i].energy << " " << fret << " - Gradient: " << rmsGrad << " " << nsteps << endl;

		*/



		double rmsGrad,old;
		old=localopt(ref_inp[i].r, ref_inp[i].nAtoms, ref_inp[i].nO, alpha, 1e-6, 1000, ener_shift[i], rmsGrad);
		if(isnan(ener_shift[i]) || rmsGrad>10.0 ) 
		{ 
			// holy cow!!!
			//ener_shift[i]+=99999;		
			return 9999999; // big huh?
		}
		//cout << ref_inp[i].nAtoms << " atoms: " << old << " " << ener_shift[i] << " " << rmsGrad << endl;

	//	ener_shift[i] = fret; // energy(ref_inp[i].r, ref_inp[i].nAtoms, ref_inp[i].nO, alpha);
	}

        //if (print)  cout << "Energy shift values: " << ener_shift[0] << " " << ener_shift[1] <<  endl;
	//cout << endl;
	if (print) cout << "@$ inputs#    nAtoms    Oss2Energy    Org.Energy    Oss2BindingEnergy    Org.BindingEnergy" << endl;
        for (int i=0; i<inputs.size(); i++)
        {
                double ener = energy(inputs[i].r, inputs[i].nAtoms, inputs[i].nO, alpha);
		double rmsEtmp, bindingE1, bindingE2;

		for(int j=0; j<3; j++)
		{
			if (inputs[i].nAtoms - 3*inputs[i].nO == ref_inp[j].nAtoms - 3*ref_inp[j].nO)
			{
		                bindingE1 = ener - (inputs[i].nO-1)*ener_shift[0] - ener_shift[j];
        		        bindingE2 = inputs[i].energy - (inputs[i].nO-1)*ref_inp[0].energy - ref_inp[j].energy;

                		rmsEtmp = pow(bindingE1 - bindingE2, 2);

		                if (print) printf("@$ %5d     %5d     %10.6lf     %10.6lf    %16.10lf     %16.10lf\n", i, inputs[i].nAtoms, ener, inputs[i].energy, bindingE1, bindingE2);

                		rmsEner += rmsEtmp/inputs.size();
        
				break;
			}
		}
	}


	return sqrt(rmsEner);

}

void FOss2Parameterization::evaluateAndPrint(double* x)
{
	double alpha[100];

        int i, j, k;
        int c = 0;

	cout << "==========================================\n";

        cout << "@# Alpha: ";
        for(i=0; i<nparam; i++)
        {
                if (!fixed[i])
                {
                        alpha[i] = x[c];
                        c++;
                }
                else
                {
                        alpha[i] = oss2best[i];
                }
                printf("%.10lf ", alpha[i]);
        }
        cout << endl;

	double totalrms = 0;

	for(i=0; i<inputs.size(); i++)
	{
		double tmp = rmsEnergy(alpha, inputs[i], true);
		totalrms += weights[i] * tmp;
		printf("@# Dataset%d:weight*rms=%.6lf*%.6lf = %.6lf\n", i, weights[i], tmp, weights[i]*tmp );
	}

        printf("@# Total rms: %.8lf\n", totalrms);
	cout << "==========================================\n";

}

int FOss2Parameterization::getNParam()
{
	return nparam;
}

void FOss2Parameterization::loadReference(const char* filename, input& ref_input)
{
	if (strlen(filename) == 0) return;

	cout << "Loading reference from file " << filename << endl;
        ifstream fi(filename);

        input tmp;
	
        fi >> tmp.nAtoms;
        fi >> tmp.energy;
        tmp.nO = 1;
        cout << "Natoms: " << tmp.nAtoms << " - Energy: " << tmp.energy << endl;
        cout << "Atoms: " << endl;
        tmp.r = new double* [tmp.nAtoms];

	char ch;
        for(int i=0; i<tmp.nAtoms; i++)
        {
                tmp.r[i] = new double[3];
        	fi >> ch >> tmp.r[i][0] >> tmp.r[i][1] >> tmp.r[i][2];
                cout << ch << " " << tmp.r[i][0] << " " << tmp.r[i][1] << " " << tmp.r[i][2] << endl;
        }

	ref_input = tmp;

        fi.close();
}

void FOss2Parameterization::loadInputs(const char* filename)
{
	if (strlen(filename) == 0) return;

	ifstream fi(filename);
	vector<input> data;

        cout << "Loading data from " << filename << endl;
	int a;
	string ss;
        while (fi >> a)
        {
          //      cout << "Input #" << data.size() << "...";
                input inp;
                //fscanf(fi, "%d\n", &inp.nAtoms);
                inp.nAtoms = a;

        //        cout << inp.nAtoms << " atoms ";
                char atom;

                fi >> inp.energy; // -764.383003508 energy
		getline(fi, ss);

                if (inp.energy >= 0)
                {
            //            cout << "invalid input, ignored." << endl;

                        for (int i=0; i<inp.nAtoms; i++)
                        {
				getline(fi, ss);
				istringstream iss(ss);
                                iss >> atom >> inp.r[i][0] >> inp.r[i][1] >> inp.r[i][2] >> inp.g[i][0] >> inp.g[i][1] >> inp.g[i][2];
                        }

                        continue;
                }

                inp.r = new double*[inp.nAtoms];
                inp.g = new double*[inp.nAtoms];

       		inp.nO = 0;
                double tmp;
                for (int i=0; i<inp.nAtoms; i++)
                {
                        inp.r[i] = new double[3];
                        inp.g[i] = new double[3];
			getline(fi, ss);
			istringstream iss(ss);
                        iss >> atom >> inp.r[i][0] >> inp.r[i][1] >> inp.r[i][2] >> inp.g[i][0] >> inp.g[i][1] >> inp.g[i][2];

                        if (atom == 'O') inp.nO ++;
                }

              //  cout << "valid." << endl;

                data.push_back(inp);
        }

	inputs.push_back(data);
	weights.resize(inputs.size());

        cout << "Data set #" << inputs.size() << ": " << data.size() << " inputs..." << endl;
        fi.close();
}

void FOss2Parameterization::setWeight(int pos, double w)
{
	weights[pos] = w;
}

double FOss2Parameterization::evaluate_( double const *x)
{
	neval++;
	double alpha[100];

	int i, j, k;
	int c = 0;

//	cout << "Alpha: " << endl;
	for(i=0; i<nparam; i++)
	{
		if (!fixed[i])
		{
			alpha[i] = x[c];
			c++;
		}
		else
		{
			alpha[i] = oss2best[i];
		}
//		cout << alpha[i] << " " << fixed[i] << endl ;
	}
//	cout << endl;

	double totalrms = 0;

        for(i=0; i<inputs.size(); i++)
        {
                double tmp = rmsEnergy(alpha, inputs[i], false);
                totalrms += weights[i] * tmp;
        }

        cout << "Eval " << neval << ": " << totalrms << endl;

	return totalrms;
}

double* FOss2Parameterization::getBest()
{
	double * res = new double[ndim];

	int c=0;
	for(int i=0; i<nparam; i++)
	{
		if (!fixed[i])
		{
			res[c] = oss2best[i];
			c++;
		}
	}

	return res;
}

void FOss2Parameterization::printBestAlpha(double* x, const char* fname)
{
	FILE* fo = fopen(fname, "w");
	double res[100];

        int c=0;
	
	//fprintf(fo, "%10d %15.8lf\n", neval, bestEvaluation());
        for(int i=0; i<nparam; i++)
        {
                if (!fixed[i])
                {
                        res[i] = x[c];
                        c++;
                }
		else res[i] = oss2best[i];

		fprintf(fo, "%15.8lf\t%15.8lf\t%15.8lf\t%d\n", res[i], lb[i], ub[i], fixed[i]);
        }
	fprintf(fo, "%10d %15.8lf best\n", neval, bestEvaluation());
	fclose(fo);
}

void FOss2Parameterization::inv(double A[NMAX][NMAX], double Ainv[NMAX][NMAX], int N)
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

double FOss2Parameterization::energy(double** x, int nAtom, int nO, double alpha[])
{
    // all kind of distance: dOO, dHH, dOH
    
	int i, j, k, l;

	RP(i, nAtom) q[i] = ((i<nO)? -2 : 1);

    RP(i, nAtom)
    {
        FR(j, i+1, nAtom-1)
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
	for(i=nO;i<nAtom;i++)
		for(j=i+1;j<nAtom;j++); 
				//VHH += exp(-50.0*(r[i][j]-1));

    // <formula 9>
    double h5tmp = pow(1.0-p_h5,2) / (pow(1.0-p_h5,2) + p_h5*p_h5);

    double VOH = 0;

    RP(i, nO)
    {
        FR(j, nO, nAtom-1)
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
		FR(j, nO, nAtom-1)
		{
			FR(k, j+1, nAtom-1)
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
	
	RP(i, nAtom)
	{
		FR(j, i+1, nAtom-1)
		{
			Vq += q[i]*q[j] / r[i][j];
		}
	}
	// </part 1: V coulomb>

	// <part 2: ScdOH, ScdOO, SddOO>
	// formula 6
	RP(i, nO)
	{
		FR(j, nO, nAtom-1)
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

		FR(j, nO, nAtom-1)
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

void FOss2Parameterization::energy_grad(double** x, int nAtom, int nO, double alpha[], double& energy, double grad[][3])
{
	int i, j, k, l, t;
	double pt1, pt2, pt3, dp, dpk;

    // clear the grad matrix
	RP(i, nAtom) RP(j, 3) grad[i][j] = 0;

	// charges
	RP(i, nAtom) q[i] = ((i<nO)? -2 : 1);

	// distance vectors, euclidian distances and their derivations
    RP(i, nAtom)
    {
        FR(j, i+1, nAtom-1)
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
	
	RP(i, nAtom)
	{
		FR(j, i+1, nAtom-1)
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
		FR(j, i+1, nAtom-1)
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

		FR(j, nO, nAtom-1)
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

	
	
	RP(i, nAtom)
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
	RP(i, nAtom) RP(j, 3) grad[i][j] *= RESCALE;
	// </part 7: Vel>

// </formula 1>

    // <formula 9: 2-body interaction OH>
    double h5tmp = SQR(1.0-p_h5) / (SQR(1.0-p_h5) + p_h5*p_h5);

    double VOH = 0;
	

    RP(i, nO)
    {
        FR(j, nO, nAtom-1)
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
	
	for(i=nO;i<nAtom;i++)
		for(j=i+1;j<nAtom;j++){ //H-H
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
		FR(j, nO, nAtom-1)
		{
			FR(k, j+1, nAtom-1)
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

	/*
	// results
	cout << "VOO: " << VOO << " VOH: " << VOH << " VHOH: " << VHOH << endl;
	cout << "Vq: " << Vq << endl;
	cout << "Self-dipole: " << Vsd << endl;
	cout << "Dipole-dipole interaction: " << Vdd << endl;
	cout << "Charge-dipole interaction: " << Vcd << endl;
	cout << "Vel: " << Vel << endl;
	cout << "Energy: " << Ener << endl;
	*/	
}

double FOss2Parameterization::func(double* x, int nAtom, int nO, double alpha[])
{
	int i, j;
//	double** _x = new double*[nAtom];
	double _x[100][3];

	int cnt = 0;
	RP(i, nAtom)
	{
//		_x[i] = new double[3];
		RP(j, 3) _x[i][j] = x[cnt++];
	}

	double res = oss2energyRed(_x, nAtom, nO, alpha);

//	RP(i, nAtom) delete _x[i];

//	delete _x;

	return res;
}

double FOss2Parameterization::dfunc(double* x, int nAtom, int nO, double alpha[], double* grad)
{
	double res;
	int i, j;

	// doing type casting here - I hate this
//	double** _x = new double*[nAtom];
	double _x[100][3];

        int cnt = 0;

	double _grad[10][3];

        RP(i, nAtom)
        {
//		_x[i] = new double[3];
                RP(j, 3) _x[i][j] = x[cnt++];
        }
	

	//energy_grad(_x, nAtom, nO, alpha, res, _grad);

	res = oss2EFRed(_x, nAtom, nO, alpha, _grad);

	cnt = 0;
	RP(i, nAtom)
        {
//		delete _x[i];
                RP(j, 3) grad[cnt++] = _grad[i][j];
        }

//	delete _x;

	return res;
}


void FOss2Parameterization::lnsrch(double xold[], double fold, double g[], double p[], double x[],double *fret, double stpmax, int *check, int nAtom, int nO, double alpha[])
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam;
	
	// dimension
	int n = nAtom * 3;
	
	*check=0;
	for (sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);

	if (sum > stpmax)
		for(i=0;i<n;i++) p[i] *= stpmax/sum; //Scale if attempted step is too big.
	for (slope=0.0,i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) {};//cout<<"Roundoff problem in lnsrch.\n";
	test=0.0; //Compute min.

	for(i=0;i<n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;  //Always try full Newton step first.
	for (;;) { //Start of iteration loop.
		
		cout << "Shit here" << endl;
		RP(i, n) cout << x[i] << " " << xold[i] << " " << p[i] << " " << alam << endl;
		
		for(i=0;i<n;i++) x[i]=xold[i]+alam*p[i];

		// <modified by Quang Huy>
		// special function call
		*fret=func(x, nAtom, nO, alpha);

		cout << "Res: " << *fret << endl;

		// debug here
	//	if (*fret != *fret)
		

		if(*fret != *fret)	exit(98753016);
		


	//	cout << alam << " " << alamin << " " << *fret << " " << fold << endl;
		if (alam < alamin) { //Convergence on x. For zero finding,the calling program should verify the convergence.
			for(i=0;i<n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*fret <= fold+ALF*alam*slope) return; //Sufficient function decrease.
		else { //Backtrack.
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*fret-fold-slope)); //First time.
			else {// Subsequent backtracks.
				rhs1 = *fret-fold-alam*slope;
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
				if (tmplam > 0.5*alam)
					tmplam=0.5*alam; 
			}
		}
		alam2=alam;
		f2 = *fret;
		alam=FMAX(tmplam,0.1*alam); 
	}     
    
}


int FOss2Parameterization::dfpmin(double p[], double gtol, int iterMax, double *fret,double *rmsGrad, int nAtom, int nO, double alpha[])
{

	// dimension 
	int n = nAtom * 3;



	int check,i,its,j;
	double fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,test;	

	double dg[n+1];
	double g[n+1];
	double hdg[n+1];
	double hessin[n+1][n+1];
	double pnew[n+1];
	double xi[n+1];
	
	//Calculate starting function value and gradient,
	fp = dfunc(p, nAtom, nO, alpha, g);

	for (i=0;i<n;i++) { //and initialize the inverse Hessian to the unit matrix. 
		for (j=0;j<n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i]; //Initial line direction.
		sum += p[i]*p[i];
	}

	stpmax=STPMX;//*n;

	cout << "DFPMin" << endl;
	for (its=0;its<iterMax;its++) { //Main loop over the iterations.

	//	cout << its << ": " << fp << endl;
		lnsrch(p,fp,g,xi,pnew,fret,stpmax,&check, nAtom, nO, alpha);

		//The new function evaluation occurs in lnsrch; save the function value in fp for the
		//next line search. It is usually safe to ignore the value of check.
		fp = *fret;
		for(i=0;i<n;i++) {
			xi[i]=pnew[i]-p[i]; //Update the line direction,
			p[i]=pnew[i];  // and the current point.
		}	
			
		/*
		test=0.0; //Test for convergence on x.
		for(i=0;i<n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}	
		if (test < TOLX)  return its; */
		
		
		for(i=0;i<n;i++) dg[i]=g[i]; //Save the old gradient,
		dfunc(p, nAtom, nO, alpha, g); //and get the new gradient.
				
		test=0.0; //Test for convergence on zero gradient.		
		for(i=0;i<n;i++) {
			test+=g[i]*g[i];
		}
		test=sqrt(test/n);
		if (test < gtol) break;		
	
		for(i=0;i<n;i++) dg[i]=g[i]-dg[i]; //Compute difference of gradients,
		for(i=0;i<n;i++) { //and difference times current matrix.
			hdg[i]=0.0;
			for(j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0; //Calculate dot products for the denominators.
		for(i=0;i<n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) { //Skip update if fac not sufficiently positive.
			fac=1.0/fac;
			fad=1.0/fae;
    //The vector that makes BFGS different from DFP:
			for(i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for(i=0;i<n;i++) { //The BFGS updating formula:
				for (j=i;j<n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j] -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin[j][i]=hessin[i][j];
				}
			}
		}
		for(i=0;i<n;i++) { //Now calculate the next direction to go,
			xi[i]=0.0;
			for(j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}// and go back for another iteration.
	*rmsGrad=0;
	for(i=0;i<n;i++)	*rmsGrad+=SQR(g[i]);	
	*rmsGrad=sqrt((*rmsGrad)/n);	
		   
	return its;//if overrun	

}
