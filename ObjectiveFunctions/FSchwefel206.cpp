//===========================================================================
/*!
 *  \file FSchwefel206.cpp
 *
 *  \brief Ackley benchmark function
 *
 *  \author  Ong Yew Soon
 *  \date    2007-01-08
 *
 *  \par Copyright (c) 2006:
 *      Nanyang Technological University
 *      42 Nanyang Avenue
 *      Singapore 639798
 *      eMail: asysong@ntu.edu.sg
 *      www:   http:/www.ntu.edu.sg/home/asysong/ <BR>
 *      <BR>
 *
 *  \par Project:
 *      Memetic Algorithm
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  This file is part of MAP. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */
//===========================================================================
#include "FSchwefel206.h"
#include "../Rng/GlobalRng.h"

FSchwefel206::FSchwefel206( int nDimensions, const char* filename  ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
	loadData(filename);
}

FSchwefel206::FSchwefel206( int nDimensions, double lowerBound, double upperBound, const char* filename  ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
	loadData(filename);
}

FSchwefel206::FSchwefel206( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds, const char* filename  ):ObjectiveFunction( -1, nDimensions, lowerBounds, upperBounds)
{
	loadData(filename);
}

FSchwefel206::~FSchwefel206()
{
}


void FSchwefel206::loadData(const char* filename)
{
	int i, j;
	ifstream f(filename);

	double o[100];
	for(i=0; i<100; i++)
		f >> o[i];

	// After loading the data, we set o[1..D/4] = -100 o[3D/4..D] = 100
	for(i=0; i<ceil(1.0*nDim/4); i++) o[i] = -100;
	for(i=floor(3.0*nDim/4)-1; i<nDim; i++) o[i] = 100;


	for(i=0; i<100; i++)
	{		
		for(j=0; j<100; j++)
		{
			f >> A[i][j];
		}
	}


	for(i=0; i<nDim; i++)
	{		
		B[i] = 0;
		for(j=0; j<nDim; j++)
		{
			B[i] += A[i][j]*o[j];
		}
	}

	cout << "Evaluate o: " << evaluate_(o) << endl;
}
	
double FSchwefel206::evaluate_( vector<double>& x )
{
	ObjectiveFunction::evaluate_( x );

	int i,j;

	double res = INT_MIN;

	
	for( i=0; i< nDim; i++)
	{
		double tmp = 0.;
		for(j=0; j<nDim; j++)
		{
			tmp+=A[i][j]*x[j];
		}		
		tmp = fabs(tmp-B[i]);
		if (res < tmp) res = tmp;
	}	

	return res;
}

vector<double> FSchwefel206::gradient_( vector<double>& x )
{
	//<QuangHuy>	
	return finiteDifference(x);
	//</QuangHuy>
}
