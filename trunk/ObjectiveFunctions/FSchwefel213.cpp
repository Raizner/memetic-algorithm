//===========================================================================
/*!
 *  \file FSchwefel213.cpp
 *
 *  \brief Schwefel's problem 2.13
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
#include "FSchwefel213.h"
#include "../Rng/GlobalRng.h"

FSchwefel213::FSchwefel213( int nDimensions, const char* filename  ):ObjectiveFunction( -1, nDimensions, -100, 100)
{
	loadData(filename);
}

FSchwefel213::FSchwefel213( int nDimensions, double lowerBound, double upperBound, const char* filename  ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
	loadData(filename);
}

FSchwefel213::FSchwefel213( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds, const char* filename  ):ObjectiveFunction( -1, nDimensions, lowerBound, upperBound)
{
	loadData(filename);
}

FSchwefel213::~FSchwefel213()
{
}


void FSchwefel213::loadData(const char* filename)
{
	int i, j;
	ifstream f(filename);

	for(i=0; i<100; i++)
	{		
		for(j=0; j<100; j++)
		{
			f >> A[i][j];			
		}
	}

	for(i=0; i<100; i++)
	{		
		for(j=0; j<100; j++)
		{
			f >> B[i][j];			
		}
	}

	for(i=0; i<100; i++)
	{
			f >> alpha[i];
	}

	f.close();
}

double FSchwefel213::evaluate_( vector<double>& x )
{
	int i,j;
	double *x_;

	ObjectiveFunction::evaluate_( x );

	x_ = new double[ nDim ];

	for( i = 0; i < nDim; i++ ) {
		x_[ i ] = x[ i ];
	}

	transform(x_);

	double res = 0.0;
	
	for( i=0; i< nDim; i++)
	{
		double sum1 = 0.0;
        double sum2 = 0.0;
        for (j=0; j<nDim; j++)
        {
            sum1 += A[i][j]*sin(alpha[j]) + B[i][j]*cos(alpha[j]);
            sum2 += A[i][j]*sin(x[j]) + B[i][j]*cos(x[j]);
        }
        res += pow((sum1-sum2),2.0);
	}

	delete x_;

	return res;
}

void FSchwefel213::gradient_( vector<double>& x )
{
	ObjectiveFunction::gradient_(x,grad);
	//<QuangHuy>	
	finiteDifference(x, grad);
	//</QuangHuy>
}
