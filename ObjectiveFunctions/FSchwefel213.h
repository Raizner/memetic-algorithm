//===========================================================================
/*!
 *  \file FSchwefel213.h
 *
 *  \brief Schwefel's problem 2.13
 *
 *  \author  Ong Yew Soon
 *  \date    2006-06-19
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

/******************************************************************************
Function #5: Schewefel's problem 2.13
******************************************************************************/

#ifndef _FSchwefel213_H
#define _FSchwefel213_H

#include "ObjectiveFunction.h"

class FSchwefel213:public ObjectiveFunction
{
public:
	FSchwefel213( int nDimensions, const char* filename );
	FSchwefel213( int nDimensions, double lowerBound, double upperBound, const char* filename );
	FSchwefel213( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds, const char* filename );
	~FSchwefel213();

// <identifier>
	virtual const string equation()
	{
		return "";
	}
	virtual const string classname()
	{
		return "Schwefel with noise";
	}
// </identifier>

private:
	virtual double evaluate_( vector<double>& x );
	virtual void gradient_( vector<double>& x );	
	virtual void loadData(const char* filename);

private:
	// private data
	double A[100][100];
	double B[100][100];
	double alpha[100];
};

#endif
