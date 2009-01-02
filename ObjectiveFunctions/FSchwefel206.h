//===========================================================================
/*!
 *  \file FSchwefel206.h
 *
 *  \brief Schewefel's problem 2.6
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
Function #5: Schewefel's problem 2.6
******************************************************************************/

#ifndef _FSchwefel206_H
#define _FSchwefel206_H

#include "ObjectiveFunction.h"

class FSchwefel206:public ObjectiveFunction
{
public:
	FSchwefel206( int nDimensions, const char* filename );
	FSchwefel206( int nDimensions, double lowerBound, double upperBound, const char* filename );
	FSchwefel206( int nDimensions, vector<double>& lowerBounds, vector<double>& upperBounds, const char* filename );
	~FSchwefel206();

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
	virtual vector<double> gradient_( vector<double>& x );	
	virtual void loadData(const char* filename);

private:
	// private data
	double A[100][100];
	double B[100];
};

#endif
