//===========================================================================
/*!
 *  \file Cauchy.cpp
 *
 *  \brief Implements methods for class Cauchy, that simulates a 
 *         "standard %Cauchy distribution".
 *
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Copyright (c) 1995,1998:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>
 *
 *  \par Project:
 *      Rng
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Cauchy.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:20 $
 *
 *  \par Changes:
 *      $Log: Cauchy.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:20  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:02  shark-admin
 *      INI Administration
 *
 *      Revision 1.2  2002/05/16 13:48:17  rudi
 *      doxygen commands added.
 *
 *
 *  This file is part of Rng. This library is free software;
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
 *
 */
//===========================================================================

#include "ZNconfig.h"

#include <cmath>
//#include "Rng/Cauchy.h"
#include "Cauchy.h"


//========================================================================
/*!
 *  \brief Returns a %Cauchy random number.
 *
 *  This method performs the \em inverse \em transformation of the
 *  original uniformally distributed random numbers of the interval
 *  (0,1) created by the used pseudo random number generator to
 *  the type of the %Cauchy distribution. <br>
 *  Therefore a method based on "Numerical Recipes in C, p. 220"
 *  is used.
 *
 *  \return a %Cauchy random number
 *
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
double Cauchy::operator( )( )
{
    double x;

    // tan( 0.5 * Pi ) is not defined !!!
    do
        x = rng( );
    while( fabs( x - 0.5 ) < 1e-30 );

    return tan( x * znPiC );
}


//========================================================================
/*!
 *  \brief Returns the probability for the occurrence of random number
 *          "x". 
 *
 *  \param x the random number for which the occurrence probability
 *           is returned
 *  \return the probability for the occurrence of random number \em x
 *
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
double Cauchy::p( const double& x ) const
{
    return 1. / znPiC * ( 1 + x*x );
}





