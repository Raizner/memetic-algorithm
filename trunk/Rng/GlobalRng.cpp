//===========================================================================
/*!
 *  \file GlobalRng.cpp
 *
 *  \brief Implements a method for class Rng that subsumes several often 
 *         used random number generators.
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
 *      $RCSfile: GlobalRng.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:20 $
 *
 *  \par Changes:
 *      $Log: GlobalRng.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:20  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:02  shark-admin
 *      INI Administration
 *
 *      Revision 1.3  2002/05/16 13:48:30  rudi
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


//#include "Rng/GlobalRng.h"
#include "GlobalRng.h"

//! Instance of class Bernoulli:
Bernoulli       Rng::coinToss;

//! Instance of class DiscreteUniform:
DiscreteUniform Rng::discrete;

//! Instance of class Uniform:
Uniform         Rng::uni;

//! Instance of class Normal:
Normal          Rng::gauss;

//! Instance of class Cauchy:
Cauchy          Rng::cauchy;

//! Instance of class Geometric:
Geometric       Rng::geom;

//! Instance of class DiffGeometric:
DiffGeometric   Rng::diffGeom;

//! Instance of class Poisson:
Poisson         Rng::poisson;


//========================================================================
/*!
 *  \brief Sets the seed for all random number generators to "s".
 *
 *  \param s the seed for the random number generators in this class 
 *  \return none
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
 *  \sa RandomVar::seed
 *
 */
void Rng::seed( long s )
{
    coinToss.seed( s );
    discrete.seed( s );
    uni     .seed( s );
    gauss   .seed( s );
    cauchy  .seed( s );
    geom    .seed( s );
    diffGeom.seed( s );
    poisson .seed( s );
}




