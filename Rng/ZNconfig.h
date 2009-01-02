/* ZNconfig.h */
/* ======================================================================
 *
 *  File    :  ZNconfig.h
 *  Created :  2001-08-17
 *
 *  Copyright (c) 1995-2001 Martin Kreutz
 *  EMail - Martin.Kreutz@zn-ag.de
 *
 *  ZN Vision Technologies AG
 *  Universitaetsstrasse 160
 *  44801 Bochum
 *  Germany
 *  
 *  Last update: 2001-08-17
 *
 * ----------------------------------------------------------------------
 *
 *  This file is part of the Shark distribution. Shark is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  Shark is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * ======================================================================
 */


/* 
 *
 * As the new Shark-version (1.0.8) was arranged at the ZN 
 * (Zentrum fuer Neuroinformatik, see above),
 * it was necessary to respect the general directory structure used
 * there. Originally "ZNconfig.h" includes several files with
 * definitions for all software projects developed at the ZN.
 * Because these files are not for public use we created this "dummy"
 * file with the same name and all the definitions necessary to
 * compile the Shark-distribution.
 *
 */ 


#ifndef ZNCONFIG_H
#define ZNCONFIG_H


// Compiler-dependent definitions:

#ifdef _WIN32

#    include <limits.h>
#    include <float.h>

#    ifdef _MSC_VER
#        ifndef NOMINMAX
#	     define NOMINMAX  // prevent Microsoft from defining min and max 
                              // in windef.h
#	 endif

	 // disable warning C4250: inherits 'ClassName::memberName' via 
         // dominance
#	 pragma warning(disable: 4250)

	 // disable warning C4786: debugging symbol length > 255 character
#	 pragma warning(disable: 4786)

#    endif		

#    ifndef __MIN_MAX__
#        define __MIN_MAX__
#if _MSC_VER < 1300 // because .NET has owned macros with the same names and causes an error
         namespace std 
         {
             //
             // undefine macros min and max to avoid conflicts with template 
	     // names
             //
#            ifdef min
#                undef min
#            endif

#            ifdef max
#                undef max
#            endif

             //
             // template functions 'min' and 'max' are not defined for _WIN32
             // due to name conflicts
             //
             template < class T > inline T min( T a, T b ) 
             { 
                 return a < b ? a : b; 
             }

             template < class T > inline T max( T a, T b ) 
             { 
                 return a > b ? a : b; 
             }
         }
#endif // _MSC_VER 
#    endif

#else

#    ifdef __GNUC__
#    endif

#    include<limits.h>
#    include<float.h>

#endif


// Compiler-independent definitions:

// Type definitions:

#ifdef _MSC_VER
    // Floating point value with double precision
    typedef double ZNdouble;
#else
#ifdef __GNUC__ 
    // Floating point value with double precision
    typedef double ZNdouble;
#else
#ifdef __QNX__
    // Floating point value with double precision
    typedef double ZNdouble;
#else
#   error "cannot determine the system and/or compiler type"
#endif
#endif
#endif


// Definition of constants:

// Contains the value of e
const ZNdouble znEC         = 2.7182818284590452354;    

// Contains the value of pi
const ZNdouble znPiC        = 3.14159265358979323846;   

// Contains the value of 2*pi
const ZNdouble znPi2C       = 6.28318530717958647692;

//! Contains the value of sqrt(2)
const ZNdouble znSqrt2C     = 1.41421356237309504880;

//! Contains the value of 1/sqrt(2)
const ZNdouble znSqrt1_2C   = 0.70710678118654752440;

// Contains the value of sqrt(pi)
const ZNdouble znSqrtPiC    = 1.7724538509055160273;

// Contains the value of sqrt(2*pi)
const ZNdouble znSqrtPi2C   = 2.50662827463100050242;

//! Contains the value of 2/sqrt(pi)
const ZNdouble zn2_SqrtPiC  = 1.12837916709551257390;

// Contains the value of sqrt(e)
const ZNdouble znSqrtEC     = 1.6487212707001281468;

//! Max. possible value of ZNdouble
const ZNdouble  ZNdouble_MAX = DBL_MAX;


#endif // ZNCONFIG_H





