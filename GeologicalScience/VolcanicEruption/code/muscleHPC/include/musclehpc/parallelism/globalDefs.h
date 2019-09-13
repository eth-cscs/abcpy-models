/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Global Definitions -- header file.
 */
#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

#ifdef PLB_MPI_PARALLEL
#include "mpi.h"
#endif

#include <limits>
#include <string>
#include <cstddef>


namespace plb {

// We take the definition that the base type of T is just T, except
// in special cases. For example, for complex numbers, the base type
// of Complex<U> is U. PlbTraits is correspondingly overloaded in
// plbComplex.h.
template<typename T>
struct PlbTraits {
    typedef T BaseType;
};

/// Integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using plint instead of int ensures
 *  64-bit compatibility of the code.
 */
#ifdef PLB_BGP
typedef long long int plint;
#else
typedef ptrdiff_t plint;
#endif

/// Unsigned integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using fluplint instead of unsigned plint 
 *  ensures 64-bit compatibility of the code.
 **/
#ifdef PLB_BGP
typedef unsigned long long int pluint;
#else
typedef size_t pluint;
#endif

/// Define the endianness
/** The endianess is little-endian by default.
 *  On the Blue Gene/P, the endianness is
 *  big-endian.
 */
#ifdef PLB_BGP
#define PLB_BIG_ENDIAN
#endif

/// Unsigned integer type for tracking global ids of individual objects.
typedef unsigned long id_t;

/// Enumeration type that sets the precision for further use.
/// Single precision (float): FLT
/// Double precision (double): DBL
/// Extended precision (long double): LDBL
/// "Infinite" precision: INF
enum Precision { FLT, DBL, LDBL, INF };

template<typename T>
inline Precision floatingPointPrecision()
{
    if (sizeof(T) == sizeof(float)) {
        return FLT;
    } else if (sizeof(T) == sizeof(long double)) {
        return LDBL;
    }

    return DBL;
}

template<typename T>
inline T getEpsilon(Precision precision)
{
    T epsilon;

    switch (precision) {
    case FLT:
        epsilon = std::numeric_limits<float>::epsilon();
        break;
    case DBL:
        epsilon = std::numeric_limits<double>::epsilon();
        break;
    case LDBL:
        epsilon = std::numeric_limits<long double>::epsilon();
        break;
    case INF: default:
        epsilon = std::numeric_limits<T>::min();
        break;
    }

    T coef = 10.0; // hack for better results

    return (coef * epsilon);
}

// Version that works also for integral types, and always refers to the
// type defined by the template argument. This means that if T is a floating point
// type, then "getEpsilon<T>()" is the same as "getEplsilon<T>(floatingPointPrecision<T>())".
// The function "getEpsilon<T>(Precision)" works for floating point types only, and can
// be used with two different precisions (e.g. T is double, but Precision is float). This
// version works with one type (T) and T can be also integral (in such a case, epsilon is
// zero).
template<typename T>
inline T getEpsilon()
{
    T epsilon = std::numeric_limits<T>::epsilon();
    T coef = 10.0; // hack for better results
    return(coef * epsilon);
}

}  // namespace plb

#endif
