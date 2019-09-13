/**
* @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>

* MUSCLE-HPC communication module
* Copyright (C) 2016  University of Geneva, Switzerland
*
* MUSCLE-HPC is free software: you can redistribute it and/or
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

#ifndef MAPPERUTIL_H
#define MAPPERUTIL_H

#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


using namespace std;

namespace unige_pasc{

//------ Added by Mohamed------------
/// Integer type for Tetra (taken from palabos)
/** On some architectures, this type is larger
 *  than int. Using plint instead of int ensures
 *  64-bit compatibility of the code.
 */
#ifdef ARC_BGP
    typedef long long int plint;
#else
    typedef ptrdiff_t plint;
#endif

/// Unsigned integer type Tetra (taken from palabos)
/** On some architectures, this type is larger
 *  than int. Using fluplint instead of unsigned plint
 *  ensures 64-bit compatibility of the code.
 **/
#ifdef ARC_BGP
    typedef unsigned long long int pluint;
#else
    typedef size_t pluint;
#endif

/**
 * @brief CantorPairingFunction created a uniq tag for the intercommunicator between two MPI sub-groups
 * @param k1 the color of the first MPI sub-group
 * @param k2 the color of the second MPI sub-group
 * @return
 */
inline int CantorPairingFunction(int handlerLeaderColor, int receiverLeaderColor){
    int min,max;
    if (handlerLeaderColor > receiverLeaderColor){
        max=handlerLeaderColor;
        min=receiverLeaderColor;
    }else{
        min=handlerLeaderColor;
        max=receiverLeaderColor;
    }
    return 0.5*(min+max)*(min+max+1)+max;
}


}
#endif // MAPPERUTIL_H
