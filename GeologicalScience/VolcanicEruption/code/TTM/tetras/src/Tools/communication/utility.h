#ifndef UTILITY_H
#define UTILITY_H

/**
* @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>
* @version 1.0
* @section LICENSE

* MAPPER communication module
* Copyright (C) 2015  University of Geneva, Switzerland
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/asio.hpp>

#ifdef USE_PARALLEL_MPI
    #include <mpi.h>
#endif
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/gregorian/greg_serialize.hpp>

#include "../../piaf/include/Simulator/SimulatorTypes.hpp"



using namespace std;
using namespace piaf;


//------------------------------------

namespace unige_pasc{

/*************************************** SerializerParticles  *************************************************/

/**
 * @brief The SerializerParticles class used to serialize/deserilaize  std::vector<BoundaryParticle>
 */
class SerializerParticles
{
    friend class boost::serialization::access;
    std::vector<BoundaryParticle> & boundaryParticles;
    bool & isLocallyConverged;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & isLocallyConverged;
        ar & boundaryParticles;

    }
public:
    SerializerParticles(std::vector<BoundaryParticle> & boundaryParticles, bool & isLocallyConverged):
        boundaryParticles(boundaryParticles), isLocallyConverged(isLocallyConverged){

    }
    std::vector<BoundaryParticle>  &getboundaryParticles ()const{
        return this->boundaryParticles;
    }
};

/*************************************** util  *************************************************/

class MapperUtil
{
public:
    MapperUtil();


    template <typename T>
    static  string  stringify(T & number);

    template <typename T>
    static void unstringify(std::string str, T& value);

    static void unstringify_pluint(std::string str, pluint& value);
#ifdef USE_PARALLEL_MPI
    static void broadCast(std::string& message, int rank, int root, MPI_Comm comm);
#endif

    static time_t to_time_t(const boost::posix_time::ptime& date );
    /// Round to next plint value

    template<typename T>
    static pluint roundTo_uint(T value) {
        return value>0 ?
                   static_cast<pluint>(value+(T)0.5) :
                   static_cast<pluint>(value-(T)0.5);
    }

    /**
     * converts to big endiant representation when copying bytes if BIG_ENDIAN_SERIALIZER is defined .
     * This is useful when distributed computation is required: e.g. between a BIG-Endian arch
     * and LITTLE-ENDIAN arch.
     * The default exchanging format should be BIG_ENDIAN.
     * @author: mohamed ben belgacem
     * @date:18-10-2015
     */
    static void* memcpyEndian(void* dest, const void* src, size_t count) ;
    static void serialize(std::vector< BoundaryParticle >  & boundaryParticles, bool isConvergedSimulation, vector<char> & toSend, int cpt=0);
    static void deserialize(std::vector< BoundaryParticle > & newBoundaryParticles, bool & isConvergedSimulation,  vector<char> & vect, int cpt=0);

    static void memcopyIntoVector(BoundaryParticle & p, vector<char> & buffer, size_t &iData);
    static void memcopyIntoParticle(BoundaryParticle &p, vector<char> & buffer, size_t &iData);

    static void serialize(bool isConverged, vector<char> & toSend);
    static void deserialize(bool &isConverged, vector<char> const& toSend);


//static void serializeHeader(vector<char> & toSend, const unige_pasc::INFO_SMART_MAPPER &info_header);
    //static void unserializeHeader(vector<char>  & receivedVect, unige_pasc::INFO_SMART_MAPPER  & info_header);

};

/**************************************** membuf ************************************************/

/**
 * @brief The membuf struct used to retrieve a vector<char> from a boost::asio::streambuf object since 'setg' is protected.
 */
struct membuf : public std::streambuf
{
    membuf(char* begin, char* end) {
        this->setg(begin, begin, end);
    }
};


}


#endif // UTILITY_H
