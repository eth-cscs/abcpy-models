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


#ifndef MAPPERCOMMON_H
#define MAPPERCOMMON_H

#include "musclehpc/utils/util.h"
#include "musclehpc/utils/mmsfperf.h"

//#include <boost/date_time/posix_time/time_serialize.hpp>
//#include <boost/date_time/local_time/local_time.hpp>
//#include <boost/date_time/gregorian/greg_serialize.hpp>



using namespace std;

namespace unige_pasc{


/*********************************************** SUBMODEL ********************************************/
class Submodel_A{

   /* int argc;
    char** argv;*/

public:
    virtual void simulate(MMSFPerf<double> * mmsfprofiler=nullptr);

    // argv and argc per submodel
   /* int getArgc() const{ return this->argc;}
    char ** getArgv() const{ return this->argv;}
    void setArgc(int argc){this->argc=argc;}
    void setArgv(char ** argv){this->argv=argv;}*/

protected:
    Submodel_A(){}
    virtual ~Submodel_A(){}
    /**
     * @brief F_init initalizes the computation (F_init operator in MML).
     */
    virtual void F_init()=0;

    /**
     * @brief computeIteration runs one iteration of the model (operator S in MML)
     * @return
     */
    virtual void S()=0;

    /**
     * @brief U updates boundary domain of a  givensubmodel (Operator U in MML)
     */
    virtual void U()=0;

    /**
     * @brief Oi makes an intermediate observation (Operator Of in MML)
     */
    virtual void Oi()=0;

    /**
     * @brief Of makes a final observation (Operator Of in MML)
     */
    virtual void Of()=0;

    /**
     * @brief isConverged tells to stop or not the computation
     */
    virtual bool isConverged()=0;

};


/************************************** MapperType **************************************************/

typedef   enum  MapperType { SUBMODEL , SMART_MAPPER , MAPPER , FILTER} MapperType;

class INFO_MAPPER{

public:
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
//        ar & this->type;
//        ar & this->timeStamp;
//    }

   //vars
   MapperType type;
   time_t timeStamp; //when computation started
};
////////////////////////

class INFO_SMART_MAPPER: public INFO_MAPPER{


public:
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
//        ar & boost::serialization::base_object<INFO_MAPPER>(*this);
//        ar & this->dx;
//        ar & this->dt;
//        ar & this->sectionId;

//    }
    //vars


    double dt,dx;
    int sectionId;
};

/****************************************************************************************************/
class MapperHeader{

public :

protected:

    MapperHeader();
    virtual ~MapperHeader();
    virtual void synchronizeCommunication()=0; // since it requires connector
    virtual void serializeHeader(vector<char> &toSend, const INFO_SMART_MAPPER &info_header);
    virtual void unserializeHeader(vector<char>  & receivedVect, INFO_SMART_MAPPER  & info_header);


    //vars
    INFO_SMART_MAPPER info;
};





}


#endif // MAPPERCOMMON_H
