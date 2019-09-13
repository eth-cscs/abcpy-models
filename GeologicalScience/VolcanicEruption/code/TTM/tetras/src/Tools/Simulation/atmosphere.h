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


#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include  "utility.h"
#include  "musclehpc/utils/mmsflogger.h"
#include  "musclehpc/mapper/mappercommon.h"
#ifdef MUSCLE
#include  "mapper.h"
#endif

#include <memory>
#ifdef PROFILING
#include  "musclehpc/utils/profiling.cpp"
#endif
#include "musclehpc/conduits/conduit.h"

#include <netcdf.h>

using namespace std;
//using namespace netCDF;
//using namespace netCDF::exceptions;


using namespace std;
/*************************** Coord3D ********************************************/

struct Coord3D{

protected:
    double x,y,z;
public:
    Coord3D(){
    }
    Coord3D(double x, double y, double z): x(x), y(y), z(z){
    }

    Coord3D(const Coord3D & coord){
        this->x=coord.getX();
        this->y=coord.getY();
        this->z=coord.getZ();
    }

    //affectation
    Coord3D & operator=(const Coord3D & coord){
        this->x=coord.getX();
        this->y=coord.getY();
        this->z=coord.getZ();
        return *(this);
    }

    virtual ~Coord3D(){}
    double getX()const {return this->x;}
    double getY()const {return this->y;}
    double getZ()const {return this->z;}

    void serialize(vector<char> & buffer){
        buffer.resize(sizeof(double)*3) ;
        char * offset= &buffer[0];
        memcpy(offset, &x, sizeof(x)); // "Serialize"
        offset+=sizeof(x);
        memcpy(offset, &y, sizeof(y)); // "Serialize"
        offset+=sizeof(y);
        memcpy(offset, &z, sizeof(z)); // "Serialize"
    }

    char * unSerialize(vector<char> & buffer, char * offset){
        memcpy(&x, offset, sizeof(x)); // "unSerialize"
        offset+=sizeof(x);
        memcpy(&y, offset, sizeof(y)); // "unSerialize"
        offset+=sizeof(y);
        memcpy(&z, offset, sizeof(z)); // "unSerialize"
        offset+=sizeof(z);
        return offset;
    }


    string toString(){
        stringstream ss;
        ss<<"("<<this->x <<","<<this->y<<","<<this->z<<")";
        return ss.str();
    }

};

/*************************** Coord4D ********************************************/
struct CoordGeo{

protected:
    double lattitude, longitude;
public:
    CoordGeo(){}
    CoordGeo(double lattitude, double longitude): lattitude(lattitude), longitude(longitude){}
    virtual ~CoordGeo(){}

    CoordGeo(const CoordGeo & coord){
        this->lattitude=coord.getLattitude();
        this->longitude=coord.getLongitude();
    }


    //affectation
    CoordGeo & operator=(const CoordGeo & coord){
        this->lattitude=coord.getLattitude();
        this->longitude=coord.getLongitude();
        return *(this);
    }


    double getLattitude() const {return this->lattitude;}
    double getLongitude() const {return this->longitude;}

    void serialize(vector<char> & buffer){
        buffer.resize(sizeof(double)*2) ;
        char * offset= &buffer[0];
        memcpy(offset, &lattitude, sizeof(lattitude)); // "Serialize"
        offset+=sizeof(lattitude);
        memcpy(offset, &longitude, sizeof(longitude)); // "Serialize"
        offset+=sizeof(longitude);
    }

    char*  unSerialize(vector<char> & buffer, char* offset){
        assert(buffer.size() >= sizeof(double)*2);
        memcpy(&lattitude, offset, sizeof(lattitude)); // "unSerialize"
        offset+=sizeof(lattitude);
        memcpy(&longitude, offset, sizeof(longitude)); // "unSerialize"
        offset+=sizeof(longitude);
        return offset;
    }



    string toString(){
        stringstream ss;
        ss<<"("<<lattitude <<","<<longitude<<")";
        return ss.str();
    }

};

/*************************** AtmosphereDataRequest ********************************************/

class AtmosphereDataRequest{

private:
    Coord3D craterPosition, domainSize, domainPosition;
    CoordGeo craterGeoPosition;

    string date;
    double time;

public:
    AtmosphereDataRequest(){}
    AtmosphereDataRequest(Coord3D const &craterPosition,
                          Coord3D const &domainSize,
                          Coord3D const &domainPosition,
                          CoordGeo const &craterGeoPosition,
                          string date, double time):
        craterPosition(craterPosition),
        domainSize(domainSize),
        domainPosition(domainPosition),
        craterGeoPosition(craterGeoPosition),
        date(date), time(time){

    }
    virtual ~AtmosphereDataRequest(){}
    //setter
    void setCraterePosition(Coord3D const craterPosition){
        this->craterPosition = craterPosition;
    }
    void setDomainSize(Coord3D const domainSize){
        this->domainSize = domainSize;
    }
    void setDomainPosition(Coord3D const domainPosition){
        this->domainPosition = domainPosition;
    }
    void setCratereGeoPosition(CoordGeo const craterGeoPosition){
        this->craterGeoPosition = craterGeoPosition;
    }
    void setTime(double time){this->time=time;}
    void setDate(string date){this->date=date;}

    //getter
    Coord3D getCraterePosition() const{
        return this->craterPosition;
    }
    Coord3D getDomainSize()const{
        return this->domainSize ;
    }
    Coord3D getDomainPosition() const{
        return this->domainPosition;
    }
    CoordGeo getGeoCraterePosition() const{
        return this->craterGeoPosition;
    }
    double getTime() const{ return this->time;}
    string getDate() const {return  this->date;}

    // convert the class attributes to vector<char>
    void serialize(vector<char> & buffer){

        vector<char> buf_craterPosition, buf_domainSize, buf_domainPosition, buf_craterGeoPosition, buf_timeDate;

        this->craterPosition.serialize(buf_craterPosition);
        this->domainSize.serialize(buf_domainSize);
        this->domainPosition.serialize(buf_domainPosition);
        this->craterGeoPosition.serialize(buf_craterGeoPosition);

        //serialize date + time
        buf_timeDate.resize(sizeof(double)+ date.size());
        memcpy(&buf_timeDate[0], &time, sizeof(time)); // "Serialize"
        memcpy(&buf_timeDate[0]+sizeof(time), date.c_str(), date.size()); // "Serialize"


        //copy class attributes in one buffer
        buffer.insert( buffer.end(), buf_craterPosition.begin(), buf_craterPosition.end() );
        buffer.insert( buffer.end(), buf_domainSize.begin(), buf_domainSize.end() );
        buffer.insert( buffer.end(), buf_domainPosition.begin(), buf_domainPosition.end() );
        buffer.insert( buffer.end(), buf_craterGeoPosition.begin(), buf_craterGeoPosition.end() );
        buffer.insert( buffer.end(), buf_timeDate.begin(), buf_timeDate.end() );
        //free tmp vectors
        vector<char>().swap(buf_craterPosition);
        vector<char>().swap(buf_domainSize);
        vector<char>().swap(buf_domainPosition);
        vector<char>().swap(buf_craterGeoPosition);
        vector<char>().swap(buf_timeDate);
    }

    // convert the class attributes to vector<char>
    void unSerialize(vector<char> & buffer){
        char *offset =&buffer[0];
        this->unSerialize(buffer, offset);
    }

    // convert the class attributes to vector<char>
    void unSerialize(vector<char> & buffer, char *offset ){
        // important: keep the same order as in serialize !!
        offset= this->craterPosition.unSerialize(buffer, offset);
        offset= this->domainSize.unSerialize(buffer, offset);
        offset= this->domainPosition.unSerialize(buffer, offset);
        offset= this->craterGeoPosition.unSerialize(buffer,offset);
        memcpy(&time, offset, sizeof(time)); // "Serialize"
        offset+=sizeof(time);
        int length= &buffer[buffer.size() -1] - offset + 1;
        this->date= std::string(offset, length);
    }

    string toString (){
        stringstream ss;
        ss<<"   craterPosition       : " <<this->craterPosition.toString()<<"\n";
        ss<<"   getDomainSize        : " <<this->domainSize.toString()<<"\n";
        ss<<"   getDomainPosition    : " <<this->domainPosition.toString()<<"\n";
        ss<<"   getGeoCraterePosition: " <<this->craterGeoPosition.toString()<<"\n";
        ss<<"   time                 :"<<this->time<<"\n";
        ss<<"   datet                :"<<this->date<<"\n";
        return ss.str();
    }

};

/*************************** AtmosphereDataResponse ********************************************/

class AtmosphereDataResponse{
    double dx, dy, dz;
    vector<double>  data;

public:
    AtmosphereDataResponse(){}
    AtmosphereDataResponse(double dx, double dy, double dz, vector<double> & data): dx(dx), dy(dy), dz(dz), data(data){
    }
    virtual ~AtmosphereDataResponse(){}
    vector<double> const & getData() const {return this->data;}
    double getDx() const {return this->dx;}
    double getDy()const {return this->dy;}
    double getDz()const {return this->dz;}

    void serialize(vector<char> & buffer){
        buffer.resize(sizeof(double)*(3 + data.size())) ;
        char * offset= &buffer[0];
        memcpy(offset, &dx, sizeof(dx)); // "Serialize"
        offset+=sizeof(dx);
        memcpy(offset, &dy, sizeof(dy)); // "Serialize"
        offset+=sizeof(dy);
        memcpy(offset, &dz, sizeof(dz)); // "Serialize"
        offset+=sizeof(dz);
        //buffer.insert( buffer.end(), data.begin(), data.end() );
        char * offset_data= (char*) &data[0];
        for (size_t i =0 ; i< data.size(); i++){
            memcpy(offset, offset_data, sizeof(double)); // "Serialize"
            offset+=sizeof(double);
            offset_data+=sizeof(double);
        }
    }

    void unSerialize(vector<char> & buffer){

        assert(buffer.size() >= sizeof(double)*3);
        char * offset=&buffer[0];
        memcpy(&dx, offset, sizeof(dx)); // "unSerialize"
        offset+=sizeof(dx);
        memcpy(&dy, offset, sizeof(dy)); // "unSerialize"
        offset+=sizeof(dy);
        memcpy(&dz, offset, sizeof(dz)); // "unSerialize"
        offset+=sizeof(dz);
        data.clear();
        int length = ((&buffer[buffer.size()-1])-offset +1) ;
        data.resize(length / sizeof(double));
        memcpy(&data[0], offset, length); // "unSerialize"
    }

    string toString(){
        stringstream ss;
        ss<<"("<<getDx()<<","<<getDy()<<","<<getDz()<<")"<<"- data(";
        for(size_t i=0; i< getData().size(); i++){
            ss<<getData().at(i)<<" ";
        }
        ss<<")"<<endl;
        return ss.str();
    }
};


/*************************** Atmosphere ********************************************/

//To be implemented by Pierre

class AtmosphereSubModel: public Submodel_A{

public:
    virtual ~AtmosphereSubModel(){}

protected:
    AtmosphereSubModel(double dx, double dt): dx(dx), dt(dt), cpt(0){

    }// protected -> can not instantiated

    /**
     * @brief F_init initalizes the computation (F_init operator in MML).
     */
    virtual void F_init(){}

    /**
     * @brief computeIteration runs one iteration of the model (operator S in MML)
     * @return
     */
    virtual void S(){
        cpt++;
    }

    //virtual void U(); Not inplemented here

    /**
     * @brief Oi makes an intermediate observation (Operator Of in MML)
     */
    virtual void Oi(){}
    /**
     * @brief Of makes a final observation (Operator Of in MML)
     */
    virtual void Of(){}

    /**
     * @brief isConverged tells to stop or not the computation
     */
    virtual bool isConverged(){
        return (cpt >= 1000);
    }

protected:
    double dx, dt;
    int cpt;
    shared_ptr<plb_ofstream> logger;
};

/************************************** DistributedPlume **************************************************/
class DistributedAtmosphereSubModel:  public AtmosphereSubModel, public MapperHeader{

public:
    virtual ~DistributedAtmosphereSubModel();
protected:
    DistributedAtmosphereSubModel(double dx, double dt);

    virtual void synchronizeCommunication();
    virtual void S();
    virtual void U();
    virtual void F_init();
    virtual void Of();
    virtual bool isConverged();

  AtmosphereDataResponse createAtmosphereDataResponse(AtmosphereDataRequest &atmRequest);

protected:
    //abstract functions
    virtual void send(vector<char> const& data)=0;
    virtual void receive(vector<char> & data)=0;
    virtual int getSubmodel_ID()=0;
    virtual void end()=0;

protected:
    //bool isRemoteConverged;//tells whether the remote submodel is converged
    MpiManager * manager;
#ifdef PROFILING
    unige_pasc::Profiler<double> profiler;
    int nbrItrToSaveStatsInFile;
#endif
    bool isRemoteConverged;
};


/************************************** LightDistributedAtmosphereSubModel **************************************************/

class LightDistributedAtmosphereSubModel: public DistributedAtmosphereSubModel, public AsynchronousRemoteMpiKernel{


public:
    LightDistributedAtmosphereSubModel(double dx, double dt);
    virtual ~LightDistributedAtmosphereSubModel();

    virtual void mainLoop();

protected:
    // implement abstract functions
    virtual void send(vector<char> const& data);
    virtual void receive(vector<char> & data);
    virtual int getSubmodel_ID();
    virtual void end();
    virtual void F_init();


protected:
    string id;
    ConduitEntrance<char> * f_out;
    ConduitExit<char> * f_in;
};


/***************************************************************************************************************/

class LightDistributedAtmosphereRequester: public DistributedAtmosphereSubModel, public AsynchronousRemoteMpiKernel{


public:
    LightDistributedAtmosphereRequester(double dx, double dt);
    virtual ~LightDistributedAtmosphereRequester();

    virtual void mainLoop();
    virtual void U();

protected:
    // implement abstract functions
    virtual void send(vector<char> const& data);
    virtual void receive(vector<char> & data);
    virtual int getSubmodel_ID();
    virtual void end();
    virtual bool isConverged();
    virtual void F_init();

protected:
    string id;
    ConduitEntrance<char> * f_out;
    ConduitExit<char> * f_in;
};











#endif // ATMOSPHERE_H
