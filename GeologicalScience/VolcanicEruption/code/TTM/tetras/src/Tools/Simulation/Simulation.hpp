/*
TEphra TRAnsport Simulator (tetras)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef Simulation_H_
#define Simulation_H_

#include  "musclehpc/utils/mmsflogger.h"


#include "../ToolsTypes.hpp"

#include <iostream>
#include <stdexcept>
#include <math.h>

#include "../FileManager/VolcanoMPIFileManager.hpp"
#include "../../Models/VolcanicTephra.hpp"
#include "../EruptionParameters/EruptionParameters.hpp"
#include "../EruptionParameters/FakeEruptionParameters.hpp"

#include "piaf.hpp"

//============ added by mohamed ===========================
#include  "../communication/utility.h"
#include  "../communication/connectors/connectorcommon.h"
#include  "musclehpc/mapper/mappercommon.h"
#ifdef MUSCLE
    #include  "../communication/mapper.h"
#endif

#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/asio.hpp>
#include <memory>

#include <functional>

#ifdef PROFILING
    #include  "../../../../communication/utils/profiling.cpp"
#endif

#include  "atmosphere.h"
#include "../../parseargs.hpp"

#include <random>

#include "eventlogger.hpp"


//=========================================================


using namespace unige_pasc;
using namespace std;
using namespace piaf;



/*********************************** Simulation_A *****************************************************/
class Plume : public Submodel_A{

public:

    virtual ~Plume();
    //virtual void simulate();
    MPITopology* getMpiTopology();
    // to be re-implemented for class inheriting from MPIKernel
    virtual void createVolcanTopology();

protected:
    Plume( params p, shared_ptr<MPITopology> topology, double U0, double L0 ); // protected -> can not instantiated
    Plume( int argc , char ** argv, MPI_Comm globalCommunicator, double U0, double L0 ); // protected -> can not instantiated

    virtual void InitCommon();

    virtual void InitSimulationReal();

    virtual void InitSimulationTest();

    virtual int InjectParticlesAtCrater();

    virtual int InjectParticlesAlongPlume();

    virtual int InjectParticlesTest();

    virtual void  initLogger(string name);

    std::function<int(void)> InjectParticles;
    std::function<void(void)> InitSimulation;

    //virtual void WimSimulation();

    /**
     * @brief F_init initalizes the computation (F_init operator in MML).
     */
    virtual void F_init();

    /**
     * @brief computeIteration runs one iteration of the model (operator S in MML)
     * @return
     */
    virtual void S();

    //virtual void U(); Not inplemented here

    /**
     * @brief Oi makes an intermediate observation (Operator Of in MML)
     */
    virtual void Oi();
    /**
     * @brief Of makes a final observation (Operator Of in MML)
     */
    virtual void Of();

    /**
     * @brief isConverged tells to stop or not the computation
     */
    virtual bool isConverged();

    /**
     * @brief serializes serializes a vector collection to a vector of char. It resize the toSend vector anc copy data into it.
     * @param boundaryParticles
     * @param toSend vector<char>
     */
    virtual void serialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & toSend);

    /**
     * @brief deserialize constructs a std::vector< BoundaryParticle > from a given vector<char>
     * @param boundaryParticles
     * @param received vector<char>
     */
    virtual void deserialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & received);

    virtual void verifyDomainContainsParticles();


    std::vector< std::shared_ptr< SpeedFunctor > > createWimSpeeds( EruptionParameters ep );

    std::vector< std::shared_ptr< SpeedFunctor > > createWoodsSpeeds( EruptionParameters ep, bool isoDif );

    int computeParticleIncrement(int i);

protected:

    MPIParticleTracker* tracker_;

    int argc;
    char **argv;
    MPI_Comm globalCommunicator;

    std::random_device rd_;
    std::mt19937 mt_;

    params p_;
    shared_ptr<MPITopology> topology_;
    MPIAbstractEulerianDomain* dom_;
    MPIParticleRepository* repository_;
    MPIGridTerrain* terrain_;
    double volatile currentTime_;

    shared_ptr<VolcanoMPIFileManager> fm;
    shared_ptr<EruptionParameters> eruptionParameters;
    // EruptionParameters eruptionParameters;

    int  totalParticlePerFamily;
    std::vector<double> scalingFactors;
    std::vector<int> totalParticles;
    int totalParticle;
    MPIAbstractSimulator *sim;

    // Variables used to perform time measurments
    boost::posix_time::ptime timeStart;
    boost::posix_time::ptime timeEnd;
    boost::posix_time::time_duration timeDif;
    int cpt;
    bool isProcess;
    int frequenceCheckConvergence; // after itrToCheckConvergence iteration check whether they are still particles in the domain
#ifdef PROFILING
    unige_pasc::Profiler<double> profiler;
    int nbrItrToSaveStatsInFile;
#endif
    int previous_particules_deposees;
    int injectedParticles;
    int sum_injectedParticles, sum_outsideParticles;

    //
    //Logger logger;
    shared_ptr<plb_ofstream> logger;

    double U0_;
    double L0_;

};

/*********************************** Continent *****************************************************/
class Continent  : public Plume{ // a modified version of Plume ?
// TO complete by Pierre
public:

    virtual ~Continent();

    virtual void verifyDomainContainsParticles();

protected:
    Continent(params p, shared_ptr<MPITopology> topology, double U0, double L0 ); // protected -> can not instantiated
    Continent(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
   //

   /**
    * @brief F_init initalizes the computation (F_init operator in MML).
    */
   virtual void F_init();

   /**
    * @brief computeIteration runs one iteration of the model (operator S in MML)
    * @return
    */
   virtual void S();

   //virtual void U(); Not inplemented here

   /**
    * @brief Oi makes an intermediate observation (Operator Of in MML)
    */
   virtual void Oi();
   /**
    * @brief Of makes a final observation (Operator Of in MML)
    */
   virtual void Of();

   /**
    * @brief isConverged tells to stop or not the computation
    */
};

/***************************************** VolanicConnector ***********************************************/

class VolanicConnector{
public:
    virtual ~VolanicConnector(){}
    virtual void send(vector<char> const& data)=0;
    virtual void receive(vector<char> & data)=0;
    virtual int getSubmodel_ID()=0;
    virtual void end()=0;


};

/***************************************** LocalConnector ***********************************************/

class LocalVolcanicConnector :public VolanicConnector, public Connector<LocalVolcanicConnector>
{
  friend class Connector<LocalVolcanicConnector>;

public:
   virtual void send( vector<char> const &buffer);
   virtual void receive( vector<char> & buffer);
   virtual  void init ();
   virtual int getSubmodel_ID();
   virtual void end();
private:

    LocalVolcanicConnector();
    LocalVolcanicConnector(const LocalVolcanicConnector&){}
    ~LocalVolcanicConnector();
    bool isInitialized;

};
///


/**************************************** LocalPlume ************************************************/
class LocalPlume : public Plume {
public:
    LocalPlume(params p, shared_ptr<MPITopology> topology, double U0, double L0 );
    LocalPlume(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~LocalPlume();
    void mainLoop();
protected:
    void U();

    //mohamed: monolithic implementation of the abstract method
private:
    LocalVolcanicConnector * connector;

};

/**************************************** LocalContinent ************************************************/
class LocalContinent : public Continent {
public:
    LocalContinent( params p, shared_ptr<MPITopology> topology, double U0, double L0 );
    LocalContinent( int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0 );
    virtual ~LocalContinent();
    void mainLoop();
protected:
    void U();
    //mohamed: monolithic implementation of the abstract method
private:
    LocalVolcanicConnector * connector;
};

/************************************** DistributedPlume **************************************************/
class DistributedPlume: public Plume , public MapperHeader{

public:
        virtual ~DistributedPlume();
protected:
    DistributedPlume(params p, shared_ptr<MPITopology> topology, double U0, double L0);
    DistributedPlume(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual bool isConverged();// redfine
    virtual void synchronizeCommunication();
    virtual void U();
    virtual void F_init();
    virtual void Of();

protected:
    void serialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & toSend);
    void deserialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & received);
    //abstract functions
    virtual void send(vector<char> const& data)=0;
    virtual void receive(vector<char> & data)=0;
    virtual void sendConvergence()=0;
    virtual bool receiveConvergence()=0;
    virtual int getSubmodel_ID()=0;
    virtual void end()=0;

protected:
    bool isRemoteConverged;//tells whether the remote submodel is converged
};

/************************************** DistributedPlume **************************************************/
class DistributedAtmospheredPlume: public DistributedPlume{

public:
        virtual ~DistributedAtmospheredPlume();
protected:
    DistributedAtmospheredPlume(params p, shared_ptr<MPITopology> topology, double U0, double L0);
    DistributedAtmospheredPlume( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

    virtual void synchronizeCommunication();
    virtual void U();
    void sendAtmosphereDataRequest();
    void receiveAtmosphereDataResponse(AtmosphereDataResponse & response);

    double lastTime;

protected:
    virtual void sendToAtmosphere(vector<char> const& data)=0;
    virtual void receiveFromAtmosphere(vector<char> & data)=0;

private:
    AtmosphereDataRequest createAtmosphereDataRequest();

};

/************************************** DistributedContinent **************************************************/
class DistributedContinent: public Continent , public MapperHeader{

public:
    DistributedContinent( params p, shared_ptr<MPITopology> topology, double U0, double L0);
    DistributedContinent( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~DistributedContinent();

protected:
    virtual bool isConverged();// redfine
    void synchronizeCommunication();
    void U();
    void F_init();
    void Of();

protected:
    void serialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & toSend);
    void deserialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & received);
    //abstract functions
    virtual void send(vector<char> const& data)=0;
    virtual void receive(vector<char> & data)=0;
    virtual void sendConvergence()=0;
    virtual bool receiveConvergence()=0;
    virtual int getSubmodel_ID()=0;
    virtual void end()=0;

protected:
    bool isRemoteConverged;//tells whether the remote submodel is converged

};

/************************************** DistributedAtmospheredContinent **************************************************/
class DistributedAtmospheredContinent: public DistributedContinent{

public:
        virtual ~DistributedAtmospheredContinent();
protected:
    DistributedAtmospheredContinent(params p, shared_ptr<MPITopology> topology, double U0, double L0);
    DistributedAtmospheredContinent(int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

    virtual void synchronizeCommunication();
    virtual void U();
    void sendAtmosphereDataRequest();
    void receiveAtmosphereDataResponse(AtmosphereDataResponse & response);

    double lastTime;

protected:
    virtual void sendToAtmosphere(vector<char> const& data)=0;
    virtual void receiveFromAtmosphere(vector<char> & data)=0;

private:
    AtmosphereDataRequest createAtmosphereDataRequest();

};

/************************************** LightDistributedPlume **************************************************/
class LightDistributedPlume: public DistributedPlume, public AsynchronousRemoteMpiKernel{

public:

    virtual void mainLoop();
    virtual void createVolcanTopology();
    LightDistributedPlume( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~LightDistributedPlume();

protected:
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    virtual void sendConvergence();
    virtual bool receiveConvergence();
    int  getSubmodel_ID();
    void end();

protected:

    string id;
    ConduitEntrance<char> * f_out;
    ConduitExit<char> * f_in;

};

/************************************** DistributeAtmosphereddPlume **************************************************/
class LightDistributedAtmospheredPlume: public DistributedAtmospheredPlume, public AsynchronousRemoteMpiKernel{

public:

    virtual void mainLoop();
    virtual void createVolcanTopology();
    LightDistributedAtmospheredPlume( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~LightDistributedAtmospheredPlume();

protected:
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    virtual void sendConvergence();
    virtual bool receiveConvergence();
    void sendToAtmosphere(vector<char> const& data);
    void receiveFromAtmosphere(vector<char> & data);
    int  getSubmodel_ID();
    void end();


protected:

    string id;
    ConduitEntrance<char> * f_out, *f_out_atm;
    ConduitExit<char> * f_in, *f_in_atm;

};
/************************************** LightDistributedContinent *********************************************/

class LightDistributedContinent: public DistributedContinent, public AsynchronousRemoteMpiKernel{
public:
    virtual void mainLoop();
    virtual void createVolcanTopology();
    LightDistributedContinent( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~LightDistributedContinent();

protected:
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    virtual void sendConvergence();
    virtual bool receiveConvergence();
    int  getSubmodel_ID();
    void end();


protected:
    string id;
    ConduitEntrance<char> * f_out;
    ConduitExit<char> * f_in;
};

/************************************** DistributeAtmosphereddContinent **************************************************/
class LightDistributedAtmospheredContinent: public DistributedAtmospheredContinent, public AsynchronousRemoteMpiKernel{

public:

    virtual void mainLoop();

    LightDistributedAtmospheredContinent( int  argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
    virtual ~LightDistributedAtmospheredContinent();

protected:
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    virtual void sendConvergence();
    virtual bool receiveConvergence();
    void sendToAtmosphere(vector<char> const& data);
    void receiveFromAtmosphere(vector<char> & data);
    int  getSubmodel_ID();
    void end();


protected:

    string id;
    ConduitEntrance<char> * f_out, *f_out_atm;
    ConduitExit<char> * f_in, *f_in_atm;

};

/*************************************************************************************************/
class TolopogyCreator{

public:
    static TolopogyCreator& getInstance() {
        static TolopogyCreator    instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }

    shared_ptr<MPITopology> createTopology(params & p, int argc, char **argv, MPI_Comm globalCommunicator);
private:
   //riend shared_ptr<MPITopology> createTopology(params & p, int argc, char **argv);
    TolopogyCreator();



};

/*inline  shared_ptr<MPITopology>createTopology(params & p, int argc, char **argv){
    static TolopogyCreator instance;
    return instance.getTopology(p, argc, argv);
}*/

/************************************** FactoryLocalPlume *********************************************/
class FactoryLocalPlume{

public:
        virtual ~FactoryLocalPlume();
        virtual LocalPlume *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
protected:
    //shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

public:
    FactoryLocalPlume();
};

/************************************** FactoryLocalContinent *********************************************/
class FactoryLocalContinent{

public:
        virtual ~FactoryLocalContinent();
        virtual LocalContinent *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);
protected:
   // shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

public:
    FactoryLocalContinent();
};

/************************************** FactorDistributedPlume *********************************************/
class FactoryDistributedPlume{

public:
        virtual ~FactoryDistributedPlume();
        virtual DistributedPlume *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)=0;
protected:
    //shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

protected:
    FactoryDistributedPlume();
};

/////
class FactoryDistributedAtmospheredPlume{

public:
        virtual ~FactoryDistributedAtmospheredPlume();
        virtual DistributedAtmospheredPlume *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)=0;
protected:
    //shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

protected:
    FactoryDistributedAtmospheredPlume();
};

/************************************** FactorDistributedPlume *********************************************/
class FactoryLightDistributedPlume: public FactoryDistributedPlume{

public:
    FactoryLightDistributedPlume();
    virtual ~FactoryLightDistributedPlume();
    virtual DistributedPlume *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

};

/************************************** FactoryLightDistributedAmospheredPlume *********************************************/
class FactoryLightDistributedAtmospheredPlume: public FactoryDistributedAtmospheredPlume{

public:
    FactoryLightDistributedAtmospheredPlume();
    virtual ~FactoryLightDistributedAtmospheredPlume();
    virtual DistributedAtmospheredPlume *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

};
/************************************** FactorDistributedContinent *********************************************/
class FactoryDistributedContinent{

public:
        virtual ~FactoryDistributedContinent();
        virtual DistributedContinent *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)=0;
protected:
    //shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

protected:
    FactoryDistributedContinent();
};
/////
class FactoryDistributedAtmospheredContinent{

public:
        virtual ~FactoryDistributedAtmospheredContinent();
        virtual DistributedAtmospheredContinent *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)=0;
protected:
    //shared_ptr<MPITopology> getTopology(params &p, int argc, char ** argv);

protected:
    FactoryDistributedAtmospheredContinent();
};


/************************************** FactorDistributedContinent *********************************************/
class FactoryLightDistributedContinent: public FactoryDistributedContinent{

public:
    FactoryLightDistributedContinent();
    virtual ~FactoryLightDistributedContinent();
    virtual DistributedContinent *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

};
/************************************** FactoryLightDistributedAmospheredContinent *********************************************/
class FactoryLightDistributedAtmospheredContinent: public FactoryDistributedAtmospheredContinent{

public:
    FactoryLightDistributedAtmospheredContinent();
    virtual ~FactoryLightDistributedAtmospheredContinent();
    virtual DistributedAtmospheredContinent *newInstance(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0);

};
#endif /* Simulation_H_ */
