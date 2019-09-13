#ifndef DISTRIBUTEDSIMULATION_H
#define DISTRIBUTEDSIMULATION_H

#include   "../communication/musclemapper.h"
#include  "Simulation.hpp"
#include  "atmosphere.h"
using namespace unige_pasc;


/************************************** MuscleDistributedPlume **************************************************/
class MuscleDistributedPlume: public DistributedPlume{

public:
    //MuscleDistributedPlume(params p, shared_ptr<MPITopology> topology, int * argc, char *** argv);
    MuscleDistributedPlume(int  *argc, char *** argv);
    virtual ~MuscleDistributedPlume();

protected:
    virtual void F_init();
    virtual void send(vector<char> const& data);
    virtual void receive(vector<char> & data);
    virtual int  getSubmodel_ID();
    virtual void end();
    void sendConvergence();
    bool receiveConvergence();


protected:
    //used for muscle
    TetrasConnector * connector;
};


/************************************** MuscleDistributedPlumeAtmosphered **************************************************/
class MuscleDistributedAtmospheredPlume:  public DistributedAtmospheredPlume{

public:
    MuscleDistributedAtmospheredPlume(int *argc, char ***argv);
    virtual ~MuscleDistributedAtmospheredPlume();

protected:
   virtual void F_init();
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    void sendToAtmosphere(vector<char> const& data);
    void receiveFromAtmosphere(vector<char> & data);
    int  getSubmodel_ID();
    void end();
    void sendConvergence();
    bool receiveConvergence();


protected:
    //used for muscle
    TetrasDoubleConnector * connector; // connects Plume + (Continent+ Atmosphere)
};


/************************************** MuscleDistributedAtmosphere **************************************************/
class MuscleDistributedAtmosphereSubModel: public DistributedAtmosphereSubModel{

public:
    MuscleDistributedAtmosphereSubModel(double dx, double dt, int * argc, char *** argv, int rank);
    virtual ~MuscleDistributedAtmosphereSubModel();

protected:
    // implement abstract functions
    virtual void send(vector<char> const& data);
    virtual void receive(vector<char> & data);
    virtual int getSubmodel_ID();
    virtual void end();
     virtual void F_init();

protected:
    //used for muscle
    TetrasConnector * connector;
};


/************************************** MuscleDistributedContinent *********************************************/

class MuscleDistributedContinent: public DistributedContinent{
public:
    MuscleDistributedContinent(int *argc, char ***argv);
    virtual ~MuscleDistributedContinent();

protected:
    void F_init();
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    int  getSubmodel_ID();
    void end();
    void sendConvergence();
    bool receiveConvergence();


    //used for muscle
    TetrasConnector * connector;

};

/************************************** MuscleDistributedContinentAtmosphered **************************************************/
class MuscleDistributedAtmospheredContinent:  public DistributedAtmospheredContinent{

public:
    MuscleDistributedAtmospheredContinent(int *argc, char ***argv);
    virtual ~MuscleDistributedAtmospheredContinent();

protected:
    virtual void F_init();
    void send(vector<char> const& data);
    void receive(vector<char> & data);
    void sendToAtmosphere(vector<char> const& data);
    void receiveFromAtmosphere(vector<char> & data);
    int  getSubmodel_ID();
    void end();
    void sendConvergence();
    bool receiveConvergence();


protected:
    //used for muscle
    TetrasDoubleConnector * connector; // connects Continent + (Continent+ Atmosphere)
};

/************************************** FactorMuscleDistributedContinent *********************************************/
class FactorMuscleDistributedPlume: public FactoryDistributedPlume{

public:
    FactorMuscleDistributedPlume();
    virtual ~FactorMuscleDistributedPlume();
    virtual DistributedPlume *newInstance(int argc, char **argv);

};

/************************************** FactorMuscleDistributedContinent *********************************************/
class FactoryMuscleDistributedAtmospheredPlume: public FactoryDistributedPlume{

public:
    FactoryMuscleDistributedAtmospheredPlume();
    virtual ~FactoryMuscleDistributedAtmospheredPlume();
    virtual DistributedAtmospheredPlume *newInstance(int argc, char **argv);

};


/************************************** FactorMuscleDistributedContinent *********************************************/
class FactorMuscleDistributedContinent: public FactoryDistributedContinent{

public:
    FactorMuscleDistributedContinent();
    virtual ~FactorMuscleDistributedContinent();
    virtual DistributedContinent *newInstance(int argc, char ** argv);

};

/************************************** FactorMuscleDistributedContinent *********************************************/
class FactoryMuscleDistributedAtmospheredContinent: public FactoryDistributedContinent{

public:
    FactoryMuscleDistributedAtmospheredContinent();
    virtual ~FactoryMuscleDistributedAtmospheredContinent();
    virtual DistributedAtmospheredContinent *newInstance(int argc, char ** argv);

};



#endif // DISTRIBUTEDSIMULATION_H
