#ifndef CONNECTORCOMMON_H
#define CONNECTORCOMMON_H

#include <iostream>
#include <string>

#include <vector>


using namespace std;


namespace unige_pasc{

/******************************************** Connector ********************************************/

template <class T>
class Connector
{

public:
    static T& getInstance(){
        static T    instance = T(); // Guaranteed to be destroyed.Instantiated on first use.
        return instance;
    }

private:
    T& operator= (const T&){}
};

/************************************** BasicConnectorInterface **************************************************/
class BasicConnectorInterface{
public:
    virtual void send(vector<char> const& data, string conduit_out)=0;
    virtual void receive(vector<char> & data, std::string conduit_in)=0;
   virtual void setCommandLineArgs( int *argc, char *** argv)=0;
    virtual void init (int rank)=0;
    virtual void end()=0;
    virtual string getName()=0;
    virtual int getMPIRank() const=0;
    virtual pluint getSubmodel_ID()=0;
    virtual pluint getsubmodels_NBR()=0;
    virtual ~BasicConnectorInterface(){}
protected:


};

/************************************** BasicSmartConnectorInterface **************************************************/
class BasicSmartConnectorInterface{
public:


    virtual void sendToCoarse(vector<char> const&buffer)=0;
    virtual void receiveFromCoarse( vector<char> & buffer)=0;
    virtual void sendToFine(vector<char> const&buffer)=0;
    virtual void receiveFromFine( vector<char> & buffer)=0;
    virtual void receiveFromFin1(vector<char>  &buffer)=0;
    virtual void receiveFromFin2( vector<char> & buffer)=0;
    virtual void setUpCoarseAndFineConduits(double dt1, double dt2)=0;

    virtual ~BasicSmartConnectorInterface(){}
protected:

};


}


#endif // CONNECTORCOMMON_H
