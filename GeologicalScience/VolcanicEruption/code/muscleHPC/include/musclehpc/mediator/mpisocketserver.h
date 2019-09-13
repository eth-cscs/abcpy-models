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

#ifndef MPISOCKETSERVER_H
#define MPISOCKETSERVER_H


#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/mediator/relayer.h"

// for socket C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <string.h> // for memset
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <atomic>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include <sstream>      // std::stringstream
#include <cstring>




using namespace std;
//using namespace unige_pasc;


/************************************ Manager ************************************************/
class Manager
{
public:
    Manager(int kernelsSize,
            map<string, shared_ptr<MpiManager> > & interMpiManagersMap,
                map<string, int> &mapColorWhereToRun,
                shared_ptr<MpiManager> &localMpiManager,
                shared_ptr<MpiManager> &globalMpiManager,
            string forwarder,
            bool useTCP=false);

    virtual ~Manager();
    void run();

private:
   // ColorInfo * infoColor;
    int kernelsSize;
   // vector<ColorInfo *> &colorsInfoMap;
    map<string, shared_ptr<MpiManager> > interMpiManagersMap;
    map<string, int> &mapColorWhereToRun;
    shared_ptr<MpiManager> &localMpiManager;
    shared_ptr<MpiManager> &globalMpiManager;
   int const MPI_TAG_INIT_MANAGER_URL = 50;
   string forwarder;
   bool useTCP;
};


/************************************ KernelConnectionInfo *************************************/
enum class HANDSHAKE_TAG: int{ CONNECT=1, LISTEN=2, END=3 }; // scoped

//************************************** InterCommUtils *************************************************/

class InterCommUtils{


private:
    InterCommUtils();

public:
    ~InterCommUtils();
    void broadCastStringToRemote(string  message, int root, shared_ptr<MpiManager> &m);
    void broadCastStringFromRemote(string & message, int root, shared_ptr<MpiManager> &m);
    friend InterCommUtils& getInterCommUtils();
};

inline InterCommUtils& getInterCommUtils() {
    static InterCommUtils instance;
    return instance;
}


/************************************ KernelConnectionInfo *************************************/
class KernelConnectionInfo{

protected:
    string name;
    string mpiUrl;
    string myRelayerUrl;
    string myForwarderUrl;
    int color; // just to enumerate kernels
    typedef typename std::map<string, shared_ptr<KernelConnectionInfo> >::iterator it_map_comm;


public:
    KernelConnectionInfo();
    virtual ~KernelConnectionInfo();

    void setName(string name);
    void setMyMPIUrl(string mpiUrl);
    void setMyRelayerUrl(string mpiUrl);
    void setMyForwarderUrl(string fwdUrl);
    void setcolor(int color);
    string getName() const;
    string getMyMPIUrl() const;
    string getMyRelayerUrl() const;
    string getMyForwarderUrl() const;
    int getcolor() const;

    void receiveInfo();
    void sendInfo();
    virtual void broadCastMyInfoToRemote( int root);
    virtual void broadCastMyInfoFromRemote( int root);

    virtual void broadCastMapToRemote(map<string, shared_ptr<KernelConnectionInfo> >  & remoteKernelsMap, int root);
    virtual void broadCastMapFromRemote(map<string, shared_ptr<KernelConnectionInfo> > & remoteKernelsMap, int root);
    virtual void broadCastInfoToRemote(shared_ptr<KernelConnectionInfo>  & info, int root);
    virtual void broadCastInfoFromRemote(shared_ptr<KernelConnectionInfo>  & info, int root);

    virtual void moveInterMpiManager(shared_ptr<MpiManager> & interMpiManager)=0;
    virtual shared_ptr<MpiManager> getInterMpiManager()=0;
    virtual shared_ptr<KernelConnectionInfo>  createKernelConnectionInfo();



    virtual bool send(string & message)=0;
    virtual bool receive( string & message)=0;
    virtual bool send(int & message)=0;
    virtual bool receive( int & message)=0;

    virtual void broadCastStringToRemote(string & message, int root)=0;
    virtual void broadCastStringFromRemote(string & message, int root)=0;
    virtual void broadCastIntToRemote(int & length, int root)=0;
    virtual void broadCastIntFromRemote(int & length, int root)=0;
};

///
class KernelMPIConnectionInfo: public KernelConnectionInfo{
private:

    shared_ptr<MpiManager> interMpiManager;
    int conduit_MPI_Tag;
    int receiverRank;
    int senderRank;



public:

    KernelMPIConnectionInfo();
    KernelMPIConnectionInfo(shared_ptr<MpiManager>& interMpiManager);
    virtual ~KernelMPIConnectionInfo();


   virtual void moveInterMpiManager(shared_ptr<MpiManager> & interMpiManager);
   virtual shared_ptr<MpiManager> getInterMpiManager();
  // virtual shared_ptr<KernelConnectionInfo>  createKernelConnectionInfo();

    virtual bool send(string & message);
    virtual bool receive( string & message);
    virtual bool send(int & value);
    virtual bool receive( int & value);

    virtual void broadCastStringToRemote(string & message, int root);
    virtual void broadCastStringFromRemote(string & message, int root);
    virtual void broadCastIntToRemote(int & length, int root);
    virtual void broadCastIntFromRemote(int & length, int root);
};

/////
class KernelTCPConnectionInfo: public KernelConnectionInfo{

private:

    int sockfd;
    shared_ptr<MpiManager> localMpiManager;

public:

    KernelTCPConnectionInfo();
    KernelTCPConnectionInfo(int sockfd, shared_ptr<MpiManager> &localMpiManager);
    virtual ~KernelTCPConnectionInfo();


   virtual void moveInterMpiManager(shared_ptr<MpiManager> &interMpiManager);
   virtual shared_ptr<MpiManager> getInterMpiManager();
   //virtual shared_ptr<KernelConnectionInfo> createKernelConnectionInfo();

    virtual bool send(string & message);
    virtual bool receive( string & message);
    virtual bool send(int & value);
    virtual bool receive( int & value);

    virtual void broadCastStringToRemote(string & message, int root);
    virtual void broadCastStringFromRemote(string & message, int root);
    virtual void broadCastIntToRemote(int & length, int root);
    virtual void broadCastIntFromRemote(int & length, int root);
};

//********************************** MpiSocketServer *******************************************
class SocketServer
{
protected:

    map<string, shared_ptr<KernelConnectionInfo>  > remoteKernelsMap; // <kernelID, KernelConnectionInfo* >
    vector<string> remoteKernelsIDs;

    char port_name[MPI_MAX_PORT_NAME];
    typedef typename std::map<string, shared_ptr<KernelConnectionInfo> >::iterator it_map_comm_string_ConnInfo;
    typedef typename std::map<int, int>::iterator it_map_comm_int_int;

    bool amImanager;
    int kernelNumbers;
    shared_ptr<KernelConnectionInfo>  myLocalInfo;

public:
    SocketServer();
    virtual ~SocketServer();

    void setIamManager(bool kernelNumbers);
    void setKernelNumbers(int  isManager);
    map<string, shared_ptr<KernelConnectionInfo> > & getRemoteKernelsMap();
    shared_ptr<KernelConnectionInfo>  getMyLocalInfo();

    virtual void initConnection()=0;
    virtual void handShake()=0;
    virtual void disconnect()=0;

    virtual shared_ptr<MpiManager>& getMpiManager() =0;
    virtual void setMyRemoteKernels(vector<string> myRemoteKerelsIDs)=0;


protected:
    virtual int openMPISocket()=0;
    virtual void closeMPISocket()=0;
    virtual shared_ptr<KernelConnectionInfo> acceptConnection()=0;

};

//********************************** SocketManager *******************************************

class SocketManager: public SocketServer{

public:
    virtual ~SocketManager();
    virtual void initConnection();
    virtual void handShake();
    virtual void disconnect();
    virtual void setMyRemoteKernels(vector<string> myRemoteKerelsIDs);
protected:
       SocketManager(shared_ptr<MpiManager> & localMpiManager, string forwarder="");
       virtual void acceptFromAllKernels( vector<shared_ptr<KernelConnectionInfo> > & queuedConnection, int numberOfConnections)=0;
       virtual void sendRemoteOrder(HANDSHAKE_TAG order, shared_ptr<KernelConnectionInfo> &clientinfo);
       virtual void informAllmyCamaradesKernelsOfMyURL();
       virtual void handShakeWithRelayers();

protected:
    map<string, int> mapColorWhereToRun;
    shared_ptr<MpiManager> myMpiManager;
    map<string, shared_ptr<MpiManager> > sameBinaryInterMpiManager;
    Endpoint forwarder;

};

//********************************** MpiSocketManager *******************************************
class MpiSocketManager: public SocketManager{

public:
    MpiSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, int> mapColorWhereToRun);
    MpiSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, shared_ptr<MpiManager> > interMpiManagers, map<string, int> mapColorWhereToRun);
    virtual ~MpiSocketManager();
    virtual shared_ptr<MpiManager>& getMpiManager() ;

protected:
    virtual void acceptFromAllKernels(vector<shared_ptr<KernelConnectionInfo> > & queuedConnection, int numberOfConnections);
    virtual int openMPISocket();
    virtual shared_ptr<KernelConnectionInfo> acceptConnection();
    virtual void closeMPISocket();
};

//********************************** TCPSocketManager *******************************************
class TCPSocketManager: public SocketManager{

    //int sock;
    string tcpUrl;
    int sockfd;
    int portno;
    struct sockaddr_in serv_addr, cli_addr;

public:
    TCPSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, int> mapColorWhereToRun, string forwarder);
    TCPSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, shared_ptr<MpiManager> > interMpiManagers, map<string, int> mapColorWhereToRun, string forwarder);
    virtual ~TCPSocketManager();
    virtual shared_ptr<MpiManager>& getMpiManager() ;
protected:
    virtual void acceptFromAllKernels( vector<shared_ptr<KernelConnectionInfo> > & queuedConnection, int numberOfConnections);
    virtual int openMPISocket();
    virtual shared_ptr<KernelConnectionInfo> acceptConnection();
    virtual void closeMPISocket();
};

//********************************** SocketKernel *******************************************

class SocketKernel: public SocketServer{

public:
    virtual ~SocketKernel();
    virtual void initConnection();
    virtual void handShake();
    virtual void disconnect();
    virtual void setMyRemoteKernels(vector<string> myRemoteKerelsIDs);
    shared_ptr<RelayerCommunicator> & getRelayerCommunicator();
    virtual shared_ptr<MpiManager>& getMpiManager() ;

protected:
       SocketKernel(shared_ptr<MpiManager> globalMpiManager,
                    shared_ptr<MpiManager> & localMpiManager,
                    shared_ptr<MpiManager> interMpiManager,
                    shared_ptr<RelayerCommunicator> & relayerCommunicator,
                    string myKernelName, string forwarderAddress, string managerPortUrl="");


       virtual shared_ptr<KernelConnectionInfo> acceptConnection();
       virtual HANDSHAKE_TAG getRemoteOrder();
       virtual void connectWithMyRelayer(string &relayerUrl, string &ForwarderUrl);
       virtual void acceptConnectionfromRemoteKernel();
       virtual void connectToMyRemoteKernels();
       virtual void checkoutManagerPortUrl();

       virtual void closeMPISocket();
       virtual int openMPISocket();

       virtual shared_ptr<KernelConnectionInfo> connectToManager();

protected:


        vector<string> myRemoteKerelsIDs;

        shared_ptr<MpiManager> globalMpiManager;
        shared_ptr<MpiManager> myMpiManager;
        shared_ptr<MpiManager> interMpiManager;
        shared_ptr<RelayerCommunicator> relayerCommunicator;
        string myKernelName;
        Endpoint forwarderAddress;
        string managerPortUrl;

private:
       Endpoint managerFWD;
};

//********************************** MpiSocketKernel *******************************************

class MpiSocketKernel: public SocketKernel{

public:
    MpiSocketKernel(shared_ptr<MpiManager> globalMpiManager, shared_ptr<MpiManager> localMpiManager, string myKernelName, string managerPortUrl, shared_ptr<RelayerCommunicator> &relayerCommunicator, string forwarderAddress);
    MpiSocketKernel(shared_ptr<MpiManager> globalMpiManager, shared_ptr<MpiManager> localMpiManager, string myKernelName, shared_ptr<MpiManager> interMpiManager, shared_ptr<RelayerCommunicator> & relayerCommunicator,string forwarderAddress);
    virtual ~MpiSocketKernel();

protected:

   // shared_ptr<KernelConnectionInfo> connectToManager();
};

//********************************** TCPSocketKernel *******************************************

class TCPSocketKernel: public SocketKernel{

    int sockfd;
public:
    TCPSocketKernel(shared_ptr<MpiManager> globalMpiManager, shared_ptr<MpiManager> localMpiManager, shared_ptr<RelayerCommunicator> &relayerCommunicator, string myKernelName, string managerPortUrl);
    virtual ~TCPSocketKernel();

protected:
   // shared_ptr<KernelConnectionInfo> connectToManager();
};



#endif // MPISOCKETSERVER_H
