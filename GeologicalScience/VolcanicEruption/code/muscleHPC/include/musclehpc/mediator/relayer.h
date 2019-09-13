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

#ifndef MANAGER_H
#define MANAGER_H

#include <cstdint>
#include <iostream>
#include <assert.h>
#include <atomic>
#include <memory>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <signal.h>
#include <exception>      // std::exception

#include <sys/ioctl.h>
#include <stdio.h>

#include "musclehpc/mediator/networking.h"
#include "musclehpc/utils/optionparser.h"
#include "musclehpc/parallelism/mpiManager.h"





/************************************* RelayerCommunicator **********************************************/

class RelayerCommunicator{

public:
    virtual ~RelayerCommunicator(){}

    virtual shared_ptr<MpiManager> getInterMpiManager()=0;


    virtual bool sendInt(int value);
    virtual bool receiveInt(int & value);
    virtual bool sendString(string message);
    virtual bool receiveString(string & message);
    virtual void sendMessage(Message & m, int header=-1)=0;
    virtual void receiveMessage(Message & m, int header=-1)=0;


    virtual bool bCastIntTo(int  value, int root);
    virtual bool bCastIntFrom(int & value , int root);
    virtual bool bCastStringTo(string  message, int root);
    virtual bool bCastStringFrom(string & message , int root);
    virtual bool bCastUrlTo(Endpoint & relayerURL, int root);
    virtual bool bCastUrlFrom(Endpoint & relayerURL , int root);


    //absrtract
    virtual bool sendRaw( char* data, int length, int value)=0;
    virtual bool send(vector<char> & buffer, int header=-1)=0;
    virtual bool receive(vector<char> & buffer, int header=-1)=0;
    virtual bool bCastTo(vector<char> & buffer, int root)=0;
    virtual bool bCastFrom(vector<char> & buffer, int root)=0;
    virtual bool disconnect()=0;

};

/************************************* MPI_Client_RelayerCommunicator **********************************************/

class MPI_Client_RelayerCommunicator: public RelayerCommunicator{

private:
    shared_ptr<MpiManager> myMpiManager;
    shared_ptr<MpiManager> interMpiManager;
    string managerPortUrl;
    int myRank;

    //char remotePortName[MPI_MAX_PORT_NAME];
public:
    //MPI_Client_RelayerCommunicator(string managerPortUrl, shared_ptr<MpiManager> & myMpiManager);
    MPI_Client_RelayerCommunicator(shared_ptr<MpiManager> &myMpiManager, shared_ptr<MpiManager>& interMpiManager);

    virtual ~MPI_Client_RelayerCommunicator();

    virtual shared_ptr<MpiManager> getInterMpiManager();
    //virtual bool initConnection();

    virtual void sendMessage(Message & m, int header=-1);
    virtual void receiveMessage(Message & m, int header=-1);

    virtual bool sendRaw( char* data, int length, int value);
    virtual bool send(vector<char>  &buffer, int header=-1);
    virtual bool receive(vector<char> & buffer, int header=-1);
    virtual bool bCastTo(vector<char> & buffer, int root);
    virtual bool bCastFrom(vector<char> & buffer, int root);
    virtual bool disconnect();
};
/************************************* TCP_Client_RelayerCommunicator **********************************************/

class TCP_Client_RelayerCommunicator: public RelayerCommunicator{

private:
    int sockfd;
    string managerPortUrl;
    shared_ptr<MpiManager> myMpiManager;
public:
    TCP_Client_RelayerCommunicator(string managerPortUrl, shared_ptr<MpiManager> & myMpiManager);
    TCP_Client_RelayerCommunicator(int  sockfd, shared_ptr<MpiManager> & myMpiManager);

    virtual ~TCP_Client_RelayerCommunicator();

    virtual shared_ptr<MpiManager> getInterMpiManager();
   // virtual bool initConnection();

    virtual void sendMessage(Message & m, int header=-1);
    virtual void receiveMessage(Message & m, int header=-1);

    virtual bool sendRaw( char* data, int length, int value);
    virtual bool send(vector<char> &buffer, int header=-1);
    virtual bool receive(vector<char> & buffer, int header=-1);
    virtual bool bCastTo(vector<char> & buffer, int root);
    virtual bool bCastFrom(vector<char> & buffer, int root);
    virtual bool disconnect();

    //
    void unsetSocket();
};


/*********************************************************************************************/
/************************************** Relayer **********************************************/
/*********************************************************************************************/


struct MsgContainer{

    MsgContainer():isHeaderSent(false), leftToreceive(0){}
    bool isHeaderSent; // mpi rank src
    int leftToreceive;  // how much left to complete receiving a message from mpi_src_rank
    Message m;
};



class Relayer{

public:
    virtual ~Relayer(){}
    //virtual bool connectWithMyKernels()=0;
    virtual void forwardMessage(shared_ptr<MsgContainer> & msgCtn, shared_ptr<RelayerCommunicator> & relayer)=0;
 virtual shared_ptr<MsgContainer> receiveLocalMessage()=0;
    virtual void run()=0;
    virtual int getBackEndCore()const=0;
    virtual int getFrontEndCore()const=0;


protected:
    virtual void waitForCommingActionOnFrontEnd()=0;
    virtual bool closefrontEndSocket()=0;
};

/************************************* MPIRelayer **********************************************/
// MPI backend + TCP fontend relayer
class MPIRelayer: public Relayer{

private:
    shared_ptr<MpiManager>  globalMpiManager;
    shared_ptr<MpiManager>  myMpiManager;

    int backendRank, frontendRank;
    int frontEnd_sockfd; //my TCP sockfd
    Endpoint relayerEndPoint;
    fd_set active_fd_set;
    int max_sd;
    Endpoint forwarderAddress;

    int manager_fd; // sockfs to connect to remote manager

    string managerPortUrl;
    vector<Relayer> remoteRelayers;
    map<string, shared_ptr<RelayerCommunicator> > localKernelMap;//<kernelName, RelayerCommunicator>
    map<string, int > remoteKernelSocketsMap;//<kernelName, socketId>
    map<string, shared_ptr<RelayerCommunicator>  > remoteRelayerMap;//<relyerUrl, TCPRelayerCommunicator>
    map<string, int> mapColorWhereToRun; //<kernelname, color>
    shared_ptr<RelayerCommunicator> managerRelayerCommunicator;
    MPI_Request request;
    map<int, shared_ptr<MsgContainer> > messagesContainerMap; //< mpi_rank_src, msgcontainer>

public:
    MPIRelayer(shared_ptr<MpiManager>  globalMpiManager, shared_ptr<MpiManager>  &myMpiManager,
               map<string, shared_ptr<RelayerCommunicator> > & localKernelMap,
               map<string, int> &mapColorWhereToRun, string forwarderAddress);
    virtual void run();
    //virtual bool connectWithMyKernels();
    virtual void forwardMessage(shared_ptr<MsgContainer> &msgCtn, shared_ptr<RelayerCommunicator>& relayer);
    virtual shared_ptr<MsgContainer> receiveLocalMessage();

    virtual int getBackEndCore()const;
    virtual int getFrontEndCore()const;

protected:
    virtual bool receiveRemoteMessage(shared_ptr<RelayerCommunicator> &relayerCommunicator, Message & m);
    virtual void waitForCommingActionOnFrontEnd();
    virtual bool closefrontEndSocket();
    virtual void receiveRelayersList(map<string, Endpoint> &relayersMap, map<string, Endpoint> &fwdMap);

    virtual int connectToForwarder(string const &relayerEndpoint);

    bool isRelayerBackend();
    bool isRelayerFrontEnd();

};

#endif // MANAGER_H
