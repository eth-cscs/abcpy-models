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

#ifndef NETWORKING_H
#define NETWORKING_H


#include <atomic>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include <sstream>      // std::stringstream
#include <cstring>
#include "musclehpc/parallelism/mpiManager.h"

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
#include <netdb.h>
#include <ifaddrs.h>
#include<signal.h>
#include <linux/if_link.h>

#include <errno.h>
#include <limits.h>
#include <fcntl.h>
#include <cstdlib>
#include <netinet/tcp.h>
#include <thread>         // std::thread
#include <mutex>          // std::mutex
using namespace std;

/************************************* enum **********************************************/

enum class MPI_CONDUIT_OP {sendE ,send, iSend, rSend, iSendRequestFree, receiveE, receive,
                           iRecv, sendRecv, sendToMaster,
                           gather, scatterV, gatherV, bCast,
                           bCastString, bCastThroughMaster,
                           reduce, reduceVect, allReduceVect, reduceAndBcast};


enum class MPI_CONDUIT_DATATYPE {CHAR, INT,  LONG, FLOAT, DOUBLE, LONG_LONG, UNSIGNED_LONG_LONG, LONG_DOUBLE};


/************************************* RelayerURL **********************************************/

struct SocketURL{

    SocketURL(string socketurl);
    SocketURL(string ipAddress, int port): ipAddress(ipAddress), port(port){}
    string ipAddress;
    int port;
};

/************************************* Endpoint **********************************************/

struct Endpoint{
    Endpoint();

    Endpoint(string rawEndpoints);

    void parse();
    void parse(string rawEndpoints);
    string rawEndpoints; //IP1;IP2;...;IP3;port;grank
    vector<SocketURL> urls;
    int port;
    int grank; // global MPI rank
};

/************************************* Message **********************************************/

struct Message{

    Message():senderCoreId(0), receiverCoreId(0), operation(0), datatype(0),  MPI_tag(0), endMessage(0){} // 1 means END_MESSAGE
    string senderName;
    string receiverName;
    int senderCoreId, receiverCoreId;
    int operation;
    int datatype;
    int MPI_tag;
    vector<char> data;
    int endMessage;

    void toBytes(vector<char> & buffer);



    void fromByte(vector<char> const  & buffer);
};


/************************************* Networking **********************************************/
class Messaging{

public:

    static Messaging& getInstance() {
        static Messaging    instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }

    virtual ~Messaging();
private:
    Messaging();
    // C++ 11
    // =======
    // We can use the better technique of deleting the methods
    // we don't want.
    Messaging(Messaging const&) ;
    void operator=(Messaging const&)  ;

public:


    int sendMessageHeader(shared_ptr<MpiManager> & mpiManager, string senderName, int senderCoreId, string receiverName,
                          int receiverCoreId, int operation, int dataType, int MPI_tag, int endMessage, int dest);
    void receiveMessageHeader(shared_ptr<MpiManager> & mpiManager, Message &m, int sender, int relayerGlobalTag);

    void send(shared_ptr<MpiManager> & mpiManager, Message & m, int dest);
    void send(shared_ptr<MpiManager> & mpiManager, string senderName, int senderCoreId, string receiverName,
              int receiverCoreId, int operation, int dataType, int MPI_tag, int endMessage, char * data, int count, int dest);
    void iSend(shared_ptr<MpiManager> &mpiManager, string senderName, int senderCoreId, string receiverName, int receiverCoreId,
                         int operation, int dataType, int MPI_tag, int endMessage, char *data, int count, int dest, MPI_Request *request);

    void receive(shared_ptr<MpiManager> & mpiManager, Message &m, int sender, int relayerGlobalTag);
    void iRecv(shared_ptr<MpiManager> & mpiManager, Message &m, int sender, int relayerGlobalTag, MPI_Request *request);



};


/************************************* Networking **********************************************/

class Networking{


public:
    static Networking& getInstance() {
        static Networking    instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }
private:
    Networking();
    // C++ 11
    // =======
    // We can use the better technique of deleting the methods
    // we don't want.
    Networking(Networking const&) ;
    void operator=(Networking const&)  ;
    int chunk_size;
    mutex m_mutex; // to avoid multithreading open socket susing the same port simultaneously

public:

    /**
     * @brief getHostIpAddresses loop over all the network interfaces and retrieve IP addresses.
     * @param ipAdresses vector of IPV4 and IPv6
     */
    void getHostIpAddresses(vector<string> & ipAdresses, bool isIPV4=true, bool isIPV6=false);
    /**
     * @brief connectToFrontEndSocket establish a TCP socket client connection
     * @param frontEndUrl IP:port
     * @return socket descriptor file.
     */
    int connectToFrontEndSocket(SocketURL & frontEndUrl);
    /**
     * @brief connectToFrontEndSocket establish a TCP socket client connection
     * @param frontEndUrl IP1;IP2;...;IPn;port
     * @return socket descriptor file.
     */
    int connectToFrontEndRelayer(Endpoint & myrelayer);
    /**
     * @brief connectToFrontEndForwarder
     * @param forwarder
     * @return
     */
    int connectToFrontEndForwarder(Endpoint & forwarder, string name, string remoteEndpoint, string correspondingForwarding="");

    /**
     * @brief openfrontEndSocket open a TCP Asynchrone Socket server
     * @param frontEnd_TCPUrl IP:port. I empty choose randomly a port and bind to all interfaces + set up frontEnd_TCPUrl parameter.
     * @return socket file descriptor
     */
    int  openfrontEndSocket(Endpoint & endPoint, int rank=0);

    /**
     * @brief openfrontendMPI open an MPI based socket
     * @param globalCommunicator
     * @return
     */
    string openfrontendMPI(MPI_Comm globalCommunicator);
    /**
     * @brief canOpenMPIPort checks if it is possible to open a MPI socket
     * @param globalCommunicator
     * @return
     */
    bool canOpenMPIPort(MPI_Comm globalCommunicator);

   void setnonblocking(int sock);
   void setBlocking(int sock);
    /*
    ** sock_read
    **
    ** Attempt to read COUNT bytes from socket SOCK into BUFFER.
    ** Returns number of bytes read, or -1 if an error occurs.
    ** Will read less than the specified count *only* if the peer sends
    ** EOF
    */

    int sock_read(int sock, void *buffer, int count, int header=-1);

    /*
    ** sock_write
    **
    ** Attempt to write COUNT bytes to socket SOCK from BUFFER.
    ** Returns number of bytes written, or -1 if an error occurs.
    ** Return value will always be either -1 or COUNT
    */

    int sock_write(int sock, const void *buffer, int count, int header=-1);

};

#endif // NETWORKING_H
