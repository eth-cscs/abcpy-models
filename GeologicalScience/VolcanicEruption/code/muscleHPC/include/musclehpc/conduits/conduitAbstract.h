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


#ifndef CONDUITABSTRACT_H
#define CONDUITABSTRACT_H


#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/mediator/networking.h"
using namespace std;
namespace unige_pasc{

/*** FOEWARD DECLARATION **/

class Kernel;
class MpiKernel;
class MMSF;



/************************************ VirtualMPICommunicator<T> ************************************************/
template<typename T>
class VirtualMPICommunicator{

protected:
    VirtualMPICommunicator(){}

public:
    virtual ~VirtualMPICommunicator(){}

    /// Returns the number of processes
    virtual int getSize()=0;
    /// Returns the process ID
    virtual int getRank() const=0;
    /// Returns process ID of main processor
    virtual int bossId() const=0;
    /// Tells whether current processor is main processor
    virtual bool isMainProcessor() const=0;
    /// Returns universal MPI-time in seconds
    virtual double getTime() const=0;
    /// Returns the global communicator for this program or library instance.
    virtual MPI_Comm getGlobalCommunicator()=0;
    /// Synchronizes the processes
    virtual void barrier()=0;
    /// Send and receive data between two partners
    virtual void sendRecv( T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = -1 )=0;


    // --- collective ---

    /// Implementation code for Gather
    virtual void gather(T* sendBuf, int sendCount, T* recvBuf, int recvCounts, int root)=0;

    /// Scatter data from one processor over multiple processors
    // virtual void scatterV( T *sendBuf, T *recvBuf, int* sendCounts, int root = 0 )=0;
    // virtual void scatterV( T *sendBuf, int* sendCounts, int* displs, T* recvBuf, int recvCount, int root);


    /// Gather data from multiple processors to one processor
    virtual void gatherV(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs, int root=0)=0;

    /// Broadcast data from one processor to remote processors
    virtual void bCast( T* sendBuf, int sendCount, int root = 0 )=0;

    /// Broadcast data from one processor from  remote processors
    ///virtual void bCastFrom( T* sendBuf, int sendCount, int root = 0 )=0;


    /// Reduction operation toward one processor
    virtual  void reduce( T sendVal, T& recvVal, MPI_Op op, int root = 0 )=0;

    /// Element-per-element reduction of a vector of data
    virtual void reduceVect( std::vector<T>& sendVal, std::vector<T>& recvVal, MPI_Op op, int root = 0 )=0;

    /// Inplace element-per-element reduction of a vector of data; result
    ///   available on all MPI threads.
    /*virtual void allReduceVect( std::vector<T>& sendRecvVal, MPI_Op op )=0;

    /// Reduction operation, followed by a broadcast
    virtual  void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0 )=0;*/

    /// Complete a non-blocking MPI operation
    virtual  void wait(MPI_Request* request, MPI_Status* status)=0;
};

/************************************ conduitCommon ************************************************/
class ConduitCommon
{
public:
    ConduitCommon(){}
    virtual ~ConduitCommon(){}

    /**
     * @brief getName
     * @return gets the conduit name.
     */
    virtual std::string getName()=0;

    /**
     * @brief isConduitSafe checks if the sender ans receiver kernels can exchange safelfy data through this confuit.
     *         For example, if the conduit is of type PTP, both of sender and receiver should run on the same number of core.
     * @return true means conduit is safe, false otherwise.
     */
    virtual bool isConduitSafe()=0;

    /**
     * @brief getHandlerKernel gets the active submodel instance which handles the conduit.
     * @returnKernel *
     */
    virtual Kernel * getHandlerKernel()const=0;
    /**
     * @brief getRemoteKernel gets the inactive remote subnmodel
     * @return
     */
    virtual Kernel * getRemoteKernel()const=0;
    /**
     * @brief setHandler  sets the active submodel instance which handles the conduit.
     * @param kernelHoldingTheConduit
     */
    virtual void setHandler(Kernel * kernelHoldingTheConduit)=0;
    /**
     * @brief setTag
     * @param conduitTag
     */
    virtual void setTag (int conduitTag)=0;
    /**
     * @brief getConduitTag
     * @return
     */
    virtual int getConduitTag()=0;
    /**
     * @brief setRelayerURL
     * @param relayerUrl
     */
    virtual void setRelayerURL(shared_ptr<Endpoint> & relayerUrl)=0;
    /**
     * @brief getHandlerInerCommunicatortManager
     * @return
     */
    virtual shared_ptr<MpiManager> getHandlerInerCommunicatortManager()=0;

protected:

};



/************************************ ConduitEntrance<E> ************************************************/
template<typename E>
class ConduitEntrance:  public virtual ConduitCommon, public virtual VirtualMPICommunicator<E>{
public:
    virtual ~ConduitEntrance(){}
    /**
     * @brief send data into a conduit. this is not a blocking method.
     * @param data
     */
    virtual bool send(E * data, int count)=0;

    /**
     * @brief stop propagates a stop signal to a conduit.
     */
    virtual  void stop()=0;

    //--- SEND ---
    virtual  void send( E *buf, int count, int dest, int tag = -1 )=0;

    /// Sends data at *buf, non blocking
    virtual  void iSend( E *buf, int count, int dest, MPI_Request* request, int tag = -1 )=0;

    /// Sends data at *buf, assuming that receiver is ready.
    virtual void rSend( E *buf, int count, int dest, int tag = -1 )=0;

    /// Sends data at *buf, non blocking and request free
    virtual  void iSendRequestFree( E *buf, int count, int dest, int tag = -1 )=0;

    /// Sends data to master processor
    virtual  void sendToMaster( E* sendBuf, int sendCount, bool iAmRoot )=0;


};

/************************************ conduitExit<F> ************************************************/
template<typename F>
class ConduitExit: public virtual  ConduitCommon, public virtual VirtualMPICommunicator<F>{
public:
    virtual ~ConduitExit(){}
    /**
     * @brief receive data from a conduit. this is a blocking method.
     * @param data
     */
    virtual F * receive(int & count)=0;

    //-- RECV ---
    /// Receives data at *buf, blocking
    virtual  void receive( F *buf, int count, int source, int tag = -1 )=0;

    /// Receives data at *buf, non blocking
    virtual void iRecv( F *buf, int count, int source, MPI_Request* request, int tag = -1 )=0;

};

/************************************ Conduit<E,F> ************************************************/
/**
 * Each conduit has a sender and a receiver of type Kernel.
 */
template<typename E>
class Conduit: public virtual ConduitEntrance <E>, public virtual ConduitExit <E>{
public:
    virtual ~Conduit(){}//<- important to avoid  memory leak
    virtual void setHandler(Kernel * kernelHoldingTheConduit)=0;
    virtual Kernel * getHandlerKernel()const=0;
    virtual Kernel * getRemoteKernel()const=0;
    virtual  void registerSender(Kernel * s)=0;
    virtual void registerReceiver(Kernel * s)=0;
    virtual void unregister(Kernel * s)=0;
    virtual void initialize(string id)=0;
    virtual bool send(E * data, int count)=0;
    virtual void stop()=0;
    virtual E * receive(int & count)=0;
    //  virtual void receive(E * data , int & count)=0;
    virtual string getName()=0;
    virtual void setTag (int conduitTag)=0;


    /**
     * @brief isConduitSafe checks if the sender ans receiver kernels can exchange safelfy data through this confuit.
     *         For example, if the conduit is of type PTP, both of sender and receiver should run on the same number of core.
     * @return true means conduit is safe, false otherwise.
     */
    virtual bool isConduitSafe()=0;

};

}// end namespace
#endif // CONDUITABSTRACT_H
