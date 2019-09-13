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


#ifndef CONDUIT_H
#define CONDUIT_H

#include <iostream>
#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/conduits/kernels.h"



using namespace std;
using namespace unige_pasc;

namespace unige_pasc{

/************************************ BasicConduit <E> ************************************************/
template<typename E>
class BasicConduit: public Conduit <E>{

protected:
    string ID; //conduit identifier
    Kernel * sender;
    Kernel * receiver;
    std::atomic<bool> stopped;

public:
    BasicConduit(){
        stopped=false;
        sender=0;
        receiver=0;
    }
    virtual ~BasicConduit(){
    }

    /**
     * @brief setHandler set the active kernel handling the conduit since the latter is shared between two kernels.
     * @param kernelHoldingTheConduit
     */
    virtual void setHandler(Kernel * kernelHoldingTheConduit)=0;

    /**
     * @brief sends the data through the conduit
     * @param data pointer to data (E)
     * @param count size of data
     * @return true if ok, false otherwise.
     */
    virtual bool send(E * data, int count)=0;
    /**
     * @brief receive receives a data from conduit.
     * @param count size of data to receive
     * @return pointer to received data of type E
     */
    virtual E * receive( int & count)=0;

    /**
     * @brief isConduitSafe
     * @return
     */
    virtual bool isConduitSafe()=0;
    /**
     * @brief getHandlerKernel
     * @return the kernel hamdling the conduit
     */
    virtual Kernel * getHandlerKernel()const=0;
    /**
     * @brief getRemoteKernel
     * @return the remote kernel instance associated with the handler kernel on the same conduit.
     */
    virtual Kernel * getRemoteKernel()const=0;
    /**
     * @brief setTag each conduit has a uniq tag. This tag is used to tag data sent/received through this conduit
     * @param conduitTag
     */
    virtual void setTag (int conduitTag)=0;

    /**
     * @brief initialize set the conduit name
     * @param id
     */
    virtual void initialize(string id) {
        this->ID = id;
    }

    /**
     * @brief getName get the conduit name;
     * @return string
     */
    virtual string getName(){
        return ID;
    }


    /**
     * @brief registerSender assigns the sender kernel associated to the current conduit
     * @param k
     */
    virtual void registerSender(Kernel * k){
        if(!sender) {
            sender = k;
        } else {
            cout<<"The sender kernel "<< k->getName() << " is already registered in conduit " << ID <<endl;
        }
    }
    /**
     * @brief registerReceiver assigns the receiver kernel associated to the current conduit
     * @param k
     */
    virtual void registerReceiver(Kernel * k){
        if(!receiver) {
            receiver = k;
        } else {
            cout<<"The receiver kernel is already registered in conduit " << ID <<endl;
        }
    }
    virtual void unregister(Kernel * k){
       // cout<<"Conduit " << ID << ": UnRegistering kernel " << k->getName() <<endl;
        if( k == sender ) {
            sender = 0;
        } else if( k == receiver ) {
            receiver = 0;
        } else {
            cout<<"Kernel " << k->getName() << " was never registred in conduit " << k->getName()<<endl;
        }

    }
    virtual void stop() {
        cout<<"Conduit " << ID << ": Received stop signal." <<endl;
        stopped=true;
    }


};

/************************************ MPI_PTP_Conduit <E> ************************************************/
/**
 * MPI_Point_To_Point_Conduit class connects two MPI submodels. Each MPI process of rank k can send/receive data
 * with the MPI process of the same rank k of the corresponding remote MPI process.
 */

template<typename E>
class MPI_PTP_Conduit: public BasicConduit <E>{
protected:
    MpiKernel * handler; // the submodel handling the current conduit instance.
    MpiKernel * remote; // the remote submodel to which the conduit is connected also. Note that
                        // this is only a inactive instance and is only used to get useful information.
    shared_ptr<MpiManager>  handlerInerCommtManager; // an MPI manager handling the Inter-communicaton communicator between two MPI subgroups.

    shared_ptr<Endpoint> relayerUrl; // used for MTO couunication when direct communication between two MPI programs is not possible.
    bool useRelayer; // true: a relayer will be used, False: direct MPI communication is used
    MPI_CONDUIT_DATATYPE conduitDataType;

    bool initialized; // true: conduit is initialized, false: not initialized

    int uniqInterCommTag;// a tag used as a reference to get the interCommunicator created between two MPI sub-groups.
    int conduit_MPI_Tag; //When two submodels S1 and S2 ahve several conduits, they share the same intercommunicator between the correspoding submodels.
                         // Each conduit uses an conduit_MPI_Tag to send/recences messages
                         // through the same shared interCommunicator.



public:
    friend class MpiKernel;
    MPI_PTP_Conduit(): BasicConduit<E>(), initialized(false) {
        conduit_MPI_Tag=99; // default value, it can be changed.
        handler=0;
        remote=0;
        useRelayer=false;
        conduitDataType =getTemplateDataType(); // conduit data types (specified in the couling file as argument):  CHAR, INT, FLOAT, DOUBLE, ...
    }

    virtual ~MPI_PTP_Conduit(){
    }


    virtual void setRelayerURL(shared_ptr<Endpoint> & relayerUrl){
        this->relayerUrl=relayerUrl;
        useRelayer=true;
    }


   MPI_CONDUIT_DATATYPE  getTemplateDataType(){

        MPI_CONDUIT_DATATYPE datatype = MPI_CONDUIT_DATATYPE::CHAR;

        if (std::is_same<E, char>::value){
            datatype = MPI_CONDUIT_DATATYPE::CHAR;
        } else if (std::is_same<E, int>::value){
            datatype = MPI_CONDUIT_DATATYPE::INT;
        }else if (std::is_same<E, long>::value){
            datatype = MPI_CONDUIT_DATATYPE::LONG;
        }else if (std::is_same<E, double>::value){
            datatype = MPI_CONDUIT_DATATYPE::DOUBLE;
        }else if (std::is_same<E, float>::value){
            datatype = MPI_CONDUIT_DATATYPE::FLOAT;
        }else if (std::is_same<E, long long>::value){
            datatype = MPI_CONDUIT_DATATYPE::LONG_LONG;
        }else if (std::is_same<E, unsigned long long>::value){
            datatype = MPI_CONDUIT_DATATYPE::UNSIGNED_LONG_LONG;
        }else if (std::is_same<E, long double>::value){
            datatype = MPI_CONDUIT_DATATYPE::LONG_DOUBLE;
        }
        return datatype;
    }


    virtual bool send(E * data, int count){

       assert(this->handler->getGlobalMpiManager());
        if(!useRelayer){
            shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
            shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
            //assert(m!=0);
            m->send<int>(&count, 1, senderMpimanager->getRank(), conduit_MPI_Tag);
            m->send<E>(data, count, senderMpimanager->getRank(),conduit_MPI_Tag);

        }else{
            int dest= this->relayerUrl->grank;
            string senderName = this->handler->getName();
            int senderCoreId= this->handler->getMpiManager()->getRank();
            string receiverName= this->remote->getName();
            int receiverCoreId= senderCoreId;
            int operation = static_cast<int>(MPI_CONDUIT_OP::sendE);
            int dataType= static_cast<int>(this->conduitDataType);
            int countData =  count * sizeof(E);
            int endMessage=0;
            Messaging::getInstance().send(this->handler->getGlobalMpiManager(),
                 senderName,  senderCoreId,  receiverName,  receiverCoreId,  operation,  dataType, conduit_MPI_Tag, endMessage, (char*) data, countData, dest);
        }

        return true;
    }



    virtual E * receive( int & count){
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!useRelayer){
           //assert(m!=0);
           m->receive<int>(&count, 1 , this->handler->getMpiManager()->getRank(),  conduit_MPI_Tag);
           E * e = new E[count];
           m->receive<E>(e, count , this->handler->getMpiManager()->getRank(),  conduit_MPI_Tag);
           return e;
       }else{
           m->receive<int>(&count, 1 , MPI_ANY_SOURCE,  conduit_MPI_Tag);
           E * e = new E[count];
           m->receive<E>(e, count , MPI_ANY_SOURCE,  conduit_MPI_Tag);
           return e;
       }
    }


    virtual void setHandler(Kernel * kernelHoldingTheConduit){

        if (this->sender == kernelHoldingTheConduit){
            this->handler= dynamic_cast<  MpiKernel *> (this->sender);
            this->remote= dynamic_cast<  MpiKernel *> (this->receiver);
        }else if (this->receiver == kernelHoldingTheConduit){
            this->handler= dynamic_cast<  MpiKernel *> (this->receiver);
            this->remote= dynamic_cast<  MpiKernel *> (this->sender);
        }else{
            //std::cerr <<" -> Error: the kernel : "<<kernelHoldingTheConduit->getName()<<" is not associated with the conduit:"<<std::endl;
            return;
        }

    }

    virtual Kernel * getHandlerKernel()const{
        return (Kernel *) this->handler;
    }
    virtual Kernel * getRemoteKernel()const{
        return (Kernel *)this->remote;
    }


    /**
     * @brief isConduitSafe checks if the sender and receiver kernels can exchange safelfy data through this confuit.
     *         For example, if the conduit is of type PTP, both of sender and receiver should run on the same number of core.
     * @return true means conduit is safe, false otherwise.
     */
    virtual bool isConduitSafe(){
        return (this->sender->getCoresOverWhichTorun() == this->receiver->getCoresOverWhichTorun());
    }

    /**
     * @brief setMPIComm_Tag assigns an MPI tag for the current conduit to be used in MPI operation. This is
     * needed since several conduits (between the same two MPI sub-groups) can share the same InterComm env.
     * @param MPI_TAG_COUNTER
     */
    virtual void setTag(int MPI_TAG_COUNTER){
        this->conduit_MPI_Tag=MPI_TAG_COUNTER;
    }

    virtual int getConduitTag(){
        return this->conduit_MPI_Tag;
    }

   virtual shared_ptr<MpiManager> getHandlerInerCommunicatortManager(){
        if(!this->initialized){

            assert(handler);
            assert(remote);
            int handlerLeaderColor = this->handler->getColorWhereItRuns();
            int remoteLeaderColor = this->remote->getColorWhereItRuns();

            uniqInterCommTag=CantorPairingFunction(handlerLeaderColor, remoteLeaderColor);

            this->handlerInerCommtManager=this->handler->existInterCommunicator(uniqInterCommTag);
            assert(handlerInerCommtManager);
        }
        return handlerInerCommtManager;
    }

    // implement virtualMpiConnector methods
    // --- main ---
    /// Returns the number of processes
    virtual int getSize(){
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        return senderMpimanager->getSize();
    }
   /// Returns the process ID
    virtual int getRank() const{
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        return senderMpimanager->getRank();
    }

    /// Returns process ID of main processor
    virtual int bossId() const {
       shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        return senderMpimanager->bossId();
    }
    /// Tells whether current processor is main processor
    virtual bool isMainProcessor() const{
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        return senderMpimanager->isMainProcessor();
    }
    /// Returns universal MPI-time in seconds
    virtual double getTime() const{
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        return senderMpimanager->getTime();
    }
    /// Synchronizes the processes
    virtual void barrier(){
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        senderMpimanager->barrier();
    }

    /// Complete a non-blocking MPI operation
    virtual  void wait(MPI_Request* request, MPI_Status* status){
        shared_ptr<MpiManager> senderMpimanager = this->handler->getMpiManager();
        senderMpimanager->wait(request, status );
    }


    /// Returns the global communicator for this program or library instance.
    virtual MPI_Comm getGlobalCommunicator() {
         shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
         return m->getGlobalCommunicator();
    }


    /// Send and receive data between two partners
    virtual void sendRecv( E *sendBuf, E *recvBuf, int count, int dest, int source, int tag = -1 ){
        if(tag==-1) tag=conduit_MPI_Tag;
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->sendRecv(sendBuf, recvBuf, count, dest, source, tag);
    }

    // --- collective ---

    /// Implementation code for Gather
    virtual void gather (E* sendBuf, int sendCount, E* recvBuf, int recvCounts, int root){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
         if(!this->useRelayer){
            m->gather(sendBuf, sendCount, recvBuf, recvCounts, root);
        }else{

            this->handler->getMpiManager()->gather(sendBuf, sendCount, recvBuf, recvCounts, root);

            if(this->handler->getMpiManager()->isMainProcessor()){
                int receiverCoreId= 0;
                this->forwardMessage(MPI_CONDUIT_OP::gather, recvCounts, recvBuf, receiverCoreId, conduit_MPI_Tag);
            }
        }
    }

    /// Scatter data from one processor over multiple processors
    /*virtual void scatterV( E *sendBuf, E *recvBuf, int* sendCounts, int root = 0 ){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->scatterV(sendBuf, recvBuf, sendCounts, root);
    }*/

   /// Gather data from multiple processors to one processor
   virtual void gatherV(E* sendBuf, int sendCount, E* recvBuf, int* recvCounts, int* displs, int root=0){
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!this->useRelayer){
           m->gatherV(sendBuf, sendCount, recvBuf, recvCounts, displs, root);
       }else{
           this->handler->getMpiManager()->gatherV(sendBuf, sendCount, recvBuf, recvCounts, displs, root);
           if(this->handler->getMpiManager()->isMainProcessor()){
               int receiverCoreId= 0;
               this->forwardMessage(MPI_CONDUIT_OP::gatherV, (*recvCounts) , recvBuf, receiverCoreId, conduit_MPI_Tag);
           }
       }
   }

   virtual void bCast( E* sendBuf, int sendCount, int root = 0 ){

       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!this->useRelayer){
           m->bCast(sendBuf, sendCount, root);
       }else {

           int val = (root == MPI_PROC_NULL)? 1: 0;
           int global_sum;
           this->handler->getMpiManager()->reduce(val, global_sum, MPI_SUM);
           this->handler->getMpiManager()->bCast(&global_sum, 1);
           bool isBcastTo=(global_sum == (this->handler->getMpiManager()->getSize() -1))? true: false;
           if(this->handler->getMpiManager()->getSize()==1){
                if(root== MPI_ROOT)   // in case group contains only 1 proc
                    isBcastTo=true;
                else
                    isBcastTo=false;
           }

           if(this->handler->getMpiManager()->getRank()  == 0) {
               if(isBcastTo){ // bcast To
                   int receiverCoreId= 0;
                   this->forwardMessage(MPI_CONDUIT_OP::bCast, sendCount, sendBuf, receiverCoreId, conduit_MPI_Tag);
            }else{//bcast From
                   int tag=conduit_MPI_Tag;
                   int source =0;
                   if(useRelayer){source = MPI_ANY_SOURCE;}
                   MPI_Request request;
                   m->iRecv<E>(sendBuf, sendCount, source, &request, tag);
                   MPI_Status status;
                   MPI_Wait(&request, &status);
               }
           }
           if(!isBcastTo){
               this->handler->getMpiManager()->bCast(sendBuf, sendCount, root);
           }
       }
   }

   /// Broadcast data from one processor from multiple processors
   /*virtual void bCastFrom( E* sendBuf, int sendCount, int root = 0 ){

       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!this->useRelayer){
           m->bCast(sendBuf, sendCount, root);
       }else{
           if(this->handler->getMpiManager()->getRank() == 0) {
               int tag=conduit_MPI_Tag;
               int source =0;
               if(useRelayer){source = MPI_ANY_SOURCE;}
               MPI_Request request;
               m->iRecv<E>(sendBuf, sendCount, source, &request, tag);
               MPI_Status status;
               MPI_Wait(&request, &status);
           }
           this->handler->getMpiManager()->bCast(sendBuf, sendCount, root);
       }
   }*/



   /// Reduction operation toward one processor
   virtual  void reduce( E sendVal, E& recvVal, MPI_Op op, int root = 0 ){
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!this->useRelayer){
           m->reduce(sendVal, recvVal, op, root);
       }else{
           this->handler->getMpiManager()->reduce(sendVal, recvVal, op, root);

           if(this->handler->getMpiManager()->isMainProcessor()){
               int receiverCoreId= 0;
               this->forwardMessage(MPI_CONDUIT_OP::reduce, 1, &recvVal, receiverCoreId, conduit_MPI_Tag);
           }
       }
   }

    /// Element-per-element reduction of a vector of data
    virtual void reduceVect( std::vector<E>& sendVal, std::vector<E>& recvVal, MPI_Op op, int root = 0 ){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        if(!this->useRelayer){
            m->reduceVect(sendVal, recvVal, op, root);
        }else{
            this->handler->getMpiManager()->reduceVect(sendVal, recvVal, op, root);

            if(this->handler->getMpiManager()->isMainProcessor()){
                int receiverCoreId= 0;
                this->forwardMessage(MPI_CONDUIT_OP::reduceVect, recvVal.size(), (E *) recvVal.data(), receiverCoreId, conduit_MPI_Tag);
            }
        }
    }

    /// Inplace element-per-element reduction of a vector of data; result
    ///   available on all MPI threads.
    /*virtual void allReduceVect( std::vector<E>& sendRecvVal, MPI_Op op ){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->allReduceVect(sendRecvVal, op);
    }

    /// Reduction operation, followed by a broadcast
    virtual  void reduceAndBcast(E& reductVal, MPI_Op op, int root = 0 ){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->reduceAndBcast(reductVal, op, root);
    }*/


   //--- SEND ---
   virtual  void send( E *buf, int count, int dest, int tag = -1 ){
      if(tag==-1) tag=conduit_MPI_Tag;
       if(!useRelayer){
           shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
           m->send(buf, count, dest, tag);
       }else{
           int receiverCoreId= dest;
           this->forwardMessage(MPI_CONDUIT_OP::send,  count, buf, receiverCoreId, conduit_MPI_Tag);
       }
   }

   /// Sends data at *buf, non blocking
   virtual  void iSend( E *buf, int count, int dest, MPI_Request* request, int tag = -1 ){
       if(tag==-1) tag=conduit_MPI_Tag;
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!useRelayer){
           m->iSend(buf, count, dest, request, tag);
       }else{
           int receiverCoreId= dest;
           this->forwardMessage(MPI_CONDUIT_OP::iSend, count, buf, receiverCoreId, conduit_MPI_Tag, request);
       }
    }

   /// Sends data at *buf, assuming that receiver is ready.
   virtual void rSend( E *buf, int count, int dest, int tag = -1 ){
       if(tag==-1) tag=conduit_MPI_Tag;
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       if(!useRelayer){
           m->rSend(buf, count, dest, tag);
       }else{
           int receiverCoreId= dest;
           this->forwardMessage(MPI_CONDUIT_OP::rSend, count, buf, receiverCoreId, conduit_MPI_Tag);
        }
   }

    /// Sends data at *buf, non blocking and request free
    virtual  void iSendRequestFree( E *buf, int count, int dest, int tag = -1 ){
        if(tag==-1) tag=conduit_MPI_Tag;
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->iSendRequestFree(buf, count, dest, tag);
    }

    /// Sends data to master processor
    virtual  void sendToMaster( E* sendBuf, int sendCount, bool iAmRoot ){
        shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
        m->sendToMaster(sendBuf, sendCount, iAmRoot);
    }

   //-- RECV ---
   /// Receives data at *buf, blocking
   virtual  void receive( E *buf, int count, int source, int tag = -1 ){

       if(tag==-1) tag=conduit_MPI_Tag;
       if(useRelayer){source = MPI_ANY_SOURCE;}
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       m->receive<E>(buf, count, source, tag);
   }

   /// Receives data at *buf, non blocking
   virtual void iRecv( E *buf, int count, int source, MPI_Request* request, int tag = -1 ){
       if(tag==-1) tag=conduit_MPI_Tag;
       if(useRelayer){source = MPI_ANY_SOURCE;}
       shared_ptr<MpiManager> m= this->getHandlerInerCommunicatortManager();
       m->iRecv<E>(buf, count, source, request, tag);
   }

private:

    virtual size_t hashConduit_MPI_TAG( string id1, string id2){
        assert (!id1.empty() && ! !id2.empty());

        // suppose that two Mpi sub-group has several conduits between them. They will share the same interMpiManager but different mpi TAG.
        // since the name of the conduit is uniq, we compute the tag based on the name of the conduit using hash function
        std::hash<std::string> hash_fn;
        size_t int_hash_1 = hash_fn(id1);
        size_t int_hash_2 = hash_fn(id2);
        return int_hash_1+int_hash_2;
    }



   void forwardMessage(MPI_CONDUIT_OP op, int count, E * buf, int receiverCoreId, int tag, MPI_Request * request=0,  bool isEndMessage=false){
       int relayerdest= this->relayerUrl->grank;
       string senderName = this->handler->getName();
       int senderCoreId= this->handler->getMpiManager()->getRank();
       string receiverName= this->remote->getName();
       int operation = static_cast<int>(op);
       int dataType= static_cast<int>(this->conduitDataType);
       int countData =  count * sizeof(E);
       int endMessage=(isEndMessage)? 1: 0;
       if(op != MPI_CONDUIT_OP::iSend){
           Messaging::getInstance().send(this->handler->getGlobalMpiManager(), senderName,  senderCoreId,  receiverName,
                         receiverCoreId,  operation,  dataType, tag, endMessage, (char*) buf, countData, relayerdest);
       }else{
           Messaging::getInstance().iSend(this->handler->getGlobalMpiManager(), senderName,  senderCoreId,  receiverName,
                         receiverCoreId,  operation,  dataType, tag, endMessage, (char*) buf, countData, relayerdest, request);
       }
   }
};


/************************************ MPI_MTM_Conduit <E> ************************************************/
template<typename E>
class MPI_MTM_Conduit: public MPI_PTP_Conduit<E>{

public:
    MPI_MTM_Conduit(): MPI_PTP_Conduit<E>(){
    }
    virtual ~MPI_MTM_Conduit(){
    }

    virtual bool send(E * data, int count){
        if(this->handler->getMpiManager()->isMainProcessor()){ // only current sub group leader
            return MPI_PTP_Conduit<E>::send(data, count);
        }else{
            return 0;
        }
    }
    virtual E * receive( int & count){
        if(this->handler->getMpiManager()->isMainProcessor()){ // only current sub group leader
            return MPI_PTP_Conduit<E>::receive(count);
        }else{
            return 0;
        }
    }
    virtual bool isConduitSafe(){
        return true;
    }
};

/************************************ MPI_Relayer_Conduit <E> ************************************************/


/************************************ Remote_MPI_PTP_Conduit <E> ************************************************/
// not yet finished
/*template<typename E>
class Remote_MPI_PTP_Conduit: public MPI_PTP_Conduit <E>{

public:
    //friend class MpiKernel;
    Remote_MPI_PTP_Conduit(): MPI_PTP_Conduit<E>() {
    }

    virtual ~Remote_MPI_PTP_Conduit(){
    }

};*/

/************************************ BasicConduit <E> ************************************************/

template<class E, template<class> class C>
class ConduitFactory{

    /**
     * Produces a new conduit instance.
     * @return A newly instanciated conduit.
     */
public:
    virtual ~ConduitFactory(){}
    virtual C<E> * newInstance()=0;
};

/************************************ BasicConduitFactory <E> ************************************************/

/*template<typename E>
class BasicConduitFactory: public ConduitFactory <E, BasicConduit> {

public:
    BasicConduit<E> * newInstance() {
        return new BasicConduit<E>();
    }
};*/

template<typename E>
class BasicConduitFactory{
public:
    BasicConduit<E> * newInstance() {
        return new BasicConduit<E>();
    }
};

template<typename E>
class MPI_MTM_ConduitFactory{
public:
    MPI_MTM_Conduit<E> * newInstance() {
        return new MPI_MTM_Conduit<E>();
    }
};

template<typename E>
class MPI_PTP_ConduitFactory{
public:
    MPI_PTP_Conduit<E> * newInstance() {
        return new MPI_PTP_Conduit<E>();
    }
};




}//end namespace

#endif // CONDUIT_H
