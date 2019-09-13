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


#ifndef KERNELS_H
#define KERNELS_H

#include <algorithm>    // std::sort
#include <memory>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf


#include "musclehpc/utils/util.h"
#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/conduits/conduitAbstract.h"
#include "musclehpc/mapper/mappercommon.h"
#include "musclehpc/mediator/mpisocketserver.h"
#include "musclehpc/utils/mmsfperf.h"

using namespace std;
using namespace unige_pasc;

namespace unige_pasc{
/************************************ Port ************************************************/
/**
 * @brief The Port class is the MMSP PORT NAME.
 */
class Port{
private:
    string kernelName;
    string portName;
    std::string delimiter = ".";
public:
    /**
     * @brief Port
     * @param str MMSF coupling information. Syntax KERNEL_NAME:PORT_NAME
     */
    Port(string str){
        this->parse(str);
    }
    string getKernelName()const {return kernelName;}
    string getPortName()const {return portName;}
private:
    /**
     * @brief parse
     * @param s Syntax KERNEL_NAME:PORT_NAME
     */
    void parse(string s){
        //reg expr of s:  (\\w+).(\\w+)
        size_t pos = 0;
        // serach kernel name
        if((pos = s.find(delimiter)) != std::string::npos){
            kernelName = s.substr(0, pos);
            portName = s.substr(pos + delimiter.length(), s.length()-1);
        }
    }
};

/************************************ KernelCommon ************************************************/

class KernelCommon{

    int kernelArgc;
    char ** kernelArgv;
public:
    virtual ~KernelCommon(){}

    // argv and argc per submodel
    int getKernelArgc() const{ return this->kernelArgc;}
    char ** getKernelArgv() const{ return this->kernelArgv;}
    void setKernelArgc(int argc){this->kernelArgc=argc;}
    void setKernelArgv(char ** argv){this->kernelArgv=argv;}

    virtual void initialize(string getName)=0;
    /**
     * @brief ID getKernel Identifier
     * @return string name of the kernel
     */
    virtual string getName()=0;
    /**
     * @brief mainLoop runs the main kernel in a loop
     */
    virtual void mainLoop()=0;
    /**
     * @brief registerFout assign a port to an exit conduit.
     * @param fout the port name
     * @param conduit a pointer to a conduit
     */
    virtual void registerFout(string fout, ConduitCommon * conduit) = 0;
    /**
     * @brief registerFin assign a port to an entry conduit.
     * @param fin the port name
     * @param conduit a pointer to a conduit
     */
    virtual void registerFin(string fin, ConduitCommon * conduit)=0 ;

};

/************************************ Kernel ************************************************/

class Kernel: public KernelCommon{
protected:
    map <string, ConduitCommon *> f_out; // mapping between an exit  conduit and its name
    map <string, ConduitCommon *> f_in;  // mapping between an entry conduit and its name
    int colorWhereiTRun; // this is  MPI group Id where the kernel should run
    typedef typename std::map<string, ConduitCommon * >::iterator it_map_condtuis;
    int coresOverWhichTorun; // number of MPI cores allocated to this kernels
    int globalLeaderRank;
    string id;
    // constructor
    Kernel():colorWhereiTRun(0),coresOverWhichTorun(0){
    }

public:

    virtual ~Kernel(){
        f_out.clear();
        f_in.clear();
    }
    virtual void initialize(string iD){
        this->id=iD;
    }
    virtual string getName(){
        return this->id;
    }

    virtual void mainLoop()=0;


    void setCoresOverWhichTorun(int cores){
        this->coresOverWhichTorun=cores;
    }
    int getCoresOverWhichTorun()const{
        return this->coresOverWhichTorun;
    }
    // mpi force to run in mpi subgroup of color  colorWhereiTRun
    /**
     * @brief setColorWhereToRun force the kernel to run in the MPI subgroup of the id colorWhereiTRun
     * @param colorWhereiTRun the MPI subgroup id
     */
    virtual void setColorWhereToRun(int colorWhereiTRun){
        this->colorWhereiTRun=colorWhereiTRun;
    }
    virtual int getColorWhereItRuns() const{
        return this->colorWhereiTRun;
    }

    virtual void setGlobalLeaderRank(int globalLeaderRank){
        this->globalLeaderRank=globalLeaderRank;
    }
    virtual int getGlobalLeaderRank() const{
        return this->globalLeaderRank;
    }

    virtual void registerFout(string fout, ConduitCommon * conduit){
        //Conduit<E>* conduit = dynamic_cast<Conduit<E>*>(conduitCommon);
        it_map_condtuis it = f_out.find(fout);
        if( it!= f_out.end()) {
            cout<<" A f_out was already declared with ID: " << it->first <<endl;
        }else{
            //cout<<"kernel "<< ID()<<"  insert f_out declared with ID: " << fout << " with conduit " << conduit->getName() <<endl;
            f_out.insert ( std::pair<string, ConduitCommon* >(fout, conduit ));
        }
    }
    virtual void registerFin(string fin, ConduitCommon * conduit){
        it_map_condtuis it = f_in.find(fin);
        if( it!= f_in.end()) {
            cout<<" A F_in was already declared with ID: " << it->first <<endl;
        }else{
            // cout<<"kernel "<< ID()<<"  insert f_in declared with ID: " << fin << " with conduit " << conduit->getName() <<endl;
            f_in.insert ( std::pair<string, ConduitCommon* >(fin, conduit ));
        }
    }

    template<typename E>
    ConduitEntrance<E> * getConduitEntrance( string id ) {
        it_map_condtuis it = f_out.find(id);
        if( it!= f_out.end()) {
            ConduitCommon * conduit = it->second;
            ConduitEntrance<E> * ptr = dynamic_cast< ConduitEntrance<E> *>(conduit);
            return ptr;
        }else{
            return 0;
        }
    }

    template<typename E>
    ConduitExit<E> * getConduitExit( string id ) {
        it_map_condtuis it = f_in.find(id);
        if( it!= f_in.end()) {
            ConduitCommon * conduit = it->second;
            ConduitExit<E> * ptr = dynamic_cast< ConduitExit<E> *>(conduit);
            return ptr;
        }else{
            return 0;
        }
    }


    /**
     * @brief prepareConduits prepares conduits to be operational:
     * create MPI intercommunicator for the case of MPI_Conduit between sender and receiver kernels.
     */
    virtual void prepareConduits(){

        for(  it_map_condtuis it = f_in.begin(); it!= f_in.end(); it++) {
            ConduitCommon * conduit = it->second;
            assert(conduit);
            conduit->setHandler(this);
        }

        for(  it_map_condtuis it = f_out.begin(); it!= f_out.end(); it++) {
            ConduitCommon * conduit = it->second;
            assert(conduit);
            conduit->setHandler(this);
        }
    }
};

/************************************ MpiKernel ************************************************/

class MpiKernel: public Kernel{

    friend class  SocketServer;
protected:
    shared_ptr<MpiManager> mpiManager;
    shared_ptr<MpiManager> globalMpiManager;
    std::map<int, shared_ptr<MpiManager> > intercommMap;// <uniqInterCommTag, MPI_Comm>
    typedef typename std::map<int, shared_ptr<MpiManager> >::iterator it_map_comm_int_MpiManager;
    typedef typename std::map<string, Kernel *>::iterator it_map_comm_string_kernel;


    MpiKernel():Kernel(){
        //mpiManager=0;
    }

public:

    virtual ~MpiKernel(){
        /*for( it_map_comm_int_MpiManager it = intercommMap.begin(); it!=intercommMap.end();it++){
            MpiManager * ptr= it->second;
            if(ptr)
                delete ptr;
        }*/
        intercommMap.clear();
    }

    void setMpiManager(shared_ptr<MpiManager> & mpiManager){
        this->mpiManager=mpiManager;
    }

    void setGlobalMpiManager(shared_ptr<MpiManager> & globalMpiManager){
        this->globalMpiManager=globalMpiManager;
    }

    shared_ptr<MpiManager> getMpiManager(){
        return this->mpiManager;
    }

    /*MpiManager * getMpiManager(){
        return this->mpiManager.get();
    }*/

    shared_ptr<MpiManager> & getGlobalMpiManager(){
        return this->globalMpiManager;
    }

    /**
     * @brief addInterCommunicator assigns a MpiManager tp each MPI interCommunicator environement
     * @param interCommTag tag of the interCommunicator
     * @param interManager MpiManager pointer
     */
    void addInterCommunicator(int uniqInterCommTag, shared_ptr<MpiManager> & interManager){
        if(! existInterCommunicator( uniqInterCommTag)){
            intercommMap.insert ( std::pair<int, shared_ptr<MpiManager>  >(uniqInterCommTag, interManager ));
        }
    }

    /**
     * @brief existInterCommunicator checks whether the kernel has already assigned an MPI interCommunicator
     * @param uniqInterCommTag tag of the interCommunicator
     * @return pointer to the assigned MpiManager, 0 otherwise.
     */
    shared_ptr<MpiManager>  existInterCommunicator(int uniqInterCommTag){
        shared_ptr<MpiManager>  m;
        it_map_comm_int_MpiManager it = intercommMap.find(uniqInterCommTag);
        if(it != intercommMap.end()){
            m= it->second;
        }
        return m;
    }


    /**
     * @brief prepareConduits prepares conduits to be operational:
     * create MPI intercommunicator for the case of MPI_Conduit between sender and receiver kernels.
     */
    virtual void prepareConduits(){

        // 1) call Kernel parent prepareConduits()
        Kernel::prepareConduits();

        // 2) select kernels
        std::map<string, Kernel *> allMyRemoteKernelsMap = getAllmyRemoteKernel();

        //3) create intercommunicators using MPI_Split mechanism
        this->createInterCommunicators(allMyRemoteKernelsMap);
    }


    virtual void createInterCommunicators(std::map<string, Kernel *> & allMyRemoteKernelsMap){

        for(  it_map_comm_string_kernel it = allMyRemoteKernelsMap.begin(); it!= allMyRemoteKernelsMap.end(); it++) {

            MpiKernel *localKernel= dynamic_cast<  MpiKernel *> (this);
            MpiKernel *remoteKernel= dynamic_cast<  MpiKernel *> (it->second);

            assert(localKernel);
            assert(remoteKernel);

            if(localKernel->getName() == remoteKernel->getName())
                continue;

            int uniqInterCommTag=CantorPairingFunctionKernels(localKernel, remoteKernel);

            int remoteLeaderGlobalRank = remoteKernel->getGlobalLeaderRank();

            if(! localKernel->existInterCommunicator(uniqInterCommTag)){
                MPI_Comm intercomm ;
                assert(localKernel->getMpiManager());
                MPI_Intercomm_create(localKernel->getMpiManager()->getGlobalCommunicator(),
                                     localKernel->getMpiManager()->bossId(),
                                     MPI_COMM_WORLD,
                                     remoteLeaderGlobalRank ,
                                     uniqInterCommTag,
                                     &intercomm);
                shared_ptr<MpiManager > interMpimanager= std::make_shared<MpiManager>();
                interMpimanager->initOnlyCommunicator(intercomm);
                localKernel->addInterCommunicator(uniqInterCommTag, interMpimanager);
            }
        }
    }


    std::map<string, Kernel *>  getAllmyRemoteKernel(){

        std::map<string, Kernel *> allMyRemoteKernelsMap;

        for(  it_map_condtuis it = f_in.begin(); it!= f_in.end(); it++) {
            ConduitCommon * conduit = it->second;
            assert(conduit);
            //Kernel * k= conduit->getSenderSubmodel();
            Kernel * k= conduit->getRemoteKernel();
            if(k)
                allMyRemoteKernelsMap.insert ( std::pair<string, Kernel* >(k->getName(), k));// here it doest not allow duplicated key
        }
        for(  it_map_condtuis it = f_out.begin(); it!= f_out.end(); it++) {
            ConduitCommon * conduit = it->second;
            assert(conduit);
            //Kernel * k= conduit->getReceiverSubmodel();
            Kernel * k= conduit->getRemoteKernel();
            if(k)
                allMyRemoteKernelsMap.insert ( std::pair<string, Kernel* >(k->getName(), k));// here it doest not allow duplicated key
        }
        return allMyRemoteKernelsMap;
    }

protected:

    int CantorPairingFunctionKernels(MpiKernel * localKernel, MpiKernel *remoteKernel){
        int senderLeaderColor = localKernel->getColorWhereItRuns();
        int receiverLeaderColor = remoteKernel->getColorWhereItRuns();
        return CantorPairingFunction(senderLeaderColor, receiverLeaderColor);
    }

};



/************************************ RemoteMpiKernel ************************************************/
class RemoteMpiKernel: public MpiKernel{

private:
    shared_ptr<SocketServer> mpiServerSocket ; //= topologyConnection(serverUri,  allMyRemoteKernelsId);
    typedef typename std::map<string, shared_ptr<KernelConnectionInfo> >::iterator it_map_comm_string_ConnInfo;

public:
    RemoteMpiKernel(): MpiKernel(){
    }
    virtual ~RemoteMpiKernel(){}

    void setMpiSocketServer(shared_ptr<SocketServer> mpiServerSocket_ptr){
        mpiServerSocket= mpiServerSocket_ptr;
    }

    /**
     * @brief prepareConduits prepares conduits to be operational:
     * create MPI intercommunicator for the case of MPI_Conduit between sender and remote receiver kernels.
     */
    virtual void prepareConduits(){

        // 1) call Kernel parent prepareConduits()
        Kernel::prepareConduits();

        // 2) select kernels
        std::map<string, Kernel *> allMyRemoteKernelsMap = getAllmyRemoteKernel();

        //3) create intercommunicators using Listen/Accept mechnism
        this->createInterCommunicatorsWithListenAcceptMechanism(allMyRemoteKernelsMap);

    }

protected:
    virtual void createInterCommunicatorsWithListenAcceptMechanism(std::map<string, Kernel *> & allMyRemoteKernelsMap){

        // construct a vector
        vector<string> allMyRemoteKernelsId;
        for(  it_map_comm_string_kernel it = allMyRemoteKernelsMap.begin(); it!= allMyRemoteKernelsMap.end(); it++) {
            MpiKernel *localKernel= dynamic_cast<  MpiKernel *> (this);
            MpiKernel *remoteKernel= dynamic_cast<  MpiKernel *> (it->second);

            assert(localKernel);
            assert(remoteKernel);

            if(localKernel->getName() == remoteKernel->getName())
                continue;

            int uniqInterCommTag=CantorPairingFunctionKernels(localKernel, remoteKernel);
            if(! localKernel->existInterCommunicator(uniqInterCommTag)){
                allMyRemoteKernelsId.push_back(it->first);
            }
            //else: I have already created an intercommunication with this remote kernel since we are within the same binary
        }
        // remove duplication from vector
        std::sort( allMyRemoteKernelsId.begin(), allMyRemoteKernelsId.end() );
        allMyRemoteKernelsId.erase( std::unique( allMyRemoteKernelsId.begin(), allMyRemoteKernelsId.end() ), allMyRemoteKernelsId.end() );
        assert(allMyRemoteKernelsId.size() <= allMyRemoteKernelsMap.size());

        // start remote registration machanism
        mpiServerSocket->setMyRemoteKernels(allMyRemoteKernelsId);

        mpiServerSocket->initConnection();
        mpiServerSocket->handShake();

        map<string, shared_ptr<KernelConnectionInfo>  > & remoteKernelsMap= mpiServerSocket->getRemoteKernelsMap();
        shared_ptr<KernelConnectionInfo> myLocalInfo= mpiServerSocket->getMyLocalInfo();
        int handlerLeaderColor = myLocalInfo->getcolor() ;
        this->setColorWhereToRun(myLocalInfo->getcolor());
        this->setMpiManager(mpiServerSocket->getMpiManager());
        assert(mpiServerSocket->getMpiManager().get());

        //assert(remoteKernelsMap.size() -1 == allMyRemoteKernelsMap.size());

        for( it_map_comm_string_kernel it = allMyRemoteKernelsMap.begin(); it!= allMyRemoteKernelsMap.end();it++){
            string kid= it->first;
            RemoteMpiKernel* kmpi= dynamic_cast<RemoteMpiKernel* > (it->second);

            shared_ptr<KernelConnectionInfo> & infoConn_ptr=remoteKernelsMap.at(kid);
            kmpi->setColorWhereToRun(infoConn_ptr->getcolor());
            int receiverLeaderColor =  infoConn_ptr->getcolor() ;

            shared_ptr<MpiManager> interManager=infoConn_ptr->getInterMpiManager();
            int uniqInterCommTag;

            if(interManager){
                uniqInterCommTag =CantorPairingFunction(handlerLeaderColor, receiverLeaderColor);
            }else{ // the kernel is: already in the same binary as me (OR) has a pb to open a MPI_PORT.

                MpiKernel *localKernel= dynamic_cast<  MpiKernel *> (this);
                MpiKernel *remoteKernel= dynamic_cast<  MpiKernel *> (it->second);
                assert(localKernel);
                assert(remoteKernel);
                if(localKernel->getName() == remoteKernel->getName())
                    continue;
                uniqInterCommTag=CantorPairingFunctionKernels(localKernel, remoteKernel);
                if(! localKernel->existInterCommunicator(uniqInterCommTag)){
                    /// so use the relayer
                   std::shared_ptr<SocketKernel> derived = std::dynamic_pointer_cast<SocketKernel> (mpiServerSocket);
                   assert(derived);
                   interManager= derived->getRelayerCommunicator()->getInterMpiManager();
                   shared_ptr<Endpoint> edp= std::make_shared<Endpoint>(mpiServerSocket->getMyLocalInfo()->getMyRelayerUrl());
                   for(  it_map_condtuis it = f_in.begin(); it!= f_in.end(); it++) {
                       ConduitCommon * conduit = it->second;
                       if (conduit->getRemoteKernel()->getName() == remoteKernel->getName())
                            conduit->setRelayerURL(edp);

                   }
                   for(  it_map_condtuis it = f_out.begin(); it!= f_out.end(); it++) {
                       ConduitCommon * conduit = it->second;
                        if (conduit->getRemoteKernel()->getName() == remoteKernel->getName())
                            conduit->setRelayerURL(edp);
                   }

                }else{/*the kernel is: already in the same binary as me*/}
            }
            this->addInterCommunicator(uniqInterCommTag, interManager);

        }
        mpiServerSocket->disconnect();

    }

};

/************************************ EmptyRemoteMpiKernel ************************************************/

class EmptyRemoteMpiKernel: public RemoteMpiKernel{
private:
    string id;
public:
    EmptyRemoteMpiKernel(): RemoteMpiKernel(){
    }


    virtual ~EmptyRemoteMpiKernel(){}
    //virtual void initialize(string id){this->id=id;}
    //virtual string ID(){return this->id;}
    virtual void mainLoop(){
    }

};


/*****************************************************************************************************************/

class AsynchronousRemoteMpiKernel: public RemoteMpiKernel{

private:
    MPI_Request  listISendReq [2];
    MPI_Request  listIRecvReq [2];

public:
    AsynchronousRemoteMpiKernel():RemoteMpiKernel(){}
    virtual ~AsynchronousRemoteMpiKernel(){}

public:

    void iSend( ConduitEntrance<char> *f_out, char * data, int s){

        MpiKernel * handler= dynamic_cast<MpiKernel *> (f_out->getHandlerKernel());
        shared_ptr<MpiManager>  handlerMpimanager =handler->getMpiManager();
        MPI_Status  status;

        f_out->iSend((char*) &s, sizeof(s), handlerMpimanager->getRank(),&listISendReq[0]);
        MPI_Wait(&listISendReq[0], &status);
        f_out->iSend(data, s, handlerMpimanager->getRank(),&listISendReq[1]);
    }

    char * iRecv(ConduitExit<char> *f_in, int & size ){

        MpiKernel * handler= dynamic_cast<MpiKernel *> (f_in->getHandlerKernel());
        shared_ptr<MpiManager>  handlerMpimanager =handler->getMpiManager();

        MPI_Status  status;
        f_in->iRecv((char*) &size, sizeof(size), handlerMpimanager->getRank(), &listIRecvReq[0]);
        MPI_Wait(&listIRecvReq[0], &status);
        char * data = new char[size];
        f_in->iRecv(data, size, handlerMpimanager->getRank(), &listIRecvReq[1]);
        return data;
    }



    char * iRecvAndWaitAll(ConduitExit<char> *f_in, int & size ){
        char * data = iRecv(f_in, size);
        waitReceive();
        waitSend();
        return data;
    }

    void waitSend(){
        MPI_Status  status;
        MPI_Wait(&listISendReq[1], &status);
    }
    void waitReceive(){
        MPI_Status  status;
        MPI_Wait(&listIRecvReq[1], &status);
    }


};

//***********************************MpiSubmodel************************************************/
class MpiSubmodel: public RemoteMpiKernel, public Submodel_A{

public:
    MpiSubmodel(): RemoteMpiKernel(), Submodel_A(){
    }
    virtual ~MpiSubmodel(){}
    virtual void mainLoop(){
        OptionParser * parser = OptionParser::getInstance();
        shared_ptr<MMSFPerf<double> > mmsfprofiler;

        if (parser && parser->isProfiling()){
            mmsfprofiler = std::make_shared<MMSFPerf<double>>();
            mmsfprofiler->setPrefixFileName(this->getName());
            mmsfprofiler->setRank(this->getMpiManager()->getRank());
        }

        this->simulate(mmsfprofiler.get());

    }
};


}// end namespace unige_pasc

#endif // KERNELS_H
