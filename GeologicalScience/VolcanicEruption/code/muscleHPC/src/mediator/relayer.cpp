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

#include "musclehpc/mediator/relayer.h"


/*********************************************************RelayerCommunicator *************************************************/

bool RelayerCommunicator::sendInt(int value){
    vector<char> buffer;
    buffer.resize(sizeof(value));
    memcpy((char*) buffer.data(), (const char *) &value, sizeof(value));
    bool res=send(buffer);
    vector<char>().swap(buffer);
    return res;
}

bool RelayerCommunicator::receiveInt(int &value){
    vector<char> buffer;
    bool res=receive(buffer);
    memcpy((char*) &value, (const char *)buffer.data(), sizeof(value));
    vector<char>().swap(buffer);
    return res;
}

bool RelayerCommunicator::sendString(string message){
    vector<char> buffer;
    buffer.resize(message.length()+1);
    std::copy(message.c_str(), message.c_str()+message.length()+1, buffer.data());
    bool res=send(buffer);
    vector<char>().swap(buffer);
    return res;
}

bool RelayerCommunicator::receiveString(string &message){
    vector<char> buffer;
    bool res=receive(buffer);
    message=buffer.data();
    vector<char>().swap(buffer);
    return res;
}


bool RelayerCommunicator::bCastIntTo(int value, int root){
    vector<char> buffer;
    buffer.resize(sizeof(value));
    memcpy((char*) buffer.data(), (const char *) &value, sizeof(value));
    bool res=this->bCastTo(buffer, root);
    vector<char>().swap(buffer);
    return res;
}

bool RelayerCommunicator::bCastIntFrom(int &value, int root){
    vector<char> buffer;
    bool res=this->bCastFrom(buffer, root);
    memcpy((char*) &value, (const char *)buffer.data(), sizeof(value));
    vector<char>().swap(buffer);
    return res;
}


bool RelayerCommunicator::bCastStringTo(string message, int root){
    vector<char> buffer;
    buffer.resize(message.length()+1);
    std::copy(message.c_str(), message.c_str()+message.length()+1, buffer.data());
    bool res=this->bCastTo(buffer, root);
    vector<char>().swap(buffer);
    return res;
}


bool RelayerCommunicator::bCastStringFrom(string &message, int root){
    vector<char> buffer;
    bool res= this->bCastFrom(buffer, root);
    message=buffer.data();
    vector<char>().swap(buffer);
    return res;
}

bool RelayerCommunicator::bCastUrlTo(Endpoint &relayerURL, int root){
    return this->bCastStringTo(relayerURL.rawEndpoints, root);
}

bool RelayerCommunicator::bCastUrlFrom(Endpoint &relayerURL, int root){

    string raw_endpoint;
    bool res= this->bCastStringFrom(raw_endpoint, root);
    relayerURL.parse(raw_endpoint);
    return res;
}







/********************************************************* MPI_Client_RelayerCommunicator *************************************************/

MPI_Client_RelayerCommunicator::MPI_Client_RelayerCommunicator(shared_ptr<MpiManager> &myMpiManager, shared_ptr<MpiManager> &interMpiManager)
    : RelayerCommunicator(), myMpiManager(myMpiManager), interMpiManager(interMpiManager){}

MPI_Client_RelayerCommunicator::~MPI_Client_RelayerCommunicator(){}

shared_ptr<MpiManager> MPI_Client_RelayerCommunicator::getInterMpiManager(){
    return interMpiManager;
}

void MPI_Client_RelayerCommunicator::sendMessage(Message &m, int header){

    this->sendString(m.senderName);
    this->sendString(m.receiverName);
    this->sendInt(m.senderCoreId);
    this->sendInt(m.receiverCoreId);
    this->sendInt(m.operation);
    this->sendInt(m.datatype);
    this->sendInt(m.MPI_tag);
    this->sendInt(m.endMessage);
    this->send(m.data);

}

void MPI_Client_RelayerCommunicator::receiveMessage(Message &m, int header){

    this->receiveString(m.senderName);
    this->receiveString(m.receiverName);
    this->receiveInt(m.senderCoreId);
    this->receiveInt(m.receiverCoreId);
    this->receiveInt(m.operation);
    this->receiveInt(m.datatype);
    this->receiveInt(m.MPI_tag);
    this->receiveInt(m.endMessage);
    this->receive(m.data);
}

bool MPI_Client_RelayerCommunicator::sendRaw( char *data, int length, int header){
    MPI_Request request;
    MPI_Status status;
    this->interMpiManager->iSend<char>((char*)data, length, myMpiManager->getRank(), &request, 0);
    MPI_Wait (&request, &status);
    return true;
}


bool MPI_Client_RelayerCommunicator::send(vector<char> &buffer, int header){
    int length = (int) buffer.size();
    this->interMpiManager->send<int>(&length, 1, myMpiManager->getRank(), 0);
    MPI_Request request;
    MPI_Status status;
    this->interMpiManager->iSend<char>(buffer.data(), buffer.size(), myMpiManager->getRank(), &request, 0);
    MPI_Wait (&request, &status);
    return true;
}

bool MPI_Client_RelayerCommunicator::receive(vector<char> &buffer, int header){
    int length;
    this->interMpiManager->receive<int>(&length, 1 , myMpiManager->getRank(),  0);
    buffer.resize(length);
    MPI_Request request;
    MPI_Status status;
    this->interMpiManager->iRecv<char>(buffer.data(), length , myMpiManager->getRank() , &request, 0);
    MPI_Wait (&request, &status);
    return true;
}

bool MPI_Client_RelayerCommunicator::bCastTo(vector<char> &buffer, int root){
    int length = (int) buffer.size();
    this->interMpiManager->bCast<int>(&length, 1, root);
    this->interMpiManager->bCast<char>(buffer.data(), buffer.size(), root);
    return true;
}
bool MPI_Client_RelayerCommunicator::bCastFrom(vector<char> &buffer, int root){
    int length;
    this->interMpiManager->bCast<int>(&length, 1, root);
    buffer.resize(length);
    this->interMpiManager->bCast<char>(buffer.data(), buffer.size(), root);
    return true;
}

bool MPI_Client_RelayerCommunicator::disconnect(){
    //return MPI_Close_port(remotePortName);
    return true;
}
/********************************************************* TCP_Client_RelayerCommunicator *************************************************/

TCP_Client_RelayerCommunicator::TCP_Client_RelayerCommunicator(string managerPortUrl, shared_ptr<MpiManager> &myMpiManager)
    : RelayerCommunicator(), managerPortUrl(managerPortUrl), myMpiManager(myMpiManager) {

}

TCP_Client_RelayerCommunicator::TCP_Client_RelayerCommunicator(int sockfd, shared_ptr<MpiManager> &myMpiManager)
    : RelayerCommunicator(), sockfd(sockfd), myMpiManager(myMpiManager){
}


TCP_Client_RelayerCommunicator::~TCP_Client_RelayerCommunicator(){
    // don't free sockfd here !!!
}

shared_ptr<MpiManager > TCP_Client_RelayerCommunicator::getInterMpiManager(){
    return nullptr;
}

void TCP_Client_RelayerCommunicator::sendMessage(Message &m, int header){

    vector<char> buffer;
    m.toBytes(buffer);
    this->send(buffer, header);
    vector<char>().swap(buffer);

    /*this->sendString(m.senderName);
    this->sendString(m.receiverName);
    this->sendInt(m.senderCoreId);
    this->sendInt(m.receiverCoreId);
    this->sendInt(m.operation);
    this->sendInt(m.datatype);
    this->sendInt(m.MPI_tag);
    this->sendInt(m.endMessage);
    this->send(m.data);*/

}

void TCP_Client_RelayerCommunicator::receiveMessage(Message &m,  int header){


    vector<char> buffer;
    this->receive(buffer, header);
    m.fromByte(buffer);
    vector<char>().swap(buffer);

    /*this->receiveString(m.senderName);
    this->receiveString(m.receiverName);
    this->receiveInt(m.senderCoreId);
    this->receiveInt(m.receiverCoreId);
    this->receiveInt(m.operation);
    this->receiveInt(m.datatype);
    this->receiveInt(m.MPI_tag);
    this->receiveInt(m.endMessage);
    this->receive(m.data);*/
}

bool TCP_Client_RelayerCommunicator::sendRaw( char *data, int length, int header){

    Networking & net = Networking::getInstance();
    int n =net.sock_write(sockfd, data, length, header);
    return n;
}


bool TCP_Client_RelayerCommunicator::send(vector<char> &buffer, int header){

    int length = buffer.size();
    Networking & net = Networking::getInstance();
    int n1=net.sock_write(sockfd, &length, sizeof(length), header);
    int n =net.sock_write(sockfd, buffer.data(), buffer.size(), header);


    return n && n1;
}

bool TCP_Client_RelayerCommunicator::receive(vector<char> &buffer, int header){

    Networking & net = Networking::getInstance();
    int n=0, n1=0;
    int value;
    n = net.sock_read(sockfd,&value,sizeof(value), header);//read(sockfd,&value,sizeof(value));
    buffer.resize(value);
    n1 = net.sock_read(sockfd, (char*) &buffer[0], value, header);//read(sockfd, (char*) &buffer[0], value);

    return (n&&n1)? true:false;
}
bool TCP_Client_RelayerCommunicator::bCastTo(vector<char> &buffer, int root){
    if(myMpiManager->getRank() == root){
        return this->send(buffer);
    }
    return true;
}

bool TCP_Client_RelayerCommunicator::bCastFrom(vector<char> &buffer, int root){

    int length ;
    if(myMpiManager->getRank() == root){
        this->receive(buffer);
        length = (int) buffer.size();
    }
    myMpiManager->bCast<int>(&length, 1, root);
    if(myMpiManager->getRank() != root){
        buffer.resize(length);
    }
    myMpiManager->bCast<char>(buffer.data(), length, root);
    return true;
}

bool TCP_Client_RelayerCommunicator::disconnect(){
    if(sockfd >0)
        return close(sockfd);
    return true;
}

void TCP_Client_RelayerCommunicator::unsetSocket(){
    this->sockfd = -1;
}


/********************************************************* MPIRelayer *************************************************/

MPIRelayer::MPIRelayer(shared_ptr<MpiManager>  globalMpiManager, shared_ptr<MpiManager> &myMpiManager,
                       map<string, shared_ptr<RelayerCommunicator> > & localKernelMap,
                       map<string, int> & mapColorWhereToRun,string forwarderAddress )
    : globalMpiManager(globalMpiManager), myMpiManager(myMpiManager),
      max_sd(0), localKernelMap(localKernelMap), mapColorWhereToRun(mapColorWhereToRun){

    if(this->myMpiManager->getSize()==1){
        this->backendRank=this->myMpiManager->bossId();
        this->frontendRank=backendRank;
    }else{
        this->backendRank=this->myMpiManager->bossId();
        this->frontendRank=1;
    }
    // init relayerEndPoint
    /*int grankBackEnd;
    if(this->isRelayerBackend()){
        grankBackEnd= this->globalMpiManager->getRank();
    }*/

    request=MPI_REQUEST_NULL;
    this->forwarderAddress.parse(forwarderAddress);

}

bool MPIRelayer::isRelayerBackend(){
    if(this->myMpiManager->getSize()==1 || this->myMpiManager->getRank() == this->backendRank){
        return true;
    }else{
        return false;
    }
}

bool MPIRelayer::isRelayerFrontEnd(){
    if(this->myMpiManager->getSize()==1 || this->myMpiManager->getRank()== this->frontendRank){
        return true;
    }else{
        return false;
    }
}

int MPIRelayer::getBackEndCore() const {return this->backendRank;}

int MPIRelayer::getFrontEndCore() const{return this->frontendRank;}

void MPIRelayer::run(){
    //core0: receives from kernels and forwared to remote.
    //Core1: listen on frontend socket and forward to local kernels

    int grankBackEnd;
    if(this->isRelayerBackend()){
        grankBackEnd= this->globalMpiManager->getRank();
    }
    this->myMpiManager->bCast(&grankBackEnd, 1, this->getBackEndCore());

    // (2) open FrontEnd Socket server
    if(this->isRelayerFrontEnd()){
        frontEnd_sockfd=Networking::getInstance().openfrontEndSocket(relayerEndPoint, grankBackEnd);
        max_sd = frontEnd_sockfd;
        FD_ZERO (&active_fd_set);
        FD_SET (frontEnd_sockfd, &active_fd_set);
    }
    this->myMpiManager->bCast(relayerEndPoint.rawEndpoints, this->getFrontEndCore());
    relayerEndPoint.parse();

    // (3) send them my TCP frond end URL
    for (map<string, shared_ptr<RelayerCommunicator> >::iterator it= localKernelMap.begin(); it!=localKernelMap.end(); ++it){
        //string kernelName = it->first;
        shared_ptr<RelayerCommunicator> & relayerCommunicator = it->second;
        int root= (this->myMpiManager->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
        //broadCastToRemote my relay endpoint
        vector<char> buffer;
        buffer.resize(relayerEndPoint.rawEndpoints.length()+1);
        std::copy(relayerEndPoint.rawEndpoints.c_str(), relayerEndPoint.rawEndpoints.c_str()+relayerEndPoint.rawEndpoints.length()+1, &buffer[0]);
        relayerCommunicator->bCastTo(buffer, root);
        // broscast my forwarder url
        buffer.resize(forwarderAddress.rawEndpoints.length()+1);
        std::copy(forwarderAddress.rawEndpoints.c_str(), forwarderAddress.rawEndpoints.c_str()+forwarderAddress.rawEndpoints.length()+1, &buffer[0]);
        relayerCommunicator->bCastTo(buffer, root);

    }

    map<string, Endpoint> relayersMap;//<kernelName, rawRelayerEndpoint>
    map<string, Endpoint> fwdMap;//<kernelName, rawForwarderEndpoint>
    int len=-1;
    //receive relayers list
    if(this->isRelayerFrontEnd()){
        receiveRelayersList(relayersMap, fwdMap);
        len= relayersMap.size();
    }

    //bcast the managerurl
    this->myMpiManager->bCast(this->managerPortUrl, this->getFrontEndCore());
    // breacast relayers list to backend core
    this->myMpiManager->bCast(&len, 1, this->getFrontEndCore());
    map<string, Endpoint>::iterator it = relayersMap.begin();
    for(int i=0; i< len; i++){
        string kernelName;
        string rel;
        string fwd;
        if(this->isRelayerFrontEnd()){
            kernelName=it->first;
            Endpoint & edp =it->second;
            rel=edp.rawEndpoints;
            fwd=fwdMap[kernelName].rawEndpoints;
            it++;
        }
        this->myMpiManager->bCast(kernelName, this->getFrontEndCore());
        this->myMpiManager->bCast(rel, this->getFrontEndCore());
        this->myMpiManager->bCast(fwd, this->getFrontEndCore());
        if(! this->isRelayerFrontEnd()){
            Endpoint edp(rel);
            Endpoint fwd_edp(fwd);
            relayersMap.insert ( std::pair<string, Endpoint>(kernelName, edp));
            fwdMap.insert ( std::pair<string, Endpoint>(kernelName, fwd_edp));
            cout<<kernelName<<" <-rel-> "<<edp.rawEndpoints<< " <-fwd-> "<< fwd<<endl;
        }
    }

    vector<int> sockfds;
    // backend connects to all the relayers
    if(this->isRelayerBackend()){

        for (map<string, Endpoint>::iterator it = relayersMap.begin(); it!=relayersMap.end(); it++){
            string kernelName=it->first;
            Endpoint& relayer =it->second;
            if(relayerEndPoint.rawEndpoints != relayer.rawEndpoints){//avoid all kernels that are within me
                map<string, shared_ptr<RelayerCommunicator> >::iterator itfind = remoteRelayerMap.find(kernelName);
                if(itfind == remoteRelayerMap.end()){
                    // try to connect to the relayer
                    int sock= Networking::getInstance().connectToFrontEndRelayer(relayer);
                    if(sock<0){
                        string fwd=fwdMap[kernelName].rawEndpoints;
                        Endpoint fwdEdp(fwd);
                        if(!fwd.empty()){
                            sock= Networking::getInstance().connectToFrontEndForwarder(fwdEdp, "Relay",relayer.rawEndpoints, fwdEdp.rawEndpoints);
                        }
                        if(sock<0)
                            sock= Networking::getInstance().connectToFrontEndForwarder(this->forwarderAddress, "Relay",relayer.rawEndpoints, fwdEdp.rawEndpoints );
                    }
                    shared_ptr<RelayerCommunicator> rel= std::make_shared<TCP_Client_RelayerCommunicator>(sock, this->myMpiManager);
                    remoteRelayerMap.insert ( std::pair<string, shared_ptr<RelayerCommunicator> >(kernelName, rel));
                    sockfds.push_back(sock);
                   // cout<<"}}} establish connection with relayer:"<<relayer.rawEndpoints<<" for kernel "<<kernelName<<endl;
                }
            }
        }
    }

    this->myMpiManager->barrier();

    // inform the manager that I have connected to other relayers
    if(this->isRelayerFrontEnd()){
        managerRelayerCommunicator->sendInt(1);
    }

    // continue
    if(this->isRelayerFrontEnd()){
        try{
            // receive from relayer + forward to local kernels
            this->waitForCommingActionOnFrontEnd();
            cout<<"closing all sockets..."<<endl;
            this->closefrontEndSocket();
            cout<<"## Frontend stopped."<<endl;
        }catch(std::exception &e){
            cout<<"Frontend exp: "<<e.what()<<endl;
        }
    }

    size_t nbr_stop=0;

    if(this->isRelayerBackend()){
        bool isProcess=true;


        while(isProcess){

            shared_ptr<MsgContainer> msgCtn= receiveLocalMessage();
            Message &m = msgCtn->m;
            if (m.endMessage ==1){ // stop from a kernel
                nbr_stop++;
                if(nbr_stop == localKernelMap.size()){ // all kernels sent stop signal
                    isProcess=false;
                }
            }else{
                map<string, shared_ptr<RelayerCommunicator> >::iterator it = remoteRelayerMap.find(m.receiverName);
                if (it != remoteRelayerMap.end()){
                    shared_ptr<RelayerCommunicator> & rel =it->second;
                    forwardMessage(msgCtn, rel );
                }else{
                    cout<<"\t [Failed]: "<< m.senderName <<": No "<<m.receiverName<<" found in remoteRelayerMap"<<endl;
                }
            }
        }
        // close all client connections
        for (auto sock: sockfds){
            close (sock);
        }
        int sock=Networking::getInstance().connectToFrontEndRelayer(relayerEndPoint);
        shared_ptr<RelayerCommunicator> local_rel_comm= std::make_shared<TCP_Client_RelayerCommunicator>(sock, this->myMpiManager);
        Message localMessage;
        localMessage.endMessage=1;
        int header=1;
        local_rel_comm->sendMessage(localMessage,header);
        local_rel_comm->disconnect();
        for (map<string, shared_ptr<RelayerCommunicator> >::iterator it= remoteRelayerMap.begin(); it != remoteRelayerMap.end(); it++)
            it->second->disconnect();
        cout<<"## BackEnd stopped."<<endl;

    }
}

//***


bool MPIRelayer::receiveRemoteMessage(shared_ptr<RelayerCommunicator> &relayerCommunicator, Message &m){
    int header=0;
    relayerCommunicator->receiveMessage(m, header);
    if(m.endMessage ==1) return false;

    //cout<<"\t Message: "<<m.senderName<<"["<< m.senderCoreId<<"] -->"<<m.receiverName<<"["<< m.receiverCoreId<<"] : data size="<<m.data.size()<<endl;
    map<string, shared_ptr<RelayerCommunicator> >::iterator it = localKernelMap.find(m.receiverName);
    if (it != localKernelMap.end()){
        shared_ptr<RelayerCommunicator> & rel =it->second;
        shared_ptr<MpiManager>  interMpiManager = rel->getInterMpiManager();

        if(interMpiManager){//MPI relayer
            if (m.operation == static_cast<int>(MPI_CONDUIT_OP::sendE)
                    || m.operation == static_cast<int>(MPI_CONDUIT_OP::send)
                    || m.operation == static_cast<int>(MPI_CONDUIT_OP::iSend) ){

                if (m.operation == static_cast<int>(MPI_CONDUIT_OP::sendE)){
                    int count = m.data.size();
                    interMpiManager->send<int>(&count, 1,  m.receiverCoreId, m.MPI_tag);
                }

                MPI_Status status;
                interMpiManager->iSend<char>(m.data.data(), m.data.size(), m.receiverCoreId, &this->request, m.MPI_tag);
                if( m.operation != static_cast<int>(MPI_CONDUIT_OP::iSend))
                    MPI_Wait (&request, &status);
            }else if (m.operation == static_cast<int>(MPI_CONDUIT_OP::bCast) || m.operation == static_cast<int>(MPI_CONDUIT_OP::bCastString)){
                if(m.operation == static_cast<int>(MPI_CONDUIT_OP::bCastString)){
                    int count = m.data.size();
                    interMpiManager->send<int>(&count, 1,  m.receiverCoreId, m.MPI_tag);
                }
                MPI_Status status;
                interMpiManager->iSend<char>(m.data.data(), m.data.size(), 0, &this->request, m.MPI_tag);
                MPI_Wait (&request, &status);
            } else if  ( (m.operation == static_cast<int>(MPI_CONDUIT_OP::gather))
                         || (m.operation == static_cast<int>(MPI_CONDUIT_OP::reduce))){
                MPI_Status status;
                interMpiManager->iSend<char>(m.data.data(), m.data.size(), 0, &this->request, m.MPI_tag);
                MPI_Wait (&request, &status);
        }
            // add here the rest of operations: such as bCast ... etc
        }else{ // TCPRelayer
            rel->sendMessage(m);
        }
    }else{
        cout<<"\t [Failed]: No "<<m.receiverName<<" found in localKernelMap"<<endl;
    }
    return true;
}

/*void MPIRelayer::forwardMessage(shared_ptr<MsgContainer> &msgCtn, shared_ptr<RelayerCommunicator> & relayer){

    Message & m = msgCtn->m;
    map<string, int>::iterator it=this->mapColorWhereToRun.find(m.receiverName);
    assert(it != mapColorWhereToRun.end());
    int destColor= it->second;
    if(! msgCtn->isHeaderSent){

        int headerLen = sizeof(int)+ m.senderName.length()+1
                +sizeof(int)+m.receiverName.size()+1
                +sizeof(m.senderCoreId)
                +sizeof(m.receiverCoreId)
                +sizeof(m.operation)
                +sizeof(m.datatype)
                +sizeof(m.MPI_tag)
                +sizeof(m.endMessage)
                +sizeof(int); // size of message data vect

        vector<char> buffer;
        buffer.resize(sizeof(int)+headerLen);// size of raw vect to send

        int originalBuffersize= headerLen + msgCtn->leftToreceive;

        int iData=0;
        memcpy((char*) &buffer[iData], (const char*)&originalBuffersize, sizeof(originalBuffersize));
        iData+=sizeof(originalBuffersize);

        int ss= m.senderName.length()+1;
        memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)m.senderName.c_str(), ss);
        iData+=ss;

        ss= m.receiverName.length()+1;
        memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)m.receiverName.c_str(), ss);
        iData+=ss;

        memcpy((char*) &buffer[iData], (const char*)&m.senderCoreId, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)&m.receiverCoreId, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)&m.operation, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)&m.datatype, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)&m.MPI_tag, sizeof(int));
        iData+=sizeof(int);
        memcpy((char*) &buffer[iData], (const char*)&m.endMessage, sizeof(int));
        iData+=sizeof(int);
        ss= msgCtn->leftToreceive;
        memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));

        relayer->sendRaw((char*) buffer.data(), buffer.size(), destColor);
        vector<char>().swap(buffer);

        msgCtn->isHeaderSent=true;

    }else{
        relayer->sendRaw((char*) msgCtn->m.data.data(), msgCtn->m.data.size(), destColor);
        //relayer->send(msgCtn->m.data, destColor);
        vector<char>().swap(msgCtn->m.data);
    }

    //relayer->sendMessage(m, destColor);
}*/

void MPIRelayer::forwardMessage(shared_ptr<MsgContainer> &msgCtn, shared_ptr<RelayerCommunicator> & relayer){
    map<string, int>::iterator it=this->mapColorWhereToRun.find(msgCtn->m.receiverName);
    assert(it != mapColorWhereToRun.end());
    int destColor= it->second;

    if(msgCtn->leftToreceive==0 && msgCtn->isHeaderSent){
        relayer->sendMessage(msgCtn->m, destColor);
    }
}

/*shared_ptr<MsgContainer>  MPIRelayer::receiveLocalMessage(){

    MPI_Status  status_pb;
    MPI_Probe(MPI_ANY_SOURCE, MPI_RELAYER_GLOBAL_TAG, this->globalMpiManager->getGlobalCommunicator(), &status_pb);
    int mpi_src =status_pb.MPI_SOURCE;

    map<int, shared_ptr<MsgContainer> >::iterator it =this->messagesContainerMap.find(mpi_src);
    shared_ptr<MsgContainer> msgCtn;
    if(it == this->messagesContainerMap.end()){
        msgCtn = std::make_shared<MsgContainer>();
        this->messagesContainerMap.insert(std::pair<int, shared_ptr<MsgContainer> > (mpi_src, msgCtn));
    }else{
        msgCtn= it->second;
    }
    assert(msgCtn);

    if(msgCtn->leftToreceive==0){
        // receive new Message header
        Messaging::getInstance().receiveMessageHeader(this->globalMpiManager, msgCtn->m, mpi_src, MPI_RELAYER_GLOBAL_TAG);
        // recev data size
        this->globalMpiManager->receive<int>(& msgCtn->leftToreceive, 1, mpi_src, MPI_RELAYER_GLOBAL_TAG);
        msgCtn->isHeaderSent=false;
    }else{
        // receive message data progressively
        msgCtn->isHeaderSent=true;
        int len=4096;
        if(msgCtn->leftToreceive < len)
            len= msgCtn->leftToreceive;
        msgCtn->m.data.resize(len);
        MPI_Request request;
        this->globalMpiManager->iRecv<char>(msgCtn->m.data.data(), len, mpi_src, &request, MPI_RELAYER_GLOBAL_TAG);
        MPI_Status status;
        MPI_Wait (&request, &status);
        msgCtn->leftToreceive -=len;
    }
    return msgCtn;
}*/

shared_ptr<MsgContainer>  MPIRelayer::receiveLocalMessage(){

    MPI_Status  status_pb;
    MPI_Probe(MPI_ANY_SOURCE, MPI_RELAYER_GLOBAL_TAG, this->globalMpiManager->getGlobalCommunicator(), &status_pb);

    int mpi_src =status_pb.MPI_SOURCE;
    shared_ptr<MsgContainer> msgCtn = std::make_shared<MsgContainer>();
    // receive new Message header
    Messaging::getInstance().receiveMessageHeader(this->globalMpiManager, msgCtn->m, mpi_src, MPI_RELAYER_GLOBAL_TAG);
    // recev data size
    this->globalMpiManager->receive<int>(& msgCtn->leftToreceive, 1, mpi_src, MPI_RELAYER_GLOBAL_TAG);
    msgCtn->isHeaderSent=false;


    // receive message data progressively
    msgCtn->isHeaderSent=true;
    int len= msgCtn->leftToreceive;
    msgCtn->m.data.resize(len);
    MPI_Request request;
    this->globalMpiManager->iRecv<char>(msgCtn->m.data.data(), len, mpi_src, &request, MPI_RELAYER_GLOBAL_TAG);
    MPI_Status status;
    MPI_Wait (&request, &status);
    msgCtn->leftToreceive -=len;

    //Messaging::getInstance().receive(this->globalMpiManager, msgCtn->m, mpi_src, MPI_RELAYER_GLOBAL_TAG);
    return msgCtn;
}

void MPIRelayer::receiveRelayersList( map<string, Endpoint> &relayersMap, map<string, Endpoint> &fwdMap){

    int client_len;
    struct sockaddr_in client_address;
    int result;
    fd_set readfds;
    FD_ZERO (&readfds);

    /// get relayers list from  the Manager
    int fd;
    int nread;
    //testfds = readfds;
    readfds = active_fd_set;
    cout<<"server waiting"<<endl;
    result = select(max_sd+1, &readfds, (fd_set *)0, (fd_set *)0, (struct timeval *) 0);
    if(result < 1) {
        cerr<<"error on select"<<endl;
        exit(1);
    }

    for(fd = 0; fd <= max_sd; fd++) {
        if(FD_ISSET(fd,&readfds)) {
            if(fd == frontEnd_sockfd) {
                client_len = sizeof(client_address);
                int client_sockfd = accept(frontEnd_sockfd, (struct sockaddr *)&client_address, (socklen_t*)&client_len);
                FD_SET(client_sockfd, &active_fd_set);
                //highest file descriptor number, need it for the select function
                if(client_sockfd > max_sd){
                    max_sd = client_sockfd;
                }
                managerRelayerCommunicator= std::make_shared<TCP_Client_RelayerCommunicator>(client_sockfd, this->myMpiManager);
                // recv the manager URL
                managerRelayerCommunicator->receiveString(this->managerPortUrl);
                //recv number of relayers
                int numberofrelayers=-1;
                managerRelayerCommunicator->receiveInt(numberofrelayers);
                // send the relayers list
                for(int i=0; i< numberofrelayers; i++){
                    //broadCastToRemote a string
                    string kernelName;
                    managerRelayerCommunicator->receiveString(kernelName);
                    string relayerUrl;
                    managerRelayerCommunicator->receiveString(relayerUrl);
                    string fwdUrl;
                    managerRelayerCommunicator->receiveString(fwdUrl);
                    //managerRelayerCommunicator->disconnect();
                    Endpoint rel_endpoint(relayerUrl);
                    Endpoint fwd_endpoint(fwdUrl);
                    relayersMap.insert ( std::pair<string, Endpoint>(kernelName, rel_endpoint));
                    fwdMap.insert ( std::pair<string, Endpoint>(kernelName, fwd_endpoint));

                }
            }else{
                ioctl(fd, FIONREAD, &nread);
                if(nread == 0) {
                    close(fd);
                    FD_CLR(fd, &active_fd_set);
                    //cout<<"removing client on fd "<< fd<<endl;
                }else{
                    // nothing here
                }
            }
        }
    }
}

int MPIRelayer::connectToForwarder(const string &relayerEndpoint){

    int sockfd=Networking::getInstance().connectToFrontEndRelayer(this->forwarderAddress);
    if(sockfd <0)
        return sockfd;

    string relayerName="Relayer_"+relayerEndPoint.rawEndpoints;
    int needManager=1;
    int len=relayerName.size()+1;
    int len1=relayerEndpoint.size()+1;
    int buffsize=sizeof(needManager)*3
            +len+len1;

    vector<char>  buffer;
    buffer.resize(buffsize);
    int idData=0;
    memcpy((char*) &buffer[idData], (const char*)&needManager, sizeof(needManager));
    idData+=sizeof(needManager);
    memcpy((char*) &buffer[idData], (const char*)&len, sizeof(len));
    idData+=sizeof(len);
    std::copy(relayerName.c_str(), relayerName.c_str()+len, &buffer[idData]);
    idData+=len;
    memcpy((char*) &buffer[idData], (const char*)&len1, sizeof(len1));
    idData+=sizeof(len1);
    std::copy(relayerEndpoint.c_str(), relayerEndpoint.c_str()+len1, &buffer[idData]);

    Networking::getInstance().sock_write(sockfd, buffer.data(), buffer.size());
    vector<char>().swap(buffer);

    return sockfd;
}


void MPIRelayer::waitForCommingActionOnFrontEnd(){

    int client_sockfd;
    int client_len;
    struct sockaddr_in client_address;
    int result;
    fd_set readfds;
    FD_ZERO (&readfds);
    bool isProcess=true;
    vector<int> descriptors;
    while(isProcess) {
        int fd;
        int nread;
        readfds = active_fd_set;
        result = select(max_sd+1, &readfds, (fd_set *)0, (fd_set *)0, (struct timeval *) 0);
        if(result < 1) {
            cerr<<"error on select"<<endl;
            exit(1);
        }

        for(fd = 0; fd <= max_sd; fd++) {
            if(FD_ISSET(fd,&readfds)) {
                if(fd == frontEnd_sockfd) {

                    client_len = sizeof(client_address);
                    client_sockfd = accept(frontEnd_sockfd, (struct sockaddr *)&client_address, (socklen_t*)&client_len);
                    FD_SET(client_sockfd, &active_fd_set);
                    //highest file descriptor number, need it for the select function
                    if(client_sockfd > max_sd){
                        max_sd = client_sockfd;
                    }
                    descriptors.push_back(client_sockfd);
                }else{
                    ioctl(fd, FIONREAD, &nread);
                    if(nread == 0) {
                        close(fd);
                        FD_CLR(fd, &active_fd_set);
                        descriptors.erase(std::remove(descriptors.begin(), descriptors.end(), fd), descriptors.end());
                    }else{

                        shared_ptr<RelayerCommunicator> relayerCommunicator= std::make_shared<TCP_Client_RelayerCommunicator>(fd, this->myMpiManager);
                        Message m;
                        bool isEnd = ! receiveRemoteMessage(relayerCommunicator, m);
                        if(isEnd) {
                            isProcess=false;
                            close (fd);
                            for(auto remotefd: descriptors){
                                close (remotefd);
                            }
                        }
                    }
                }
            }
        }
    }
    // wait the last iSend if it exists

    if (request != MPI_REQUEST_NULL){
        MPI_Status status;
        MPI_Wait (&request, &status);
    }
    //cout<<"waitForCommingActionOnFrontEnd ->end()"<<endl;
}




bool MPIRelayer::closefrontEndSocket(){
    managerRelayerCommunicator->disconnect();

    int fd;
    for(fd = 0; fd < max_sd; fd++) {
        if(FD_ISSET(fd,&active_fd_set)) {
            close(fd);
            FD_CLR(fd, &active_fd_set);
        }
    }
    close(frontEnd_sockfd);
    return true;
}
