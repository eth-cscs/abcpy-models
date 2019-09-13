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

#include "musclehpc/mediator/mpisocketserver.h"

/********************************************************* Manager *************************************************/

Manager::Manager(int kernelsSize,
                 map<string, shared_ptr<MpiManager> > &interMpiManagersMap,
                 map<string, int> &mapColorWhereToRun,
                 shared_ptr<MpiManager> &localMpiManager,
                 shared_ptr<MpiManager> &globalMpiManager, string forwarder, bool useTCP):
    kernelsSize(kernelsSize),
    interMpiManagersMap(interMpiManagersMap),
    mapColorWhereToRun(mapColorWhereToRun),
    localMpiManager(localMpiManager),
    globalMpiManager(globalMpiManager),
    forwarder(forwarder),
    useTCP(useTCP){

    if(!this->useTCP && !Networking::getInstance().canOpenMPIPort(localMpiManager->getGlobalCommunicator())){
        this->useTCP=true;
    }
}

Manager::~Manager(){}

void Manager::run(){

    try{
        shared_ptr<SocketServer> server;
        if(!useTCP)
            server= std::move(std::make_shared<MpiSocketManager>(localMpiManager, interMpiManagersMap, mapColorWhereToRun));
        else
            server= std::move(std::make_shared<TCPSocketManager>(localMpiManager, interMpiManagersMap, mapColorWhereToRun, this->forwarder));

        server->setKernelNumbers(kernelsSize);

        server->initConnection();
        server->handShake();
        server->disconnect();
    }catch (std::exception & e){
        cout<<"--> exception "<<e.what()<<endl;
    }
    cout<<"# Manager ended."<<endl;
}

//********************************************************* SocketServer ********************************
SocketServer::SocketServer()
{
    //myLocalInfo=0;
    memset(this->port_name, 0, MPI_MAX_PORT_NAME);
}
SocketServer::~SocketServer()
{

}

void SocketServer::setKernelNumbers(int kernelNumbers){
    this->kernelNumbers=kernelNumbers;
}
void SocketServer::setIamManager(bool isManager){
    this->amImanager=isManager;
}

map<string, shared_ptr<KernelConnectionInfo> > &SocketServer::getRemoteKernelsMap(){
    return this->remoteKernelsMap;
}

shared_ptr<KernelConnectionInfo> SocketServer::getMyLocalInfo(){
    return this->myLocalInfo;
}

//**********************************************************************************************/
//************************************ SocketManager *******************************************/
//**********************************************************************************************/

SocketManager::SocketManager(shared_ptr<MpiManager> &localMpiManager, string forwarder):myMpiManager(localMpiManager) {
    if(!forwarder.empty())
        this->forwarder.parse(forwarder);
}

SocketManager::~SocketManager(){

}

void SocketManager::initConnection(){

    bool debug=false;
    int open_res =this->openMPISocket();
    if(open_res != MPI_SUCCESS) return ;
    cout<<"Manager available at: \""<< port_name<<"\""<<endl;
    this->informAllmyCamaradesKernelsOfMyURL();

    if (debug && this->getMpiManager()->isMainProcessor()) cout <<"wainting client connection:"<<endl;
    vector<shared_ptr<KernelConnectionInfo> > queuedConnection;

    this->acceptFromAllKernels(queuedConnection, this->kernelNumbers);

    if (debug && this->getMpiManager()->isMainProcessor()) cout <<"all kernels are connected"<<endl;
    // int remoteKernelColor=0;
    for (size_t i=0; i<queuedConnection.size();i++){
        shared_ptr<KernelConnectionInfo> &clientinfo = queuedConnection.at(i);
        /*int root=0;
        /clientinfo->broadCastMyInfoFromRemote(root);*/
        clientinfo->setcolor(this->mapColorWhereToRun[clientinfo->getName()]);
        if (debug && this->getMpiManager()->isMainProcessor()) {
            string kenpt=(!clientinfo->getMyMPIUrl().empty())? clientinfo->getMyMPIUrl():clientinfo->getMyRelayerUrl();
            cout <<"received: (kname="<< clientinfo->getName() <<"), (Url="<< kenpt<<")"<<endl;
        }
        this->remoteKernelsMap.insert ( std::pair<string, shared_ptr<KernelConnectionInfo> >(clientinfo->getName(), clientinfo));
        //remoteKernelColor++;
    }
    queuedConnection.clear();

    if (debug && this->getMpiManager()->isMainProcessor()) cout <<"Broadcast KernelsMapInfo to all remote kernels "<<endl;
    int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;

    for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        assert(ptr!=0);
        ptr->broadCastMapToRemote(remoteKernelsMap, root);
    }
}



void SocketManager::handShake(){
    bool debug=false;
    /// ------------- HandShake with Relayers ------------------
    this->handShakeWithRelayers();

    /// ------------- HandShake with kernels ------------------
    if ( Networking::getInstance().canOpenMPIPort(myMpiManager->getGlobalCommunicator())){
        try{
            for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
                /// choose a kernel
                shared_ptr<KernelConnectionInfo> & client= it->second;
                /// tell him to connect to its remote kernels
                sendRemoteOrder(HANDSHAKE_TAG::CONNECT, client);
                if (debug && this->getMpiManager()->isMainProcessor())
                    cout<<"receive the list of requested kernels to which "<< client->getName() <<" wants to connect"<<endl;
                int root= 0;
                vector<string> requestedkids;
                int length=0 ;
                /// I receive from the kernel its remote kernels number
                client->broadCastIntFromRemote(length,root);
                if (debug && this->getMpiManager()->isMainProcessor())
                    cout<< client->getName()<<"  wants to connect to "<<length<<" kernels" <<endl;
                /// I receive from the kernel all its remote kernels string names
                for (int i=0; i< length; i++){
                    string kid;
                    client->broadCastStringFromRemote(kid,root);
                    //getInterCommUtils().broadCastStringFromRemote(kid,root, client->getInterMpiManager());
                    requestedkids.push_back(kid);
                }
                if (debug && this->getMpiManager()->isMainProcessor()){
                    cout<<"     Here is the kernsl name: ";
                    for (int i=0; i< length; i++){
                        cout<< requestedkids.at(i)<<"  ";
                    }
                    cout<<endl;
                }
                /// inform the requested remote kernels to listen
                for (size_t i=0; i< requestedkids.size(); i++){
                    shared_ptr<KernelConnectionInfo> &requestedclient= this->remoteKernelsMap[requestedkids.at(i)];
                    string remoteMpiUrlport=requestedclient->getMyMPIUrl();
                    if(!remoteMpiUrlport.empty()){// remote kernel has opened an MPI port
                        sendRemoteOrder(HANDSHAKE_TAG::LISTEN, requestedclient);
                    }
                }

            }// end while
            if (debug && this->getMpiManager()->isMainProcessor()) cout<<" ending handShake ..." <<endl;
            for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
                shared_ptr<KernelConnectionInfo> & client= it->second;
                sendRemoteOrder(HANDSHAKE_TAG::END, client);

            }
            if (debug && this->getMpiManager()->isMainProcessor())
                cout<<"# HandShake ended." <<endl;

        }catch( const std::exception& e ) { // reference to the base of a polymorphic object
            std::cout << "[manager] "<<e.what(); // information from error printed
        }
    }


}


void SocketManager::handShakeWithRelayers(){

    vector<string> relayerUrls;
    map<string, string> relayerUrlsFwd; //<relayer, fwd>
    for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        assert(ptr!=0);
        if(! ptr->getMyRelayerUrl().empty())
            relayerUrls.push_back(ptr->getMyRelayerUrl());
        relayerUrlsFwd.insert ( std::pair<string, string>(ptr->getMyRelayerUrl(), ptr->getMyForwarderUrl()));
    }

    if(relayerUrls.empty()) return;

    // remove duplicated relayers
    std::sort( relayerUrls.begin(), relayerUrls.end() );
    relayerUrls.erase( std::unique( relayerUrls.begin(), relayerUrls.end() ), relayerUrls.end() );
    // connection to relayers
    map<string, shared_ptr<RelayerCommunicator> > map_relayer_sockfd; //<relayer,sockfd>
    for(auto r: relayerUrls){
        int sock=-1;
        assert (! r.empty());
        if(this->getMpiManager()->isMainProcessor()){
            Endpoint edp(r);
            sock= Networking::getInstance().connectToFrontEndRelayer(edp);
            if(sock<0){
                assert(! relayerUrlsFwd[r].empty());
                Endpoint fwd(relayerUrlsFwd[r]);
                if(! fwd.rawEndpoints.empty())
                    sock=Networking::getInstance().connectToFrontEndForwarder(fwd, "Manager",edp.rawEndpoints, fwd.rawEndpoints);
                if(sock < 0)
                    sock=Networking::getInstance().connectToFrontEndForwarder(this->forwarder, "Manager",edp.rawEndpoints, fwd.rawEndpoints);
            }
            assert(sock>=0);
        }
        shared_ptr<RelayerCommunicator> rel= std::make_shared<TCP_Client_RelayerCommunicator>(sock, this->getMpiManager());
        map_relayer_sockfd.insert ( std::pair<string, shared_ptr<RelayerCommunicator> >(r, rel));

    }
    // send a map to each relayer: the managerURL + <kernelName,relayer> map
    for (map<string, shared_ptr<RelayerCommunicator> >::iterator it= map_relayer_sockfd.begin(); it!=map_relayer_sockfd.end(); ++it){
        //string relayerURL = it->first;
        shared_ptr<RelayerCommunicator> & relayerCommunicator = it->second;
        string mngUrl(port_name);
        if(this->getMpiManager()->isMainProcessor()){
            //send Manager url
            relayerCommunicator->sendString(mngUrl);
            //send number of relayers
            int len=remoteKernelsMap.size();
            relayerCommunicator->sendInt(len);
        }
        // send the the <kernelName,relayer> map
        for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
            shared_ptr<KernelConnectionInfo> & ptr= it->second;
            assert(ptr!=0);
            if(this->getMpiManager()->isMainProcessor()){
                relayerCommunicator->sendString(ptr->getName());
                relayerCommunicator->sendString(ptr->getMyRelayerUrl());
                relayerCommunicator->sendString(ptr->getMyForwarderUrl());
            }
        }
    }

    // receive from them value 1.
    for (map<string, shared_ptr<RelayerCommunicator> >::iterator it= map_relayer_sockfd.begin(); it!=map_relayer_sockfd.end(); ++it){
        shared_ptr<RelayerCommunicator> & relayerCommunicator = it->second;
        //string kernelName = it->first;
        if(this->getMpiManager()->isMainProcessor()){
            int akk;
            relayerCommunicator->receiveInt(akk);
            relayerCommunicator->disconnect();
        }
    }

}

void SocketManager::disconnect(){

    try{
        remoteKernelsMap.clear();
        this->closeMPISocket();
        cout<<" Registration finished" <<endl;
    }catch( const std::exception& e ) { // reference to the base of a polymorphic object
        std::cout << e.what(); // information from length_error printed
    }
}

void SocketManager::setMyRemoteKernels(vector<string> myRemoteKerelsIDs){}

void SocketManager::sendRemoteOrder(HANDSHAKE_TAG order, shared_ptr<KernelConnectionInfo> &clientinfo){

    int root= (this->myMpiManager->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
    int myorder= static_cast<int>(order);
    clientinfo->broadCastIntToRemote(myorder, root);
    //clientinfo->getInterMpiManager()->bCast<int>(&myorder,1, root);
}

void SocketManager::informAllmyCamaradesKernelsOfMyURL(){
    int root= (this->myMpiManager->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;

    for(map<string, shared_ptr<MpiManager> >::iterator it = sameBinaryInterMpiManager.begin(); it != sameBinaryInterMpiManager.end(); it++) {
        //for(size_t i=0; i< sameBinaryInterMpiManager.size(); i++){
        shared_ptr<MpiManager> & inter = it->second; //sameBinaryInterMpiManager.at(i);
        //broadCastToRemote
        string managerPortUrl = this->port_name;
        int length = (int) managerPortUrl.size();
        inter->bCast(&length, 1, root);
        char* buffer = new char[length+1];
        std::copy(managerPortUrl.c_str(), managerPortUrl.c_str()+length+1, buffer);
        inter->bCast(buffer, length+1, root);
        delete [] buffer;

    }

}

//********************************** MpiSocketManager *******************************************

MpiSocketManager::MpiSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, int> mapColorWhereToRun):
    SocketManager(localMpiManager){

    this->mapColorWhereToRun.insert(mapColorWhereToRun.begin(), mapColorWhereToRun.end());

}

MpiSocketManager::MpiSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, shared_ptr<MpiManager> > interMpiManagers, map<string, int> mapColorWhereToRun):
    SocketManager(localMpiManager){
    this->sameBinaryInterMpiManager.insert(interMpiManagers.begin(), interMpiManagers.end());
    this->mapColorWhereToRun.insert(mapColorWhereToRun.begin(), mapColorWhereToRun.end());
}


shared_ptr<MpiManager> & MpiSocketManager::getMpiManager() {
    return this->myMpiManager;
}

MpiSocketManager::~MpiSocketManager(){

}

int  MpiSocketManager::openMPISocket(){
    int res=MPI_ERR_UNKNOWN;
    string endpoint = Networking::getInstance().openfrontendMPI(myMpiManager->getGlobalCommunicator());
    if(!endpoint.empty()){
        strncpy (port_name, endpoint.c_str(), endpoint.length()+1);
        res=MPI_SUCCESS;
    }
    return res;
}

shared_ptr<KernelConnectionInfo> MpiSocketManager::acceptConnection(){
    MPI_Comm client_interComm;
    MPI_Comm_set_errhandler(myMpiManager->getGlobalCommunicator(), MPI_ERRORS_RETURN);
    int res_connect=MPI_Comm_accept(port_name, MPI_INFO_NULL, myMpiManager->bossId(), myMpiManager->getGlobalCommunicator(),  &client_interComm);
    MPI_Comm_set_errhandler(this->myMpiManager->getGlobalCommunicator(), MPI_ERRORS_ARE_FATAL);
    //if( debug && this->getMpiManager()->isMainProcessor()) cout<<" get connection "<< currenConnectionsNumber<<endl;

    if( res_connect  != MPI_SUCCESS )  return 0;

    shared_ptr<MpiManager> interMpiManager= std::make_shared<MpiManager>();
    interMpiManager->initOnlyCommunicator(client_interComm);
    shared_ptr<KernelConnectionInfo>  clientinfo= make_shared<KernelMPIConnectionInfo>(interMpiManager);
    return clientinfo;
}

void MpiSocketManager::acceptFromAllKernels(vector<shared_ptr<KernelConnectionInfo> > &queuedConnection, int numberOfConnections){

    if (this->myMpiManager->isMainProcessor()) cout <<"wainting client connection:"<<endl;
    int currenConnectionsNumber=0;
    while ( currenConnectionsNumber < numberOfConnections) {//remoteKernelsIDs.size()) {

        shared_ptr<KernelConnectionInfo>  clientinfo= acceptConnection();
        assert(clientinfo);
        int root=0;
        clientinfo->broadCastMyInfoFromRemote(root);
        bool isKernelAlreadyConnected=false;
        for(shared_ptr<KernelConnectionInfo>  & cl: queuedConnection){
            if(clientinfo->getName() == cl->getName()){
                cerr<<"[Warning] "<<clientinfo->getName()<<" is already connected -> will be ignored !"<<endl;
                isKernelAlreadyConnected=true;
                break;
            }
        }
        if(! isKernelAlreadyConnected){
            queuedConnection.push_back(clientinfo);
            currenConnectionsNumber++;
        }
    }//end while accept
}

void MpiSocketManager::closeMPISocket(){
    if(strlen(port_name)>0)
        MPI_Close_port(port_name);
}


//********************************** TCPSocketManager *******************************************

TCPSocketManager::TCPSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, int> mapColorWhereToRun, string forwarder):
    SocketManager(localMpiManager, forwarder){
    this->mapColorWhereToRun.insert(mapColorWhereToRun.begin(), mapColorWhereToRun.end());
}

TCPSocketManager::TCPSocketManager(shared_ptr<MpiManager> localMpiManager, map<string, shared_ptr<MpiManager> > interMpiManagers,
                                   map<string, int> mapColorWhereToRun, string forwarder):
    SocketManager(localMpiManager, forwarder){
    this->sameBinaryInterMpiManager.insert(interMpiManagers.begin(), interMpiManagers.end());
    this->mapColorWhereToRun.insert(mapColorWhereToRun.begin(), mapColorWhereToRun.end());
}


shared_ptr<MpiManager> &TCPSocketManager::getMpiManager() {
    return this->myMpiManager;
}

TCPSocketManager::~TCPSocketManager(){
}


int TCPSocketManager::openMPISocket(){

    int res_open=MPI_ERR_UNKNOWN;

    if (this->myMpiManager->isMainProcessor()){ // only on root
        Endpoint frontend("0.0.0.0:5001");
        sockfd= Networking::getInstance().openfrontEndSocket(frontend);
        string suffix;
        if (!this->forwarder.rawEndpoints.empty()){
            suffix= "/"+this->forwarder.rawEndpoints;
        }
        string portname = frontend.rawEndpoints+suffix;
        strncpy (this->port_name, portname.c_str(), portname.length()+1);
        res_open=(sockfd > 0 )? MPI_SUCCESS: MPI_ERR_UNKNOWN;
    }

    return res_open;
}


void TCPSocketManager::acceptFromAllKernels(vector<shared_ptr<KernelConnectionInfo> > &queuedConnection, int numberOfConnections){

    //bool debug= false;
    if (this->myMpiManager->isMainProcessor()) {

        int client_len;
        struct sockaddr_in client_address;
        int result;
        fd_set readfds, active_fd_set;
        FD_ZERO (&readfds);
        FD_ZERO (&active_fd_set);
        FD_SET (sockfd, &active_fd_set);
        int max_sd =sockfd;

        int currenConnectionsNumber=0;
        while ( currenConnectionsNumber < numberOfConnections) {//remoteKernelsIDs.size()) {

            int fd;
            readfds = active_fd_set;
            result = select(max_sd+1, &readfds, (fd_set *)0, (fd_set *)0, (struct timeval *) 0);
            if(result < 1) {
                cerr<<"error on select"<<endl;
                exit(1);
            }

            for(fd = 0; fd <= max_sd; fd++) {
                if(FD_ISSET(fd,&readfds)) {
                    if(fd == sockfd) {
                        client_len = sizeof(client_address);
                        int client_sockfd = accept(sockfd, (struct sockaddr *)&client_address, (socklen_t*)&client_len);
                        if(client_sockfd <= 0) continue;

                        shared_ptr<KernelConnectionInfo> clientinfo= std::make_shared <KernelTCPConnectionInfo>(client_sockfd, myMpiManager);
                        assert(clientinfo);
                        int root=0;
                        clientinfo->broadCastMyInfoFromRemote(root);
                        bool isKernelAlreadyConnected=false;
                        for(shared_ptr<KernelConnectionInfo>  & cl: queuedConnection){
                            if(clientinfo->getName() == cl->getName()){
                                cerr<<"[Warning] "<<clientinfo->getName()<<" is already connected -> will be ignored !"<<endl;
                                isKernelAlreadyConnected=true;
                                clientinfo.reset();
                                break;
                            }
                        }
                        if(! isKernelAlreadyConnected){
                            queuedConnection.push_back(clientinfo);
                            currenConnectionsNumber++;
                            FD_SET(client_sockfd, &active_fd_set);
                            if(client_sockfd > max_sd){
                                max_sd = client_sockfd;
                            }
                        }
                        if( currenConnectionsNumber >= numberOfConnections)
                            return;
                    }
                }
            }

        }//end while accept
    }
}

shared_ptr<KernelConnectionInfo> TCPSocketManager::acceptConnection(){


    if (this->myMpiManager->isMainProcessor()){ // only on root
        struct sockaddr  cli_addr;
        int newsockfd, clilen;
        clilen = sizeof(cli_addr);
        newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr,  (socklen_t*)&clilen);
        if (newsockfd < 0) {
            cerr<<"ERROR on accept"<<endl;
            return 0;
        }
        shared_ptr<KernelConnectionInfo> clientinfo= std::make_shared <KernelTCPConnectionInfo>(newsockfd, myMpiManager);
        return clientinfo;
    }
    return 0;
}

void TCPSocketManager::closeMPISocket(){
    if(sockfd>0)
        close(sockfd);
}


//**********************************************************************************************/
//************************************ SocketKernel *******************************************/
//**********************************************************************************************/

SocketKernel::SocketKernel(shared_ptr<MpiManager> globalMpiManager,
                           shared_ptr<MpiManager> & localMpiManager,
                           shared_ptr<MpiManager> interMpiManager,
                           shared_ptr<RelayerCommunicator> &relayerCommunicator,
                           string myKernelName, string forwarderAddress,
                           string managerPortUrl):
    globalMpiManager(globalMpiManager), myMpiManager(localMpiManager),
    interMpiManager(interMpiManager), relayerCommunicator(relayerCommunicator),
    myKernelName(myKernelName), managerPortUrl(managerPortUrl){
    this->forwarderAddress.parse(forwarderAddress);

}


SocketKernel::~SocketKernel()
{
}

void SocketKernel::initConnection(){

    bool debug=false;
    bool isBossId=this->getMpiManager()->isMainProcessor();
    string myMPIUri;
    //-------------------- connect to my relayer if it exits --------------------
    string relayerUrl, fwdUrl;
    connectWithMyRelayer(relayerUrl, fwdUrl);
    //--------------------
    if ( debug && isBossId) cout<<"--- openMPISocket --------\n";
    int res_open=this->openMPISocket();
    if(res_open == MPI_SUCCESS){
        myMPIUri=this->port_name;
    }
    if (debug && isBossId) cout<<"--- checkoutManagerPortUrl --------\n";
    this->checkoutManagerPortUrl();
    if (debug && isBossId) cout<<"--- connectToManager: ---------\n";

    myLocalInfo = this->connectToManager(); //new KernelMPIConnectionInfo(inter);
    assert(myLocalInfo);
    myLocalInfo->setName(this->myKernelName);
    myLocalInfo->setMyMPIUrl(myMPIUri);
    myLocalInfo->setMyRelayerUrl(relayerUrl);
    myLocalInfo->setMyForwarderUrl(fwdUrl);
    int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
    myLocalInfo->broadCastMyInfoToRemote(root);

    if (debug && isBossId) cout<<"--- receive from server manager all URL of connected kernels ---------------\n";
    int remoteRoot=0;
    myLocalInfo->broadCastMapFromRemote( this->remoteKernelsMap, remoteRoot);

    if (debug && isBossId) cout<<"--- display received data --------------------------------------------------\n";
    for( it_map_comm_string_ConnInfo it = this->remoteKernelsMap.begin(); it!= this->remoteKernelsMap.end();it++){
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        if (this->getMpiManager()->isMainProcessor()) {
            cout<<"["<< myLocalInfo->getName()<< "]: "<<ptr->getName()<<" located at: \n\t"<<ptr->getMyMPIUrl()<<"\n\t"<<ptr->getMyRelayerUrl()<<endl;
        }
    }
    if (debug && isBossId) cout<<"--- find my assigned color -------------------------------------------------\n";
    it_map_comm_string_ConnInfo it = remoteKernelsMap.find(myLocalInfo->getName());
    if( it!= remoteKernelsMap.end()) {
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        myLocalInfo->setcolor(ptr->getcolor());
    }else if (this->getMpiManager()->isMainProcessor()) {
        cout<<"-> my color not found "<<endl;
    }
    if (this->getMpiManager()->isMainProcessor()) {
        cout<<"-> my color =  "<<myLocalInfo->getcolor()<<endl;
    }
}

void SocketKernel::handShake(){

    if (! Networking::getInstance().canOpenMPIPort(myMpiManager->getGlobalCommunicator())) return ;
    bool isDone=false;
    while (! isDone){
        HANDSHAKE_TAG order=this->getRemoteOrder();
        int remoteOreder= static_cast<int>(order);
        if(order == HANDSHAKE_TAG::LISTEN ){ // accept connect from one kernel
            acceptConnectionfromRemoteKernel();
        }else if (order == HANDSHAKE_TAG::CONNECT ){ // connet to all my remote kernels
            connectToMyRemoteKernels();
        }else if (order == HANDSHAKE_TAG::END){ // handshake ended
            isDone=true;
        }else{ //unknow order
            if (this->getMpiManager()->isMainProcessor()) cout<<" Unknown handshake order: "<< remoteOreder <<endl;
            exit (EXIT_FAILURE);
        }
    }//end while

}

void SocketKernel::acceptConnectionfromRemoteKernel(){
    //accept
    MPI_Comm client_interComm;
    MPI_Comm_set_errhandler(myMpiManager->getGlobalCommunicator(), MPI_ERRORS_RETURN);
    int res_connect=MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, this->getMpiManager()->getGlobalCommunicator(),  &client_interComm);
    MPI_Comm_set_errhandler(this->myMpiManager->getGlobalCommunicator(), MPI_ERRORS_ARE_FATAL);
    if( res_connect  != MPI_SUCCESS )  return ;

    //if (this->getMpiManager()->isMainProcessor()) cout<<"YYYY Acepting MPI YYYY"<<endl;
    shared_ptr<MpiManager> interMpiManager = std::make_shared<MpiManager>();
    interMpiManager->initOnlyCommunicator(client_interComm);

    //receive
    string kid;
    int root=0;
    getInterCommUtils().broadCastStringFromRemote(kid, root, interMpiManager);

    // find the remote kernel who connected to me
    it_map_comm_string_ConnInfo it = remoteKernelsMap.find(kid);
    if( it!= remoteKernelsMap.end()) {
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        ptr->moveInterMpiManager(interMpiManager);
        if (this->getMpiManager()->isMainProcessor()) cout<<myLocalInfo->getName()<<": got handshake from: "<< kid<<endl;

    }else if (this->getMpiManager()->isMainProcessor()) {
        cout<<"-> sender not found "<<endl;
    }
}

void SocketKernel::connectToMyRemoteKernels(){

    bool debug=false;
    // bool isROOT=this->getMpiManager()->isMainProcessor();

    int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
    // send to MANAGER my remote kernel list: only those that ihave not contact them

    vector<shared_ptr<KernelConnectionInfo> > notYetContacted;
    for (size_t i=0; i< myRemoteKerelsIDs.size(); i++){
        string kid= myRemoteKerelsIDs.at(i);
        shared_ptr<KernelConnectionInfo> & remoteKernel=remoteKernelsMap[kid];
        if(! remoteKernel->getInterMpiManager() && !remoteKernel->getMyMPIUrl().empty()){
            notYetContacted.push_back(remoteKernel);
        }else{
            // if (this->getMpiManager()->isMainProcessor()) cout << this->getMyLocalInfo()->getId()<< ": I am altready connected to "<< remoteKernel->getId()<<endl;
        }
    }

    int length=notYetContacted.size();
    if (debug && this->getMpiManager()->isMainProcessor()) cout<<myLocalInfo->getName()<<" wants to connect to " << length <<" kernels"<<endl;
    //send size of kernel to connect to the manager
    //myLocalInfo->getInterMpiManager()->bCast<int>(&length, 1, root);
    myLocalInfo->broadCastIntToRemote(length, root);
    // send to manager  their names

    stringstream ss;
    ss<<this->getMyLocalInfo()->getName()<< "kernels are [";
    for (size_t i=0; i< notYetContacted.size(); i++){
        shared_ptr<KernelConnectionInfo> & remoteKernel=notYetContacted.at(i);
        assert(remoteKernel && (! remoteKernel->getInterMpiManager()));
        //getInterCommUtils().broadCastStringToRemote(remoteKernel->getId(), root,  myLocalInfo->getInterMpiManager());
        string idRemote=remoteKernel->getName();
        ss<<idRemote<<" ";
        myLocalInfo->broadCastStringToRemote(idRemote, root);
    }
    ss<<"]"<<endl;
    if (debug && this->getMpiManager()->isMainProcessor()) cout<<ss.str();
    // wait until all kernels are notified <-- not needed
    //HANDSHAKE_TAG order=(HANDSHAKE_TAG)this->getRemoteOrder();

    // connect to my remote kernels
    for (size_t i=0; i< notYetContacted.size(); i++){
        // string kid= myRemoteKerelsIDs.at(i);
        shared_ptr<KernelConnectionInfo> & remoteKernel=notYetContacted.at(i);

        //connect
        MPI_Comm inter_comm;
        char remotePortName[MPI_MAX_PORT_NAME];
        strcpy( remotePortName, remoteKernel->getMyMPIUrl().c_str() );

        MPI_Comm_set_errhandler(myMpiManager->getGlobalCommunicator(), MPI_ERRORS_RETURN);
        int res_connect=MPI_Comm_connect( remotePortName, MPI_INFO_NULL, 0, this->getMpiManager()->getGlobalCommunicator(),  &inter_comm );
        MPI_Comm_set_errhandler(this->myMpiManager->getGlobalCommunicator(), MPI_ERRORS_ARE_FATAL);
        if(res_connect == MPI_SUCCESS){

            shared_ptr<MpiManager> interMpiManager = std::make_shared<MpiManager>();
            interMpiManager->initOnlyCommunicator(inter_comm);
            //send my kernelName to remote kernel
            string mykernelName=this->myLocalInfo->getName();
            getInterCommUtils().broadCastStringToRemote(mykernelName, root, interMpiManager);
            // assign the new created intercomm
            remoteKernel->moveInterMpiManager(interMpiManager);
            if (this->getMpiManager()->isMainProcessor()){
                cout<< myLocalInfo->getName() << ": handshake to: "<< remoteKernel->getName()<<endl;
            }
        }else{
            throw std::runtime_error("Fatal Error on MPI_Comm_connect(...) in function connectToMyRemoteKernels()");
        }

    }
}

void SocketKernel::checkoutManagerPortUrl(){
    if(managerPortUrl.empty() && this->interMpiManager){
        int root= 0;
        //broadCastFromRemote manager url
        int length ;
        this->interMpiManager->bCast(&length, 1, root);
        char* buffer = new char[length+1];
        this->interMpiManager->bCast(buffer, length+1, root);
        managerPortUrl = buffer;
        delete [] buffer;
    }

    assert(! managerPortUrl.empty());
}

void SocketKernel::connectWithMyRelayer(string & relayerUrl, string &forwarderUrl){

    if(this->relayerCommunicator){
        //receive the frontend TCP URL of the relayer
        vector<char> buffer;
        int root=0;
        this->relayerCommunicator->bCastFrom( buffer, root);
        relayerUrl = buffer.data();
        vector<char>().swap(buffer);
        this->relayerCommunicator->bCastFrom( buffer, root);
        forwarderUrl = buffer.data();
        vector<char>().swap(buffer);
    }

}


shared_ptr<KernelConnectionInfo> SocketKernel::acceptConnection(){
    MPI_Comm client_interComm;
    MPI_Comm_set_errhandler(myMpiManager->getGlobalCommunicator(), MPI_ERRORS_RETURN);
    int res_connect=MPI_Comm_accept(port_name, MPI_INFO_NULL, this->getMpiManager()->bossId(), this->getMpiManager()->getGlobalCommunicator(),  &client_interComm);
    MPI_Comm_set_errhandler(this->myMpiManager->getGlobalCommunicator(), MPI_ERRORS_ARE_FATAL);
    if( res_connect  != MPI_SUCCESS )  return 0;
    shared_ptr<MpiManager> interMpiManager= std::make_shared<MpiManager>();
    interMpiManager->initOnlyCommunicator(client_interComm);
    shared_ptr<KernelConnectionInfo>  clientinfo= std::make_shared <KernelMPIConnectionInfo>(interMpiManager);
    return clientinfo;

}

void SocketKernel::disconnect(){

    closeMPISocket();
}

void SocketKernel::setMyRemoteKernels(vector<string> myRemoteKerelsIDs){
    this->myRemoteKerelsIDs.swap(myRemoteKerelsIDs);
}


HANDSHAKE_TAG SocketKernel::getRemoteOrder(){
    int root= 0;
    int order;
    //myLocalInfo->getInterMpiManager()->bCast<int>(&order,1, root);
    myLocalInfo->broadCastIntFromRemote(order,root);
    //HANDSHAKE_TAG remoteOrder=static_cast< typename std::underlying_type<HANDSHAKE_TAG>::type >(order);
    HANDSHAKE_TAG remoteOrder=static_cast<HANDSHAKE_TAG>(order);
    return remoteOrder;
}

shared_ptr<MpiManager> &SocketKernel::getMpiManager() {
    return this->myMpiManager;
}

shared_ptr<RelayerCommunicator> &SocketKernel::getRelayerCommunicator(){
    return this->relayerCommunicator;
}


void SocketKernel::closeMPISocket(){
    if(strlen(port_name) > 0){
        MPI_Close_port(port_name);
    }
}

int SocketKernel::openMPISocket(){

    int res=MPI_ERR_UNKNOWN;
    string endpoint = Networking::getInstance().openfrontendMPI(myMpiManager->getGlobalCommunicator());
    if(!endpoint.empty()){
        strncpy (port_name, endpoint.c_str(), endpoint.length()+1);
        res=MPI_SUCCESS;
    }
    return res;
}

shared_ptr<KernelConnectionInfo> SocketKernel::connectToManager(){


    std::vector<std::string> elems;
    char delim=';';
    std::stringstream ss(managerPortUrl);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    bool contains_tcp_prefix=false; //true for MPI
    for (auto e: elems){
        if (e.find("tcp://") == 0){//yes it starts starts with
            contains_tcp_prefix=true;
            break;
        }
    }

    shared_ptr<KernelConnectionInfo>  info;
    if(contains_tcp_prefix){ //========= MPI connect to manager ========
        char remotePortName[MPI_MAX_PORT_NAME];
        strcpy( remotePortName, managerPortUrl.c_str() );/* assume server's name is cmd-line arg */
        /* find port associated with a service name , collective for all clients */
        //MPI_Lookup_name("test server", MPI_INFO_NULL, port_name);
        assert(this->getMpiManager().get());
        MPI_Comm serverManager;// MANAGER
        //if ( debug && isBossId) cout<<"server Socket "<<managerPortUrl<<endl;
        MPI_Comm_set_errhandler(myMpiManager->getGlobalCommunicator(), MPI_ERRORS_RETURN);
        int res_connect=MPI_Comm_connect( remotePortName, MPI_INFO_NULL, this->getMpiManager()->bossId(), this->getMpiManager()->getGlobalCommunicator(),  &serverManager );
        MPI_Comm_set_errhandler(this->myMpiManager->getGlobalCommunicator(), MPI_ERRORS_ARE_FATAL);

        if(res_connect != MPI_SUCCESS) return 0;

        shared_ptr<MpiManager> inter = std::make_shared<MpiManager>();
        inter->initOnlyCommunicator(serverManager);
        info= std::make_shared<KernelMPIConnectionInfo>(inter);
    }else{ //========= TCP connect to manager ========
        int sock=-1;
        if(this->myMpiManager->isMainProcessor()){

            // managerPortUrl snytax: ip1;ip2;...;port;core/ip1;ip2;...;port;core
            std::vector<std::string> elems;
            std::size_t found=this->managerPortUrl.find("/");

            if(found != std::string::npos){
                std::stringstream ss(managerPortUrl);
                std::string item;
                while (std::getline(ss, item, '/')) {
                    elems.push_back(item);
                }
                managerPortUrl= elems.at(0);
                string fwd= elems.at(1);
                this->managerFWD.parse(fwd);
            }
            Endpoint mgr (managerPortUrl);
            sock= Networking::getInstance().connectToFrontEndRelayer(mgr);
            if(sock<0){
                sock=Networking::getInstance().connectToFrontEndForwarder(this->managerFWD, myKernelName, mgr.rawEndpoints, this->managerFWD.rawEndpoints );
                if(sock < 0)
                    sock=Networking::getInstance().connectToFrontEndForwarder(this->forwarderAddress, myKernelName, mgr.rawEndpoints, this->managerFWD.rawEndpoints);
            }
        }
        info= std::make_shared<KernelTCPConnectionInfo>(sock, this->myMpiManager);
    }
    return info;
}

//********************************** MpiSocketKernel ************************************************

MpiSocketKernel::MpiSocketKernel(shared_ptr<MpiManager> globalMpiManager, shared_ptr<MpiManager> localMpiManager,
                                 string myKernelName, string managerPortUrl,
                                 shared_ptr<RelayerCommunicator> &relayerCommunicator, string forwarderAddress):
    SocketKernel(globalMpiManager, localMpiManager, 0, relayerCommunicator, myKernelName, forwarderAddress, managerPortUrl){

}

MpiSocketKernel::MpiSocketKernel(shared_ptr<MpiManager> globalMpiManager, shared_ptr<MpiManager>localMpiManager,
                                 string myKernelName, shared_ptr<MpiManager> interMpiManager,
                                 shared_ptr<RelayerCommunicator> &relayerCommunicator, string forwarderAddress):
    SocketKernel(globalMpiManager, localMpiManager, interMpiManager, relayerCommunicator, myKernelName, forwarderAddress){
}

MpiSocketKernel::~MpiSocketKernel(){

}





//********************************** SocketKernel ************************************************

TCPSocketKernel::TCPSocketKernel(shared_ptr<MpiManager> globalMpiManager,shared_ptr<MpiManager> localMpiManager, shared_ptr<RelayerCommunicator> &relayerCommunicator, string myKernelName, string managerPortUrl):
    SocketKernel(globalMpiManager, localMpiManager, 0, relayerCommunicator, myKernelName, managerPortUrl),  sockfd(-1){
}

TCPSocketKernel::~TCPSocketKernel(){
}







//***************************************************** InterCommUtils ************************************
InterCommUtils::InterCommUtils(){}

InterCommUtils::~InterCommUtils(){}

void InterCommUtils::broadCastStringToRemote(string message, int root, shared_ptr<MpiManager> &m){

    int length = (int) message.size();
    m->bCast(&length, 1, root);
    char* buffer = new char[length+1];
    //if (m->getRank()==root) {
    std::copy(message.c_str(), message.c_str()+length+1, buffer);
    //}
    m->bCast(buffer, length+1, root);
    delete [] buffer;
}

void InterCommUtils::broadCastStringFromRemote(string &message, int root, shared_ptr<MpiManager> & m){

    int length ;
    m->bCast(&length, 1, root);
    char* buffer = new char[length+1];
    m->bCast(buffer, length+1, root);
    message = buffer;
    delete [] buffer;
}

//***************************************************** KernelConnectionInfo_A ************************************
KernelConnectionInfo::KernelConnectionInfo():color(0){}
KernelConnectionInfo::~KernelConnectionInfo(){}

void  KernelConnectionInfo::setName(string id){this->name=id;}

void  KernelConnectionInfo::setMyMPIUrl(string uri){this->mpiUrl=uri;}

void KernelConnectionInfo::setMyRelayerUrl(string uri){
    this->myRelayerUrl=uri;
}

void KernelConnectionInfo::setMyForwarderUrl(string fwdUrl){
    this->myForwarderUrl=fwdUrl;
}

void  KernelConnectionInfo::setcolor(int color){this->color=color;}


string  KernelConnectionInfo::getName() const{
    return this->name;
}

string  KernelConnectionInfo::getMyMPIUrl() const{
    return this->mpiUrl;
}

string KernelConnectionInfo::getMyRelayerUrl() const{
    return this->myRelayerUrl;
}

string KernelConnectionInfo::getMyForwarderUrl() const{
    return this->myForwarderUrl;
}

int  KernelConnectionInfo::getcolor() const{
    return this->color;
}

void KernelConnectionInfo::receiveInfo(){
    this->receive(this->name);
    this->receive(this->mpiUrl);
    this->receive(this->myRelayerUrl);
    this->receive(this->myForwarderUrl);
    this->receive(this->color);
    //interMpiManager->receive<int>(&this->color, 1 , senderRank,  conduit_MPI_Tag);

}

void KernelConnectionInfo::sendInfo(){
    this->send(this->name);
    this->send(this->mpiUrl);
    this->send(this->myRelayerUrl);
    this->send(this->myForwarderUrl);
    this->send(this->color);
    //interMpiManager->send<int>(&this->color, 1, receiverRank, conduit_MPI_Tag);
}


void  KernelConnectionInfo::broadCastMyInfoFromRemote(int root){
    this->broadCastStringFromRemote(this->name, root);
    this->broadCastStringFromRemote(this->mpiUrl, root);
    this->broadCastStringFromRemote(this->myRelayerUrl, root);
    this->broadCastStringFromRemote(this->myForwarderUrl, root);
    this->broadCastIntFromRemote(this->color, root);
}

void  KernelConnectionInfo::broadCastMyInfoToRemote(int root){
    this->broadCastStringToRemote(this->name, root);
    this->broadCastStringToRemote(this->mpiUrl, root);
    this->broadCastStringToRemote(this->myRelayerUrl, root);
    this->broadCastStringToRemote(this->myForwarderUrl, root);
    this->broadCastIntToRemote(this->color, root);
}

void  KernelConnectionInfo::broadCastMapFromRemote(map<string, shared_ptr<KernelConnectionInfo> > &remoteKernelsMap, int root){
    int length;
    //interMpiManager->bCast<int>(&length, 1, root);
    broadCastIntFromRemote(length, root);
    for (int i=0; i< length; i++){
        shared_ptr<KernelConnectionInfo>  info = this->createKernelConnectionInfo();
        this->broadCastInfoFromRemote(info, root);
        remoteKernelsMap.insert ( std::pair<string,  shared_ptr<KernelConnectionInfo>  >(info->getName(), info));
    }
}

void  KernelConnectionInfo::broadCastMapToRemote(map<string, shared_ptr<KernelConnectionInfo> > &remoteKernelsMap, int root){
    //int length=remoteKernelsMap.size()-1; // not broascast to itself
    int length=remoteKernelsMap.size();
    //interMpiManager->bCast<int>(&length, 1, root);
    this->broadCastIntToRemote(length, root);

    for( it_map_comm it = remoteKernelsMap.begin(); it!= remoteKernelsMap.end();it++){
        //string kernelID=it->first;
        shared_ptr<KernelConnectionInfo> & ptr= it->second;
        //if(kernelID != this->id){
        // remote must know its color
        this->broadCastInfoToRemote(ptr, root);
        // }
    }
}

void  KernelConnectionInfo::broadCastInfoFromRemote(shared_ptr<KernelConnectionInfo> &info, int root){
    this->broadCastStringFromRemote(info->name, root);
    this->broadCastStringFromRemote(info->mpiUrl, root);
    this->broadCastStringFromRemote(info->myRelayerUrl, root);
    this->broadCastStringFromRemote(info->myForwarderUrl, root);
    this->broadCastIntFromRemote(info->color, root);
}

void  KernelConnectionInfo::broadCastInfoToRemote(shared_ptr<KernelConnectionInfo> &info, int root){
    this->broadCastStringToRemote(info->name, root);
    this->broadCastStringToRemote(info->mpiUrl, root);
    this->broadCastStringToRemote(info->myRelayerUrl, root);
    this->broadCastStringToRemote(info->myForwarderUrl, root);
    this->broadCastIntToRemote(info->color, root);
}



shared_ptr<KernelConnectionInfo> KernelConnectionInfo::createKernelConnectionInfo()
{
    return std::make_shared<KernelMPIConnectionInfo>();
}

//***************************************************** KernelConnectionInfo ************************************


KernelMPIConnectionInfo::KernelMPIConnectionInfo():
    KernelConnectionInfo(), conduit_MPI_Tag(5), receiverRank(0), senderRank(0){}


KernelMPIConnectionInfo::KernelMPIConnectionInfo(shared_ptr<MpiManager> &interMpiManager):
    KernelConnectionInfo(), interMpiManager(interMpiManager),conduit_MPI_Tag(5), receiverRank(0), senderRank(0){}

KernelMPIConnectionInfo::~KernelMPIConnectionInfo(){
}


void KernelMPIConnectionInfo::moveInterMpiManager(shared_ptr<MpiManager> &interMpiManager){
    this->interMpiManager=interMpiManager;
}

shared_ptr<MpiManager> KernelMPIConnectionInfo::getInterMpiManager(){
    return this->interMpiManager;
}



bool KernelMPIConnectionInfo::send(string & message){

    int length = (int) message.size();
    this->interMpiManager->send<int>(&length, 1, receiverRank, conduit_MPI_Tag);
    // bCast(&length, 1, root);
    char* buffer = new char[length+1];
    std::copy(message.c_str(), message.c_str()+length+1, buffer);
    //bCast(buffer, length+1, root);
    this->interMpiManager->send<char>(buffer, length+1, receiverRank, conduit_MPI_Tag);
    delete [] buffer;
    return true;
}

bool KernelMPIConnectionInfo::send(int & value){
    this->interMpiManager->send<int>(&value, 1, receiverRank, conduit_MPI_Tag);
    return true;
}
bool KernelMPIConnectionInfo::receive(int & value){
    interMpiManager->receive<int>(&value, 1 , senderRank,  conduit_MPI_Tag);
    return true;
}


bool KernelMPIConnectionInfo::receive( string & message){
    int length;
    this->interMpiManager->receive<int>(&length, 1 , senderRank,  conduit_MPI_Tag);
    char* buffer = new char[length+1];
    this->interMpiManager->receive<char>(buffer, length+1 , senderRank,  conduit_MPI_Tag);

    message = buffer;
    delete [] buffer;
    return true;
}

void KernelMPIConnectionInfo::broadCastStringToRemote(string &message, int root){

    int length = (int) message.size();
    interMpiManager->bCast<int>(&length, 1, root);
    char* buffer = new char[length+1];
    std::copy(message.c_str(), message.c_str()+length+1, buffer);
    interMpiManager->bCast<char>(buffer, length+1, root);
    delete [] buffer;

}

void KernelMPIConnectionInfo::broadCastStringFromRemote(string &message, int root){

    int length ;
    interMpiManager->bCast<int>(&length, 1, root);
    char* buffer = new char[length+1];
    interMpiManager->bCast<char>(buffer, length+1, root);
    message = buffer;
    delete [] buffer;

}

void KernelMPIConnectionInfo::broadCastIntToRemote(int &length, int root)
{
    interMpiManager->bCast(&length, 1, root);
}

void KernelMPIConnectionInfo::broadCastIntFromRemote(int &length, int root)
{
    interMpiManager->bCast(&length, 1, root);
}


//***************************************************** KernelTCPConnectionInfo ************************************

KernelTCPConnectionInfo::KernelTCPConnectionInfo():
    KernelConnectionInfo(){}


KernelTCPConnectionInfo::KernelTCPConnectionInfo(int  sockfd, shared_ptr<MpiManager> &localMpiManager):
    KernelConnectionInfo(), sockfd(sockfd), localMpiManager(localMpiManager){}

KernelTCPConnectionInfo::~KernelTCPConnectionInfo(){
    if(sockfd>0)
        close(sockfd);
}


void KernelTCPConnectionInfo::moveInterMpiManager(shared_ptr<MpiManager> &interMpiManager){
    //this->interMpiManager=interMpiManager;
    throw -1;
}

shared_ptr<MpiManager> KernelTCPConnectionInfo::getInterMpiManager(){
    //return this->interMpiManager;
    throw -2;
    shared_ptr<MpiManager> p;
    return p;
}



bool KernelTCPConnectionInfo::send(string & message){

    int n=0, n1=0;
    if(localMpiManager->isMainProcessor()){
        int length = (int) message.size();
        n = Networking::getInstance().sock_write(sockfd, &length, sizeof(length));
        //this->interMpiManager->send<int>(&length, 1, receiverRank, conduit_MPI_Tag);
        char* buffer = new char[length+1];
        std::copy(message.c_str(), message.c_str()+length+1, buffer);
        n1 = Networking::getInstance().sock_write(sockfd, buffer, length+1);
        //this->interMpiManager->send<char>(buffer, length+1, receiverRank, conduit_MPI_Tag);
        delete [] buffer;
    }
    return (n && n1)? true: false;

}

bool KernelTCPConnectionInfo::send(int & value){
    //this->interMpiManager->send<int>(&value, 1, receiverRank, conduit_MPI_Tag);
    int n=0;
    if(localMpiManager->isMainProcessor()){
        n = Networking::getInstance().sock_write(sockfd, &value, sizeof(value));
    }
    return (n)? true:false;
}
bool KernelTCPConnectionInfo::receive(int & value){
    // interMpiManager->receive<int>(&value, 1 , senderRank,  conduit_MPI_Tag);
    int n=0;
    if(localMpiManager->isMainProcessor()){

        n =  Networking::getInstance().sock_read(sockfd,&value,sizeof(value));//read(sockfd,&value,sizeof(value));
    }
    return (n)? true:false;
}


bool KernelTCPConnectionInfo::receive( string & message){

    int n=0, n1=0;
    if(localMpiManager->isMainProcessor()){
        int length;
        // this->interMpiManager->receive<int>(&length, 1 , senderRank,  conduit_MPI_Tag);
        n = Networking::getInstance().sock_read(sockfd, &length, sizeof(length));
        char* buffer = new char[length+1];
        //this->interMpiManager->receive<char>(buffer, length+1 , senderRank,  conduit_MPI_Tag);
        n1 = Networking::getInstance().sock_read(sockfd, buffer, length+1);
        message = buffer;
        delete [] buffer;
    }
    return (n1 && n)? true:false;
}

void KernelTCPConnectionInfo::broadCastStringToRemote(string &message, int root){
    if(localMpiManager->isMainProcessor()){
        this->send(message);
    }

}

void KernelTCPConnectionInfo::broadCastStringFromRemote(string &message, int root){
    int length ;
    if(localMpiManager->isMainProcessor()){
        this->receive(message);
        length = (int) message.size();
    }
    localMpiManager->bCast<int>(&length, 1);
    char* buffer = new char[length+1];
    if (localMpiManager->isMainProcessor()) {
        std::copy(message.c_str(), message.c_str()+length+1, buffer);
    }
    localMpiManager->bCast<char>(buffer, length+1);
    if (! localMpiManager->isMainProcessor()) {
        message = buffer;
    }
    delete [] buffer;
}

void KernelTCPConnectionInfo::broadCastIntToRemote(int &length, int root)
{
    if(localMpiManager->isMainProcessor()){
        this->send(length);
    }
    //nterMpiManager->bCast(&length, 1, root);
}

void KernelTCPConnectionInfo::broadCastIntFromRemote(int &length, int root)
{
    //interMpiManager->bCast(&length, 1, root);
    if(localMpiManager->isMainProcessor()){
        this->receive(length);
    }
    localMpiManager->bCast<int>(&length, 1);
}
