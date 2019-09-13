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

#ifndef MMSF_HPP
#define MMSF_HPP

#include "musclehpc/conduits/mmsf.h"

using namespace std;
namespace unige_pasc{

/************************************************ ColorInfo **************************************************/

ColorInfo::ColorInfo(Kernel *k, int color, int rankInf, int rankSup, bool isManager):
    kernel(k), color(color), rankInf(rankInf),  rankSup(rankSup), manager(isManager), obsolet(false), relayer(false){
}

Kernel *ColorInfo::getKernel() const {
    return this->kernel;
}

int ColorInfo::getColor() const{
    return this->color;
}

bool ColorInfo::isManager() const{
    return this->manager;
}

int ColorInfo::getGlobalLeaderRank() const{
    return this->rankInf;
}

bool ColorInfo::isObsolet(){
    if(!this->manager && !this->relayer){
        return (kernel==0);
    }
    return false;
}

bool ColorInfo::doesBelongToMe(int processRank){
    if( (processRank >= rankInf) && (processRank <= rankSup)){
        return true;
    }
    return false;
}

bool ColorInfo::isRelayer(){
    return this->relayer;
}

void ColorInfo::setRelayer(bool value){
    this->relayer=value;
}

/************************************ MMSF ************************************************/

MMSF::MMSF( ) {

}

MMSF::~MMSF(){
    for ( it_map_condtuis it = conduits.begin(); it != conduits.end(); it++ )
    {
        ConduitCommon* ptr=it->second ;
        if(ptr)
            delete ptr;
    }

    conduits.clear();
    addedKernels.clear();
}


void MMSF::addConduit( string id, ConduitCommon* conduit) {
    it_map_condtuis it = conduits.find(id);
    if( it!= conduits.end()) {
        cout<<" A conduit was already declared with ID: " << id <<endl;
    }
    conduits.insert ( std::pair<string,ConduitCommon* >(id, conduit));
}



template<typename E>
Binder<E> * MMSF::connect( string fromStr ) {
    Binder<E> * b= new Binder<E>(this);
    return (b->From(fromStr));
}


template<typename E>
void MMSF::connect( string fromStr, string toStr, Conduit<E> * conduit, string conduitID ) {
    conduit->initialize(conduitID );
    Kernel *fromKernel = assignConduit( fromStr, conduit, true);
    conduit->registerSender( fromKernel );
    // Attaching the to kernel
    Kernel *toKernel = assignConduit( toStr, conduit, false);
    conduit->registerReceiver( toKernel );
    this->addConduit(conduitID, conduit);
}

template<typename E>
Kernel * MMSF::assignConduit(string portString, Conduit<E> * conduit, bool isIn){
    Port * port = new Port( portString );
    Kernel * kernel = getKernel(port->getKernelName());
    assert(kernel!=0);
    if(isIn)
        kernel->registerFout(port->getPortName(),  conduit); // F_out is used to send into it
    else
        kernel->registerFin(port->getPortName(), conduit); // F_in is used to receive from
    return kernel;
}


Kernel * MMSF::getKernel( string id ) {
    it_map_kernel it = addedKernels.find(id);
    if( it!= addedKernels.end()) {
        return it->second;
    }else{
        return 0;
    }
}

map<string, Kernel*> const & MMSF::getListKernels()const {
    return this->addedKernels;
}

void MMSF::compute(){

    int cpt=0;
    for (it_map_kernel it = this->addedKernels.begin(); it!= this->addedKernels.end(); ++it){
        Kernel * k = it->second;
        k->setColorWhereToRun(cpt);
        cpt++;
        k->mainLoop();
    }
}

// protected
Kernel * MMSF::getkernelOfColor(int color){
    it_map_kernel it = addedKernels.begin();
    int cpt=0;
    while (cpt < color && it!=addedKernels.end()){
        it++;
        cpt++;
    }
    if(it==addedKernels.end()){
        return 0;
    }else{
        Kernel * k = it->second;
        return k;
    }
}


/************************************ MMSF_MPI ************************************************/

MMSF_MPI::MMSF_MPI(int argc, char** argv, MPI_Comm globalCommunicator): MMSF(){

    shared_ptr<MpiManager> interMpiManager ( make_shared<MpiManager>() );
    this->globalMpiManager= std::move(interMpiManager);
    this->globalMpiManager->init(&argc, &argv, globalCommunicator); //<-- MPI Init here

    //shared_ptr<OptionParser> parser (make_shared<OptionParser>(argc, argv) );
    //this->optionsParser= std::move(parser);

    this->optionsParser = OptionParser::getInstance(argc, argv);
    // if '-h' option display help and quit program
    if( ! this->optionsParser->getDisplayHelp().empty()){
        if (this->globalMpiManager->isMainProcessor()){
            cout<<this->optionsParser->getDisplayHelp();
        }
        this->globalMpiManager.reset();
        exit(EXIT_SUCCESS); // <- may be we sould throw an exception instead of using exit ??
    }

}

MMSF_MPI::MMSF_MPI(int argc, char** argv): MMSF_MPI(argc, argv, MPI_COMM_WORLD){}

MMSF_MPI::~MMSF_MPI(){
    mappingCoresColor.clear();
    for (size_t i=0; i<this->colorsInfoMap.size(); i++){
        ColorInfo * ptr= this->colorsInfoMap.at(i);
        if(ptr){
            delete ptr;
            ptr=0;
        }
    }
    this->colorsInfoMap.clear();
}

void MMSF_MPI::addKernel( Kernel * k, string kernelId, int cores) {


    it_map_kernel it = addedKernels.find(kernelId);
    if( it!= addedKernels.end()) {
        cout<<" A kernel was already declared with ID: " << kernelId <<endl;
    }

    vector<string> kernelsToRun = this->optionsParser->getkernelsToRun();

    bool ishere= (this->optionsParser->allowAllKernelToRun())?true:(std::find(kernelsToRun.begin(), kernelsToRun.end(), kernelId) != kernelsToRun.end());
    if(!ishere){
        if(this->globalMpiManager->isMainProcessor()) cout<<"\n=== [Warning]: -> The kernel "<<kernelId<<" will not run here ==="<<endl;
        k->setCoresOverWhichTorun(0); // will not run it
    }else{
        k->setCoresOverWhichTorun(cores);
    }
    k->initialize(kernelId);
    addedKernels.insert ( std::pair<string,Kernel * >(kernelId, k));
    k->setKernelArgc(this->optionsParser->getArgc());
    k->setKernelArgv(this->optionsParser->getArgv());
}

void MMSF_MPI::addKernel( Kernel * k, string kernelId) {
    map<string, int> const &mapcores= this->optionsParser->getRequestedCores();
    it_map_cores itmc = mapcores.find(kernelId);
    int reqCores=this->globalMpiManager->getSize();
    if( itmc!= mapcores.end()) {
        reqCores=itmc->second;
    }
    this->addKernel(k, kernelId, reqCores);
}



void MMSF_MPI::addKernel(string kernelId) {
    EmptyRemoteMpiKernel * k = new EmptyRemoteMpiKernel();

    it_map_kernel it = addedKernels.find(kernelId);
    if( it!= addedKernels.end()) {
        cout<<" A kernel was already declared with ID: " << kernelId <<endl;
    }
    k->setCoresOverWhichTorun(0);
    k->initialize(kernelId);
    addedKernels.insert ( std::pair<string,Kernel * >(kernelId, k));
}

bool MMSF_MPI::loadCxACoupling (){

    vector<CouplingInfo*> vect;
    vector<CommandLineInfo*>  vectCmdInfo;
    bool isok= this->optionsParser->getCouplingFromCxaFile(vect, vectCmdInfo);

    if(!isok){
        return false;
    }else{
        //for(size_t i=0; i<vect.size();i++){
        //    CouplingInfo* info = vect.at(i);
        for(std::vector<CouplingInfo*>::iterator it = vect.begin(); it != vect.end(); ++it) {
            CouplingInfo* info = *it;
            cxaAllKernels.push_back(info->getSenderName());
            cxaAllKernels.push_back(info->getReceiverName());
            string k1=info->getSenderName()+"."+info->getSenderPort();
            string k2=info->getReceiverName()+"."+info->getReceiverPort();
            string cd=info->getConduitName();
            string dataType=info->getDataType();

            if(dataType=="char"){
                this->connect<char>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="int"){
                this->connect<int>(k1)->To(k2)->with_MTM_MPI(cd);
            }/*else if(dataType=="bool"){
                this->connect<bool>(k1)->To(k2)->with_MTM_MPI(cd);
            }*/else if(dataType=="longlong"){
                this->connect<long long>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="unsignedlonglong"){
                this->connect<unsigned long long>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="long"){
                this->connect<long>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="float"){
                this->connect<float>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="double"){
                this->connect<double>(k1)->To(k2)->with_MTM_MPI(cd);
            }else if(dataType=="longdouble"){
                this->connect<long double>(k1)->To(k2)->with_MTM_MPI(cd);
            }
            delete info;
        }

        // remove duplicated kernels name in cxaAllKernels
        std::sort( cxaAllKernels.begin(), cxaAllKernels.end() );
        cxaAllKernels.erase( std::unique( cxaAllKernels.begin(), cxaAllKernels.end() ), cxaAllKernels.end() );

        vect.clear();

        for(vector<CommandLineInfo*>::iterator it= vectCmdInfo.begin() ; it != vectCmdInfo.end();it++){
            CommandLineInfo *cinfo = *it;

            it_map_kernel mapit = addedKernels.find(cinfo->getKernelName());
            if( mapit!= addedKernels.end()) {
                Kernel * k= mapit->second;
                k->setKernelArgc(cinfo->getargc());
                k->setKernelArgv(cinfo->getargv());

                //cout<<cinfo->getKernelName()<<" ++++++++> "<<k->getKernelArgv()[0] <<" "<<k->getKernelArgv()[1]<<endl;
            }else{
                //cout<<cinfo->getKernelName()<<" ---------> "<<cinfo->getargv()[0] <<" "<<cinfo->getargv()[1]<<endl;
            }
        }
    }
    return true;
}

int MMSF_MPI::getCoresNumber(){
    return this->globalMpiManager->getSize();
}


int MMSF_MPI::getGlobalRank(){
    return this->globalMpiManager->getRank();
}

void MMSF_MPI::compute(bool useGlobalManager){

    

    bool isManager=this->optionsParser->isManager();

    //bool isManager; // true: means Socket/connect ativated and manger is running here
    // false: manger is running here and no socket/connect mechnism (only for MMSF_MPI)

    assert(this->globalMpiManager.get()!=0);

    


    int  requestedKernelTotalCores;
    bool isok1= computeTotalRequestedCores(requestedKernelTotalCores, isManager);
    this->populateColorsInfoMap(isManager);
    this->perapreConduitsTag();
    bool isok2 = verifyPTPConduitsAndCores();
    if(!isok2){
        //fix this after
        isok2=true;
    }

    int grank=this->globalMpiManager->getRank();

    ColorInfo * infoColor= this->getColorInfoByRank(grank);
    bool isok3 =(infoColor)? true: false;

    int color=infoColor->getColor();

    MPI_Comm comm_work;
    MPI_Comm_split(globalMpiManager->getGlobalCommunicator(), color, grank, &comm_work);
    shared_ptr<MpiManager> localMpiManager(std::make_shared<MpiManager>());
    localMpiManager->init(comm_work, color);

    

    if(!(isok1 && isok2 &&isok3 && !infoColor->isObsolet())){
        if(infoColor->isObsolet())
            cout<<"-> Process of rank="<<grank<<" is not used"<<endl;
        //cout<<"error isok1="<<isok1<<", isok2="<<isok2<< ", isok3="<<isok3<<endl;
    }else{
        launch(localMpiManager,isManager);//default
    }

    // globalMpiManager->getRank();

}

//private:



bool MMSF_MPI::verifyPTPConduitsAndCores(){
    for (it_map_condtuis it = conduits.begin(); it!= conduits.end(); ++it){

        ConduitCommon * cd = it->second;
        if (! cd->isConduitSafe()){
            if(this->globalMpiManager->isMainProcessor()){
                cout<<"****************************************************************************************************\n"
                   <<"[Warning]: MPI_PTP_Conduit is used -> Kernels should have the same number of core."  <<"\n"
                  <<"****************************************************************************************************"<<endl;
            }
            return false;
        }
    }
    return true;
}

bool MMSF_MPI::computeTotalRequestedCores(int & requestedTotalCores, bool isManager){

    string info="";
    int addedcore=0;
    if(isManager){
        info="+ Manager ";
        addedcore=1;
    }

    if(! Networking::getInstance().canOpenMPIPort(this->globalMpiManager->getGlobalCommunicator())
            && !this->optionsParser->getkernelsToRun().empty()){
        addedcore+=2; // relay needs 2 cores
        info+="+ Relay";
    }

    requestedTotalCores=0;
    for (it_map_kernel it = addedKernels.begin(); it!= addedKernels.end(); ++it){
        Kernel * k = it->second;
        int cores= k->getCoresOverWhichTorun();
        requestedTotalCores+=cores;
    }

    int availablCores=this->getCoresNumber();

    if( (requestedTotalCores+addedcore) <= availablCores){
        if(this->globalMpiManager->isMainProcessor()){
            cout<<"\n[INFO]: Total available cores         = "<<getCoresNumber()<<"\n"
               <<"           Total requested kernels cores "<<info <<" = "<<requestedTotalCores+addedcore<<endl;
        }
    }else {
        if(this->globalMpiManager->isMainProcessor()){
            cout<<"\n[Info]:  total available cores                = "<<getCoresNumber()<<"\n"
               <<"         sum of total requested kernels cores "<<info <<" = "<<requestedTotalCores+addedcore<<"\n"
              <<"****************************************************************************************************\n"
             <<"[Error]: Sum of total requested cores for all kernels should not exceed the total available cores\n"
            <<"****************************************************************************************************"<<endl;
        }
        return false;
    }

    return true;
}

void MMSF_MPI::populateColorsInfoMap(bool isManager){

    int color=0;
    int rankInf, rankSup;
    rankInf=0;
    rankSup=0;



    // manager rank should have a ColorInfo with no kernel;
    if(isManager){
        ColorInfo * info = new ColorInfo(0, color, rankInf, rankSup, isManager);
        this->colorsInfoMap.push_back(info);
        color=1; //because color 0 is reservred for the Manager
        rankInf=1;//because rank 0 is reservred for the Manager
        rankSup=1;
    }

    if(color==0) color++; // color of kernels starts from 1

    // assign a color to each kernel
    for (it_map_kernel it = addedKernels.begin(); it!= addedKernels.end(); ++it){
        Kernel * k = it->second;
        if(k) mapColorWhereToRun.insert(pair<string,int>(k->getName(), color++));
    }

    // loop over kernels;

    for (it_map_kernel it = addedKernels.begin(); it!= addedKernels.end(); ++it){
        Kernel * k = it->second;
        int rqcores= k->getCoresOverWhichTorun();
        if(rqcores>0){//only kernel with correct cores number
            rankSup= rankInf+ k->getCoresOverWhichTorun()-1;
            k->setColorWhereToRun(mapColorWhereToRun[k->getName()]);
            ColorInfo * info = new ColorInfo(k, k->getColorWhereItRuns(), rankInf, rankSup);
            k->setGlobalLeaderRank(rankInf);
            this->colorsInfoMap.push_back(info);
            rankInf= rankSup+1;
        }
    }

    if(! Networking::getInstance().canOpenMPIPort(this->globalMpiManager->getGlobalCommunicator())){
        //----- the relayer ----------
        rankSup= rankInf+1; // only 2 cores
        ColorInfo * info = new ColorInfo(0, color, rankInf, rankSup);
        info->setRelayer(true);
        this->colorsInfoMap.push_back(info);
        //---------------------------
        rankInf= rankSup+1;
        color++;
    }

    // obselet rank should have a ColorInfo with no kernel;
    if(rankInf<= (this->globalMpiManager->getSize() -1)){
        rankSup= this->globalMpiManager->getSize() - 1;
        ColorInfo * info = new ColorInfo(0, color, rankInf, rankSup);
        this->colorsInfoMap.push_back(info);
    }
}

ColorInfo * MMSF_MPI::getColorInfoByRank( int rank){
    for (size_t i=0; i< this->colorsInfoMap.size(); i++){
        ColorInfo * info = this->colorsInfoMap.at(i);
        if (info->doesBelongToMe(rank)){
            return info;
        }
    }
    return 0;
}

ColorInfo * MMSF_MPI::getRelayer(){
    for (size_t i=0; i< this->colorsInfoMap.size(); i++){
        ColorInfo * info = this->colorsInfoMap.at(i);
        if (info->isRelayer()){
            return info;
        }
    }
    return 0;
}


void MMSF_MPI::perapreConduitsTag(){
    vector<string> conduitsNames;
    for (it_map_condtuis it = conduits.begin(); it!= conduits.end(); ++it){
        ConduitCommon * cd = it->second;
        conduitsNames.push_back(cd->getName());
    }
    // remove duplication from vector
    std::sort( conduitsNames.begin(), conduitsNames.end() );
    conduitsNames.erase( std::unique( conduitsNames.begin(), conduitsNames.end() ), conduitsNames.end() );
    int tag=10;
    for (size_t i=0; i<conduitsNames.size(); i++){
        this->conduitsTagMap.insert ( std::pair<string, int >(conduitsNames.at(i), tag));
        ConduitCommon* tmp=conduits.at(conduitsNames.at(i));
        tmp->setTag(tag);
        tag++;
    }
}

int MMSF_MPI::get_InterComm_MPI_Unig_Tag(ConduitCommon * conduit){
    /*int uptag;
    int flag;
    MPI_Comm_get_attr(globalMpiManager->getGlobalCommunicator(), MPI_TAG_UB, &uptag, &flag);
    if(flag){
        std::hash<std::string> hash_fn;
        size_t int_hash_1 = hash_fn(conduit->getName());
        return int_hash_1 % uptag;
    }*/
    return conduitsTagMap.at(conduit->getName());
    //return this->current_MPI_tag++;
}


//shared_ptr<OptionParser> & MMSF_MPI::getOptionParser(){
//    return this->optionsParser;
//}


int MMSF_MPI::connectSubmodelsWithinsameBinary(shared_ptr<MpiManager> &localMpiManager){

    vector<string> activeSubmodels;
    for (size_t i=0; i<this->colorsInfoMap.size();i++){
        ColorInfo * info= colorsInfoMap.at(i);
        Kernel * k= info->getKernel();
        if(k){
            activeSubmodels.push_back(k->getName());
        }
    }


    int grank=this->globalMpiManager->getRank();
    ColorInfo * info= this->getColorInfoByRank(grank);
    Kernel * k= info->getKernel();
    if(k){

        RemoteMpiKernel* kmpi= dynamic_cast<RemoteMpiKernel* > (k);
        kmpi->setMpiManager(localMpiManager);
        // 1) call Kernel parent prepareConduits()
        kmpi->Kernel::prepareConduits();

        std::map<string, Kernel *> allMyRemoteKernelsMap=kmpi->getAllmyRemoteKernel();
        std::map<string, Kernel *> allSameKernelsMap;

        // filter kernels and consider only those running within the same binary
        for(  std::map<string, Kernel *>::iterator it = allMyRemoteKernelsMap.begin(); it!= allMyRemoteKernelsMap.end(); it++) {
            string remoteKernelName = it->first;
            if(std::find(activeSubmodels.begin(), activeSubmodels.end(), remoteKernelName) != activeSubmodels.end()) {
                /* remoteKernelName is in the same binary as k */
                allSameKernelsMap.insert ( std::pair<string, Kernel* >(remoteKernelName, it->second));// here it doest not allow duplicated key
            }
        }

        kmpi->createInterCommunicators(allSameKernelsMap);
    }
    return activeSubmodels.size();

}

void MMSF_MPI::launch(shared_ptr<MpiManager> &localMpiManager, bool isManager){

    int activeSubmodels=this->connectSubmodelsWithinsameBinary(localMpiManager);

    int grank=this->globalMpiManager->getRank();

    ColorInfo * infoColor= this->getColorInfoByRank(grank);
    bool allowAllKernelsToRun=this->optionsParser->allowAllKernelToRun();
            //|| (activeSubmodels == this->cxaAllKernels.size());

    shared_ptr<SocketServer> server;
    if(!allowAllKernelsToRun ){// here we use listen/connect mechanism

        if(infoColor && infoColor->isRelayer() && activeSubmodels>0){
            try{
                cout<<"-> Process of rank="<<grank<<" is a <<relayer>>"<<endl;
                map<string, shared_ptr<RelayerCommunicator> > rmap= this->createRelayerKernelsIntercommunicators(localMpiManager);
                MPIRelayer relayer (this->globalMpiManager, localMpiManager, rmap, mapColorWhereToRun, this->optionsParser->getForwarderFrontEnd());
                relayer.run();
                //cout<<"## Relayer ended."<<endl;
            }
            catch (std::exception &e){
                cout<<"<<<Error>>> -> handler: "<<e.what()<<endl;
            }
            return;
        }

        if(infoColor->isManager()){
            map<string, shared_ptr<MpiManager> > interMpiManagersMap= this->createManagerKernelsInterCommunicators(localMpiManager);
            Manager manager (this->addedKernels.size(), interMpiManagersMap, mapColorWhereToRun, localMpiManager, globalMpiManager, this->optionsParser->getForwarderFrontEnd(), this->optionsParser->useManagerTCP());
            manager.run();
           // cout<<"## manager ended."<<endl;
            return;
        }

        if(infoColor->getKernel()){

            ColorInfo * relayer=this->getRelayer();
            shared_ptr<RelayerCommunicator> relayerCommunicator;
            if(relayer){
                shared_ptr<MpiManager> interMpiRelayer= this->createInterCommunicator(infoColor, relayer, localMpiManager);
                relayerCommunicator=std::make_shared<MPI_Client_RelayerCommunicator>( localMpiManager, interMpiRelayer);
            }

            if(isManager){
                ColorInfo * mng= this->getColorInfoByRank(0);
                shared_ptr<MpiManager> interMpiManager  = this->createInterCommunicator(infoColor, mng, localMpiManager);
                server=std::make_shared<MpiSocketKernel>(globalMpiManager, localMpiManager, infoColor->getKernel()->getName(), std::move(interMpiManager), relayerCommunicator, this->optionsParser->getForwarderFrontEnd());
            } else {// manager not running here
                server=std::make_shared<MpiSocketKernel>(globalMpiManager, localMpiManager, infoColor->getKernel()->getName(), this->optionsParser->getManagerUrl(), relayerCommunicator, this->optionsParser->getForwarderFrontEnd());
            }
            RemoteMpiKernel* kmpi= dynamic_cast<RemoteMpiKernel* > (infoColor->getKernel());
            kmpi->setMpiSocketServer(server);
        }
    } // end  listen/connect mechanism

    if(infoColor->getKernel()){
        try{
            RemoteMpiKernel* kmpi= dynamic_cast<RemoteMpiKernel* > (infoColor->getKernel());
            kmpi->setGlobalMpiManager(globalMpiManager);
            if(allowAllKernelsToRun){
                kmpi->setMpiManager(localMpiManager);
                kmpi->MpiKernel::prepareConduits();// using split mechanism
            }else{
                kmpi->prepareConduits(); // using listen/connect mechanism
            }
            if(kmpi->getMpiManager()->isMainProcessor()) cout<< kmpi->getName()<< " runs on "<<kmpi->getMpiManager()->getSize()<< " cores." <<endl;
            kmpi->mainLoop();
            // == wait all perocesses== .
        }catch (std::exception &e){
            cout<<"<<<Kernel mainloop Error>>> -> "<<e.what()<<endl;
        }
        try{
            localMpiManager->barrier();
            // == inform relayer that I finished ==
            if(!allowAllKernelsToRun && !server->getMyLocalInfo()->getMyRelayerUrl().empty()){
                // send end Message to relyer

                Endpoint edp(server->getMyLocalInfo()->getMyRelayerUrl());
                int dest= edp.grank;
                if(localMpiManager->isMainProcessor()){
                    Message m;
                    m.endMessage=1;
                    Messaging::getInstance().send(globalMpiManager,m , dest);
                }
            }
        }catch (std::exception &e){
            cout<<"<<<Ending Ralay Error>>> -> "<<e.what()<<endl;
        }
    }
}


map<string, shared_ptr<MpiManager> >  MMSF_MPI::createManagerKernelsInterCommunicators( shared_ptr<MpiManager> &localMpiManager){

    int grank=this->globalMpiManager->getRank();
    ColorInfo * infoColor= this->getColorInfoByRank(grank);
    map<string, shared_ptr<MpiManager> > interMpiManagersMap;
    for (size_t i=0; i<this->colorsInfoMap.size();i++){
        ColorInfo * remoteInfo= colorsInfoMap.at(i);
        if(remoteInfo->getKernel()){
            int uniqInterCommTag=CantorPairingFunction(infoColor->getColor(), remoteInfo->getColor())+MPI_TAG_INIT_MANAGER_URL;
            MPI_Comm intercomm;
            //int remoteLeaderGlobalRank = ;
            MPI_Intercomm_create(localMpiManager->getGlobalCommunicator(),
                                 localMpiManager->bossId(),
                                 globalMpiManager->getGlobalCommunicator(),
                                 remoteInfo->getGlobalLeaderRank() ,
                                 uniqInterCommTag,
                                 &intercomm);

            shared_ptr<MpiManager>  interMpiManager= std::make_shared<MpiManager>();
            interMpiManager->initOnlyCommunicator(intercomm);
            interMpiManagersMap.insert( pair<string, shared_ptr<MpiManager> > (remoteInfo->getKernel()->getName(), interMpiManager));

        }
    } //end for
    return interMpiManagersMap;
}

map<string, shared_ptr<RelayerCommunicator> >  MMSF_MPI::createRelayerKernelsIntercommunicators(shared_ptr<MpiManager> &localMpiManager){
    int grank=this->globalMpiManager->getRank();
    ColorInfo * infoColor= this->getColorInfoByRank(grank);
    map<string, shared_ptr<RelayerCommunicator> > localKernelMap;

    for (size_t i=0; i<this->colorsInfoMap.size();i++){
        ColorInfo * remoteInfo= colorsInfoMap.at(i);
        if(remoteInfo->getKernel()){
            int uniqInterCommTag=CantorPairingFunction(infoColor->getColor(), remoteInfo->getColor())+MPI_TAG_INIT_MANAGER_URL;
            MPI_Comm intercomm;
            int remoteLeaderGlobalRank = remoteInfo->getGlobalLeaderRank();
            MPI_Intercomm_create(localMpiManager->getGlobalCommunicator(),
                                 localMpiManager->bossId(),
                                 globalMpiManager->getGlobalCommunicator(),
                                 remoteLeaderGlobalRank,
                                 uniqInterCommTag,
                                 &intercomm);

            shared_ptr<MpiManager>  interMpiManager=  std::make_shared<MpiManager>();
            interMpiManager->initOnlyCommunicator(intercomm);
            shared_ptr<RelayerCommunicator> rm= std::make_shared<MPI_Client_RelayerCommunicator>(localMpiManager, interMpiManager);
            localKernelMap.insert( pair<string, shared_ptr<RelayerCommunicator> > (remoteInfo->getKernel()->getName(), rm));
        }
    }
    return localKernelMap;
}


shared_ptr<MpiManager> MMSF_MPI::createInterCommunicator(ColorInfo* localInfo, ColorInfo* remoteInfo, shared_ptr<MpiManager> &localMpiManager){

    int uniqInterCommTag=CantorPairingFunction(localInfo->getColor(), remoteInfo->getColor())+ MPI_TAG_INIT_MANAGER_URL;

    MPI_Comm intercomm;
    MPI_Intercomm_create(localMpiManager->getGlobalCommunicator(), localMpiManager->bossId(),
                         globalMpiManager->getGlobalCommunicator(), remoteInfo->getGlobalLeaderRank(),
                         uniqInterCommTag, &intercomm);

    shared_ptr<MpiManager> interMpiManager ( make_shared<MpiManager>() );
    interMpiManager->initOnlyCommunicator(intercomm);
    return interMpiManager;
}






}//end namespace
#endif
