/**
* @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>
* @version 1.0
* @section LICENSE

* MAPPER communication module
* Copyright (C) 2015  University of Geneva, Switzerland
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mapper.h"
#include <unistd.h>

namespace unige_pasc{


/**************************************** SmartMapper ************************************************/

TransportInterpolation::TransportInterpolation()
    : isProcess(true), it_factor(1), cpt(0), isAllSubmodelsConverged(false), lastFimeTimeStamp(0.0)
{


    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.type=SMART_MAPPER;
    this->info.dt=0.0;
    this->info.dx=0.0;

}

TransportInterpolation::~TransportInterpolation(){}


void TransportInterpolation::synchronizeCommunication(){

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);

    //requires only on

    vector<char> toSend;
    this->serializeHeader(toSend, this->info);

    vector<char> toReceiveF1in, toReceiveF2in;
    this->receiveFromFin1(toReceiveF1in);
    this->receiveFromFin2(toReceiveF2in);
    INFO_SMART_MAPPER I_f1in, I_f2in;
    this->unserializeHeader(toReceiveF1in, I_f1in);
    this->unserializeHeader(toReceiveF2in, I_f2in);
    this->setUpCoarseAndFineConduits(I_f1in.dt, I_f2in.dt);

    this->sendToCoarse(toSend);//HERE: send toSend vector through MUSCLE
    this->sendToFine(toSend);//HERE: send toSend vector through MUSCLE
    double q=  (I_f1in.dt > I_f2in.dt)? (I_f1in.dt/I_f2in.dt):(I_f2in.dt/I_f1in.dt);
    this->it_factor=MapperUtil::roundTo_uint<pluint>(q);


    logger->getWriter()<<"Interpoler -> Detecting Coarse and Fine Submodels"<<"\n";
    if(it_factor==1){
        coarseId=I_f1in.sectionId;
        fineId  =I_f2in.sectionId;
        this->coarseDt=I_f1in.dt;
        this->fineDt=I_f2in.dt;
        logger->getWriter()<<"Interpoler -> Both of connected submodels have the same dt="<< I_f1in.dt<<"\n";
    }else{
        if (I_f1in.dt > I_f2in.dt){
            logger->getWriter()<<"-> Coarse submodel has dt= "<< I_f1in.dt<<"\n";
            logger->getWriter()<<"-> Fine   submodel has dt= "<< I_f2in.dt<<"\n";
            coarseId=I_f1in.sectionId;
            fineId  =I_f2in.sectionId;
            this->coarseDt=I_f1in.dt;
            this->fineDt=I_f2in.dt;

        }else{
            logger->getWriter()<<"-> Coarse submodel has dt= "<< I_f2in.dt<<"\n";
            logger->getWriter()<<"-> Fine   submodel has dt= "<< I_f1in.dt<<"\n";
            coarseId=I_f2in.sectionId;
            fineId  =I_f1in.sectionId;
            this->coarseDt=I_f2in.dt;
            this->fineDt=I_f1in.dt;
        }
    }
    logger->getWriter().flush();
}

bool TransportInterpolation::isConverged(){
    return (! this->isProcess);
}


void TransportInterpolation::F_init(){
    this->cpt=0;
    this->info.sectionId=this->getSubmodel_ID();
    string kname= this->getName()+".log";
    // plb_ofstream  pcout(kname.c_str());
    this->logger= std::move(std::make_shared<plb_ofstream>(this->getMPIRank(), kname.c_str()));
    this->synchronizeCommunication();

}
void TransportInterpolation::Oi(){

}
void TransportInterpolation::S(){

    try{
        vector<char> fromCoarse, fromFine;

        //bool debug=true;
        bool debug=false;
        this->lastFimeTimeStamp+=this->fineDt;
        treatCoarse=false;
        if ((this->cpt%this->it_factor) == 0){treatCoarse=true;}

        //=========== receive ============================
        if (debug && this->getMPIRank()==0) logger->getWriter()<<"Iteration ["<<cpt<<"]:\n\treceive from submodel["<<fineId<<"]";
        this->isFineConverged = this->receiveFineConvergence();
        if(treatCoarse){ this->isCoarseConverged = this->receiveCoarseConvergence();}

        this->receiveFromFine(fromFine);
        if(treatCoarse){
            if (debug && this->getMPIRank()==0) logger->getWriter()<<"+ submodel["<<coarseId<<"]";
            this->receiveFromCoarse(fromCoarse);
        }else{
            if (debug && this->getMPIRank()==0) logger->getWriter()<<"             ";
        }
        if (debug && this->getMPIRank()==0)  logger->getWriter()<<"\t: [OK]"<<"\n";
        //=========== check convergence ==================
        string statusMessage;
        isAllSubmodelsConverged= this->checkConvergence(statusMessage, fromFine, fromCoarse );
        if (debug && this->getMPIRank()==0) logger->getWriter()<<"\tcheck global convergence\t\t: "<<statusMessage<<"\n";
        //=========== interpolation ======================
        vector<char> ToCoarse, ToFine;

        if (debug && this->getMPIRank()==0) logger->getWriter()<<"\tinterpolate";
        this->interpolate(ToCoarse, ToFine);
        if (debug && this->getMPIRank()==0) logger->getWriter()<<"\t\t\t\t: [OK]"<<"\n";
        //==================send==========================
        if (debug && this->getMPIRank()==0) logger->getWriter()<<"\tsend to submodel["<<fineId<<"]";
        this->sendToFine(ToFine);
        if(treatCoarse || isAllSubmodelsConverged){
            if (debug && this->getMPIRank()==0) logger->getWriter()<<" + submodel["<<coarseId<<"]";
            this->sendToCoarse(ToCoarse);
        }else{
            if (debug && this->getMPIRank()==0)  logger->getWriter()<<"             ";
        }
        if (debug && this->getMPIRank()==0) logger->getWriter()<<"\t: [OK]"<<"\n";
        this->cpt++;
        //==================clear=========================
        vector<char>().swap(fromCoarse);
        vector<char>().swap(fromFine);
        //fromCoarse.clear();
        //fromFine.clear();
        logger->getWriter().flush();
    }catch (const std::exception &exc)
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr <<"#######> mapper: "<< exc.what()<<endl;
        exit (EXIT_FAILURE);
    }
}

void TransportInterpolation::U(){
    if (isAllSubmodelsConverged){
        logger->getWriter()<<"will terminte ..."<<"\n";
        this->isProcess=false; //stop
    }


}
void TransportInterpolation::Of(){
    this->end();
    this->logger->getWriter().close();
}

bool TransportInterpolation::checkConvergence(string & statusMessage, vector<char> & fromFine, vector<char> & fromCoarse){

    bool converged=false;

    vector<BoundaryParticle> fromCoarseParticles, fromFineParticles;
    // check if both of submodels converged
    if(treatCoarse) {
        bool tmp = isCoarseConverged;
        MapperUtil::deserialize(fromCoarseParticles, isCoarseConverged, fromCoarse );
        isCoarseConverged=tmp;
        this->particlesToFine.insert(this->particlesToFine.end(),  fromCoarseParticles.begin(), fromCoarseParticles.end());
        fromCoarseParticles.clear();
    }

    bool tmp_fine = isFineConverged;
    MapperUtil::deserialize(fromFineParticles, isFineConverged, fromFine );
    isFineConverged=tmp_fine;

    this->particlesToCoarse.insert(particlesToCoarse.end(), fromFineParticles.begin(), fromFineParticles.end());
    fromFineParticles.clear();

    // if (this->getMPIRank()==0) logger->getWriter()<<"\tparticlesToCoarse size= "<< particlesToCoarse.size() << "\n";
    // if (this->getMPIRank()==0) logger->getWriter()<<"\tparticlesToFine size= "<< particlesToFine.size() << "\n";
    // if (this->getMPIRank()==0) logger->getWriter()<<"\tisFineConverged="<<isFineConverged<< " isCoarseConverged="<<isCoarseConverged<<"\n";

    //------ check convergence -------------
    statusMessage="NOT REACHED";

    //cout<<"Iteration ["<<cpt<<"]:" <<"isFineConverged= "<< isFineConverged <<", isCoarseConverged= "<< isCoarseConverged << endl;
    if(isCoarseConverged && isFineConverged ){
        statusMessage="convergence REACHED on both submodels -> I will inform them !";
        //here we need to inform both of submodel to not send to me anymore.
        converged=true;
    }else{
        std::ostringstream o;
        if (isCoarseConverged){
            o << coarseId;
            statusMessage="temporarily convergence is detected only on the submodel of [id="+ o.str() +"].";
        }
        if (isFineConverged){
            o << fineId;
            statusMessage="temporarily convergence is detected only on the submodel of  [id="+ o.str() +"].";
        }
    }

    return converged;
}

void TransportInterpolation::interpolate(vector<char> & ToCoarse, vector<char> & ToFine){


    vector<char>().swap(ToFine);
    vector<char>().swap(ToCoarse);
    //----------------------------------------
    //----- interpolate ----------------------

    if(treatCoarse){
        // prepare the send to coarse all accumulated from fine

        for (uint i=0; i< particlesToCoarse.size();++i){
            BoundaryParticle &p =  particlesToCoarse.at(i);
            Double3 newPosition = p.getLagrangianSpeed()*(this->lastFimeTimeStamp - p.getTimeStamp());
            p.setDisplacement(newPosition);
            p.setTimeStamp(this->lastFimeTimeStamp);
        }
        MapperUtil::serialize(particlesToCoarse,this->isAllSubmodelsConverged, ToCoarse);
        particlesToCoarse.clear();
    }

    MapperUtil::serialize(particlesToFine,this->isAllSubmodelsConverged, ToFine);
    particlesToFine.clear();

}

/**********************************************************************************************/
LightTransportInterpolation::LightTransportInterpolation():
    TransportInterpolation(), AsynchronousRemoteMpiKernel(){

    this->fin1= string("f1_in");
    this->fin2= string("f2_in");
    this->fout1= string("f1_out");
    this->fout2= string("f2_out");

}
LightTransportInterpolation::~LightTransportInterpolation(){}


void LightTransportInterpolation::initialize(string id){this->id=id;}
string LightTransportInterpolation::getName(){return this->id;}

void LightTransportInterpolation::mainLoop(){
    //string kname(this->getKernelArgv()[0]);// 1st argument is the kernel name
    ////string kname= this->getName()+".log";
    ////this->logger= std::move(std::make_shared<plb_ofstream>(this->getMPIRank(), kname.c_str()));

    std::cout << "Hey from LightTransportInterpolation::mainLoop(), my rank is " << getMPIRank() << std::endl;

    this->simulate();
}

int LightTransportInterpolation::getSubmodel_ID(){
    return 1;
}

void LightTransportInterpolation::end(){

}
int LightTransportInterpolation::getMPIRank(){
    return this->getMpiManager()->getRank();
}


void LightTransportInterpolation::receiveFromFin1( vector<char> & data){
    /*ConduitExit<char> * tmp = this->getConduitExit<char>(this->fin1);
    int count;
    char * buffer= tmp->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);*/


    ConduitExit<char> * tmp = this->getConduitExit<char>(this->fin1);
    int count;
    char * buffer= this->iRecv(tmp, count);
    this->waitReceive();
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;

}
void LightTransportInterpolation::receiveFromFin2( vector<char> & data){
    /*    ConduitExit<char> * tmp = this->getConduitExit<char>(this->fin2);
    int count;
    char * buffer= tmp->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);*/
    ConduitExit<char> * tmp = this->getConduitExit<char>(this->fin2);
    int count;
    char * buffer= this->iRecv(tmp, count);
    this->waitReceive();
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

void LightTransportInterpolation::sendToCoarse(vector<char> const&data){
    //coarse_out->send((char*)data.data(), data.size());
    this->iSend(coarse_out, (char*)data.data(), data.size());
    this->waitSend();


}
void LightTransportInterpolation::receiveFromCoarse( vector<char> & data){
    /* int count;
    char * buffer= this->coarse_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);*/

    int count;
    char * buffer= this->iRecv(coarse_in, count);
    this->waitReceive();
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;

}

void LightTransportInterpolation::sendToFine(vector<char> const&data){
    //this->fine_out->send((char*)data.data(), data.size());
    this->iSend(fine_out, (char*)data.data(), data.size());
    this->waitSend();
}


void LightTransportInterpolation::receiveFromFine( vector<char> & data){
    /*int count;
    char * buffer= this->fine_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);*/
    int count;
    char * buffer= this->iRecv(fine_in, count);
    this->waitReceive();
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

bool LightTransportInterpolation::receiveFineConvergence()
{
    int count;
    char * buffer= this->fine_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete [] buffer;
    bool isConverged;
    MapperUtil::deserialize(isConverged, vect);
    vector<char>().swap(vect);
    return isConverged;
}
bool LightTransportInterpolation::receiveCoarseConvergence()
{
    int count;
    char * buffer= this->coarse_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete [] buffer;
    bool isConverged;
    MapperUtil::deserialize(isConverged, vect);
    vector<char>().swap(vect);
    return isConverged;
}

void LightTransportInterpolation::setUpCoarseAndFineConduits(double dt1, double dt2){
    if (dt1 < dt2){
        this->fine_in=this->getConduitExit<char>(this->fin1);
        this->coarse_in=this->getConduitExit<char>(this->fin2);
        this->fine_out=this->getConduitEntrance<char>(this->fout1);
        this->coarse_out=this->getConduitEntrance<char>(this->fout2);
    }else {
        this->fine_in=this->getConduitExit<char>(this->fin2);
        this->coarse_in=this->getConduitExit<char>(this->fin1);
        this->fine_out=this->getConduitEntrance<char>(this->fout2);
        this->coarse_out=this->getConduitEntrance<char>(this->fout1);
    }
}



}// namespace

////
