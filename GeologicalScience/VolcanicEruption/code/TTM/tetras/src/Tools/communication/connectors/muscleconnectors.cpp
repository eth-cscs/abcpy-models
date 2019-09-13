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

#include "muscleconnectors.h"

namespace unige_pasc{

/************************************** BasicMuscleConnector **************************************************/

BasicMuscleConnector::BasicMuscleConnector():isInitialized (false){}

BasicMuscleConnector::~BasicMuscleConnector(){}


void BasicMuscleConnector::setCommandLineArgs( int *argc, char *** argv){
    this->argc=*argc;
    this->argv=*argv;
}

void BasicMuscleConnector::init (int rank){

    if (this->isInitialized) return;
    this->rankMPI=rank;
    this->isInitialized=true;

    if( rankMPI == 0 )
    {
        muscle::env::init(&argc,&argv);// MUSCLE calls are only permitted in the rank 0 process
        muscle::logger::info("the temporary path is %s", muscle::env::get_tmp_path().c_str());
        //MBB: here  read parmeter from the CxA file
        this->name=(muscle::cxa::kernel_name());
    }
}

void BasicMuscleConnector::send(vector<char> const & buffer, std::string conduit_out){

    if(this->rankMPI == 0){
        uint64_t ss=(uint64_t) (buffer.size());
        muscle::env::send(conduit_out, &ss, 1, MUSCLE_INT64);
        muscle::env::send(conduit_out, &buffer[0], ss, MUSCLE_RAW);
    }
}

void  BasicMuscleConnector::receive(vector<char> & buffer, std::string conduit_in){

    if(this->rankMPI == 0){
        size_t len=1;
        uint64_t ss;
        muscle::env::receive(conduit_in, (void*) &ss,  len, MUSCLE_INT64);
        len = ss;
        buffer.resize(len);
        muscle::env::receive(conduit_in, (void*) &buffer[0], len, MUSCLE_RAW);
    }
}

void  BasicMuscleConnector::end(){
    if(this->rankMPI == 0){
        muscle::env::finalize();
    }
}

string BasicMuscleConnector::getName(){return this->name;}

int  BasicMuscleConnector::getMPIRank() const{
    return this->rankMPI;
}
/*********************************** TetrasConnector *****************************************************/

TetrasConnector::TetrasConnector(): BasicMuscleConnector()
{}

TetrasConnector::~TetrasConnector()
{}

void TetrasConnector::init (int rank){

    //init MUSCLE
    BasicMuscleConnector::init(rank);

    std::string index_str, number_str;
    if( this->rankMPI == 0 )
    {
        cout<<"----------------------------------------------------------------"<<endl;
        //MBB: here  read parmeter from the CxA file
        index_str=muscle::cxa::get_property(this->name+":index");
        number_str=muscle::cxa::get_property("global:sectionsNumber");
        cout <<"[submodel name] : "<<this->name<<endl;
        cout <<"[submodel idx]  : "<<index_str<<endl;
        cout <<"[submodels nbr ]: "<<number_str<<endl;
        cout<<"----------------------------------------------------------------"<<endl;
        MapperUtil::unstringify_pluint(index_str,  this->submodel_ID);
        MapperUtil::unstringify_pluint(number_str, this->submodels_NBR);
        this->conduit_fin=string("fin");
        this->conduit_fout=string("fout");
    }

}

void TetrasConnector::send(vector<char> const & buffer){
    BasicMuscleConnector::send(buffer, this->conduit_fout);
}


void  TetrasConnector::receive(vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->conduit_fin);
}

void  TetrasConnector::end(){

    if( this->rankMPI == 0 ){
        cout<<"-> synchronize ending ...";
        try{
            uint64_t ss=1;
            vector<char> tofinalize;
            tofinalize.resize(sizeof(ss)); //juste send something to synchronize ending
            std::memcpy(&tofinalize[0],&ss, sizeof(ss));
            this->send(tofinalize);
            this->receive(tofinalize);
            BasicMuscleConnector::end();
            cout<<" [OK]"<<endl;
        }
        catch (int e){
            cout<<" [Failed]: "<<e<<endl;
        }
    }

}


/*********************************** PlumeConnector *****************************************************/

TetrasDoubleConnector::TetrasDoubleConnector(): BasicMuscleConnector(){}

TetrasDoubleConnector::~TetrasDoubleConnector(){}
void TetrasDoubleConnector::init (int rank){

    //init MUSCLE
    BasicMuscleConnector::init(rank);

    std::string index_str, number_str;
    if( this->rankMPI == 0 )
    {
        cout<<"----------------------------------------------------------------"<<endl;
        //MBB: here  read parmeter from the CxA file
        index_str=muscle::cxa::get_property(this->name+":index");
        number_str=muscle::cxa::get_property("global:sectionsNumber");
        cout <<"[submodel name] : "<<this->name<<endl;
        cout <<"[submodel idx]  : "<<index_str<<endl;
        cout <<"[submodels nbr ]: "<<number_str<<endl;
        cout<<"----------------------------------------------------------------"<<endl;
        MapperUtil::unstringify_pluint(index_str,  this->submodel_ID);
        MapperUtil::unstringify_pluint(number_str, this->submodels_NBR);
        this->conduit_fin=string("fin");
        this->conduit_fout=string("fout");
        this->fin_atm=string("fin_atm"); // port for the atmosphere model
        this->fout_atm=string("fout_atm");
    }

}

void TetrasDoubleConnector::sendToAtmosphere(vector<char> const & buffer){
    BasicMuscleConnector::send(buffer,  this->fout_atm);
}

void  TetrasDoubleConnector::receiveFromAtmosphere(vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->fin_atm);
}

void TetrasDoubleConnector::sendToContinent(vector<char> const & buffer){
    BasicMuscleConnector::send(buffer, this->conduit_fout);
}

void  TetrasDoubleConnector::receiveFromContinent(vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->conduit_fin);
}

void  TetrasDoubleConnector::end(){

    if( this->rankMPI == 0 ){
        cout<<"-> synchronize ending ...";
        try{
            uint64_t ss=1;
            vector<char> tofinalize, tmp;
            tofinalize.resize(sizeof(ss)); //juste send something to synchronize ending
            std::memcpy(&tofinalize[0],&ss, sizeof(ss));
            this->sendToContinent(tofinalize);
            this->receiveFromContinent(tmp);
            vector<char>().swap(tmp);
            this->sendToAtmosphere(tofinalize);
            this->receiveFromAtmosphere(tmp);
            BasicMuscleConnector::end();
            cout<<" [OK]"<<endl;
        }
        catch (int e){
            cout<<" [Failed]: "<<e<<endl;
        }
    }

}


/*********************************** InterpolationConnector *****************************************************/

TransportInterpolationConnector::TransportInterpolationConnector(): BasicMuscleConnector()
{}

TransportInterpolationConnector::~TransportInterpolationConnector()
{}

void TransportInterpolationConnector::init (int rank){

    //init MUSCLE
    BasicMuscleConnector::init(rank);

    std::string index_str, number_str;
    if( this->rankMPI == 0 )
    {
        cout<<"----------------------------------------------------------------"<<endl;
        //MBB: here  read parmeter from the CxA file
        index_str=muscle::cxa::get_property(this->name+":index");
        number_str=muscle::cxa::get_property("global:sectionsNumber");
        cout <<"[submodel name] : "<<this->name<<endl;
        cout <<"[submodel idx]  : "<<index_str<<endl;
        cout <<"[submodels nbr ]: "<<number_str<<endl;
        cout<<"----------------------------------------------------------------"<<endl;

        MapperUtil::unstringify_pluint(index_str,  this->submodel_ID);
        MapperUtil::unstringify_pluint(number_str, this->submodels_NBR);

        this->fin1= string("fin1");
        this->fin2= string("fin2");
        this->fout1= string("fout1");
        this->fout2= string("fout2");
    }

}

pluint TransportInterpolationConnector::getSubmodel_ID(){return this->submodel_ID;}

pluint TransportInterpolationConnector::getsubmodels_NBR(){return this->submodels_NBR;}

void TransportInterpolationConnector::receiveFromFin1( vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->fin1);
}
void TransportInterpolationConnector::receiveFromFin2( vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->fin2);
}

void TransportInterpolationConnector::sendToCoarse(vector<char> const&buffer){
    BasicMuscleConnector::send(buffer, this->coarse_out);
}
void TransportInterpolationConnector::receiveFromCoarse( vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->coarse_in);
}
void TransportInterpolationConnector::sendToFine(vector<char> const&buffer){
    BasicMuscleConnector::send(buffer, this->fine_out);
}
void TransportInterpolationConnector::receiveFromFine( vector<char> & buffer){
    BasicMuscleConnector::receive(buffer, this->fine_in);
}

void TransportInterpolationConnector:: setUpCoarseAndFineConduits(double dt1, double dt2){
    if (dt1 < dt2){
        this->fine_in=this->fin1;
        this->fine_out=this->fout1;
        this->coarse_in=this->fin2;
        this->coarse_out=this->fout2;
    }else {
        this->coarse_in=this->fin1;
        this->coarse_out=this->fout1;
        this->fine_in=this->fin2;
        this->fine_out=this->fout2;
    }
}

void TransportInterpolationConnector::setCoarseIn(string conduit){this->coarse_in=conduit;}
void TransportInterpolationConnector::setCoarseOut(string conduit){this->coarse_out=conduit;}
void TransportInterpolationConnector::setFineIn(string conduit){this->fine_in=conduit;}
void TransportInterpolationConnector::setFineOut(string conduit){this->fine_out=conduit;}

void  TransportInterpolationConnector::end(){
    if( this->rankMPI == 0 ){
        cout<<"-> synchronize ending ...";
        try{uint64_t ss=1;
            vector<char> tofinalize;
            tofinalize.resize(sizeof(ss)); //juste send something to synchronize ending
            std::memcpy(&tofinalize[0],&ss, sizeof(ss));
            sendToCoarse(tofinalize);
            sendToFine(tofinalize);
            receiveFromCoarse(tofinalize);
            receiveFromFine(tofinalize);
            cout<<" [OK]"<<endl;
        }catch(int e){
            cout<<" [Failed]: "<<e<<endl;
        }
    }
    BasicMuscleConnector::end();

}
/*************************************************************************************************/

}//end name space
