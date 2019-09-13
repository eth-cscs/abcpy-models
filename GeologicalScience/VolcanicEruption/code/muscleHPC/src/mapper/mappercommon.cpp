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

#include "musclehpc/mapper/mappercommon.h"

namespace unige_pasc{


/***************************************** util ***********************************************/

void Submodel_A::simulate(MMSFPerf<double> *mmsfprofiler){

    if (mmsfprofiler) mmsfprofiler->start();
    int cpt = 0;
    const int SAVE_ITR=200;
    // SEL loop
    if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::F_init);
    this->F_init();
    if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::F_init);

    while(! this->isConverged()){

        if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::S);
        this->S();
        if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::S);

        if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::U);
        this->U();
        if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::U);

        if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::Send);
        if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::Send);
        if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::Receive);
        if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::Receive);

        if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::Oi);
        this->Oi();
        if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::Oi);
        if (cpt % SAVE_ITR  == 0 && cpt > 0){
            if (mmsfprofiler) mmsfprofiler->saveToFile();
            cpt = 0;
        }
        cpt ++;
    }

    if (mmsfprofiler) mmsfprofiler->record.registerBeginOperationStamp(SELOp::Of);
    this->Of();
    if (mmsfprofiler) mmsfprofiler->record.registerEndOperationStamp(SELOp::Of);
    if (mmsfprofiler) mmsfprofiler->end();
    if (mmsfprofiler) cout<<" Writing time = "<<mmsfprofiler->getWritingTime()<< endl;

}


/**************************************** BasicMuscleMapper ************************************************/


MapperHeader::MapperHeader(){}
MapperHeader::~MapperHeader(){}


void MapperHeader::serializeHeader(vector<char> & toSend, const INFO_SMART_MAPPER &info_header){
 toSend.resize(sizeof(INFO_SMART_MAPPER)) ;
    memcpy(&toSend[0], &info_header, sizeof(info_header)); // "Serialize"
 //util::serializeHeader(toSend,info_header);
/*
    boost::asio::streambuf write_buffer;
    std::ostream streamToSend(&write_buffer);
    boost::archive::text_oarchive oa(streamToSend);
    // write class instance to archive
    oa << info_header;

    pluint buffSize=   write_buffer.size();
    toSend.resize(buffSize);
    write_buffer.sgetn (&toSend[0],buffSize);
*/
}


void MapperHeader::unserializeHeader(vector<char>  & receivedVect, INFO_SMART_MAPPER  & info_header){
 receivedVect.resize(sizeof(INFO_SMART_MAPPER));
    memcpy(&info_header, &receivedVect[0], sizeof(info_header)); // "unSerialize"
// util::unserializeHeader(receivedVect,info_header);
/*
    membuf write_reader(&receivedVect[0],&receivedVect[0]+receivedVect.size());
    std::istream streamReceived(&write_reader);
    boost::archive::text_iarchive ia(streamReceived);
    // read class state from archive
    INFO_SMART_MAPPER info_tmp;
    ia >> info_tmp;
    info_header=info_tmp;
*/
}

//
}
