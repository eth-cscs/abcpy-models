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

#ifndef MMSFPERF_H
#define MMSFPERF_H


#include <string>       // std::string
#include <iostream>
#include <sstream>      // std::stringstream
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
using namespace std;

//*******************TIME stamp*******************
template<typename T>
T getTimeStamp(){
#ifdef USE_PARALLEL_MPI
    return MPI_Wtime();
#else
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    double endTime = (double) ts.tv_nsec * (double) 1.0e-9;
    return endTime;
#endif
}



enum SELOp {F_init, Oi, S, U, Of, Send, Receive};
//||||---------
template<typename T>
class SELTimeStamp{

public:
    SELTimeStamp(): fieldSeparator("\t"){
        computeTime=0.0;
        Finit_Time=0.0;
        Oi_Time=0.0;
        S_Time=0.0;
        U_Time=0.0;
        Of_Time=0.0;
        receive_Time=0.0;
        send_Time=0.0;
    }

    //affectation
    SELTimeStamp<T> & operator=(const SELTimeStamp<T> & rhs)
    {

        if (this != &rhs){
            if (this->stamp_begin.size()>0){
                this->stamp_begin.clear();
                this->stamp_begin.insert(rhs.stamp_begin.begin(), rhs.stamp_begin.end());
            }
            if (this->stamp_end.size()>0){
                this->stamp_end.clear();
                this->stamp_end.insert(rhs.stamp_end.begin(), rhs.stamp_end.end());
            }

            computeTime=rhs.computeTime;
            Finit_Time=rhs.Finit_Time;
            Oi_Time=rhs.Oi_Time;
            S_Time=rhs.S_Time;
            U_Time=rhs.U_Time;
            Of_Time=rhs.Of_Time;
            receive_Time=rhs.receive_Time;
            send_Time=rhs.send_Time;
        }
        return *(this);
    }


    int getIterationNumberSofar(){
        return stamp_begin.at(S).size();
    }

    void updateTotalTimeAndClear(){


        Oi_Time+=this->getOperationTime(Oi);
        S_Time+=this->getOperationTime(S);
        U_Time+=this->getOperationTime( U);
        send_Time+=this->getOperationTime( Send);
        receive_Time+=this->getOperationTime( Receive);

        stamp_begin.at(Oi).clear();
        stamp_end.at(Oi).clear();
        stamp_begin.at(S).clear();
        stamp_end.at(S).clear();
        stamp_begin.at(U).clear();
        stamp_end.at(U).clear();
        stamp_begin.at(Send).clear();
        stamp_end.at(Send).clear();
        stamp_begin.at(Receive).clear();
        stamp_end.at(Receive).clear();
    }

    void registerBeginOperationStamp(SELOp op, T stamp){
        if(stamp_begin.find(op) == stamp_begin.end()){
            vector<T> vect;
            stamp_begin.insert ( std::pair<SELOp, vector<T> >(op, vect));
        }
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        stamp_begin_op.push_back(stamp);
    }
    void registerBeginOperationStamp(SELOp op){
       registerBeginOperationStamp(op, getTimeStamp<T>());
    }

    void registerEndOperationStamp(SELOp op, T stamp){
        if(stamp_end.find(op) == stamp_end.end()){
            vector<T> vect;
            stamp_end.insert ( std::pair<SELOp, vector<T> >(op, vect));
        }
        vector<T> & stamp_end_op=stamp_end.at(op);
        stamp_end_op.push_back(stamp);

        if(op == SELOp::F_init){
            Finit_Time+=this->getOperationTime(F_init);
        }
        if(op == SELOp::Of){
            Of_Time+=this->getOperationTime(Of);
        }
    }
    void registerEndOperationStamp(SELOp op){
        registerEndOperationStamp(op, getTimeStamp<T>());
    }

    //getters
    T getOperationBeginStamp(int itr, SELOp op){
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        return stamp_begin_op.at(itr);
    }
    T getOperationEndStamp(int itr, SELOp op){
        vector<T> & stamp_end_op=stamp_end.at(op);
        return stamp_end_op.at(itr);
    }

    std::string toString(){
        std::stringstream ss;
        ss.precision(numeric_limits<double>::digits10 + 1);
        vector<T> & stamp_begin_op=stamp_begin.at(S);
        for(sz_vector it = 0; it < stamp_begin_op.size(); it++) {
            ss<<getOperationBeginStamp(it, Oi)<<fieldSeparator<<getOperationEndStamp(it, Oi)<<fieldSeparator;
            ss<<getOperationBeginStamp(it, S)<<fieldSeparator<<getOperationEndStamp(it, S)<<fieldSeparator;
            ss<<getOperationBeginStamp(it, U)<<fieldSeparator<<getOperationEndStamp(it, U)<<fieldSeparator;
            ss<<getOperationBeginStamp(it, Send)<<fieldSeparator<<getOperationEndStamp(it,Send)<<fieldSeparator;
            ss<<getOperationBeginStamp(it, Receive)<<fieldSeparator<<getOperationEndStamp(it, Receive)<<"\n";
        }
        return  ss.str();
    }

    T getOperationTime(int itr, SELOp op){
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        vector<T> & stamp_end_op=stamp_end.at(op);
        return stamp_end_op.at(itr)-stamp_begin_op.at(itr);
    }

    T getOperationTime(SELOp op){
        T timeOp=0.0;
        for(sz_vector it = 0; it < stamp_begin.at(op).size(); it++) {
            timeOp += this->getOperationTime(it, op);
        }
        return timeOp;
    }


    T getOp_Time(SELOp op){
        if(op == SELOp::F_init){
            return this->Finit_Time;
        }else if(op == SELOp::Oi){
            return this->Oi_Time;
        }else if(op == SELOp::S){
            return this->S_Time;
        }else if(op == SELOp::U){
            return this->U_Time;
        }else if(op == SELOp::Send){
            return this->send_Time;
        }else if(op == SELOp::Receive){
            return this->receive_Time;
        }else if(op == SELOp::Of){
            return this->Of_Time;
        }else{
            throw std::invalid_argument( "getOp_Time(): received invalide operation" );
        }
        return -1;
    }





private:
    typedef typename std::vector<T>::size_type sz_vector;
    map<SELOp, vector<T> >  stamp_begin;
    map<SELOp, vector<T> >  stamp_end;
    string fieldSeparator;
    T computeTime, Finit_Time, Oi_Time, S_Time, U_Time, Of_Time, receive_Time, send_Time;


};

//||||---------
template<typename T>
class MMSFPerf{

private:

    T startComputationTimeStamp;
    T endComputationTimeStamp;
    bool isStarted;
    bool isStopped;
    bool isFileCreated;
    int rank;
    string filename;
    int nbr_itr;
    std::shared_ptr<ofstream> ofile;
    T saveTime;

public:
    SELTimeStamp<T> record;
    MMSFPerf(){
        isStarted=false;
        isStopped=false;
        isFileCreated=false;
        rank=0;
        filename="submodel";
        nbr_itr=0;
        saveTime=0.0;
    }


    void setRank(int rank){
        this->rank = rank;
    }
    void setPrefixFileName(string prefix){
        filename = prefix;
    }



    void start(){
        if(! isStarted){
            startComputationTimeStamp=getTimeStamp<T>();
            isStarted=true;
            //if(rank==0)
            //        cout<<"# started timestamp "<< filename<< ": "<<std::setprecision (numeric_limits<double>::digits10 + 1)<< startComputationTimeStamp <<endl;
        }
    }

    void reStart(){
        isStarted=false;
        start();
    }

    void end(){
        endComputationTimeStamp=getTimeStamp<T>();
        isStopped=true;
        saveToFile();
    }

    T getWritingTime(){
        return this->saveTime;
    }

    void saveToFile(){
        if(rank==0){
            try{
                double t1=MPI_Wtime();
                setFileName();
                //ofstream ofile(filename.c_str(),std::ostream::app);//append mode
                cout<<"-> Saving profiling to "<<filename.c_str()<<endl;
                *ofile << getIterationDataString();
                record.updateTotalTimeAndClear();
                if(isStopped){//save average stats
                    *ofile << getAveragedStats();
                    *ofile << "#Profiling writing (overhead time) =" << this->saveTime <<endl;
                     ofile->close();
                     return;
                }
                this->saveTime += MPI_Wtime() - t1;

            } catch (std::exception & e){
                cerr<< "-> Exception occured: "<< e.what()<<endl;
            }
        }
    }

    std::string  getIterationDataString(){
        this->nbr_itr+=T(record.getIterationNumberSofar());
        return record.toString();
    }


    std::string getAveragedStats(){

        std::stringstream ss;
        ss<<"#=========================================================="<<endl;
        ss<<"#T_compute\tFinit\tOi\tS\tU\tSend\tReceive\tOf"<<endl;
        ss<<"#nbr_it="<<nbr_itr<<endl;
        ss<<"#T_Finit_b(0)\tT_Finit_e(1)\tOf_b(2)\tOf_e(3)\n";
        ss<<"#"<<std::setprecision (numeric_limits<double>::digits10 + 1)<<record.getOperationBeginStamp(0, SELOp::F_init)<<"\t"<<record.getOperationEndStamp(0, SELOp::F_init)<<"\t";
        ss<<std::setprecision (numeric_limits<double>::digits10 + 1)<<record.getOperationBeginStamp(0, SELOp::Of)<<"\t"<<record.getOperationEndStamp(0, SELOp::Of)<<"\n";
        ss<<"begin time="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<startComputationTimeStamp<<endl;
        ss<<"end time="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<endComputationTimeStamp<<endl;

        T computeTime, Finit_Time, Oi_Time, S_Time, U_Time, send_Time, receive_Time, Of_Time;
        computeTime = 0.0;
        Finit_Time = record.getOp_Time(SELOp::F_init);
        computeTime += Finit_Time;
        T tmpValue=record.getOp_Time(SELOp::Oi);
        computeTime += tmpValue;
        Oi_Time = tmpValue/nbr_itr;
        tmpValue=record.getOp_Time(SELOp::S);
        computeTime += tmpValue;
        S_Time = tmpValue/nbr_itr;
        tmpValue = record.getOp_Time(SELOp::U);
        computeTime += tmpValue;
        U_Time = tmpValue/nbr_itr;
        tmpValue = record.getOp_Time(SELOp::Send);
        computeTime += tmpValue;
        send_Time = tmpValue/nbr_itr;
        tmpValue = record.getOp_Time(SELOp::Receive);
        computeTime += tmpValue;
        receive_Time = tmpValue/nbr_itr;
        Of_Time = record.getOp_Time(SELOp::Of);
        computeTime += Of_Time;

        //T computeTime= Finit_Time + (Oi_Time + S_Time + U_Time + send_Time + receive_Time)* nbr_itr + Of_Time;

        ss<<(computeTime)<<"\t";
        ss<<(Finit_Time)<<"\t";
        ss<<(Oi_Time)<<"\t";
        ss<<(S_Time)<<"\t";
        ss<<(U_Time)<<"\t";
        ss<<(send_Time)<<"\t";
        ss<<(receive_Time)<<"\t";
        ss<<(Of_Time)<<endl;

        //convert the stream buffer into a string
        ss<<"#=========================================================="<<endl;
        return  ss.str();
    }

private:

    void setFileName(){

        if(!isFileCreated){
            std::stringstream fname;
            fname<< filename<<"_stamp_"<< std::fixed<< std::setprecision (8)<<startComputationTimeStamp<<".profiler";
            filename=fname.str();
            std::stringstream ss;

            ss<<"T_Oi_b(0)\tT_Oi_e(1)";
            ss<<"\tT_S_b(2)\tT_S_e(3)";
            ss<<"\tT_U_b(4)\tT_U_e(5)";
            ss<<"\tT_send_b(6)\tT_send_e(7)";
            ss<<"\tT_receive_b(8)\tT_receive_e(9)";
            ss<<endl;
            ofile = std::make_shared<ofstream>(filename.c_str(),std::ostream::app);
            //ofstream ofile(filename.c_str(),std::ostream::app);//append mode
            *ofile << ss.str();
            //ofile->close();
            isFileCreated=true;
        }
    }
};

#endif // MMSFPERF_H
