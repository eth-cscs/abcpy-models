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

#ifndef PROFILING_H
#define PROFILING_H

#include <string>       // std::string
#include <iostream>
#include <sstream>      // std::stringstream
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
#ifdef USE_PARALLEL_MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <iostream>
#include <fstream>
#include <memory>

using namespace std;

namespace unige_pasc {

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

/*void writeStringToFile( std::string outputFilename, std::string content  ){
   // std::vector<char> buffer(content.begin(),content.end());
    std::ofstream outputBuffer(outputFilename.c_str(), std::ios::binary|std::ios::out);
    outputBuffer<< content;
    outputBuffer.close();
}*/



// This class will register all the time stamp of the operation during one iteration
// of the 3DFS model. This includes all the 5 levels of data processors
template<typename T>
struct IterationTimeStamp {


    enum Operation { copyOverlap, preProcess, process, gather, send, receive, update, sendrec};
    typedef typename std::vector<T>::size_type sz_vector;



    IterationTimeStamp( bool useParallelGather=false, int LEVELS_NBR=1,  bool isOverlap=false)
        :LEVELS_NBR(LEVELS_NBR), useParallelGather(useParallelGather), isOverlap(isOverlap)
    {
        vector<T> stamp_begin_copyOverlap; // for each level computation to prepare the gather data for border
        vector<T> stamp_begin_preProcess; // for each level computation to prepare the gather data for border
        vector<T> stamp_begin_process; // for each level computation
        vector<T> stamp_begin_gather;// for each level gather time
        vector<T> stamp_begin_send; // for each level send time
        vector<T> stamp_begin_receive;// for each level receive time
        vector<T> stamp_begin_update;// for each level update time

        stamp_begin.insert ( std::pair<Operation,vector<T> >(copyOverlap,stamp_begin_copyOverlap) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(preProcess,stamp_begin_preProcess) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(process,stamp_begin_process) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(gather,stamp_begin_gather) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(send,stamp_begin_send) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(receive,stamp_begin_receive) );
        stamp_begin.insert ( std::pair<Operation,vector<T> >(update,stamp_begin_update) );

        vector<T> stamp_end_copyOverlap;
        vector<T> stamp_end_preProcess; // for each level computation to prepare the gather data for border
        vector<T> stamp_end_process;
        vector<T> stamp_end_gather;
        vector<T> stamp_end_send;
        vector<T> stamp_end_receive;
        vector<T> stamp_end_update;

        stamp_end.insert ( std::pair<Operation,vector<T> >(copyOverlap,stamp_end_copyOverlap) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(preProcess,stamp_end_preProcess) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(process,stamp_end_process) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(gather,stamp_end_gather) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(send,stamp_end_send) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(receive,stamp_end_receive) );
        stamp_end.insert ( std::pair<Operation,vector<T> >(update,stamp_end_update) );
    }


    // copy and affectation constructors are they needed ??
    //copy constructor
    IterationTimeStamp(const IterationTimeStamp<T> & rhs)
        :
          stamp_begin(rhs.stamp_begin),
          stamp_end(rhs.stamp_end),
          useParallelGather (rhs.useParallelGather)
    {}

    //affectation
    IterationTimeStamp<T> & operator=(const IterationTimeStamp<T> & rhs)
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
            this->useParallelGather =rhs.useParallelGather;
        }
        return *(this);
    }

    void registerBeginOperationStamp(int level, Operation op){
        registerBeginOperationStamp(level,op,unige_pasc::getTimeStamp<T>());
    }
    void registerEndOperationStamp(int level, Operation op){
        registerEndOperationStamp(level,op,unige_pasc::getTimeStamp<T>());
    }

    void registerBeginOperationStamp(int level, Operation op, T stamp){
        assert(level<LEVELS_NBR);
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        stamp_begin_op.insert(stamp_begin_op.begin()+level,stamp);
    }
    void registerEndOperationStamp(int level, Operation op, T stamp){
        assert(level<LEVELS_NBR);
        vector<T> & stamp_end_op=stamp_end.at(op);
        stamp_end_op.insert(stamp_end_op.begin()+level,stamp);
    }

    void registerBeginCopyOverlapStamp(int level){
        registerBeginOperationStamp(level,copyOverlap);
    }
    void registerEndCopyOverlapStamp(int level){
        registerEndOperationStamp(level,copyOverlap);
    }

    void registerBeginPreProcessLevelStamp(int level){
        registerBeginOperationStamp(level,preProcess);
    }
    void registerEndPreProcessLevelStamp(int level){
        registerEndOperationStamp(level,preProcess);
    }

    void registerBeginProcessLevelStamp(int level){
        registerBeginOperationStamp(level,process);
    }
    void registerBeginProcessLevelStamp(int level, T const&timeStamp){
        registerBeginOperationStamp(level,process,timeStamp);
    }
    void registerEndProcessLevelStamp(int level){
        registerEndOperationStamp(level,process);
    }
    void registerEndProcessLevelStamp(int level,T const&timeStamp){
        registerEndOperationStamp(level,process,timeStamp);
    }

    void registerBeginGatherLevelStamp(int level){
        registerBeginOperationStamp(level,gather);
    }
    void registerEndGatherLevelStamp(int level){
        registerEndOperationStamp(level,gather);
    }

    void registerBeginUpdateLevelStamp(int level){
        registerBeginOperationStamp(level,update);
    }
    void registerEndUpdateLevelStamp(int level){
        registerEndOperationStamp(level,update);
    }

    void registerSendReceiveLevelStamp(T WTimeSendBegin, T WTimeSendEnd, T WTimeReceiveBegin, T WTimeReceiveEnd, int level){
        registerBeginOperationStamp(level,send,WTimeSendBegin);
        registerEndOperationStamp(level,send, WTimeSendEnd);
        registerBeginOperationStamp(level,receive,WTimeReceiveBegin);
        registerEndOperationStamp(level,receive,WTimeReceiveEnd);
    }

    //getters
    T getOperationBeginStamp(int level, Operation op){
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        return stamp_begin_op.at(level);
    }
    T getOperationEndStamp(int level, Operation op){
        vector<T> & stamp_end_op=stamp_end.at(op);
        return stamp_end_op.at(level);
    }

    T getOperationTime(int level,Operation op){
        if(op != sendrec){
            vector<T> & stamp_begin_op=stamp_begin.at(op);
            vector<T> & stamp_end_op=stamp_end.at(op);
            return stamp_end_op.at(level)-stamp_begin_op.at(level);
        }else{
            vector<T> & stamp_begin_send=stamp_begin.at(send);
            vector<T> & stamp_end_receive=stamp_end.at(receive);
            return stamp_end_receive.at(level)-stamp_begin_send.at(level);
        }
    }

    T getOperationTime(Operation op){
        T timeOp=0.0;
        Operation op_local=op;
        if(op ==sendrec){
            op_local=send;
        }
        for(sz_vector pl = 0; pl != stamp_begin.at(op_local).size(); pl++) {
            timeOp+=this->getOperationTime(pl,op);
        }
        return timeOp;
    }

    T getOperationTimeAverage(Operation op){
        vector<T> & stamp_begin_op=stamp_begin.at(op);
        T levels=T(stamp_begin_op.size());
        return getOperationTime(op)/levels;
    }


    T getCopyTime(int level){
        return getOperationTime(level,copyOverlap);
    }
    T getGatherTime(int level){
        return getOperationTime(level,gather);
    }
    T getPreProcessingTime(int level){
        return getOperationTime(level,preProcess);
    }
    T getProcessingTime(int level){
        return getOperationTime(level,process);
    }
    T getSendReceiveTime(int level){
        return getOperationTime(level,sendrec);
    }
    T getUpdateTime(int level){
        return getOperationTime(level,update);
    }

    //util
    T computeTotalTime(int level){

        T iterationLevelTime=this->getCopyTime(level)+this->getPreProcessingTime(level);

        T effectiveComputingTime=this->getProcessingTime(level);

        if(isOverlap){
            T timeActionInParallel=this->getSendReceiveTime(level);

            if(this->useParallelGather)
                timeActionInParallel+=this->getGatherTime(level);

            if(effectiveComputingTime < timeActionInParallel)
                effectiveComputingTime=timeActionInParallel;
            if(!this->useParallelGather)
                iterationLevelTime+=this->getGatherTime(level);
        }else{// not overlap
            effectiveComputingTime+=this->getGatherTime(level)+this->getSendReceiveTime(level);
        }


        iterationLevelTime+=effectiveComputingTime+this->getUpdateTime(level);
        return iterationLevelTime;

    }

    T computeTotalTime(){
        T TotalTime=0.0;
        vector<T> & stamp_begin_op=stamp_begin.at(copyOverlap);
        for(sz_vector pl = 0; pl != stamp_begin_op.size(); pl++) {
            TotalTime+=this->computeTotalTime(pl);
        }
        return TotalTime;
    }

    std::string toString(){
        std::stringstream ss;
        ss.precision(numeric_limits<double>::digits10 + 1);
        /*  ss<<"#T_copy_b\tT_copy_e";
        ss<<"\tT_preProcess_b\tT_preProcess_e";
        ss<<"\tT_process_b\tT_process_e";
        ss<<"\tT_gather_b\tT_gather_e";
        ss<<"\tT_send_b\tT_send_e";
        ss<<"\tT_receive_b\tT_receive_e";
        ss<<"\tT_update_b\tT_update_e"<<endl;*/
        vector<T> & stamp_begin_op=stamp_begin.at(copyOverlap);
        for(sz_vector pl = 0; pl != stamp_begin_op.size(); pl++) {
            ss<<getOperationBeginStamp(pl,copyOverlap)<<"\t"<<getOperationEndStamp(pl,copyOverlap)<<"\t";
            ss<<getOperationBeginStamp(pl,preProcess)<<"\t"<<getOperationEndStamp(pl,preProcess)<<"\t";
            ss<<getOperationBeginStamp(pl,process)<<"\t"<<getOperationEndStamp(pl,process)<<"\t";
            ss<<getOperationBeginStamp(pl,gather)<<"\t"<<getOperationEndStamp(pl,gather)<<"\t";
            ss<<getOperationBeginStamp(pl,send)<<"\t"<<getOperationEndStamp(pl,send)<<"\t";
            ss<<getOperationBeginStamp(pl,receive)<<"\t"<<getOperationEndStamp(pl,receive)<<"\t";
            ss<<getOperationBeginStamp(pl,update)<<"\t"<<getOperationEndStamp(pl,update);
        }
        return  ss.str();
    }

    std::string getStats(){

        std::stringstream ss;
        ss.precision(numeric_limits<double>::digits10 + 1);
        // vector<T> & stamp_begin_op=stamp_begin.at(copyOverlap);
        //for(sz_vector pl = 0; pl != stamp_begin_op.size(); pl++) {
        ss<<this->computeTotalTime()<<"\t";
        ss<<this->getOperationTime(copyOverlap)<<"\t";
        ss<<this->getOperationTime(preProcess)<<"\t";
        ss<<this->getOperationTime(process)<<"\t";
        ss<<this->getOperationTime(gather)<<"\t";
        ss<<this->getOperationTime(send)<<"\t";
        ss<<this->getOperationTime(receive)<<"\t";
        ss<<this->getOperationTime(sendrec)<<"\t";
        ss<<this->getOperationTime(update)<<"\n";

        //}
        return  ss.str();
    }



    //number of levels in the 3DFS data processors call
    //static const int LEVELS_NBR = 5;
    int LEVELS_NBR;
    map<Operation, vector<T> >  stamp_begin;
    map<Operation, vector<T> >  stamp_end;

    bool useParallelGather;
    bool isOverlap;

};

/****************************************************************************************/
template<typename T>
struct TimeStampMuscle {
    // time stamp to measure the synchronization time of muscle coupling

    typedef typename std::vector<IterationTimeStamp<T> >::iterator it_type;

    TimeStampMuscle(){
        isStarted=false;
        isStopped=false;
        isFileCreated=false;
        computeTime=0.0;
        copyTime=0.0;
        preProcessTime=0.0;
        processTime=0.0;
        gatherTime=0.0;
        sendTime=0.0;
        receiveTime=0.0;
        sendrecTime=0.0;
        updateTime=0.0;
        nbr_itr=0.0;
        this->rank=0;
        prefix="submodel";
    }


    //copy constructor
    TimeStampMuscle(const TimeStampMuscle<T> & timeStamp)
        : //stampTimeIteraions(timeStamp.stampTimeIteraions),
          startComputationTimeStamp (timeStamp.startComputationTimeStamp),
          endComputationTimeStamp(timeStamp.endComputationTimeStamp),
          isStarted(timeStamp.isStarted),
          isStopped(timeStamp.isStopped)
    {
        computeTime=timeStamp.computeTime;
        copyTime=timeStamp.copyTime;
        preProcessTime=timeStamp.preProcessTime;
        processTime=timeStamp.processTime;
        gatherTime=timeStamp.gatherTime;
        sendTime=timeStamp.endTime;
        receiveTime=timeStamp.receiveTime;
        sendrecTime=timeStamp.sendrecTime;
        updateTime=timeStamp.updateTime;
        isFileCreated=timeStamp.isFileCreated;
        rank=timeStamp.rank;
        prefix=timeStamp.prefix;

        nbr_itr=0.0;
        if (stampTimeIteraions.size()>0){
            this->stampTimeIteraions.clear();
            this->stampTimeIteraions.resize(timeStamp.stampTimeIteraions.size());
            copy(timeStamp.stampTimeIteraions.begin(), timeStamp.stampTimeIteraions.end(), this->stampTimeIteraions.begin());
        }
    }

    //affectation
    TimeStampMuscle<T> & operator=(const TimeStampMuscle<T> & timeStamp)
    {

        if (this != &timeStamp ){
            startComputationTimeStamp =timeStamp.startComputationTimeStamp;
            endComputationTimeStamp=timeStamp.endComputationTimeStamp;
            isStarted=timeStamp.isStarted;
            isStopped=timeStamp.isStopped;
            computeTime=timeStamp.computeTime;
            copyTime=timeStamp.copyTime;
            preProcessTime=timeStamp.preProcessTime;
            processTime=timeStamp.processTime;
            gatherTime=timeStamp.gatherTime;
            sendTime=timeStamp.endTime;
            receiveTime=timeStamp.receiveTime;
            sendrecTime=timeStamp.sendrecTime;
            updateTime=timeStamp.updateTime;
            isFileCreated=timeStamp.isFileCreated;
            rank=timeStamp.rank;
            prefix=timeStamp.prefix;


            if (stampTimeIteraions.size()>0){
                this->stampTimeIteraions.clear();
                this->stampTimeIteraions.resize(timeStamp.stampTimeIteraions.size());
                copy(timeStamp.stampTimeIteraions.begin(), timeStamp.stampTimeIteraions.end(), this->stampTimeIteraions.begin());
            }
        }

        return *(this);
    }


    void setRank(int rank){
        this->rank=rank;
    }

    void start(){
        if(! isStarted){
            startComputationTimeStamp=unige_pasc::getTimeStamp<T>();
            isStarted=true;
            //setFileName();
            //writeHeader();
        }
    }

    void reStart(){
        isStarted=false;
        start();
    }

    void end(bool considerParticles=true){
        endComputationTimeStamp=unige_pasc::getTimeStamp<T>();
        isStopped=true;
        saveToFile(considerParticles);
    }



    void setFileName(){
        if(!isFileCreated){
            std::stringstream fname;
            fname<< prefix<<"_stamp_"<< std::fixed<< std::setprecision (8)<<startComputationTimeStamp<<".dat";
            filename=fname.str();
            std::stringstream ss;
            ss<<"#T_copy_b(0)\tT_copy_e(1)";
            ss<<"\tT_preProcess_b(2)\tT_preProcess_e(3)";
            ss<<"\tT_process_b(4)\tT_process_e(5)";
            ss<<"\tT_gather_b(6)\tT_gather_e(7)";
            ss<<"\tT_send_b(8)\tT_send_e(9)";
            ss<<"\tT_receive_b(10)\tT_receive_e(11)";
            ss<<"\tT_update_b(12)\tT_update_e(13)";
            ss<<"\tinjected(14)\tdeposed(15)\toutside(16)";
            ss<<"\ttotalParticlesInDomain(17)\tsum_injectedParticles(18)\ttotalParticlesInTerrain(19)\tsum_outsideParticles(20)"<<endl;


            ofstream ofile(filename.c_str(),std::ostream::app);//append mode
            ofile << ss.str();
            ofile.close();
            isFileCreated=true;
        }
    }

    void saveToFile(bool considerParticles=true){
        if(rank==0){
            setFileName();
            cout<<"-> Saving profiling to "<<filename.c_str()<<endl;
            ofstream ofile(filename.c_str(),std::ostream::app);//append mode
            ofile << registerComputedTimeIterations(considerParticles);
            if(isStopped){//save average stats
                ofile << getAveragedStats() << "\n";
            }
            ofile.close();
        }
    }

    std::string  registerComputedTimeIterations(bool considerParticles=true){
        std::stringstream ss;
        int itr=0;
        for (it_type it=stampTimeIteraions.begin(); it != stampTimeIteraions.end(); ++it){
            //get string info
            string particles="";
            if(considerParticles)
                particles=getParticlesInformation(itr);
            ss<<it->toString()<<"\t"<<particles<<"\n";
            //ccompte the rest
            computeTime+=it->computeTotalTime();
            copyTime+=it->getOperationTime(IterationTimeStamp<T>::copyOverlap);
            preProcessTime+=it->getOperationTime(IterationTimeStamp<T>::preProcess);
            processTime+=it->getOperationTime(IterationTimeStamp<T>::process);
            gatherTime+=it->getOperationTime(IterationTimeStamp<T>::gather);
            sendTime+=it->getOperationTime(IterationTimeStamp<T>::send);
            receiveTime+=it->getOperationTime(IterationTimeStamp<T>::receive);
            sendrecTime+=it->getOperationTime(IterationTimeStamp<T>::sendrec);
            updateTime+=it->getOperationTime(IterationTimeStamp<T>::update);
            itr++;

        }
        this->nbr_itr+=T(stampTimeIteraions.size());
        stampTimeIteraions.clear();
        this->clearParticlesvectors(); // <-- particles
        return ss.str();
    }


    std::string toString(){
        std::stringstream ss;
        ss<<"begin="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<startComputationTimeStamp<<endl;
        ss<<"end="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<endComputationTimeStamp<<endl;
        ss<<"#-----------------------------------------------------------------"<<endl;
        ss<<"#size of the time-stamps:"<<stampTimeIteraions.size()<<endl;
        for (it_type it=stampTimeIteraions.begin(); it != stampTimeIteraions.end(); ++it){
            //IterationTimeStamp<T> * itr=it;
            ss<<it->toString()<<endl;
        }
        //convert the stream buffer into a string
        ss<<"-----------------------------------------------------------------"<<endl;
        return  ss.str();
    }

    std::string getStats(){
        std::stringstream ss;
        ss<<"#=========================================================="<<endl;
        ss<<"#size of the time-stamps:"<<stampTimeIteraions.size()<<endl;
        ss<<"#T_total\tT_copy\tT_preProcess\tT_process\tT_gather\tT_send\tT_receive\tT_transfer\tT_update"<<endl;
        ss<<"begin="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<startComputationTimeStamp<<endl;
        ss<<"end="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<endComputationTimeStamp<<endl;
        for (it_type it=stampTimeIteraions.begin(); it != stampTimeIteraions.end(); ++it){
            //IterationTimeStamp<T> * itr=it;
            ss<<it->getStats();
        }
        //convert the stream buffer into a string
        ss<<"#=========================================================="<<endl;

        return  ss.str();
    }

    std::string getAveragedStats(){
        std::stringstream ss;
        ss<<"#=========================================================="<<endl;
        ss<<"#size of the time-stamps:"<<stampTimeIteraions.size()<<endl;
        ss<<"#T_total\tCompute\tT_copy\tT_preProcess\tT_process\tT_gather\tT_send\tT_receive\tT_transfer\tT_update"<<endl;
        ss<<"begin time="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<startComputationTimeStamp<<endl;
        ss<<"end time="<<std::setprecision (numeric_limits<double>::digits10 + 1)<<endComputationTimeStamp<<endl;
        ss<<"nbr_it="<<nbr_itr<<"\t";
        ss<<(computeTime/nbr_itr)<<"\t";
        ss<<(copyTime/nbr_itr)<<"\t";
        ss<<(preProcessTime/nbr_itr)<<"\t";
        ss<<(processTime/nbr_itr)<<"\t";
        ss<<(gatherTime/nbr_itr)<<"\t";
        ss<<(sendTime/nbr_itr)<<"\t";
        ss<<(receiveTime/nbr_itr)<<"\t";
        ss<<(sendrecTime/nbr_itr)<<"\t";
        ss<<(updateTime/nbr_itr)<<endl;
        //convert the stream buffer into a string
        ss<<"#=========================================================="<<endl;

        return  ss.str();
    }




    T getTotalIterationsTime(){
        /*T compuetTime=0.0;
        for (it_type it=stampTimeIteraions.begin(); it != stampTimeIteraions.end(); ++it){
            //IterationTimeStamp<T> * itr=it;
            compuetTime+= it->computeTotalTime();
        }
        return compuetTime;*/
        return computeTime;
    }

    T getAverageIterationTime(){
        // T numberofIterations=T(stampTimeIteraions.size());
        return getTotalIterationsTime()/nbr_itr;
    }


    void registerInjectedparticles (int injectedPar){
        injectedParticles.push_back(injectedPar);
    }
    void registerDeposedparticles (int deposedPar){
        deposedPrticules.push_back(deposedPar);
    }
    void registerOutSideparticles (int outPar){
        outsideparticles.push_back(outPar);
    }

    void registerTotalParticlesInDomain (int parts){
        totalParticlesInDomain.push_back(parts);
    }
    void registerSumInjectedParticles (int parts){
        sum_injectedParticles.push_back(parts);
    }
    void registerTotalParticlesInTerrain (int parts){
        totalParticlesInTerrain.push_back(parts);
    }
    void registerSumOutsideParticles (int parts){
        sum_outsideParticles.push_back(parts);
    }


    string getParticlesInformation(int cpt){
        std::stringstream ss;
        ss<< injectedParticles.at(cpt)<<"\t"<<deposedPrticules.at(cpt)<<"\t"<<outsideparticles.at(cpt)<<"\t";
        ss<< totalParticlesInDomain.at(cpt)<<"\t"<<sum_injectedParticles.at(cpt)<<"\t"<<totalParticlesInTerrain.at(cpt)<<"\t"<<sum_outsideParticles.at(cpt);
        return ss.str();
    }

    void clearParticlesvectors(){
        injectedParticles.clear();
        deposedPrticules.clear();
        outsideparticles.clear();

        totalParticlesInDomain.clear();
        sum_injectedParticles.clear();
        totalParticlesInTerrain.clear();
        sum_outsideParticles.clear();
    }


    vector<IterationTimeStamp<T> > stampTimeIteraions;
    T startComputationTimeStamp;
    T endComputationTimeStamp;
    bool isStarted;
    bool isStopped;
    bool isFileCreated;

    T computeTime;
    T copyTime;
    T preProcessTime;
    T processTime;
    T gatherTime;
    T sendTime;
    T receiveTime;
    T sendrecTime;
    T updateTime;
    std::string filename;

    T nbr_itr;
    int rank;
    string prefix;
    vector<int> injectedParticles, deposedPrticules , outsideparticles;
    vector<int> totalParticlesInDomain, sum_injectedParticles, totalParticlesInTerrain, sum_outsideParticles;

    //totalParticlesInDomain= Sum(injectedParticles) - (totalParticlesInTerrain+sum(particules_sorties) )
};

//***************************************************************
template<typename T>
class Profiler{
public:
    Profiler(){}
    virtual ~Profiler(){}

    void setPrefixFileName(string prefix){
        record.prefix=prefix;
    }

    void setRank(int rank){
        record.setRank(rank);
    }
    void start(){
        record.start();
    }

    void end(bool considerParticles=true){
        record.end(considerParticles);
    }

    void createStampRecording(){
        this->iterationTimeStamp.reset(new IterationTimeStamp<T>());
        // not needed
        this->iterationTimeStamp->registerBeginCopyOverlapStamp(0);
        this->iterationTimeStamp->registerEndCopyOverlapStamp(0);
        this->iterationTimeStamp->registerBeginPreProcessLevelStamp(0);
        this->iterationTimeStamp->registerEndPreProcessLevelStamp(0);


    }
    void commitStampRecording(){
        this->record.stampTimeIteraions.push_back(*iterationTimeStamp.get());
    }

    // begin ---

    void registerBeginProcessLevelStamp(){
        this->iterationTimeStamp->registerBeginProcessLevelStamp(0);
    }
    void registerBeginGatherLevelStamp(){
        this->iterationTimeStamp->registerBeginGatherLevelStamp(0);
    }
    void registerBeginUpdateLevelStamp(){
        this->iterationTimeStamp->registerBeginUpdateLevelStamp(0);
    }

    // send-rec
    void registerSendReceiveLevelStamp(T WTimeSendBegin, T WTimeSendEnd, T WTimeReceiveBegin, T WTimeReceiveEnd){
        this->iterationTimeStamp->registerSendReceiveLevelStamp(WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd, 0);
    }

    // end --
    void registerEndProcessLevelStamp(){
        this->iterationTimeStamp->registerEndProcessLevelStamp(0);
    }
    void registerEndGatherLevelStamp(){
        this->iterationTimeStamp->registerEndGatherLevelStamp(0);
    }
    void registerEndUpdateLevelStamp(){
        this->iterationTimeStamp->registerEndUpdateLevelStamp(0);
    }

    // particles
    void registerParticles (int injected, int deposed, int outside,
                            int totalParticlesInDomain, int sum_injectedParticles,
                            int totalParticlesInTerrain, int sum_outsideParticles){
        record.registerInjectedparticles(injected);
        record.registerDeposedparticles(deposed);
        record.registerOutSideparticles(outside);

        record.registerTotalParticlesInDomain (totalParticlesInDomain);
        record.registerSumInjectedParticles (sum_injectedParticles);
        record.registerTotalParticlesInTerrain (totalParticlesInTerrain);
        record.registerSumOutsideParticles (sum_outsideParticles);
    }

    void saveToFile(bool considerParticles=true){
        record.saveToFile(considerParticles);
    }

    string getAveragedStats(){
        return record.getAveragedStats();
    }


private:
    void registerBeginCopyOverlapStamp(){
        this->iterationTimeStamp->registerBeginCopyOverlapStamp(0);
    }
    void registerBeginPreProcessLevelStamp(){
        this->iterationTimeStamp->registerBeginPreProcessLevelStamp(0);
    }
    void registerEndPreProcessLevelStamp(){
        this->iterationTimeStamp->registerEndPreProcessLevelStamp(0);
    }
    void registerEndCopyOverlapStamp(){
        this->iterationTimeStamp->registerEndCopyOverlapStamp(0);
    }

    TimeStampMuscle<T> record;
    std::unique_ptr<IterationTimeStamp<T> > iterationTimeStamp;

};

//
}

#endif
