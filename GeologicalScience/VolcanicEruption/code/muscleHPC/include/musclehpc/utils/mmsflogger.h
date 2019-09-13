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

#ifndef MMSFLOGGER_H
#define MMSFLOGGER_H

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
#include <algorithm>
#include "musclehpc/parallelism/mpiManager.h"

using namespace std;


/*class global{

public:
    static global& getInstance() {
        static global    instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }

    bool isMainProcessor(){return true;}
private:
    global();
};*/


/// A buffer which reads and writes nothing.
/** This buffer is designed for use in parallel programming.
 *  It is assigned to processors which have potentially no access
 *  to the file system.
 */
class DevNullBuffer : public std::streambuf {
protected:
    virtual int_type overflow(int_type c) {
        return EOF;
    }
    virtual int_type underflow() {
        return EOF;
    }
};

// Parallel Output Streams.

struct Parallel_ostream {
    virtual ~Parallel_ostream() { }
    virtual std::ostream& getOriginalStream() =0;
};


template<typename Value>
Parallel_ostream& operator<< (Parallel_ostream& lhs, Value const& rhs) {
    lhs.getOriginalStream() << rhs;
    return lhs;
}

inline Parallel_ostream& operator<< (Parallel_ostream& lhs, std::ostream& (*op)(std::ostream&)) {
    lhs.getOriginalStream() << op;
    return lhs;
}


class plb_ofstream : public Parallel_ostream {

public:
    plb_ofstream(int rank);
    explicit plb_ofstream(int rank, const char* filename,
                          std::ostream::openmode mode = std::ostream::out | std::ostream::trunc);
    ~plb_ofstream();
    virtual std::ostream& getOriginalStream();
    virtual plb_ofstream& getWriter();

    bool is_open();
    void open(const char* filename, std::ostream::openmode mode = std::ostream::out | std::ostream::trunc);
    void close();
    void flush();
private:
    plb_ofstream(plb_ofstream const& rhs);
    plb_ofstream& operator=(plb_ofstream const& rhs);
private:
    int rank;
    DevNullBuffer devNullBuffer;
    std::ostream  devNullStream;
    std::ofstream *original;
};

/*class Parallel_referring_ostream : public Parallel_ostream {
public:
    Parallel_referring_ostream(std::ostream& original_ostream_)
        : devNullStream(&devNullBuffer),
          original_ostream(original_ostream_)
    { }
    virtual std::ostream& getOriginalStream() {
        if () {
            return original_ostream;
        }
        else {
            return devNullStream;
        }
    }
private:
    DevNullBuffer devNullBuffer;
    std::ostream  devNullStream;
    std::ostream& original_ostream;
};*/

//---------------------------------------



/*extern Parallel_referring_ostream pcout;
extern Parallel_referring_ostream pcerr;
extern Parallel_referring_ostream pclog;*/


//------------------------------------------
/*Parallel_referring_ostream pcout(std::cout);
Parallel_referring_ostream pcerr(std::cerr);
Parallel_referring_ostream pclog(std::clog);*/
//------------------------------------------
class PlumeLogger
{
public:
    plb_ofstream & pcout;
private:
    string logFileName;
    int rank;
    bool isFileCreated;
    double startComputationTimeStamp;
    ofstream * ofile;
    bool isDebug;


public:
    PlumeLogger(string logFileName, int rank, plb_ofstream & pcout, double startComputationTimeStamp =0 )
        : pcout(pcout), logFileName(logFileName), rank(rank),
          isFileCreated(false), startComputationTimeStamp(startComputationTimeStamp), isDebug(true){
        this->createLogFile();
        // pcout = new Parallel_referring_ostream(std::cout);
    }

    /*PlumeLogger(PlumeLogger const& logger) {
        this->logFileName=logger.logFileName;
        this->rank=logger.rank;
        this->isFileCreated=logger.isFileCreated;
        this->startComputationTimeStamp=logger.startComputationTimeStamp;
        this->isDebug=logger.isDebug;
        this->pcout=logger.pcout;
        this->ofile = new ofstream(logFileName.c_str(), std::ostream::app);//append mode

    }
    void operator=(PlumeLogger const& logger){
        this->logFileName=logger.logFileName;
        this->rank=logger.rank;
        this->isFileCreated=logger.isFileCreated;
        this->startComputationTimeStamp=logger.startComputationTimeStamp;
        this->isDebug=logger.isDebug;
        this->pcout=logger.pcout;
        this->ofile = new ofstream(logFileName.c_str(), std::ostream::app);//append mode
    }*/

    virtual ~PlumeLogger(){
        if(ofile){
            delete ofile;
            ofile=0;
        }
    }

    /*ofstream & getWriter(){
        return *ofile;
    }*/

    plb_ofstream & getWriter(){
        return pcout;
    }

    void createLogFile(){
        if(!isDebug) return;
        if(!isFileCreated && rank==0){
            std::stringstream fname;
            fname<< logFileName<<"_"<< std::fixed<< std::setprecision (8)<<startComputationTimeStamp<<".log";
            logFileName=fname.str();
            ofile = new ofstream(logFileName.c_str(), std::ostream::app);//append mode
            isFileCreated=true;
        }
    }

    void write(string str){
        if(!isDebug) return;
        if(rank==0){
            this->createLogFile();
            (*ofile) << str;
        }
    }

    void close(){
        if(!isDebug) return;
        if(rank==0 && ofile->is_open()){
            ofile->close();
        }
    }

};

#endif // MMSFLOGGER_H
