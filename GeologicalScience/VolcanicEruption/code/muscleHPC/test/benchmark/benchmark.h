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

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <memory>
#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/conduits/mmsf.h"


using namespace std;
using namespace unige_pasc;

class BenchMark
{

protected:
    int myid;
    MPI_Status reqstat;
    char *s_buf, *r_buf;
    double t_start , t_end;


    int MAX_MSG_SIZE, LARGE_MESSAGE_SIZE;
    int MYBUFSIZE;
    int loop, loop_large;
    int skip, skip_large;




public:
    BenchMark()
    {

        t_start = 0.0;
        t_end = 0.0;
        MAX_MSG_SIZE= 1<<26;
        LARGE_MESSAGE_SIZE=8192;
        MYBUFSIZE=MAX_MSG_SIZE;
        this->loop = 1000;
        this->skip = 10;
        this->loop_large = 100;
        this->skip_large = 10;
    }

    virtual ~BenchMark(){}


    void touch_data (void * sbuf, void * rbuf, int rank, size_t size)
    {
        if (0 == rank){
            memset(sbuf, 'a', size);
            memset(rbuf, 'b', size);
        }
    }


    int allocate_memory (char ** sbuf, char ** rbuf, int rank)
    {
        unsigned long align_size = sysconf(_SC_PAGESIZE);
        int res=posix_memalign((void**)sbuf, align_size, MYBUFSIZE);

        if (res) {
            cerr<<"EINVAL:"<<EINVAL<< "  ENOMEM:"<<ENOMEM<<endl;
            cerr<< "Error allocating host memory: "<< res<<endl;
            return 1;
        }
        res=posix_memalign((void**)rbuf, align_size, MYBUFSIZE);
        if (res) {
            cerr<<"EINVAL:"<<EINVAL<< "  ENOMEM:"<<ENOMEM<<endl;
            cerr<< "Error allocating host memory: "<< res<<endl;
            return 1;
        }
        return res;
    }


    virtual void start()=0;
    virtual void print_header(int rank)=0;

    void printHostname(){
        char hostname[1024];
        hostname[1023] = '\0';
        gethostname(hostname, 1023);
        printf("Hostname: %s\n", hostname);
        struct hostent* h;
        h = gethostbyname(hostname);
        printf("h_name: %s\n", h->h_name);
    }
};

/*********************************** PingPong *****************************************/
class PingPong: public BenchMark
{

private:
    shared_ptr <MpiManager> mpiManager;
public:
    PingPong(shared_ptr <MpiManager> mpiManager);
    virtual ~PingPong();
    virtual void start();
    virtual void print_header(int rank);
};

/*********************************** Ping *****************************************/
class Ping: public BenchMark, public RemoteMpiKernel
{
private:
    //shared_ptr <MpiManager> mpiManager;
   /////string id;
public:
    Ping();
    virtual ~Ping();

   ///// virtual void initialize(string id){this->id=id;}
   ///// virtual string getName(){return this->id;}
    virtual void mainLoop();

    virtual void start();
    virtual void print_header(int rank);
};

/*********************************** Pong *****************************************/
class Pong: public BenchMark,  public RemoteMpiKernel
{
private:
    string id;
    //shared_ptr <MpiManager> mpiManager;
public:
    Pong();
    virtual ~Pong();

   ///// virtual void initialize(string id){this->id=id;}
   ///// virtual string getName(){return this->id;}
    virtual void mainLoop();

    virtual void start();
     virtual void print_header(int rank);
};

#endif // BENCHMARK_H
