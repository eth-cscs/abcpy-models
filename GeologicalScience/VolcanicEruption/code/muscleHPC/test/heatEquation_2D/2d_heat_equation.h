/*
 * https://people.sc.fsu.edu/~jburkardt/c_src/heat_mpi/heat_mpi.html
http://www.nic.funet.fi/~magi/opinnot/mpi/
https://people.sc.fsu.edu/~jburkardt/c_src/heated_plate/heated_plate.c
https://people.sc.fsu.edu/~jburkardt/c_src/heat_mpi/heat_mpi.c
http://www.cas.usf.edu/~cconnor/parallel/2dheat/2dheat.html
http://www.cas.usf.edu/~cconnor/parallel/2dheat/2d_heat_equation.c

*/
/****************************************************************************
 * HEAT2D Example - Parallelized C Version
 * FILE: mpi_heat2D.c
 * DESCRIPTIONS: This example is based on a simplified two-dimensional heat
        * equation domain decomposition. The initial temperature is computed to be
        * high in the middle of the domain and zero at the boundaries. The
        * boundaries are held at zero throughout the simulation. During the
        * time-stepping, an array containing two domains is used; these domains
        * alternate between old data and new data.
        *
        * In this parallelized version, the grid is decomposed by the master
        * process and then distributed by rows to the worker processes. At each
        * time step, worker processes must exchange border data with neighbors,
        * because a grid point's current temperature depends upon it's previous
        * time step value plus the values of the neighboring grid points. Upon
        * completion of all time steps, the worker processes return results
        * to the master process.
        *
        * AUTHOR: Blaise Barney - adapted from D. Turner's serial version
        * CONVERTED TO MPI: George L. Gusciora (1/25/95)
        * MODIFIED BY: C. B. Connor (6/6/02)
        *
        *****************************************************************************
        * MODIFIED BY: Mohamed Ben Belgacem (13/09/2016)
        ****************************************************************************/
#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H


#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <algorithm>
#include <functional>
#include "musclehpc/conduits/mmsf.h"

#include "grid_to_bmp.hpp"

#define TIME_STEPS 500 /* number of time steps */
#define MAXWORKER 8 /* maximum number of worker tasks */
#define MINWORKER 3 /* minimum number of worker tasks */
#define BEGIN 1 /* message type */
#define NGHBOR1 2 /* message type */
#define NGHBOR2 3 /* message type */
#define NONE 0 /* indicates no neighbor */
#define DONE 4 /* message type */
#define MASTER 0 /* taskid of first process */


using namespace std;

struct Parms {
    float cx;
    float cy;
} diffusivity = {0.1, 0.1};


struct Boundary{
   // int isBoundary; //0 : is not boundary, 1: is boundary
    int boundaryPosition; // can be "0":none, "1":begin or "2":end;
    int holderRank;
    int remoteRank;
};

class Heat{

protected:
    //float u[2][NXPROB][NYPROB]; /* array for grid */
    //vector<vector<vector<float>>> u;

    int NXPROB; /* x dimension of problem grid */
    int NYPROB; /* y dimension of problem grid */
    float *** u;
    int taskid; /* this task's unique id */
    int   numworkers; /* number of worker processes */
    int   numtasks; /* number of tasks */
    int   min_number_rows,number_rows,offset,extra_rows;/* for sending rows of data */
    int   destination, source; /* to - from for message send-receive */
    int   worker_number, neighbor1,neighbor2; /* neighbor tasks */
    int   message_tag; /* for message types */
    int   nbytes; /* number of bytes received */
    int   rc,start,end; /* misc */
    int   i,ix,iy,iz,it; /* loop variables */
    Boundary  boundaryInfo;
    MPI_Status status;
    shared_ptr<MpiManager> mpiManager;
    string initial_filename, final_filename;

    void allocate();
    void saveGridInFile(string fileName);
public:

    void prepare(int argc, char** argv);

    void addHeatMpiManager(shared_ptr<MpiManager> mpiManager);

    virtual ~Heat();

    /**************************************************************************
     * subroutine update
     ****************************************************************************/
    //void solve(int start, int end, int ny, float *u1, float *u2);
    void solve();
    virtual void boundary();

    /*****************************************************************************
     * subroutine inidat
     *****************************************************************************/
    void inidat(int nx, int ny, float *** u);

    /**************************************************************************
     * subroutine prtdat
     **************************************************************************/
    void prtdat(int nx, int ny, float *** u, const char *fnam);

    void init();
    void simulate();
    void reduceOnMaster();

};


class HeatSubmodel: public Heat, public MpiSubmodel{

private:
    ConduitExit<char> *f_in;
    ConduitEntrance<char> *f_out;
    int iterationCounter;
    MPI_Request listISendReq[2], listIRecvReq[2];
    MPI_Status  status [2];
    std::stringstream ss;
    int row_boundary_send, row_boundary_recv;

public:
    HeatSubmodel(): Heat(), row_boundary_send(-2), row_boundary_recv(-2) {
    }

    virtual ~HeatSubmodel(){
    }

    virtual void F_init(){

        iterationCounter =0;
        ss<<"[id: "<<this->getName()<< "] rank="<< this->getMpiManager()->getRank()<<" -> Received: ";
        f_out=this->getConduitEntrance<char>("f_out");
        assert(f_out);
        f_in=this->getConduitExit<char>("f_in");
        assert(f_in);

        this->addHeatMpiManager(this->getMpiManager());
        this->prepare(this->getKernelArgc(), this->getKernelArgv());
        initial_filename = "initial_"+this->getName()+".dat";
        final_filename = "final_"+this->getName()+".dat";
        this->init();

        if (taskid != MASTER && boundaryInfo.holderRank == taskid){
            assert (boundaryInfo.boundaryPosition == 1 || boundaryInfo.boundaryPosition == 2);

            if (boundaryInfo.boundaryPosition == 1){
                row_boundary_send = start;
                row_boundary_recv = offset;
            }else { // (boundaryPosition == 2)
                row_boundary_send = end;
                row_boundary_recv = end+1;
            }
        }

        int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
        f_out->bCast((char*)&boundaryInfo.holderRank, sizeof(boundaryInfo.holderRank), root);//send
        f_in->bCast((char*)&boundaryInfo.remoteRank, sizeof(boundaryInfo.remoteRank), 0);//recv
    }

    virtual bool isConverged(){
        return (this->iterationCounter >= TIME_STEPS );
    }

    virtual void Oi(){
        // nothing
        /*stringstream out;
        out <<"oi_"<<this->getName()<<"_"<<iterationCounter<<".dat";
        this->reduceOnMaster();
        this->saveGridInFile(out.str());*/
    }

    virtual void S(){
        this->iterationCounter++;
        this->boundary();
        this->solve();
    }

    virtual void U(){

        if (taskid != MASTER && boundaryInfo.holderRank == taskid){
            int buf_len = NYPROB* sizeof(float);
            /// ----- non blocking sending of data -----
            f_out->iSend((char*) &buf_len, sizeof(buf_len), boundaryInfo.remoteRank, &listISendReq[0]);
            f_out->iSend((char*) &u[iz][row_boundary_send][0], buf_len, boundaryInfo.remoteRank, &listISendReq[1]);
            /// ----- non blocking receiving of data ----
            int bufRecv_len=-1;
            // receives data sizes
            f_in->iRecv((char*)&bufRecv_len, sizeof(bufRecv_len), boundaryInfo.remoteRank, &listIRecvReq[0]);
            MPI_Wait(&listIRecvReq[0], &status[0]);// wait until the iRecv of the length of the data to receive is done
            f_in->iRecv((char*) &u[iz][row_boundary_recv][0], bufRecv_len, boundaryInfo.remoteRank, &listIRecvReq[1]);
            MPI_Wait(&listIRecvReq[1], &status[1]);

            ///---- wait isend to finish ----
            MPI_Waitall(2, &listISendReq[0], &status[0]);
        }

    }

    virtual void Of(){
        this->reduceOnMaster();
        this->saveGridInFile(this->final_filename);
        //cout<<"------------------"<<endl;
        //cout<<ss.str()<<endl;
        //cout<<"------------------"<<endl;
    }


};

class Program{
private:
    int argc; char **argv;

public:
    Program( int argc, char **argv ): argc(argc), argv(argv) {
    }


    int run() {

        // do it as first statement
        MMSF_MPI * mmsf = new MMSF_MPI(argc, argv);

        // Create Kernels' Instances
        MpiSubmodel * s1 = new HeatSubmodel();
        MpiSubmodel * s2 = new HeatSubmodel();

        //register kernels
        mmsf->addKernel(s1, "S1");
        mmsf->addKernel(s2, "S2");

        if (! mmsf->loadCxACoupling()){
            cerr<<"Issue in loading coupling file"<<endl;
            return -1;
        }

        // start simulation
        mmsf->compute();

        //free kernels then the mmsf
        delete s1; delete s2;
        //delete currentKernel1;

        // do it as  last statement
        delete mmsf;
        return 0;

    }
};

/*============== MAIN ==============*/
int main(int argc, char **argv ) {

    Program prg(argc, argv);
    prg.run();
    return 0;
}


/*============== MAIN ==============*/
/*int main(int argc, char **argv ) {
    if (argc < 3){
        cerr<<"usage:"<<argv[0]<<" NXPROB NYPROB"<<endl;
        return -1;
    }

    shared_ptr<MpiManager> mpiManager= std::make_shared<MpiManager>();
    mpiManager->init(&argc, &argv);

    Heat prg;
    prg.addHeatMpiManager(mpiManager);
    prg.prepare(argc, argv);
    prg.simulate();
    return 0;
}*/

#endif
