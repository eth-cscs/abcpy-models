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

#include "musclehpc/conduits/mmsf.h"

#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <algorithm>
#include <functional>

using namespace std;

using random_bytes_engine = std::independent_bits_engine<
std::default_random_engine, CHAR_BIT, unsigned char>;

class MPIRingData_Asynchrone: public MpiSubmodel{
private:

    int iterationCounter;
    // MML ports
    ConduitEntrance<char> *f_out;
    ConduitExit<char> *f_in;

    // kernels parameters
    int argc;
    char ** argv;

    MPI_Request  listISendReq [2];
    MPI_Request  listIRecvReq [2];

public:

    MPIRingData_Asynchrone():MpiSubmodel() {
        iterationCounter=0;
    }
    virtual ~MPIRingData_Asynchrone(){}

    virtual void F_init(){

        // read argc and argv parameters form the coupling file
        this->argc=this->getKernelArgc();
        this->argv=this->getKernelArgv();
        // init parallel conduits
        f_out=this->getConduitEntrance<char>("f_out");
        f_in=this->getConduitExit<char>("f_in");
        // print message
	stringstream ss;
        ss<< this->getName() <<" rank "<< this->getMpiManager()->getRank() <<" is started --> argv: ";
	for (int i = 0; i< this->argc; i++){
		ss<<this->getKernelArgv()[i] <<" ";
	}
	cout<<ss.str()<<endl;
    }
    virtual void Oi(){ /*put an observation*/}
    virtual void Of(){
        /*put final observation*/
    }
    virtual void S(){
        iterationCounter++;
    }

    virtual bool isConverged(){
        return !(this->iterationCounter < 10 );
    }

    virtual void U(){

        const int buff_size = 100;
        // prepare a vector of char with random data
        /*random_bytes_engine rbe;
        std::vector< char> vect(buff_size);
        std::generate(begin(vect), end(vect), std::ref(rbe));*/
	std::vector<int> vect(buff_size);
	for (int i = 0; i< buff_size ; i++){
		vect[i]=i+this->iterationCounter;
	}

        // Here each MPI process will send/receive from remote process remote having the same rank
        int dest = this->getMpiManager()->getRank();

        MPI_Status  status [2];
        /// ----- non blocking sending of data -----
        int buf_len = vect.size()*sizeof(int);
        f_out->iSend((char*)&buf_len, sizeof(buf_len), dest, &listISendReq[0]);
        f_out->iSend((char*) vect.data(), buf_len, dest, &listISendReq[1]);

        /// ----- non blocking receiving of data ----
        int bufRecv_len_1=0;
        // receives data sizes
        f_in->iRecv((char*)&bufRecv_len_1, sizeof(bufRecv_len_1), dest, &listIRecvReq[0]);
        // wait until the iRecv of the length of the data to receive is done
        MPI_Waitall(1, &listIRecvReq[0], &status[0]);
        // then receives data
        vector<int> vectRecv;
        vectRecv.resize(bufRecv_len_1/sizeof(int));
        f_in->iRecv ((char*) &vectRecv[0], bufRecv_len_1, dest, &listIRecvReq[0]);
        // wait until the iRecv of data is done
        MPI_Waitall(1, &listIRecvReq[0], &status[0]);

	// ----- verify received data -----
	for (int i = 0; i< buff_size ; i++){
                assert ((vectRecv[i] == ( i+this->iterationCounter))  && "error in assertion");
        }

        /// ----- free data ------
        vector<int>().swap(vect);
        vector<int>().swap(vectRecv);

        /// ---- wait asynchrone iSend to finish
        MPI_Waitall(2, &listISendReq[0], &status[0]);
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
        MpiSubmodel * s1 = new MPIRingData_Asynchrone();
        MpiSubmodel * s2 = new MPIRingData_Asynchrone();

        //register kernels
        mmsf->addKernel(s1, "S1");
        mmsf->addKernel(s2, "S2");

        if (! mmsf->loadCxACoupling()){
            cerr<<"Issue in loading coupling file"<<endl;
            return -1;
        }

        // start simulation
        mmsf->compute(true);

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
