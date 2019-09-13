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

#define ITR (10)

class MPIRing_Asynchrone: public MpiSubmodel{
private:

    int iterationCounter;
    ConduitEntrance<char> *f1_out, *f2_out;
    ConduitExit<char> *f1_in, *f2_in;
    std::stringstream ss;
    MPI_Request  listISendReq [2];
    MPI_Request  listIRecvReq [2];

public:
    MPIRing_Asynchrone():MpiSubmodel() {
        iterationCounter=0;
    }
    virtual ~MPIRing_Asynchrone(){}

    virtual void F_init(){

        cout<< this->getName() <<" rank "<< this->getMpiManager()->getRank() <<" is started "<<endl;
        f1_out=this->getConduitEntrance<char>("f1_out");
        f2_out=this->getConduitEntrance<char>("f2_out");
        f1_in=this->getConduitExit<char>("f1_in");
        f2_in=this->getConduitExit<char>("f2_in");
        ss<<"[Name: "<<this->getName()<<"("<< this->getMpiManager()->getRank()<<") ] -> Received: ";
    }
    virtual void Oi(){ /*nothing*/}
    virtual void Of(){
        cout<<ss.str()<<endl;
    }
    virtual void S(){
        iterationCounter++;
    }

    virtual bool isConverged(){
        return !(this->iterationCounter < 10 );
    }

    virtual void U(){

        string message1, message2;
        int buf_len_1,buf_len_2;

        // prepare messages
        std::stringstream str_buff1,str_buff2;
        str_buff1<<this->getMpiManager()->getRank()<<"_f1_in_"<<iterationCounter<<" ";
        str_buff2<<this->getMpiManager()->getRank()<<"_f2_in_"<<iterationCounter<<" ";
        message1= str_buff1.str();
        message2= str_buff2.str();
        const char * buffer1 = message1.c_str();
        const char * buffer2 = message2.c_str();
        buf_len_1=(int) message1.size()+1;
        buf_len_2=(int) message2.size()+1;

        // here each rank will send/receive from the same remote rank
        int dest = this->getMpiManager()->getRank();

        MPI_Status  status [2];
        /// ----- non blocking sending of data -----
        f2_out->iSend((char*)&buf_len_2, sizeof(buf_len_2), dest, &listISendReq[0]);
        f1_out->iSend((char*)&buf_len_1, sizeof(buf_len_1), dest, &listISendReq[1]);
        //MPI_Waitall(2, &listISendReq[0], &status[0]);
        f2_out->iSend( (char*)buffer2, buf_len_2, dest, &listISendReq[0]);
        f1_out->iSend((char*) buffer1, buf_len_1, dest, &listISendReq[1]);

        /// ----- non blocking receiving of data ----
        int bufRecv_len_1=-1, bufRecv_len_2=-1;
        // receives data sizes
        f1_in->iRecv((char*)&bufRecv_len_1, sizeof(bufRecv_len_1), dest, &listIRecvReq[0]);
        f2_in->iRecv((char*)&bufRecv_len_2, sizeof(bufRecv_len_1), dest, &listIRecvReq[1]);
        // wait until the iRecv of the length of the data to receive is done
        MPI_Waitall(2, &listIRecvReq[0], &status[0]);
        // then receives data
        char * bufferRecv1 = new char[bufRecv_len_1];
        char * bufferRecv2= new char[bufRecv_len_2];
        f1_in->iRecv (bufferRecv1, bufRecv_len_1, dest, &listIRecvReq[0]);
        f2_in->iRecv (bufferRecv2, bufRecv_len_2, dest, &listIRecvReq[1]);
        // wait until the iRecv of data is done
        MPI_Waitall(2, &listIRecvReq[0], &status[0]);

        /// ----- display the received data ------
        std::string myString(bufferRecv1, bufRecv_len_1);
        std::string myString1(bufferRecv2, bufRecv_len_2);
        ss<<"("<<myString<<","<< myString1<< ") ";
        delete bufferRecv2;
        delete bufferRecv1;
        /// ---- wait asynchrone isend to finish
        // wait until the iSend is done
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
        MpiSubmodel * s1 = new MPIRing_Asynchrone();
        MpiSubmodel * s2 = new MPIRing_Asynchrone();

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
