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

#define ITR (10)

/************************************ Gate *********************************************/
class Gate: public MpiSubmodel{

private:
    int iterationCounter;
    ConduitEntrance<char> *f1_out, *f2_out;
    ConduitExit<int> *f1_in, *f2_in;
    std::stringstream ss;

public:
    Gate(){
        iterationCounter=0;
    }
    virtual ~Gate(){}
    virtual void F_init(){

        cout<< this->getName() <<" rank "<< this->getMpiManager()->getRank() <<" is started "<<endl;
        if(this->getMpiManager()->isMainProcessor()){
            f1_out=this->getConduitEntrance<char>("f1_out");
            f2_out=this->getConduitEntrance<char>("f2_out");
            f1_in=this->getConduitExit<int>("f1_in");
            f2_in=this->getConduitExit<int>("f2_in");
            ss<<"[Name: "<<this->getName()<<"("<< this->getMpiManager()->getRank()<<") ] -> Received: ";
        }
    }
    virtual void Oi(){}

    virtual void S(){ //increment
        iterationCounter++;
    }

    virtual bool isConverged(){
        return !(this->iterationCounter<= ITR );
    }

    virtual void U(){

            int e1,e2;
            f2_in->bCast(&e2, 1, 0); // bcast recv an int
            f1_in->bCast(&e1, 1, 0); // b cast recv an int
            ss<<"(f1in:"<<e1<<", f2in:"<<e2<< ") ";

            std::stringstream ss_str1,ss_str2;
            ss_str1<<"G->f1in_"<<iterationCounter<<" ";
            ss_str2<<"G->f2in_"<<iterationCounter<<" ";
            string msg1= ss_str1.str();
            string msg2= ss_str2.str();
            const char * buffer1 = msg1.c_str();
            const char * buffer2 = msg2.c_str();
            int buffer1_len,buffer2_len;
            buffer1_len=(int) msg1.size()+1;
            buffer2_len=(int) msg2.size()+1;

            int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;

            f1_out->bCast((char*) &buffer1_len, sizeof(buffer1_len), root);// bcast a char* length
            f1_out->bCast((char*) buffer1, buffer1_len, root);// bcast a char*

            f2_out->bCast((char*)&buffer2_len, sizeof(buffer2_len), root); // bcast a char* length
            f2_out->bCast((char*)buffer2, buffer2_len, root); // bcast a char*
    }

    virtual void Of(){
        cout<<"------------------"<<endl;
        cout<<ss.str()<<endl;
        cout<<"------------------"<<endl;
    }

};


/************************************ SWD2 *********************************************/
class SWD2: public MpiSubmodel{
private:
    int iterationCounter;
    ConduitEntrance<int> *f_out;
    ConduitExit<char> *f_in;
    std::stringstream ss;

public:
    virtual ~SWD2(){}

    virtual void F_init(){

        iterationCounter=0;
        ss<< this->getName() <<" rank "<< this->getMpiManager()->getRank() <<" is started "<<endl;
        f_out=this->getConduitEntrance<int>("f_out");
        f_in=this->getConduitExit<char>("f_in");
        ss<<"[Name: "<<this->getName()<<"("<< this->getMpiManager()->getRank()<<") ] -> Received: ";

    }
    virtual void Oi(){}

    virtual void S(){ //increment
        iterationCounter++;
    }

    virtual bool isConverged(){
        return !(this->iterationCounter<= ITR );
    }


    virtual void U(){
            int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;
            f_out->bCast(&iterationCounter, 1, root); // bcast send an int

            int len;
            f_in->bCast((char*) &len, sizeof(len), 0);// bcast recv an int
            char * buffer = new char[len];
            f_in->bCast((char*) buffer, len, 0);// bcast recv a char*
            std::string msg (buffer, len);
            delete buffer;
            ss<<msg<<", ";
    }

    virtual void Of(){
        cout<<"------------------"<<endl;
        cout<<ss.str()<<endl;
        cout<<"------------------"<<endl;
    }
};

/************************************ Main *********************************************/

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
        MpiSubmodel * s1 = new SWD2();
        MpiSubmodel * s2 = new SWD2();
        MpiSubmodel * gate = new Gate();

        //register kernels
        mmsf->addKernel(s1, "S1");
        mmsf->addKernel(s2, "S2");
        mmsf->addKernel(gate, "G");

        if (! mmsf->loadCxACoupling()){
            cerr<<"Issue in loading coupling file"<<endl;
            return -1;
        }

        // start simulation
        mmsf->compute(true);

        //free kernels then the mmsf
        delete s1; delete s2; delete gate;
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
/*
conduit1<int>: S1.f_out -> G.f1_in
conduit2<int>: S2.f_out -> G.f2_in
conduit3<char>: G.f1_out -> S1.f_in
conduit4<char>: G.f2_out -> S2.f_in
cores<S1>:4
cores<S2>:4
cores<G>:1
cmdline<S1>:
cmdline<S2>:
cmdline<G>:
*/
