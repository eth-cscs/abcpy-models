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


class BroadCastSample: public MpiSubmodel{
private:
    std::stringstream ss;
    ConduitExit<int> *exit;
    ConduitEntrance<int> *entrance;
    int iterationCounter;

public:
    virtual ~BroadCastSample(){}
    virtual void F_init(){
        iterationCounter =0;
        ss<<"[id: "<<this->getName()<< "] rank="<< this->getMpiManager()->getRank()<<" -> Received: ";
        entrance=this->getConduitEntrance<int>("f_out");
        assert(entrance);
        exit=this->getConduitExit<int>("f_in");
        assert(exit);
    }

    virtual bool isConverged(){
        return (this->iterationCounter >= 10 );
    }

    virtual void Oi(){
        // nothing
    }

    virtual void S(){
        this->iterationCounter++;
    }

    virtual void U(){
        int root= (this->getMpiManager()->isMainProcessor())? MPI_ROOT: MPI_PROC_NULL;

        int number=100;
        if(id=="S2") number=200;
        int data = number + this->getMpiManager()->getRank();
        entrance->bCast(&data, 1, root);
        int result;
        exit->bCast(&result, 1, 0);
        ss<<result<<" ";
    }

    virtual void Of(){
        cout<<"------------------"<<endl;
        cout<<ss.str()<<endl;
        cout<<"------------------"<<endl;
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
        MpiSubmodel *k1 = new BroadCastSample();
        MpiSubmodel *k2= new BroadCastSample();

        //register kernels
        mmsf->addKernel(k1, "S1");
        mmsf->addKernel(k2, "S2");

        // start simulation
        if (! mmsf->loadCxACoupling()){
            cerr<<"Issue in loading coupling file"<<endl;
             return -1;
        }
        mmsf->compute();

        //free kernels then the mmsf
        delete k1;
        delete k2;

        // do it as the last statement
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
