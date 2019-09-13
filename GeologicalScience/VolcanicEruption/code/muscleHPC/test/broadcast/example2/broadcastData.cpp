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

class   BroadCastData: public MpiSubmodel{
private:
    std::stringstream ss;
    ConduitExit<char> *exit;
    ConduitEntrance<char> *entrance;
    int iterationCounter;
    long totalBroadCastedSize;

public:
    virtual ~BroadCastData(){}

    virtual void F_init(){
        iterationCounter =0;
        totalBroadCastedSize=0;
        ss<<"[id: "<<this->getName()<< "] rank="<< this->getMpiManager()->getRank()<<" -> ";
        entrance=this->getConduitEntrance<char>("f_out");
        assert(entrance);
        exit=this->getConduitExit<char>("f_in");
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
        random_bytes_engine rbe;
        const int buff_size = 500;
        std::vector< char> data(buff_size);
        std::generate(begin(data), end(data), std::ref(rbe));
        std::vector< char> recvBuf(buff_size); //5 MB
        entrance->bCast(data.data(), data.size(), root);
        exit->bCast(recvBuf.data(), recvBuf.size() );
        totalBroadCastedSize+=buff_size;
    }

    virtual void Of(){
        cout<<"------------------"<<endl;
        ss<< "Total BroadCasted data size: "<< totalBroadCastedSize <<endl;
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
        MpiSubmodel *k1 = new BroadCastData();
        MpiSubmodel *k2= new BroadCastData();


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
