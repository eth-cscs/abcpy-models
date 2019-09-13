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

#ifndef BENCHMARKMUSCLE_H
#define BENCHMARKMUSCLE_H

#include "benchmark.h"
#include <cppmuscle.hpp>


using namespace muscle;
class BenchMarkMuscle: public BenchMark
{

private:
    string kernelName;
    shared_ptr <MpiManager> mpiManager;
    bool isPing;
public:
    BenchMarkMuscle(int argc, char** argv): BenchMark(){

        shared_ptr<MpiManager> tmp (std::make_shared<MpiManager>()) ;//getManagerMPI();
        tmp->init(& argc, &argv);
        mpiManager=std::move(tmp);

        env::init(&argc, &argv);
        //kernelName= cxa::get_property("kernelName");
        if(mpiManager->getRank()==0){
            this->kernelName=muscle::cxa::kernel_name();
            cout<<"kernelName="<<kernelName<<endl;
            isPing=(kernelName == "Ping")? true: false;
        }
    }

    virtual ~BenchMarkMuscle(){}

    virtual void start(){

        if (this->allocate_memory(&s_buf, &r_buf, mpiManager->getRank())) {
            /* Error allocating memory */
            cerr<<"Error allocating memory"<<endl;
            return ;
        }

        print_header(mpiManager->getRank());

        /* Latency test */
        for(int size = 0; size <= MAX_MSG_SIZE; size = (size ? size * 2 : 1)) {

            try{
                touch_data(s_buf, r_buf, mpiManager->getRank(), size);
            }catch( const std::exception& e ) { // reference to the base of a polymorphic object
                std::cout << e.what(); // information from length_error printed
            }

            if(size > LARGE_MESSAGE_SIZE) {
                loop = loop_large;
                skip = skip_large;
            }

            mpiManager->barrier();

            if(mpiManager->getRank() == 0) {
                if(isPing){
                    for(int i = 0; i < loop + skip; i++) {
                        if(i == skip) t_start = MPI_Wtime();
                        env::send("fout", s_buf, (size_t) size, MUSCLE_RAW);
                        env::receive("fin", r_buf, (size_t &)size, MUSCLE_RAW);
                    }
                }else{
                    for(int i = 0; i < loop + skip; i++) {
                        if(i == skip) t_start = MPI_Wtime();
                        env::receive("fin", r_buf, (size_t &)size, MUSCLE_RAW);
                        env::send("fout", s_buf, (size_t) size, MUSCLE_RAW);
                    }
                }

                t_end = MPI_Wtime();

                double latency = (t_end - t_start) * 1e6 / (2.0 * loop);
                cout<<size<<"\t\t\t"<<latency<<endl;
            }
        }

        free(s_buf);
        free(r_buf);

    }
    virtual void print_header(int rank){
        if(rank==0){
            cout<<"# OSU MPI Latency Test :: Muscle Ping"<<endl;
            cout<<"#Size\t\t\tLatency (us)"<<endl;
        }
    }
};



/*********************************** main *****************************************/


int main( int argc, char **argv ) {

    BenchMarkMuscle b(argc, argv);
    b.start();
    return 0;

}


#endif // BENCHMARKMUSCLE_H
