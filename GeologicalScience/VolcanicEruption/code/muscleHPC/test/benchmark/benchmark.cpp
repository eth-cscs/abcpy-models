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

#include "benchmark.h"



/*********************************** BenchMark *****************************************/


/*********************************** PingPong *****************************************/


PingPong::PingPong(shared_ptr<MpiManager> mpiManager):BenchMark(){
    this->mpiManager=mpiManager;
}

PingPong::~PingPong(){}


void PingPong::print_header(int rank){
    if(rank==0){
        cout<<"# OSU MPI Latency Test :: PingPong"<<endl;
        cout<<"#Size\t\t\tLatency (us)"<<endl;
    }
}


void PingPong::start(){

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
            for(int i = 0; i < loop + skip; i++) {
                if(i == skip) t_start = MPI_Wtime();
                mpiManager->send<char>(s_buf, size, 1, 1);
                mpiManager->receive<char>(r_buf, size, 1, 1);
            }

            t_end = MPI_Wtime();
        }else{
            for(int i = 0; i < loop + skip; i++) {
                mpiManager->receive<char>(r_buf, size, 0, 1);
                mpiManager->send<char>(s_buf, size, 0, 1);
            }

        }
        if(mpiManager->getRank() == 0) {
            double latency = (t_end - t_start) * 1e6 / (2.0 * loop);
            cout<<size<<"\t\t\t"<<latency<<endl;

        }
    }

    free(s_buf);
    free(r_buf);
}


/*********************************** Ping *****************************************/

Ping::Ping():BenchMark(), RemoteMpiKernel(){
}
Ping::~Ping(){}


void Ping::print_header(int rank){
    if(rank==0){
        cout<<"# OSU MPI Latency Test :: Ping"<<endl;
        cout<<"#Size\t\t\tLatency (us)"<<endl;
    }
}



void Ping::start(){

}

void Ping::mainLoop(){

    ConduitEntrance<char> *fout=this->getConduitEntrance<char>("fout");
    ConduitExit<char> *fin=this->getConduitExit<char>("fin");

    MPI_PTP_Conduit<char> * ffout= dynamic_cast<MPI_PTP_Conduit<char> *> (fout);
    shared_ptr<MpiManager>  m_fout = ffout->getHandlerInerCommunicatortManager();
    assert(m_fout);
    int fout_tag=ffout->getConduitTag();
    MpiKernel * fout_handler= dynamic_cast<MpiKernel * > (ffout->getHandlerKernel());


    MPI_PTP_Conduit<char> * ffin= dynamic_cast<MPI_PTP_Conduit<char> *> (fin);
    shared_ptr<MpiManager> m_fin = ffin->getHandlerInerCommunicatortManager();
    assert(m_fin);
    int fin_tag=ffin->getConduitTag();
    MpiKernel * fin_handler= dynamic_cast<MpiKernel * > (ffin->getHandlerKernel());

    if (this->allocate_memory(&s_buf, &r_buf, this->getMpiManager()->getRank())) {
        /* Error allocating memory */
        cerr<<"Error allocating memory"<<endl;
        return ;
    }

    print_header(this->getMpiManager()->getRank());

    /* Latency test */
    for(int size = 0; size <= MAX_MSG_SIZE; size = (size ? size * 2 : 1)) {

        touch_data(s_buf, r_buf, this->getMpiManager()->getRank(), size);

        if(size > LARGE_MESSAGE_SIZE) {
            loop = loop_large;
            skip = skip_large;
        }

        this->getMpiManager()->barrier();

        if(this->getMpiManager()->getRank() == 0) {

            for(int i = 0; i < loop + skip; i++) {
                if(i == skip) t_start = MPI_Wtime();
                m_fout->send<char>(s_buf, size, fout_handler->getMpiManager()->getRank(), fout_tag );
                m_fin->receive<char>(r_buf, size, fin_handler->getMpiManager()->getRank(), fin_tag);
            }

            t_end = MPI_Wtime();

            double latency = (t_end - t_start) * 1e6 / (2.0 * loop);
            cout<<size<<"\t\t\t"<<latency<<endl;
        }
    }

    free(s_buf);
    free(r_buf);

}



/*********************************** Pong *****************************************/

Pong::Pong():BenchMark(), RemoteMpiKernel(){
   // this->mpiManager=mpiManager;
}

Pong::~Pong(){}

void Pong::print_header(int rank){
    if(rank==0){
        cout<<"# OSU MPI Latency Test :: Pong"<<endl;
        cout<<"#Size\t\t\tLatency (us)"<<endl;
    }
}

void Pong::start(){

}

void Pong::mainLoop(){

    ConduitEntrance<char> *fout=this->getConduitEntrance<char>("fout");
    ConduitExit<char> *fin=this->getConduitExit<char>("fin");

    MPI_PTP_Conduit<char> * ffout= dynamic_cast<MPI_PTP_Conduit<char> *> (fout);
    shared_ptr<MpiManager> m_fout = ffout->getHandlerInerCommunicatortManager();
    int fout_tag=ffout->getConduitTag();
    MpiKernel * fout_handler= dynamic_cast<MpiKernel * > (ffout->getHandlerKernel());


    MPI_PTP_Conduit<char> * ffin= dynamic_cast<MPI_PTP_Conduit<char> *> (fin);
    shared_ptr<MpiManager>  m_fin = ffin->getHandlerInerCommunicatortManager();
    int fin_tag=ffin->getConduitTag();
    MpiKernel * fin_handler= dynamic_cast<MpiKernel * > (ffin->getHandlerKernel());

    if (this->allocate_memory(&s_buf, &r_buf, this->getMpiManager()->getRank())) {
        /* Error allocating memory */
        cerr<<"Error allocating memory"<<endl;
        return ;
    }

    print_header(this->getMpiManager()->getRank());

    /* Latency test */
    for(int size = 0; size <= MAX_MSG_SIZE; size = (size ? size * 2 : 1)) {

        touch_data(s_buf, r_buf, this->getMpiManager()->getRank(), size);

        if(size > LARGE_MESSAGE_SIZE) {
            loop = loop_large;
            skip = skip_large;
        }

        this->getMpiManager()->barrier();

        if(this->getMpiManager()->getRank() == 0) {

            for(int i = 0; i < loop + skip; i++) {
                if(i == skip) t_start = MPI_Wtime();
                m_fin->receive<char>(r_buf, size, fin_handler->getMpiManager()->getRank(), fin_tag);
                m_fout->send<char>(s_buf, size, fout_handler->getMpiManager()->getRank(), fout_tag );
            }

            t_end = MPI_Wtime();

           // double latency = (t_end - t_start) * 1e6 / (2.0 * loop);
            //cout<<size<<"\t\t\t"<<latency<<endl;
        }
    }

    free(s_buf);
    free(r_buf);

}





/*********************************** main *****************************************/


int main( int argc, char **argv ) {

    if (argc < 2){
        cerr<<"Usage: "<<argv[0]<<" Ping|Pong"<<endl;
        return -1;
    }


    std::string first_arge;
    std::vector<std::string> all_args;

    if (argc > 1) {
        first_arge = argv[1];
        all_args.assign(argv + 1, argv + argc);
    }

    string type= all_args.at(0);



    if(type == "PingPong"){
        shared_ptr<MpiManager> tmp (std::make_shared<MpiManager>()) ;//getManagerMPI();
        tmp->init(& argc, &argv);
        PingPong b(std::move(tmp));
        b.start();
        return 0;
    }



    // do it as first statement
   MMSF_MPI * mmsf = new MMSF_MPI(argc, argv);

    // Create Kernels' Instances
    Kernel *k1 = new Ping();
    Kernel *k2= new Pong();


    //register kernels
    mmsf->addKernel(k1, "Ping", 2);
    mmsf->addKernel(k2, "Pong", 2);

    //Setup entrance connections
   // mmsf->connect<char>("Ping.fout")->To("Pong.fin")->with_MTM_MPI("conduit1");
   // mmsf->connect<char>("Pong.fout")->To("Ping.fin")->with_MTM_MPI("conduit2");

     mmsf->loadCxACoupling();

    // start simulation

    mmsf->compute();

    //free kernels then the mmsf
    delete k1;
    delete k2;

    // do it as  last statement
    delete mmsf;



}
