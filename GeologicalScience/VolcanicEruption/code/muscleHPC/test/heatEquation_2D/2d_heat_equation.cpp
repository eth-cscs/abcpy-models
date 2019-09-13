

#include "2d_heat_equation.h"



void Heat::addHeatMpiManager(shared_ptr<MpiManager> mpiManager){
    this->mpiManager = mpiManager;
}

Heat::~Heat(){
    // De-Allocate memory to prevent memory leak
    if (u){
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < NXPROB; ++j){
                delete [] u[i][j];
            }
            delete [] u[i];
        }
        delete [] u;
    }
}

void Heat::prepare(int argc, char **argv){

    if (mpiManager->getSize() < 2){
         throw std::invalid_argument( "nunber of cores should be >= 2" );
    }

    initial_filename = "initial.dat";
    final_filename = "final.dat";
    NXPROB = atoi(argv[1]);
    NYPROB = atoi(argv[2]);

    boundaryInfo.holderRank= -1;
    boundaryInfo.boundaryPosition=0;
    if (argc > 3){
        this->boundaryInfo.boundaryPosition = atoi(argv[3]);
    }

    /* First, find out my taskid and how many tasks are running */
    numtasks= mpiManager->getSize();
    taskid= mpiManager->getRank();
    /*rc = MPI_Init(&argc,&argv);
        rc |= MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
        rc |= MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        if (rc != 0)
            cout<<"error initializing MPI and obtaining task ID information"<<endl;*/

    /* one of the processors is completely reserved for coordination */
    numworkers = numtasks-1;
}


void Heat::inidat(int nx, int ny, float ***u) {
    int ix, iy;
    /*for (ix = 0; ix <= nx-1; ix++)
        for (iy = 0; iy <= ny-1; iy++)
            *(u+ix*ny+iy) = (float)(nx/2* iy * (ny - iy - 1))/(4*nx*ny);*/
    /*for (ix = 0; ix <= nx-1; ix++)
            for (iy = 0; iy <= ny-1; iy++)
                u[0][ix][iy] = (float)(nx/2* iy * (ny - iy - 1))/(4*nx*ny);*/

    for (ix = 0; ix <= nx-1; ix++)
        for (iy = 0; iy <= ny-1; iy++)
            u[0][ix][iy] = (float)((ny - iy - 1))/(4*ny)*100.0;
}

void Heat::prtdat(int nx, int ny, float ***u1, const char* fnam) {
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    fprintf ( fp, "%d\n", ny ); //height
    fprintf ( fp, "%d\n", nx ); //width

    for (iy = ny-1; iy >= 0; iy--) {
        for (ix = 0; ix <= nx-1; ix++) {
            //fprintf(fp, "%d %d %6.1f\n", iy, ix, *(u1+ix*ny+iy));
            fprintf ( fp, "%6.2f ", u1[0][ix][iy] );
        }
        fputc ( '\n', fp);
    }
    fclose(fp);
    string sfname = fnam;
    ImageWriter w(sfname, sfname+".bmp");
    w.write();
}

void Heat::init(){

    this->allocate();
    /* the MASTER node subdivides the problem by grid decomposition,
      in this case by rows. The MASTER also initializes the starting values, sets the
      boundary conditions on the problem , and receives results from the
      nodes after TIME_STEPS */
    if (taskid == MASTER)
    {
        /************************* master code *******************************/
        /* Check if numworkers is within range - quit if not */
        if ((numworkers > MAXWORKER) || (numworkers < MINWORKER))
        {
            printf("MP_PROCS needs to be between %d and %d for this exercise\n",
                   MINWORKER+1,MAXWORKER+1);
           // MPI_Finalize();
        }

        /* Initialize grid */
        printf("Grid size: X= %d Y= %d Time steps= %d\n",NXPROB,NYPROB,TIME_STEPS);
        printf("Initializing grid and writing initial.dat file...\n");
        inidat(NXPROB, NYPROB,  u);
        prtdat(NXPROB, NYPROB, u, this->initial_filename.c_str());

        /* Distribute work to workers. Must first figure out how many rows to
           send and what to do with extra rows. */
        min_number_rows = NXPROB/numworkers;
        extra_rows = NXPROB%numworkers;
        offset = 0;
        for (i=1; i<=numworkers; i++)
        {

            /* The following is a particularly compact and efficient  way of distributing the
               grid. It assures that the number of rows received by each node is only up to one more
               than any other node. It can be read as:
               if (i<=extra_rows) number_rows = min_number_rows + 1
               else number_rows = min_number_rows
            */
            number_rows = (i <= extra_rows) ? min_number_rows+1 : min_number_rows;
            /* Tell each worker who its neighbors are, since they must exchange
               data with each other. */
            if (i == 1){
                neighbor1 = NONE;
                if(boundaryInfo.boundaryPosition == 1)
                    boundaryInfo.holderRank = i;
            }else{
                neighbor1 = i - 1;
            }
            if (i == numworkers){
                neighbor2 = NONE;
                if(boundaryInfo.boundaryPosition == 2)
                    boundaryInfo.holderRank = i;
            }else{
                neighbor2 = i + 1;
            }
            /* Now send startup information to each worker Note that this
               information is "tagged" as the BEGIN message*/
            destination = i;
            worker_number = i;
            message_tag = BEGIN;

            /* Send the required information to each node */
            //    void send( T *buf, int count, int dest, int tag = 0 );
            mpiManager->send<int>(&worker_number, 1, destination, message_tag);
            mpiManager->send<int>(&offset, 1, destination, message_tag);
            mpiManager->send<int>(&number_rows, 1, destination, message_tag);
            mpiManager->send<int>(&neighbor1, 1, destination, message_tag);
            mpiManager->send<int>(&neighbor2, 1, destination, message_tag);
            for (int k=0; k<number_rows; k++ ){
                mpiManager->send<float>(&u[0][offset+k][0], NYPROB, destination, message_tag);
            }
            //mpiManager->send(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, destination, message_tag, MPI_COMM_WORLD);

            /*let the world know how the problem has been divided up */
            printf("Sent to= %d offset= %d number_rows= %d neighbor1= %d neighbor2=%d\n",
                   destination,offset,number_rows,neighbor1,neighbor2);

            /* increment the offset by the number_rows so the next node will
               know where its grid begins */
            offset += number_rows;

        } /* continue doing the above for each node, i */
    } /* End init of master code */

    if (taskid != MASTER)
    {
        /************************* worker code**********************************/
        /* Initialize everything - including the borders - to zero */
        /* iz is a flag indicating one of two grids used in the analysis */
        for (iz=0; iz<2; iz++)
            for (ix=0; ix<NXPROB; ix++)
                for (iy=0; iy<NYPROB; iy++)
                    u[iz][ix][iy] = 0.0;

        /* Now receive my offset, rows, neighbors and grid partition from master
             */
        source = MASTER;
        message_tag = BEGIN;
        mpiManager->receive<int>(&worker_number, 1,  source, message_tag);
        mpiManager->receive<int>(&offset, 1,  source, message_tag);
        mpiManager->receive<int>(&number_rows, 1,  source, message_tag);
        mpiManager->receive<int>(&neighbor1, 1,  source, message_tag);
        mpiManager->receive<int>(&neighbor2, 1,  source, message_tag);
        for (int k=0; k<number_rows; k++ ){
            mpiManager->receive<float>(&u[0][offset+k][0], NYPROB,  source, message_tag);
        }

        //mpiManager->receive(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);

        /* Determine border elements. This takes into account that the
      first and last rows have fixed temperature*/

        if (offset==0)
            start=1; /* do not include row zero */
        else
            start=offset;
        if ((offset+number_rows)==NXPROB)
            end= offset + number_rows-2;  /*do not include the last row */
        else
            end = offset + number_rows-1;
    }

    // brocast the rank of the boundary process
   // MPI_Bcast(&boundaryInfo.holderRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    mpiManager->bCast<int>(&boundaryInfo.holderRank, 1,0);

    if (taskid != MASTER)
    {

        /* take a look at how the work is partiioned among processors */
        printf("worker number = %d offset= %d number_rows= %d start = %d end =%d bd_holder=%d\n",
               worker_number, offset,number_rows,start, end, boundaryInfo.holderRank);
    }

    iz = 0;
}

void Heat::simulate(){

    this->init();
    /* Begin doing TIME_STEPS iterations. Must communicate border rows with
           neighbors. If I have the first or last grid row, then I only need to
           communicate with one neighbor */

    for (it = 1; it <= TIME_STEPS; it++)
    {
        this->boundary();
        /* Now call update to update the value of grid points */
        //this->solve(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
        this->solve();

    }
    this->reduceOnMaster();
    this->saveGridInFile(this->final_filename);
}

void Heat::allocate(){
    u = new float**[2];
    for (int k = 0; k < 2; ++k) {
        u[k] = new float*[NXPROB];

        for (int i = 0; i < NXPROB; ++i)
            u[k][i] = new float[NYPROB];
    }

    /*    u.resize(2);
    for (int iz = 0; iz< 2; iz++){
        u[iz].resize(NXPROB);
        for (int i = 0; i< NXPROB; i++){
            u[iz][i].resize(NYPROB);
        }
    }
*/
}




void Heat::solve()
{

    float** u1= u[iz];
    float ** u2= u[1-iz];
    int ny = NYPROB;

    if (taskid != MASTER){
        int ix, iy;
        for (ix = start; ix <= end; ix++)
            for (iy = 1; iy <= ny-2; iy++){
                /**(u2+ix*ny+iy) = *(u1+ix*ny+iy) +
                    diffusivity.cx * (*(u1+(ix+1)*ny+iy) +
                                      *(u1+(ix-1)*ny+iy) -
                                      2.0 * *(u1+ix*ny+iy)) +
                    diffusivity.cy * (*(u1+ix*ny+iy+1) +
                                      *(u1+ix*ny+iy-1) -
                                      2.0 * *(u1+ix*ny+iy));*/

                u2[ix][iy] = u1[ix][iy] +
                        diffusivity.cx * (u1[ix+1][iy] +
                        u1[ix-1][iy] -
                        2.0 * u1[ix][iy]) +
                        diffusivity.cy * (u1[ix][iy+1] +
                        u1[ix][iy-1] -
                        2.0 * u1[ix][iy]);
            }
    }
    iz = 1 - iz;
}

void Heat::boundary(){
    if (taskid != MASTER)
    {
        if (neighbor1 != NONE)
        {
            mpiManager->send<float>(&u[iz][offset][0], NYPROB,  neighbor1, NGHBOR2);
            source = neighbor1;
            message_tag = NGHBOR1;
            mpiManager->receive<float>(&u[iz][offset-1][0], NYPROB, source, message_tag);
        }
        if (neighbor2 != NONE)
        {
            mpiManager->send<float>(&u[iz][offset+number_rows-1][0], NYPROB,  neighbor2,NGHBOR1);
            source = neighbor2;
            message_tag = NGHBOR2;
            mpiManager->receive<float>(&u[iz][offset+number_rows][0], NYPROB, source, message_tag);
        }

    }
}

void Heat::reduceOnMaster(){
    if (taskid == MASTER)
    {

        /* Now wait for results from all worker tasks */
        for (i=1; i<=numworkers; i++)
        {
            source = i;
            message_tag = DONE;
            mpiManager->receive<int>(&offset, 1, source, message_tag);
            mpiManager->receive<int>(&number_rows, 1, source, message_tag);
            for (int k=0; k<number_rows; k++ ){
                mpiManager->receive<float>(&u[0][offset+k][0], NYPROB, source, message_tag);
            }
            // mpiManager->receive(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);
        }
    }
    if (taskid != MASTER){
        /* Finally, send my portion of final results back to master */
        mpiManager->send<int>(&offset, 1,  MASTER, DONE);
        mpiManager->send<int>(&number_rows, 1,  MASTER, DONE);
        for (int k=0; k<number_rows; k++ ){
            mpiManager->send<float>(&u[iz][offset+k][0], NYPROB, MASTER, DONE);
        }
        //mpiManager->send(&u[iz][offset][0], number_rows*NYPROB, MPI_FLOAT, MASTER, DONE,MPI_COMM_WORLD);
    }
    /*gracefully exit MPI */
   // MPI_Finalize();
}

void Heat::saveGridInFile(string fileName){
    if (taskid == MASTER)
    {
        /* Write final output*/
        prtdat(NXPROB, NYPROB, u, fileName.c_str());

    }
}
