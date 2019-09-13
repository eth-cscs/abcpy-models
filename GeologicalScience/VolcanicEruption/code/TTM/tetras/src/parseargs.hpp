/*
TEphra TRAnsport Simulator (tetras)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARSEARGS_H_
#define PARSEARGS_H_

#include "piaf.hpp"
#include "version.h"


class ArgParser{

private:
    ArgParser(){}
    ~ArgParser(){}
public:

    void printUsage(string prgName){
        cout<<"Usage: "<<prgName<<" -s [SIMULATOR TYPE] -f [FILENAME]... [OPTION]..."<<endl;
        cout<<"Simulate particles in an advection field."<<endl<<endl;
        cout<<" -s [SIMULATOR TYPE]          simulator type to use (CA or EXACT)"<<endl;
        cout<<" -f [FILENAME]                the file containing the definition of the eruption"<<endl;
        cout<<" -o [FILENAME]                the file where the output is stored (hdf5), no file is written if no value specified"<<endl;
        cout<<" -npart [PARTICLE NUMBER]     the number of particles per class to insert"<<endl;
#ifdef DISPLAY
        cout<<" -d                           display, activate the 2D display"<<endl;
        cout<<" -v                           video, save the simulation in a video (2D), display mandatory"<<endl;
#endif
#ifdef CLI
        cout<<" -i                           interactive, activate the cli"<<endl;
#endif
        cout<<" -h, --help                   help, show this help"<<endl;
        cout<<" -dx [DX VALUE]               spatial integration step ( dx >= 0 )"<<endl;
        cout<<" -dt [DT VALUE]               time integration step ( dt >= 0 )"<<endl;
        cout<<" -ddt                         use the dynamic dt functionnality"<<endl;
        cout<<" -ddtValue [DDT VALUE]        use the given ddt value"<<endl;
        cout<<" -dx_t [DX TERRAIN VALUE]     spatial integration step for terrain ( dx_t >= 0, dx_t = dx if not specified )"<<endl;
        cout<<" -x [X VALUE]                 size of the domain along x axis ( x >= 0)"<<endl;
        cout<<" -y [Y VALUE]                 size of the domain along y axis ( y >= 0)"<<endl;
        cout<<" -z [X VALUE]                 size of the domain along z axis ( z >= 0)"<<endl;
        cout<<" -test                        perform some tests instead of running the simulation"<<endl;
        cout<<" -blockcyclic [VALUE]         use a block cyclic domain with given number of block per dimension for each process"<<endl;
        cout<<" -simple [DIM]                use a simple domain splitted along the given axis (x, y or z)"<<endl;
        cout<<" -determinist [seed]          perform a determinist simulation with given seed"<<endl;
        cout<<" -id                          use isotropic diffusion (in this case, horizontal diffusion of column and atmosphere is used)"<<endl;
        cout<<" -ox [X VALUE]                origin of the domain along x axis ( x >= 0)"<<endl;
        cout<<" -oy [Y VALUE]                origin of the domain along y axis ( y >= 0)"<<endl;
        cout<<" -oz [X VALUE]                origin of the domain along z axis ( z >= 0)"<<endl;
#ifdef MUSCLE
        cout<<" -muscle                      use muscle coupling. Options that follow the -muscle option will not be considered by the program)"<<endl;
#endif
        cout<<" -track                       enable tracking of each individual particle"<<endl;
        cout<<" -trackinterval               interval between two particles tracking"<<endl;
        cout<<" -trackwriteinterval          interval between two writing into file of particle tracking"<<endl;
        cout<<" -trackname                   base name for the tracking file"<<endl;
        cout<<" -tracktype                   file and type of tracking, can be hdf5 (default) or vtk"<<endl;
        cout<<" -particlescsv                writes particles into a csv file each iteration"<<endl;
        cout<<" -injectalongplume            inject particles vertically along plume instead of from the crater"<<endl;
        cout<<" -timeout [MINUTES]           set a time out after which the simulation return a result in its current state"<<endl;
    }

    void parseArgs(int argc, char **argv, params &p){

        cout << "ARGUMENTS PARSING" << endl;
        for(int i=0; i<argc; i++){
            cout << argv[i] << " ";
        }
        cout << endl;

        p.simulatorType_ = UNDEFINED_SIM;
        p.eruptionFile_ = "";
        p.outputFile_ = "";
        p.npart_ = -1;
#ifdef DISPLAY
        p.display_ = false;
        p.video_ = false;
#endif
#ifdef CLI
        p.isInteractive_ = false;
#endif
        p.dx_ = 500.0;
        p.dt_ = 0.5;
        p.ddt_ = false;
        p.ddtValue_ = 0.3;
        p.dxT_ = -1.0;
        p.domainSizeX_ = 30000.0;
        p.domainSizeY_ = 500.0;
        p.domainSizeZ_ = 50000.0;
        p.test_ = false;
        p.simpleDom_ = -1;
        p.blockCyclic_ = -1;
        p.seed_ = -1;
        p.id_ = false;
        p.origin_ = { 0.0, 0.0, 0.0 };
        p.nodif_ = false;
#ifdef MUSCLE
        p.isMuscle_=false;
#endif
        p.trackinterval_ = 0.0;
        p.trackwriteinterval_ = 100.0;
        p.track_ = false;
        p.trackFileBaseName_ = "trackpart";
        p.trackFileType_ = "hdf5";
        p.particlescsv_ = false;
        p.injectAtCrater_ = true;
        p.timeOut_ = -1.0;

        for (int i = 1; i < argc; i++) {


            /* domainSizeX_ */
            /*******************************************/
            if(strcmp(argv[i],"-x") == 0){
                try{
                    string text = argv[i+1];
                    p.domainSizeX_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter x"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* domainSizeY_ */
            /*******************************************/
            else if(strcmp(argv[i],"-y") == 0){
                try{
                    string text = argv[i+1];
                    p.domainSizeY_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter y"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* domainSizeZ_ */
            /*******************************************/
            else if(strcmp(argv[i],"-z") == 0){
                try{
                    string text = argv[i+1];
                    p.domainSizeZ_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter z"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            else if(strcmp(argv[i],"-ox") == 0){
                try{
                    string text = argv[i+1];
                    p.origin_.x_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter ox"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }

            else if(strcmp(argv[i],"-oy") == 0){
                try{
                    string text = argv[i+1];
                    p.origin_.y_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter oy"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }

            else if(strcmp(argv[i],"-oz") == 0){
                try{
                    string text = argv[i+1];
                    p.origin_.z_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter oz"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }

            /* np_ */
            /*******************************************/
            else if(strcmp(argv[i],"-npart") == 0){
                try{
                    string text = argv[i+1];
                    p.npart_ = boost::lexical_cast< int >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter np"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* dx_ */
            /*******************************************/
            else if(strcmp(argv[i],"-dx") == 0){
                try{
                    string text = argv[i+1];
                    p.dx_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter dx"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* dt_ */
            /*******************************************/
            else if(strcmp(argv[i],"-dt") == 0){
                try{
                    string text = argv[i+1];
                    p.dt_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter dt"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* nodif_ */
            /*******************************************/
            else if(strcmp(argv[i],"-nodif") == 0){
                p.nodif_ = true;
            }
            /*******************************************/

            /* ddt_ */
            /*******************************************/
            else if(strcmp(argv[i],"-ddt") == 0){
                p.ddt_ = true;
            }
            /*******************************************/

            else if(strcmp(argv[i],"-ddtValue") == 0){
                try{
                    string text = argv[i+1];
                    p.ddtValue_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter ddtValue"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }

            /* dxT_ */
            /*******************************************/
            else if(strcmp(argv[i],"-dx_t") == 0){
                try{
                    string text = argv[i+1];
                    p.dxT_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter dx_t"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* p.eruptionFile_ */
            /*******************************************/
            else if (strcmp(argv[i],"-f") == 0) {
                p.eruptionFile_ = string(argv[i+1]);
                i++;
            }
            /*******************************************/

            /* p.outputFile_ */
            /*******************************************/
            else if (strcmp(argv[i],"-o") == 0) {
                p.outputFile_ = string(argv[i+1]);
                i++;
            }
            /*******************************************/

            /* simulatorType_ */
            /*******************************************/
            else if (strcmp(argv[i],"-s") == 0) {
                if(strcmp(argv[i+1],"CA") == 0) p.simulatorType_ = CA;
                else if(strcmp(argv[i+1],"EXACT") == 0) p.simulatorType_ = EXACT;
                else{
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid simulator type."<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* simpleDom_ */
            /*******************************************/
            else if (strcmp(argv[i],"-simple") == 0) {
                if(strcmp(argv[i+1],"x") == 0) p.simpleDom_ = 0;
                else if(strcmp(argv[i+1],"y") == 0) p.simpleDom_ = 1;
                else if(strcmp(argv[i+1],"z") == 0) p.simpleDom_ = 2;
                else{
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid dimension."<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* blockCyclic_ */
            /*******************************************/
            else if (strcmp(argv[i],"-blockcyclic") == 0) {
                try{
                    string text = argv[i+1];
                    p.blockCyclic_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter blockcyclic"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* seed_ */
            /*******************************************/
            else if (strcmp(argv[i],"-determinist") == 0) {
                try{
                    string text = argv[i+1];
                    p.seed_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter determinist"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* trackinterval_ */
            /*******************************************/
            else if (strcmp(argv[i],"-trackinterval") == 0) {
                try{
                    string text = argv[i+1];
                    p.trackinterval_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "<<VERSION<<endl;
                    cout<<"Invalid value for parameter trackinterval"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* trackwriteinterval_ */
            /*******************************************/
            else if (strcmp(argv[i],"-trackwriteinterval") == 0) {
                try{
                    string text = argv[i+1];
                    p.trackwriteinterval_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "<<VERSION<<endl;
                    cout<<"Invalid value for parameter trackwriteinterval"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/


            /* trackFileBaseName_ */
            /*******************************************/
            else if (strcmp(argv[i],"-trackname") == 0) {
                try{
                    string text = argv[i+1];
                    p.trackFileBaseName_ = text;
                }
                catch(...){
                    cout<<"tetras version "<<VERSION<<endl;
                    cout<<"Invalid value for parameter trackname"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* trackFileType_ */
            /*******************************************/
            else if (strcmp(argv[i],"-tracktype") == 0) {
                try{
                    string text = argv[i+1];
                    p.trackFileType_ = text;
                    if( p.trackFileType_.compare("hdf5") != 0 && p.trackFileType_.compare("vtk") != 0 ){
                        throw 0;
                    }
                }
                catch(...){
                    cout<<"tetras version "<<VERSION<<endl;
                    cout<<"Invalid value for parameter tracktype, must be either hdf5 or vtk"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/

            /* timeOut_ */
            /*******************************************/
            else if(strcmp(argv[i],"-timeout") == 0){
                try{
                    string text = argv[i+1];
                    p.timeOut_ = boost::lexical_cast< double >( text );
                }
                catch(...){
                    cout<<"tetras version "+p.version_<<endl;
                    cout<<"Invalid value for parameter timeout"<<endl<<endl;
                    printUsage(argv[0]);
                    throw -1;
                }
                i++;
            }
            /*******************************************/



#ifdef DISPLAY
            /* display_ */
            /*******************************************/
            else if(strcmp(argv[i],"-d") == 0) p.display_ = true;
            /*******************************************/

            /* video_ */
            /*******************************************/
            else if(strcmp(argv[i],"-v") == 0) p.video_ = true;
            /*******************************************/
#endif
#ifdef CLI
            /* isInteractive_ */
            /*******************************************/
            else if(strcmp(argv[i],"-i") == 0) p.isInteractive_ = true;
            /*******************************************/
#endif
            /* test_ */
            /*******************************************/
            else if(strcmp(argv[i],"-test") == 0) p.test_ = true;
            /*******************************************/

            /* id_ */
            /*******************************************/
            else if(strcmp(argv[i],"-id") == 0) p.id_ = true;
            /*******************************************/

            /* track_ */
            /*******************************************/
            else if(strcmp(argv[i],"-track") == 0) p.track_ = true;
            /*******************************************/

            /* particlescsv_ */
            /*******************************************/
            else if(strcmp(argv[i],"-particlescsv") == 0) p.particlescsv_ = true;
            /*******************************************/

            /* injectAtCrater_ */
            /*******************************************/
            else if(strcmp(argv[i],"-injectalongplume") == 0){ 
                std::cout << "Setting injectAtCrater to false" << std::endl;
                p.injectAtCrater_ = false;
            }
            //std::cout << "Value of intectAtCrater in parseargs : " << p.injectAtCrater_ << std::endl;
            /*******************************************/

            /* help */
            /*******************************************/
            else if((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"--help") == 0)){
                cout<<"tetras version "+p.version_<<endl;
                printUsage(argv[0]);
                throw -1;
            }
#ifdef MUSCLE
            // MUSCLE arguments:
            // muscle args should be prefixed by '--'
            // so what comes after '-muscle' should not be processed by MUSCLE
            // what comes before '-muscle' should be treated by the program
            else if(strcmp(argv[i],"-muscle") == 0){
                cout<<"muscle argumenet detected"<<endl<<endl;
                p.isMuscle_=true;
                break;
            }
#endif
            /*******************************************/

            /* invalid argument */
            /*******************************************/
            else {
                //cout<<"tetras version "+p.version_<<endl;
                cout<<"Invalid argument : "<<argv[i]<<endl<<endl;
                //printUsage(argv[0]);
                //throw -1;
            }
            /*******************************************/
        }// end for parsing argv

        /* undefined simulator type */
        /*******************************************/
        if(p.simulatorType_ == UNDEFINED_SIM){
            cout<<"tetras version "+p.version_<<endl;
            cout<<"No simulator type specified"<<endl<<endl;
            printUsage(argv[0]);
            throw -1;
        }
        /*******************************************/

        /* input file undefined */
        /*******************************************/
        if(p.eruptionFile_ == ""){
            cout<<"tetras version "+p.version_<<endl;
            cout<<"No input file specified"<<endl<<endl;
            printUsage(argv[0]);
            throw -1;
        }
        /*******************************************/

        /* output file undefined */
        /*******************************************/
        if(p.outputFile_ == ""){
            cout<<"Warning, no output file specified !!!"<<endl<<endl;
        }
        /*******************************************/

        /* size incompatible with dx */
        /*******************************************/
        if(p.domainSizeX_ < p.dx_){cout<<"Domain size X < dx"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.domainSizeY_ < p.dx_){cout<<"Domain size Y < dx"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.domainSizeZ_ < p.dx_){cout<<"Domain size Z < dx"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        /*******************************************/

        /* Invalid value for dx, dt, size x y z */
        if(p.dx_ <= 0){cout<<"Invalid value for parameter dx"; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.dt_ <= 0){cout<<"Invalid value for parameter dt"; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.domainSizeX_ <= 0){cout<<"Invalid value for parameter x"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.domainSizeY_ <= 0){cout<<"Invalid value for parameter y"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
        if(p.domainSizeZ_ <= 0){cout<<"Invalid value for parameter z"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}

        if(p.dxT_ <= 0.0 ){p.dxT_ = p.dx_;}

#ifdef DISPLAY
        if(p.video_ && !p.display_){cout<<"If video is selected, display is mendatory"<<endl; cout<<"tetras version "+p.version_<<endl; printUsage(argv[0]); throw -1;}
#endif

    std::cout << "Value of injectAtCrater in parseargs : " << p.injectAtCrater_ << std::endl;
    
    }

    friend ArgParser& getArgParser();
};

inline ArgParser& getArgParser() {
    static ArgParser instance;
    return instance;
}

#endif
