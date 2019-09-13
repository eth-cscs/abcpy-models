#include "distributedsimulation.h"

//---------------------------------------------------------------------------------------------
/*MuscleDistributedPlume::MuscleDistributedPlume(params p, shared_ptr<MPITopology> topology )
    : DistributedPlume(p, topology) {
    shared_ptr<MpiManager> mpiManager(std::make_shared<MpiManager>()) ;//getManagerMPI();
    mpiManager->init(& argc, &argv);
    int rank= mpiManager->getRank();
    connector = &TetrasConnector::getInstance();
    connector->init(rank, argc, argv);
}*/

MuscleDistributedPlume::MuscleDistributedPlume( int * argc, char ***argv)
    : DistributedPlume(*argc, *argv) {

    connector = &TetrasConnector::getInstance();
    connector->setCommandLineArgs(argc, argv);
}

MuscleDistributedPlume::~MuscleDistributedPlume(){}


int MuscleDistributedPlume::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}

void MuscleDistributedPlume::F_init(){
    this->createVolcanTopology();
    connector->init(this->topology_->getRank());
    this->initLogger(connector->getName());
    DistributedPlume::F_init();

}


void MuscleDistributedPlume::send(vector<char> const& data){
    connector->send(data);
}

void MuscleDistributedPlume::receive(vector<char> & data){
    connector->receive(data);
}

void MuscleDistributedPlume::end(){
    this->connector->end();
}

void MuscleDistributedPlume::sendConvergence()
{
        bool islocalConverged= !this->isProcess;
        vector<char> buffer;
        buffer.resize(sizeof(islocalConverged));
        memcpy((void *) &buffer[0], (const char*) &islocalConverged, sizeof(islocalConverged));
        this->send(buffer);
}

bool MuscleDistributedPlume::receiveConvergence(){
        int count=0;
        vector<char> vect;
        this->receive(vect);
        count=vect.size();
        MapperUtil::deserialize(this->isRemoteConverged, vect);
        vector<char>().swap(vect);
        return (count)? true:false;
}


//---------------------------------------------------------------------------------------------
MuscleDistributedAtmosphereSubModel::MuscleDistributedAtmosphereSubModel(double dx, double dt, int *argc, char ***argv, int rank)
    : DistributedAtmosphereSubModel(dx, dt) {
    this->manager = new MpiManager();
    this->manager->init(argc, argv);
    connector = &TetrasConnector::getInstance();
    connector->setCommandLineArgs(argc, argv);
    connector->init(rank);
}

MuscleDistributedAtmosphereSubModel::~MuscleDistributedAtmosphereSubModel(){
    delete manager;
}

int MuscleDistributedAtmosphereSubModel::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}


void MuscleDistributedAtmosphereSubModel::send(vector<char> const& data){
    connector->send(data);
}

void MuscleDistributedAtmosphereSubModel::receive(vector<char> & data){
    connector->receive(data);
}
void MuscleDistributedAtmosphereSubModel::end(){
    this->connector->end();
}

void MuscleDistributedAtmosphereSubModel::F_init()
{
    DistributedAtmosphereSubModel::F_init();
    string kname= this->connector->getName()+".log";
    // plb_ofstream  pcout(kname.c_str());
    this->logger= std::move(std::make_shared<plb_ofstream>(this->manager->getRank(), kname.c_str()));
}

//---------------------------------------------------------------------------------------------
MuscleDistributedAtmospheredPlume::MuscleDistributedAtmospheredPlume(int *argc, char ***argv)
    : DistributedAtmospheredPlume(*argc,* argv) {

    connector = &TetrasDoubleConnector::getInstance();
    connector->setCommandLineArgs( argc, argv);
}

MuscleDistributedAtmospheredPlume::~MuscleDistributedAtmospheredPlume(){}

int MuscleDistributedAtmospheredPlume::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}

void MuscleDistributedAtmospheredPlume::F_init(){
    this->createVolcanTopology();
    connector->init(this->topology_->getRank());
    this->initLogger(connector->getName());
    DistributedAtmospheredPlume::F_init();

}

void MuscleDistributedAtmospheredPlume::send(vector<char> const& data){
    return this->connector->sendToContinent(data);
}

void MuscleDistributedAtmospheredPlume::receive(vector<char> & data){
    this->connector->receiveFromContinent(data);
}
void MuscleDistributedAtmospheredPlume::sendToAtmosphere(vector<char> const& data){
    this->connector->sendToAtmosphere(data);
}

void MuscleDistributedAtmospheredPlume::receiveFromAtmosphere(vector<char> & data){
    this->connector->receiveFromAtmosphere(data);
}
void MuscleDistributedAtmospheredPlume::end(){
    this->connector->end();
}

void MuscleDistributedAtmospheredPlume::sendConvergence()
{
        bool islocalConverged= !this->isProcess;
        vector<char> buffer;
        buffer.resize(sizeof(islocalConverged));
        memcpy((void *) &buffer[0], (const char*) &islocalConverged, sizeof(islocalConverged));
        this->send(buffer);
}

bool MuscleDistributedAtmospheredPlume::receiveConvergence(){
        int count=0;
        vector<char> vect;
        this->receive(vect);
        count=vect.size();
        MapperUtil::deserialize(this->isRemoteConverged, vect);
        vector<char>().swap(vect);
        return (count)? true:false;
}

//---------------------------------------------------------------------------------------------

MuscleDistributedContinent::MuscleDistributedContinent(int *argc, char ***argv)
    : DistributedContinent(*argc, *argv) {

    connector = &TetrasConnector::getInstance();
    connector->setCommandLineArgs(argc, argv);
}

MuscleDistributedContinent::~MuscleDistributedContinent(){}

int MuscleDistributedContinent::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}

void MuscleDistributedContinent::F_init(){
    this->createVolcanTopology();
    connector->init(this->topology_->getRank());
    this->initLogger(connector->getName());
    DistributedContinent::F_init();

}

void MuscleDistributedContinent::send(vector<char> const& data){
    connector->send(data);
}

void MuscleDistributedContinent::receive(vector<char> & data){
    connector->receive(data);
}
void MuscleDistributedContinent::end(){
    this->connector->end();
}

void MuscleDistributedContinent::sendConvergence()
{
        bool islocalConverged= !this->isProcess;
        vector<char> buffer;
        buffer.resize(sizeof(islocalConverged));
        memcpy((void *) &buffer[0], (const char*) &islocalConverged, sizeof(islocalConverged));
        this->send(buffer);
}

bool MuscleDistributedContinent::receiveConvergence(){
        int count=0;
        vector<char> vect;
        this->receive(vect);
        count=vect.size();
        MapperUtil::deserialize(this->isRemoteConverged, vect);
        vector<char>().swap(vect);
        return (count)? true:false;
}

//---------------------------------------------------------------------------------------------
MuscleDistributedAtmospheredContinent::MuscleDistributedAtmospheredContinent( int *argc, char ***argv)
    : DistributedAtmospheredContinent(*argc, *  argv) {

    connector = &TetrasDoubleConnector::getInstance();
    connector->setCommandLineArgs(argc, argv);
}

MuscleDistributedAtmospheredContinent::~MuscleDistributedAtmospheredContinent(){}

int MuscleDistributedAtmospheredContinent::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}

void MuscleDistributedAtmospheredContinent::F_init(){
    this->createVolcanTopology();
    connector->init(this->topology_->getRank());
    this->initLogger(connector->getName());
    DistributedAtmospheredContinent::F_init();

}
void MuscleDistributedAtmospheredContinent::send(vector<char> const& data){
    return this->connector->sendToContinent(data);
}

void MuscleDistributedAtmospheredContinent::receive(vector<char> & data){
    this->connector->receiveFromContinent(data);
}
void MuscleDistributedAtmospheredContinent::sendToAtmosphere(vector<char> const& data){
    this->connector->sendToAtmosphere(data);
}

void MuscleDistributedAtmospheredContinent::receiveFromAtmosphere(vector<char> & data){
    this->connector->receiveFromAtmosphere(data);
}
void MuscleDistributedAtmospheredContinent::end(){
    this->connector->end();
}

void MuscleDistributedAtmospheredContinent::sendConvergence()
{
        bool islocalConverged= !this->isProcess;
        vector<char> buffer;
        buffer.resize(sizeof(islocalConverged));
        memcpy((void *) &buffer[0], (const char*) &islocalConverged, sizeof(islocalConverged));
        this->send(buffer);
}

bool MuscleDistributedAtmospheredContinent::receiveConvergence(){
        int count=0;
        vector<char> vect;
        this->receive(vect);
        count=vect.size();
        MapperUtil::deserialize(this->isRemoteConverged, vect);
        vector<char>().swap(vect);
        return (count)? true:false;
}
/************************************** FactorMuscleDistributedContinent *********************************************/

FactorMuscleDistributedPlume::FactorMuscleDistributedPlume():FactoryDistributedPlume(){}
FactorMuscleDistributedPlume::~FactorMuscleDistributedPlume(){}
DistributedPlume * FactorMuscleDistributedPlume::newInstance(int argc, char **argv){
    params p;
    p.version_ = VERSION;
    ArgParser & parser= getArgParser();

    try{
        parser.parseArgs(argc, argv, p);
    }catch(exception& e){
        cout << "Standard exception: " << e.what() << endl;
    }
    bool isMuscle=false;
#ifdef MUSCLE
    isMuscle=p.isMuscle_;
#endif
    if (!isMuscle){
        LightDistributedPlume* sim = new LightDistributedPlume( argc,  argv);
        return (DistributedPlume * )sim;
    }else{//muscle
        MuscleDistributedPlume* sim  = new MuscleDistributedPlume( &argc,  &argv);
        return (DistributedPlume * ) sim;
    }

}
/************************************** FactorMuscleDistributedContinent *********************************************/

FactoryMuscleDistributedAtmospheredPlume::FactoryMuscleDistributedAtmospheredPlume():FactoryDistributedPlume(){}
FactoryMuscleDistributedAtmospheredPlume::~FactoryMuscleDistributedAtmospheredPlume(){}
DistributedAtmospheredPlume *FactoryMuscleDistributedAtmospheredPlume::newInstance(int argc, char **argv){
   params p;
   p.version_ = VERSION;
   ArgParser & parser= getArgParser();

   try{
       parser.parseArgs(argc, argv, p);
   }catch(exception& e){
       cout << "Standard exception: " << e.what() << endl;
   }
    bool isMuscle=false;
#ifdef MUSCLE
    isMuscle=p.isMuscle_;
#endif
    if (!isMuscle){
        LightDistributedAtmospheredPlume* sim = new LightDistributedAtmospheredPlume( argc,  argv);
        return (DistributedAtmospheredPlume * )sim;
    }else{//muscle
        MuscleDistributedAtmospheredPlume* sim  = new MuscleDistributedAtmospheredPlume( &argc,  &argv);
        return (DistributedAtmospheredPlume * ) sim;
    }

}

/************************************** FactorMuscleDistributedContinent *********************************************/

FactorMuscleDistributedContinent::FactorMuscleDistributedContinent():FactoryDistributedContinent(){}
FactorMuscleDistributedContinent::~FactorMuscleDistributedContinent(){}
DistributedContinent * FactorMuscleDistributedContinent::newInstance(int argc, char **argv){
    params p;
    p.version_ = VERSION;
    ArgParser & parser= getArgParser();

    try{
        parser.parseArgs(argc, argv, p);
    }catch(exception& e){
        cout << "Standard exception: " << e.what() << endl;
    }
    bool isMuscle=false;
#ifdef MUSCLE
    isMuscle=p.isMuscle_;
#endif
    if (!isMuscle){
        LightDistributedContinent* sim = new LightDistributedContinent( argc,  argv);
        return (DistributedContinent * )sim;
    }else{//muscle
        MuscleDistributedContinent* sim  = new MuscleDistributedContinent( &argc,  &argv);
        return (DistributedContinent * ) sim;
    }

}

/************************************** FactorMuscleDistributedContinent *********************************************/

FactoryMuscleDistributedAtmospheredContinent::FactoryMuscleDistributedAtmospheredContinent():FactoryDistributedContinent(){}
FactoryMuscleDistributedAtmospheredContinent::~FactoryMuscleDistributedAtmospheredContinent(){}
DistributedAtmospheredContinent *FactoryMuscleDistributedAtmospheredContinent::newInstance(int argc, char **argv){
    params p;
    p.version_ = VERSION;
    ArgParser & parser= getArgParser();

    try{
        parser.parseArgs(argc, argv, p);
    }catch(exception& e){
        cout << "Standard exception: " << e.what() << endl;
    }
    bool isMuscle=false;
#ifdef MUSCLE
    isMuscle=p.isMuscle_;
#endif
    if (!isMuscle){
        LightDistributedAtmospheredContinent* sim = new LightDistributedAtmospheredContinent( argc,  argv);
        return (DistributedAtmospheredContinent * )sim;
    }else{//muscle
        MuscleDistributedAtmospheredContinent* sim  = new MuscleDistributedAtmospheredContinent( &argc,  &argv);
        return (DistributedAtmospheredContinent * ) sim;
    }

}
