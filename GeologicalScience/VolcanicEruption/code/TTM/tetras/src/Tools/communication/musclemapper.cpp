#include "musclemapper.h"



namespace unige_pasc{
/**************************************** SmartMapper ************************************************/

MuscleTransportInterpolation::MuscleTransportInterpolation(int rank, int *argc, char ***argv):
TransportInterpolation(){
    // init the connector
    this->connector = &TransportInterpolationConnector::getInstance();
    this->connector->setCommandLineArgs( argc, argv);
    this->connector->init(rank);

}
MuscleTransportInterpolation::~MuscleTransportInterpolation(){}

int MuscleTransportInterpolation::getSubmodel_ID(){
    return this->connector->getSubmodel_ID();
}
void MuscleTransportInterpolation::end(){
    this->connector->end();
}

string MuscleTransportInterpolation::getName()
{
    return this->connector->getName();
}


void MuscleTransportInterpolation::sendToCoarse(vector<char> const&buffer){
    this->connector->sendToCoarse(buffer);

}
void MuscleTransportInterpolation::receiveFromCoarse( vector<char> & buffer){
    this->connector->receiveFromCoarse(buffer);
}
void MuscleTransportInterpolation::sendToFine(vector<char> const&buffer){
    this->connector->sendToFine(buffer);
}
void MuscleTransportInterpolation::receiveFromFine( vector<char> & buffer){
    this->connector->receiveFromFine(buffer);
}
void MuscleTransportInterpolation::receiveFromFin1(vector<char>  &buffer){
    this->connector->receiveFromFin1(buffer);
}
void MuscleTransportInterpolation::receiveFromFin2( vector<char> & buffer){
    this->connector->receiveFromFin2(buffer);
}

int MuscleTransportInterpolation::getMPIRank(){
    return this->connector->getMPIRank();
}

void MuscleTransportInterpolation::setUpCoarseAndFineConduits(double dt1, double dt2){
    this->connector->setUpCoarseAndFineConduits(dt1, dt2);
}

bool MuscleTransportInterpolation::receiveFineConvergence()
{
    std::vector<char> vect;
    this->connector->receiveFromFine(vect);
    bool isConverged;
    MapperUtil::deserialize(isConverged, vect);
    vector<char>().swap(vect);
    return isConverged;
}

bool MuscleTransportInterpolation::receiveCoarseConvergence()
{
    std::vector<char> vect;
    this->connector->receiveFromCoarse(vect);
    bool isConverged;
    MapperUtil::deserialize(isConverged, vect);
    vector<char>().swap(vect);
    return isConverged;
}



//-- ns--
}
