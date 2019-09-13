#ifndef MUSCLEMAPPER_H
#define MUSCLEMAPPER_H

#include "../connectors/muscleconnectors.h"
#include "mapper.h"
#include "../../../tetras/src/Tools/Simulation/atmosphere.h"
namespace unige_pasc{

/**************************************MuscleTransportInterpolation **************************************************/
class MuscleTransportInterpolation: public TransportInterpolation{

public:
    MuscleTransportInterpolation(int rank, int *argc, char ***argv);
    virtual ~MuscleTransportInterpolation();
protected:

    void sendToCoarse(vector<char> const&buffer);
    void receiveFromCoarse( vector<char> & buffer);
    void sendToFine(vector<char> const&buffer);
    void receiveFromFine( vector<char> & buffer);
    void receiveFromFin1(vector<char>  &buffer);
    void receiveFromFin2( vector<char> & buffer);
    int getMPIRank();
    void setUpCoarseAndFineConduits(double dt1, double dt2);
    virtual bool receiveFineConvergence();
    virtual bool receiveCoarseConvergence();
    //
    virtual int getSubmodel_ID();
    virtual void end();
    virtual string getName();

protected:
    TransportInterpolationConnector * connector;
};





//-- name space
}
#endif // MUSCLEMAPPER_H
