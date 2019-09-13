/**
* @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>
* @version 1.0
* @section LICENSE

* MAPPER communication module
* Copyright (C) 2015  University of Geneva, Switzerland
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MAPPER_H
#define MAPPER_H

#include "musclehpc/mapper/mappercommon.h"
#include "utility.h"
#include "musclehpc/utils/mmsflogger.h"
#include "musclehpc/conduits/conduit.h"


#include <boost/date_time/posix_time/time_serialize.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/gregorian/greg_serialize.hpp>
//#include <cereal/archives/binary.hpp>

using namespace boost::posix_time;
using namespace boost::gregorian;

namespace unige_pasc{



/**************************************TransportInterpolation **************************************************/
// Need to be run on one core without mpiexec
class TransportInterpolation: public Submodel_A, public MapperHeader{

public:

    virtual ~TransportInterpolation();

protected:
    TransportInterpolation();
    void synchronizeCommunication();
    virtual void F_init();
    virtual void S();
    virtual void U();
    virtual void Oi();
    virtual void Of();
    virtual bool isConverged();

protected:

    bool checkConvergence(string &status, vector<char> & fromFine, vector<char> & fromCoarse);
    void interpolate(vector<char> & ToCoarse, vector<char> & ToFine);
    //virtual
    virtual void sendToCoarse(vector<char> const&buffer)=0;
    virtual void receiveFromCoarse( vector<char> & buffer)=0;
    virtual void sendToFine(vector<char> const&buffer)=0;
    virtual bool receiveFineConvergence()=0;
    virtual bool receiveCoarseConvergence()=0;
    virtual void receiveFromFine( vector<char> & buffer)=0;
    virtual void receiveFromFin1(vector<char>  &buffer)=0;
    virtual void receiveFromFin2( vector<char> & buffer)=0;
    virtual int getMPIRank()=0;
    virtual void setUpCoarseAndFineConduits(double dt1, double dt2)=0;
    //
    virtual int getSubmodel_ID()=0;
    virtual void end()=0;
    virtual string getName()=0;


protected:
    bool isProcess;
    pluint it_factor;
    pluint cpt;
    vector<BoundaryParticle> particlesToCoarse, particlesToFine;// used to accumulate particles coming from both domains (fine and coarse if grid refinement)
    bool treatCoarse;// if true: read from Coarse and sent to it.
    bool isAllSubmodelsConverged;// true: both submodels converged, false otherwise.
    bool isCoarseConverged, isFineConverged; // convergence status per submodel :  check whether each submodel has no particle in its domain
    int coarseId, fineId;// muscle ID of each submodel
    double coarseDt, fineDt; // dt of each submodel
    double lastFimeTimeStamp; // last time stamp received from the fine submodel

    //Logger logger;
    shared_ptr<plb_ofstream> logger;
};



//******************************************** LightTransportInterpolation *************************************/

class LightTransportInterpolation: public TransportInterpolation, public AsynchronousRemoteMpiKernel{

public:
    LightTransportInterpolation();
    virtual ~LightTransportInterpolation();
protected:

    void sendToCoarse(vector<char> const&buffer);
    void receiveFromCoarse( vector<char> & buffer);
    void sendToFine(vector<char> const&buffer);
    virtual bool receiveFineConvergence();
    virtual bool receiveCoarseConvergence();
    void receiveFromFine( vector<char> & buffer);
    void receiveFromFin1(vector<char>  &buffer);
    void receiveFromFin2( vector<char> & buffer);
    int getMPIRank();
    void setUpCoarseAndFineConduits(double dt1, double dt2);
    //
    virtual int getSubmodel_ID();
    virtual void end();

    virtual void initialize(string id);
    virtual string getName();
    virtual void mainLoop();

protected:
    string fin1, fin2, fout1, fout2; // muscle consuits names.
    // string coarse_in, coarse_out, fine_in, fine_out; //conduits

    string id;
    ConduitEntrance<char> * coarse_out, *fine_out;
    ConduitExit<char> * coarse_in, *fine_in;


};



} //end name space

#endif // MAPPER_H
