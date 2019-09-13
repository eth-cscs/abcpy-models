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


#ifndef CONNECTORS_H
#define CONNECTORS_H

#include "cppmuscle.hpp"
#include "../utility.h"
#include "connectorcommon.h"


using namespace unige_pasc ;
using namespace std;

namespace unige_pasc{



/************************************** BasicMuscleConnector **************************************************/
class BasicMuscleConnector: public BasicConnectorInterface{
public:
   virtual void send(vector<char> const& data, std::string conduit_out);
   virtual void receive(vector<char> & data, std::string conduit_in);
   virtual void setCommandLineArgs (int *argc, char *** argv);
   virtual void init (int rank);
   virtual void end();
   virtual string getName();
   virtual int getMPIRank() const;

protected:
   BasicMuscleConnector();
   virtual ~BasicMuscleConnector();
   bool isInitialized;
   int argc;
   char **argv;
   int rankMPI; //MPI_Rank
   std::string name;
};

/************************************** TetrasConnector **************************************************/

class TetrasConnector : public BasicMuscleConnector, public Connector<TetrasConnector>
{
    friend class Connector<TetrasConnector>;

public:
   virtual void send(vector<char> const&buffer);
   virtual void receive( vector<char> & buffer);
   virtual void init (int rank);
   virtual pluint getSubmodel_ID(){return this->submodel_ID;}
   virtual pluint getsubmodels_NBR(){return this->submodels_NBR;}
   virtual void end();
protected:
    TetrasConnector();
    TetrasConnector(const TetrasConnector&){}
    virtual ~TetrasConnector();

//vars
    std::string conduit_fin, conduit_fout; // muscle consuits names.
    pluint submodel_ID; // unique id (number) of the submodel in the workflow (included in |1..n|)
    pluint submodels_NBR; // total number of submodels involved in the workflow

};

/********************************************PlumeConnector **********************************************/
class TetrasDoubleConnector : public BasicMuscleConnector , public Connector<TetrasDoubleConnector>
{
       friend class Connector<TetrasDoubleConnector>;


public:
    virtual void init (int rank);
    void sendToContinent(vector<char> const& data);
    void receiveFromContinent(vector<char> & data);
    void sendToAtmosphere(vector<char> const& data);
    void receiveFromAtmosphere(vector<char> & data);
    virtual pluint getSubmodel_ID(){return this->submodel_ID;}
    virtual pluint getsubmodels_NBR(){return this->submodels_NBR;}
    virtual void end();
private:
    TetrasDoubleConnector();
    TetrasDoubleConnector(const TetrasDoubleConnector&){}
    virtual ~TetrasDoubleConnector();


//vars
    std::string conduit_fin, conduit_fout; // muscle consuits names.
    pluint submodel_ID; // unique id (number) of the submodel in the workflow (included in |1..n|)
    pluint submodels_NBR; // total number of submodels involved in the workflow
    std::string fin_atm, fout_atm; // muscle conduits names.

};
/************************************** TransportInterpolationConnector **************************************************/

class TransportInterpolationConnector :  public BasicMuscleConnector, public BasicSmartConnectorInterface,  public Connector<TransportInterpolationConnector>
{
public:

    void sendToCoarse(vector<char> const&buffer);
    void receiveFromCoarse( vector<char> & buffer);
    void sendToFine(vector<char> const&buffer);
    void receiveFromFine( vector<char> & buffer);
    void receiveFromFin1(vector<char>  &buffer);
    void receiveFromFin2( vector<char> & buffer);
    void setUpCoarseAndFineConduits(double dt1, double dt2);

    void init (int rank);
    virtual pluint getSubmodel_ID();
    virtual pluint getsubmodels_NBR();
    virtual void end();

protected:


   friend class Connector<TransportInterpolationConnector>;
    TransportInterpolationConnector();
    TransportInterpolationConnector(const TransportInterpolationConnector&){}
    virtual ~TransportInterpolationConnector();
    void setCoarseIn(string conduit);
    void setCoarseOut(string conduit);
    void setFineIn(string conduit);
    void setFineOut(string conduit);


//vars
    string fin1, fin2, fout1, fout2; // muscle consuits names.
    string coarse_in, coarse_out, fine_in, fine_out; //conduits
    pluint submodel_ID; // unique id (number) of the submodel in the workflow (included in |1..n|)
    pluint submodels_NBR; // total number of submodels involved in the workflow

};

////
} //namespace
#endif // CONNECTORS_H
