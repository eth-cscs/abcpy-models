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


#ifndef MMSF_H
#define MMSF_H

#include <iostream>
#include <assert.h>
#include <atomic>
#include <memory>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream

#include "musclehpc/utils/util.h"
#include "musclehpc/utils/optionparser.h"
#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/conduits/conduit.h"
#include "musclehpc/mediator/relayer.h"

using namespace std;
using namespace unige_pasc;

namespace unige_pasc{
//-------------------------------- forward declaration -----------------------------------------

template<typename E>
class Binder;


/************************************************ ColorInfo **************************************************/
class ColorInfo{
private:
    Kernel * kernel;
    int color;
    int rankInf, rankSup;
    bool manager;
    bool obsolet;
    bool relayer;
public:
    ColorInfo(Kernel * k, int color, int rankInf, int rankSup, bool isManager=false);
    Kernel * getKernel()const;
    int getColor () const;
    bool isManager()const;
    int getGlobalLeaderRank()const;
    bool isObsolet();
    bool doesBelongToMe(int processRank);
    bool isRelayer();
    void setRelayer(bool value);

};
/************************************ MMSFCommon ************************************************/
class MMSF{

protected:
    map<string, Kernel* > addedKernels;
    map<string, ConduitCommon* > conduits;
    vector<string> cxaAllKernels;
    //string ID;
    typedef typename std::map<string, Kernel* >::iterator it_map_kernel;
    typedef typename std::map<string, ConduitCommon* >::iterator it_map_condtuis;
    typedef typename std::map<string, int>::const_iterator it_map_cores;

public:

    MMSF( );
    virtual ~MMSF();

    /**
     * Adds a kernel to the CxA.
     * @param Kernel The kernel class
     * @param id An unique name for the kernel.
     */
    virtual void addKernel( Kernel * k, string id, int cores=1) =0;

    /**
     * Adds a Remote kernel to the CxA. This kernel is Empty !
     * @param Kernel The kernel class
     * @param id An unique name for the kernel.
     */
    virtual void addKernel( string id)=0 ;

    /**
     * @brief addConduit insert a conduit
     * @param id the name of the conduit
     * @param conduit pointer to the conduit object
     */
    virtual  void addConduit( string id, ConduitCommon* conduit);
    virtual  bool loadCxACoupling ()=0;



    /**
     * @brief connect connects a  kernel port.
     * @param fromStr syntax kernel:port
     * @return Binder pointer
     */
    template<typename E>
    Binder<E> * connect( string fromStr ) ;


    /**
       * Connects two kernel ports together using a specified conduit factory.
       * @param fromStr The port of the sending kernel.
       * @param toStr The port of the receiving kernel.
       * @param conduit A conduit pointer.
       * @param conduitID An unique name for the conduit.
       */
    template<typename E>
    void connect( string fromStr, string toStr, Conduit<E> * conduit, string conduitID );

    /**
     *  @brief assignConduit assigns a conduit to an MMSF port
     *  @param isIn is true the conduit should be of type Fout, else of type Fin
     */
    template<typename E>
    Kernel * assignConduit(string portString, Conduit<E> * conduit, bool isIn);

    /**
     * @brief getKernel get pointers of the kernel having as name id.
     * @param id the name of the kernel
     * @return pointer to the kernel, 0 otherwise.
     */
    virtual Kernel * getKernel( string id );

    virtual  map<string, Kernel*> const & getListKernels()const ;
    /**
     * @brief SEL runs the Submodel Sxecution Loop compliantly to the MMSF methodology
     */
    virtual void compute();

    

protected:
    Kernel * getkernelOfColor(int color);

};



/************************************ MMSF_MPI ************************************************/

class MMSF_MPI: public MMSF{

protected:
    shared_ptr<MpiManager> globalMpiManager;
    map<int, int> mappingCoresColor;//<coresNumber, color>
    map<string,int> conduitsTagMap; // <conduitname, MPI-TAG>
    vector<ColorInfo *> colorsInfoMap;
    OptionParser * optionsParser; // parses the command line options

    // each submodel (kernel) is identified by a color (int).
    map<string, int> mapColorWhereToRun; //<kernelaName, color>
    // tag used for message between the manager Rank and the other MPI processes
    int const MPI_TAG_INIT_MANAGER_URL = 50;

public:

    MMSF_MPI(int argc, char** argv, MPI_Comm globalCommunicator);
    MMSF_MPI(int argc, char** argv);
    virtual ~MMSF_MPI();
    virtual void addKernel(Kernel * k, string kernelId, int cores);
    virtual void addKernel(Kernel * k, string kernelId);
    virtual void addKernel(string kernelId);
    virtual  bool loadCxACoupling ();

    /**
     * @brief getCoresNumber return the MPI cores number
     * @return the total cores Number
     */
    virtual int getCoresNumber();
    /**
     * @brief SEL runs the Submodel Sxecution Loop compliantly to the MMSF methodology
     */
    virtual void compute(bool useGlobalManager=false);

    /**
     * @brief get_InterComm_MPI_Unig_Tag returns a uniq MPI_TAG for the conduit
     * @return
     */
    virtual int get_InterComm_MPI_Unig_Tag(ConduitCommon *conduit);

   // shared_ptr<OptionParser> & getOptionParser();

   int getGlobalRank();
    //
protected:

    /**
     * @brief verifyPTPConduitsAndCores verifies all conduits are safe
     * @return boolean true means all conduits are safe, false otherwise.
     */
    virtual bool verifyPTPConduitsAndCores();

    /**
     * @brief perapreConduitsTag gives a uniq tag for each conduit. This tag is used to exchange data between the two MPI submodels.
     * When Two submodels S1 and S2 are interconnected with several conduits, each conduit uses a uniq tag for its messages.
     */
    virtual void perapreConduitsTag();
    /**
     * @brief computeTotalRequestedCores reads the Coupling description file and computes the rquested cores.
     * @param requestedTotalCores the computed result will be stored here
     * @param isManager Treu: the manager run in this binary on 1 core. False: manager will not run in this binary.
     * @return True if requestedTotalCores <= available cores. False othewise.
     */
    virtual bool computeTotalRequestedCores(int &requestedTotalCores, bool isManager=false);

    /**
     * @brief getColorInfoByRank gets the color to which an MPI process belongs
     * @param rank tha rank of process.
     * @return ColorInfo object
     */
    virtual ColorInfo * getColorInfoByRank( int rank);
    /**
     * @brief getRelayer return the relayer object is it exists
     * @returnColorInfo object
     */
    virtual ColorInfo * getRelayer();
    /**
     * @brief populateColorsInfoMap assigns colors to submodels in order to prepare the MPI intercomm split.
     * This populates the  colorsInfoMap map data structure
     * @param isManager True: manager will run here on 1 core. False otherwise.
     */
    virtual void populateColorsInfoMap(bool isManager=false);
    /**
     * @brief connectSubmodelsWithinsameBinary creates intercommunicators between the current submodel and the submodels running on the same binary.
     *  Here there is no need to pass through the manager.
     * @param localMpiManager the mpiManager object used by the current submodel
     * @returnnumber of active submodels running in this binary.
     */
    int connectSubmodelsWithinsameBinary(shared_ptr<MpiManager> &localMpiManager);
    /**
     * @brief launch starts the whole execution
     * @param localMpiManager the mpiManager object used by the current submodel
     * @param isManager True: manager will run here on 1 core. False otherwise.
     */
    virtual void launch(shared_ptr<MpiManager> &localMpiManager, bool isManager=false);

private:

    map<string, shared_ptr<MpiManager> > createManagerKernelsInterCommunicators(shared_ptr<MpiManager> &localMpiManager);
    map<string, shared_ptr<RelayerCommunicator> >   createRelayerKernelsIntercommunicators(shared_ptr<MpiManager> &localMpiManager);
    shared_ptr<MpiManager>  createInterCommunicator(ColorInfo* localInfo, ColorInfo* remoteInfo, shared_ptr<MpiManager> &localMpiManager);


};



/************************************ Binder ************************************************/
template<typename E>
class Binder{

protected:
    MMSF * cxa;
    Conduit<E> * conduit;
    string conduitID, from, to;

public:
    Binder( MMSF * cxa ):cxa(cxa) ,conduit(0){
    }
    virtual ~Binder(){}

    Binder<E> * From( string from ) {
        this->from = from;
        return this;
    }

    Binder<E> * To( string to ) {
        this->to = to;
        tryConnection();
        return this;
    }

    Binder<E> * with(string conduitID ) {
        this->conduitID = conduitID;
        BasicConduitFactory<E> defaultFactory ;
        conduit= defaultFactory.newInstance();
        tryConnection();
        return this;
    }

    Binder<E> * with_MTM_MPI(string conduitID ) {
        this->conduitID = conduitID;
        MPI_MTM_ConduitFactory<E> defaultFactory ;
        MPI_MTM_Conduit<E> * conduits_MTM_ptr= defaultFactory.newInstance();
        conduit=(Conduit<E> *) conduits_MTM_ptr;
        //MMSF_MPI * cxa_ptr= dynamic_cast<MMSF_MPI*>(cxa);
        //int uniq_MPI_tag= cxa_ptr->get_InterComm_MPI_Unig_Tag(conduitID);
        //conduits_MTM_ptr->setMPIComm_Tag(uniq_MPI_tag);
        tryConnection();
        return this;
    }

    Binder<E> * with_PTP_MPI(string conduitID ) {
        this->conduitID = conduitID;
        MPI_PTP_ConduitFactory<E> defaultFactory ;
        MPI_PTP_Conduit<E> * conduits_PTP_ptr= defaultFactory.newInstance();
        conduit=(Conduit<E> *) conduits_PTP_ptr;
        //MMSF_MPI * cxa_ptr= dynamic_cast<MMSF_MPI*>(cxa);
        //int uniq_MPI_tag= cxa_ptr->get_InterComm_MPI_Unig_Tag(conduitID);
        //conduits_PTP_ptr->setMPIComm_Tag(uniq_MPI_tag);
        tryConnection();
        return this;
    }


    virtual void tryConnection() {
        if ( !this->to.empty() && !conduitID.empty() && conduit) {
            assert(conduit);
            cxa->connect<E>( from, to, conduit, conduitID );
        }
    }
};


/************************************************************************************************/

}//end namespace

#endif // MMSF_H
