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

#include "atmosphere.h"


/**************************************** DistributedSimulation ************************************************/

DistributedAtmosphereSubModel::DistributedAtmosphereSubModel(double dx, double dt) : AtmosphereSubModel (dx,dt){
  boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
  this->info.timeStamp=MapperUtil::to_time_t(tm);
  this->info.type=MapperType::SUBMODEL;
  this->info.dt=dt;
  this->info.dx=dx;
  //this->info.sectionId= this->getSubmodel_ID();
  #ifdef PROFILING
  this->profiler.setPrefixFileName("DistributedAtmosphereSubModel");
  #endif
}



DistributedAtmosphereSubModel::~DistributedAtmosphereSubModel(){}

void DistributedAtmosphereSubModel::synchronizeCommunication(){

  boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
  this->info.timeStamp=MapperUtil::to_time_t(tm);

  if( manager->getRank() == 0 ) { //requires only on
    vector<char> toSend;
    this->serializeHeader(toSend, this->info);
    this->send(toSend);//HERE: send toSend vector through to remote

    vector<char> toReceive;
    this->receive(toReceive);
    INFO_SMART_MAPPER info_tmp;
    this->unserializeHeader(toReceive, info_tmp);
    isRemoteConverged=false;
  }
}

void DistributedAtmosphereSubModel::F_init(){
  // call parent F_init
  AtmosphereSubModel::F_init();
  this->info.sectionId= this->getSubmodel_ID();
  this->synchronizeCommunication();
  #ifdef PROFILING
  std::stringstream  ss;
  ss<<"DistributedAtmosphereSubModel"<<"\n";
  this->profiler.setPrefixFileName(ss.str());
  nbrItrToSaveStatsInFile=50;
  //this->profiler.setRank(mpimanager->getRank());
  this->profiler.setRank(this->manager->getRank());
  #endif
  #ifdef PROFILING
  this->profiler.start();
  #endif
}


void DistributedAtmosphereSubModel::S(){
  #ifdef PROFILING
  this->profiler.createStampRecording();
  this->profiler.registerBeginProcessLevelStamp();
  #endif
  AtmosphereSubModel::S();
  #ifdef PROFILING
  this->profiler.registerEndProcessLevelStamp();
  #endif

}

bool DistributedAtmosphereSubModel::isConverged(){
  return AtmosphereSubModel::isConverged() && this->isRemoteConverged;
}

// compute new lat, long by moving from an offset in meters
// http://gis.stackexchange.com/questions/2951/algorithm-for-offsetting-a-latitude-longitude-by-some-amount-of-meters
std::pair<double, double> dispCoord(double lat, double lon, double ofx, double ofy){
  double latrad = lat * (M_PI/180.0);
  double newLat = lat + ofy / 111111.0;
  double newLon = lon + ofx / (111111.0 * cos(latrad));
  return std::pair<double, double>(newLat, newLon);
}

AtmosphereDataResponse DistributedAtmosphereSubModel::createAtmosphereDataResponse(AtmosphereDataRequest &atmRequest){

  double ox = atmRequest.getDomainPosition().getX();
  double oy = atmRequest.getDomainPosition().getY();
  double vx = atmRequest.getCraterePosition().getX();
  double vy = atmRequest.getCraterePosition().getY();
  double sx = atmRequest.getDomainSize().getX();
  double sy = atmRequest.getDomainSize().getY();
  double vlat = atmRequest.getGeoCraterePosition().getLattitude();
  double vlon = atmRequest.getGeoCraterePosition().getLongitude();
  string date = atmRequest.getDate();

  char hostname[150];
  memset(hostname, 0, 150);
  gethostname(hostname, 150);
  int pid = getpid();
  stringstream ss1;
  ss1 << hostname << "-" << pid;
  string uid = ss1.str();

  // toward west -> longitude becomes lower
  double dxll = ox - vx;
  // toward south -> latitude becomes lower
  double dyll = oy - vy;
  double dxur = (ox + sx) - vx;
  double dyur = (oy + sy) - vy;

  std::pair<double, double> latLonLL = dispCoord( vlat, vlon, dxll, dyll );
  std::pair<double, double> latLonUR = dispCoord( vlat, vlon, dxur, dyur );

  logger->getWriter() << "vent latitude : " << vlat << "\n";
  logger->getWriter() << "vent longitude : " << vlon << "\n";
  logger->getWriter() << "dxll : " << dxll << "\n";
  logger->getWriter() << "dyll : " << dyll << "\n";
  logger->getWriter() << "lower left latitude : " << get<0>(latLonLL) << "\n";
  logger->getWriter() << "lower left longitude : " << get<1>(latLonLL) << "\n";
  logger->getWriter() << "dxur : " << dxur << "\n";
  logger->getWriter() << "dyur : " << dyur << "\n";
  logger->getWriter() << "upper right latitude : " << get<0>(latLonUR) << "\n";
  logger->getWriter() << "upper right longitude : " << get<1>(latLonUR) << "\n";
  logger->getWriter() << "Date : " << date << "\n";

  double n = get<0>(latLonUR); double s = get<0>(latLonLL);
  if(n>s){n = ceil(n); s = floor(s);} else{s = ceil(s); n = floor(n);}
  double w = get<1>(latLonLL); double e = get<1>(latLonUR);
  if(w>e){w = ceil(w); e = floor(e);} else{e = ceil(e); w = floor(w);}

  stringstream ss;
  ss << "python download_ECMWF.py -bd "<<date<<" -ed "<<date;
  ss<<" -n "<<int(n);
  ss<<" -s "<<int(s);
  ss<<" -w "<<int(w);
  ss<<" -e "<<int(e);
  ss<<" -id "<<uid;

  std::string command = ss.str();
 // std::string command = "python download_ECMWF.py -bd "+date+" -ed "+date+" -n "+std::to_string(int(n))+" -s "+std::to_string(int(s))+" -w "+std::to_string(int(w))+" -e "+std::to_string(int(e));
  logger->getWriter() << "command : " << command << "\n";
  int retval = system(command.c_str());
  logger->getWriter() << "return value : " << retval << "\n";

  std::string filename = uid+"_"+date+"_"+date+"_"+std::to_string(int(n))+"_"+std::to_string(int(s))+"_"+std::to_string(int(w))+"_"+std::to_string(int(e))+".nc";
  logger->getWriter() << "filename : " << filename << "\n";

  logger->getWriter() << "latitude " << get<0>(latLonUR) << " has been rounded to " << n << "\n";
  logger->getWriter() << "latitude " << get<0>(latLonLL) << " has been rounded to " << s << "\n";
  logger->getWriter() << "longitude " << get<1>(latLonLL) << " has been rounded to " << w << "\n";
  logger->getWriter() << "longitude " << get<1>(latLonUR) << " has been rounded to " << e << "\n";

  try{
    //NcFile dataFile(filename, NcFile::read);
    int dataFileId, retval;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &dataFileId)))
      throw(retval);
    // TODO recode with C api... add close file...

    /*

    logger->getWriter()<<"there are "<<dataFile.getVarCount()<<" variables"<<"\n";
    logger->getWriter()<<"there are "<<dataFile.getAttCount()<<" attributes"<<"\n";
    logger->getWriter()<<"there are "<<dataFile.getDimCount()<<" dimensions"<<"\n"<<"\n";

    logger->getWriter() << "variables are : " << "\n";
    for(auto dim : dataFile.getVars()){
      logger->getWriter() << get<0>(dim) << "\n";
    }

    logger->getWriter() << "\n";

    logger->getWriter() <<"attributes are : " << "\n";
    for(auto dim : dataFile.getAtts()){
      logger->getWriter() << get<0>(dim) << "\n";
    }

    logger->getWriter() << "\n";

    logger->getWriter() << "dimensions are : " << "\n";
    for(auto dim : dataFile.getDims()){
      logger->getWriter() << get<0>(dim) << "\n";
    }

    */

    /*
    auto lats = dataFile.getVar("latitude");
    auto lons = dataFile.getVar("longitude");

    logger->getWriter() << "latitude var has " << lats.getDims().size() << " dimension of size "<<  lats.getDims().at(0).getSize() <<  "\n";
    double * latitudes = new double[lats.getDims().at(0).getSize()];
    lats.getVar( latitudes );

    logger->getWriter() << "latitudes : " << "\n";
    for(int i=0; i<lats.getDims().at(0).getSize(); i++) logger->getWriter() << latitudes[i] << "\n";

    double * longitudes = new double[lons.getDims().at(0).getSize()];
    lons.getVar( longitudes );

    logger->getWriter() << "longitudes : " << "\n";
    for(int i=0; i<lons.getDims().at(0).getSize(); i++) logger->getWriter() << longitudes[i] << "\n";

    auto levs = dataFile.getVar("level");
    double * levels = new double[levs.getDims().at(0).getSize()];
    levs.getVar( levels );
    logger->getWriter() << "levels : " << "\n";
    for(int i=0; i<levs.getDims().at(0).getSize(); i++) logger->getWriter() << levels[i] << "\n";

    auto tims = dataFile.getVar("time");
    double * times = new double[tims.getDims().at(0).getSize()];
    tims.getVar( times );
    logger->getWriter() << "times : " << "\n";
    for(int i=0; i<tims.getDims().at(0).getSize(); i++) logger->getWriter() << times[i] << "\n";
    */

  //}catch(NcException& e){
  }catch(int& e){
    //e.what();
    logger->getWriter()<<"UNABLE TO OPEN ATMOSPHERE FILE : "<<filename<<" -  error code : "<< e <<"\n";
  }






  /*===================================================================================================*/
  vector<double> data;
  data.push_back(90.0+cpt);data.push_back(91.0+cpt);data.push_back(92.0+cpt);data.push_back(93.0+cpt);
  double dx=100.0, dy=200.0, dz=300.0;

  AtmosphereDataResponse atmResponse(dx, dy, dz, data);
  return atmResponse;
}

void DistributedAtmosphereSubModel::U(){


  bool debug=true;
  //== remove the particles from the domain ==

  #ifdef PROFILING
  this->profiler.registerBeginGatherLevelStamp();
  #endif

  //
  #ifdef PROFILING
  this->profiler.registerEndGatherLevelStamp();
  #endif

  //
  #ifdef PROFILING
  double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
  #endif
  stringstream ss;
  //========== receive from Plume a request ==============

  AtmosphereDataRequest atmRequest;
  if(this->manager->isMainProcessor()){
    vector<char> buffer;
    this->receive(buffer);
    char * offset= &buffer[0];
    memcpy(&isRemoteConverged, offset, sizeof(isRemoteConverged)); // "Serialize"
    offset+=sizeof(isRemoteConverged);
    atmRequest.unSerialize(buffer, offset);

    ss<<"---------------"<<"\n";
    ss<<"responser["<<this->cpt<<"]: "<<"\n";
    ss<<"remote convergence="<<this->isConverged()<<"\n";
    ss<<"<-\n"<<atmRequest.toString()<<"\n";
    logger->getWriter()<<ss.str()<<"\n";
  }
  this->manager->bCast<bool>(&isRemoteConverged, 1, this->manager->bossId() );
  //========== treat the request ==============
  // treat request


  #ifdef PROFILING
  double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
  double WTimeReceiveBegin=  WTimeSendEnd;
  #endif
  //========== prepare response ==============
  vector<char> toSend;
  if(this->manager->isMainProcessor()){
    AtmosphereDataResponse atmResponse= createAtmosphereDataResponse(atmRequest);
    atmResponse.serialize(toSend);
    ss<<"->\n"<<atmResponse.toString()<<"\n";
    ss<<"---------------"<<"\n";
    logger->getWriter() <<ss.str();
    //========== send response ==============
    this->send(toSend);
  }


  #ifdef PROFILING
  double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
  this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
  this->profiler.registerBeginUpdateLevelStamp();
  #endif
  //

  #ifdef PROFILING
  this->profiler.registerEndUpdateLevelStamp();
  #endif

  //
  #ifdef PROFILING
  this->profiler.commitStampRecording();
  //if(this->cpt% nbrItrToSaveStatsInFile==0 && this->cpt>0){
  if(cpt% nbrItrToSaveStatsInFile==0 && cpt>0){
    this->profiler.saveToFile();
  }
  #endif
  if( debug && manager->getRank() == 0 ) logger->getWriter() <<  "====================" << "\n";
}

void DistributedAtmosphereSubModel::Of(){
  // call parent
  AtmosphereSubModel::Of();
  #ifdef PROFILING
  bool considerParticles=false; // i take it from volcano :)
  this->profiler.end(considerParticles);
  if( manager->getRank() == 0 ) logger->getWriter() << this->profiler.getAveragedStats()<<"\n";
  #endif
  this->end();
}

/************************************LightDistributedAtmosphereSubModel*************************************************/
LightDistributedAtmosphereSubModel::LightDistributedAtmosphereSubModel(double dx, double dt):
DistributedAtmosphereSubModel(dx,dt), AsynchronousRemoteMpiKernel(){

}

LightDistributedAtmosphereSubModel::~LightDistributedAtmosphereSubModel(){}

int LightDistributedAtmosphereSubModel::getSubmodel_ID(){
  return 3;
}


void LightDistributedAtmosphereSubModel::send(vector<char> const& toSend){
  this->iSend(f_out, (char*)toSend.data(), (int) toSend.size());
  this->waitSend();
}

void LightDistributedAtmosphereSubModel::receive(vector<char> & data){
  int count;
  char * buffer= this->iRecv(f_in, count);
  this->waitReceive();
  std::vector<char> vect(buffer, buffer + count);
  vect.swap(data);
   delete[] buffer;
}

void LightDistributedAtmosphereSubModel::end(){}

void LightDistributedAtmosphereSubModel::F_init()
{
    DistributedAtmosphereSubModel::F_init();
    string kname= this->getName()+".log";
    // plb_ofstream  pcout(kname.c_str());
    this->logger= std::move(std::make_shared<plb_ofstream>(this->manager->getRank(), kname.c_str()));
}


void LightDistributedAtmosphereSubModel::mainLoop(){
  this->manager= this->mpiManager.get();
   // string kname = // 1st argument is the kernel name


  f_out=this->getConduitEntrance<char>("f_out");
  f_in=this->getConduitExit<char>("f_in");
  this->simulate();
}


/************************************LightDistributedAtmosphereRequester*************************************************/
LightDistributedAtmosphereRequester::LightDistributedAtmosphereRequester(double dx, double dt):
DistributedAtmosphereSubModel(dx,dt), AsynchronousRemoteMpiKernel(){

}

LightDistributedAtmosphereRequester::~LightDistributedAtmosphereRequester(){}

int LightDistributedAtmosphereRequester::getSubmodel_ID(){
  return 3;
}

void LightDistributedAtmosphereRequester::send(vector<char> const& toSend){
  this->iSend(f_out, (char*)toSend.data(), (int) toSend.size());
  this->waitSend();
}

void LightDistributedAtmosphereRequester::receive(vector<char> & data){
  int count;
  char * buffer= this->iRecv(f_in, count);
  this->waitReceive();
  std::vector<char> vect(buffer, buffer + count);
  vect.swap(data);
  delete[] buffer;
}

void LightDistributedAtmosphereRequester::end(){}
//<--redefine
bool LightDistributedAtmosphereRequester::isConverged(){
    return AtmosphereSubModel::isConverged();
}

void LightDistributedAtmosphereRequester::F_init()
{
    DistributedAtmosphereSubModel::F_init();
    string kname= this->getName()+".log";
    // plb_ofstream  pcout(kname.c_str());
    this->logger= std::move(std::make_shared<plb_ofstream>(this->manager->getRank(), kname.c_str()));
}


void LightDistributedAtmosphereRequester::mainLoop(){
  this->manager= this->mpiManager.get();
  f_out=this->getConduitEntrance<char>("f_out");
  f_in=this->getConduitExit<char>("f_in");
  this->simulate();
}


void LightDistributedAtmosphereRequester::U(){


  bool debug=true;
  //== remove the particles from the domain ==

  #ifdef PROFILING
  this->profiler.registerBeginGatherLevelStamp();
  #endif

  //
  #ifdef PROFILING
  this->profiler.registerEndGatherLevelStamp();
  #endif

  //
  #ifdef PROFILING
  double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
  #endif

  //========== prepare  a request ==============
  vector<char> buffer;
  stringstream ss;
  if(this->manager->isMainProcessor()){
    Coord3D craterPosition(1.+cpt,2.+cpt,3.+cpt);
    Coord3D domainSize (4.0, 5.0, 6.0);
    Coord3D domainPosition(7.0, 8.0, 9.0);
    CoordGeo craterGeoPosition(10.0,11.0);
    string date="08.03.2016";
    double time= 1243.2513;
    AtmosphereDataRequest atm(craterPosition, domainSize, domainPosition, craterGeoPosition, date,  time);
    vector<char> tmp;
    bool isconv= AtmosphereSubModel::isConverged();
    buffer.resize(sizeof(isconv));
    memcpy(&buffer[0], &isconv, sizeof(isconv) );
    atm.serialize(tmp);
    buffer.insert(buffer.end(), tmp.begin(), tmp.end());

    ss<<"---------------"<<"\n";
    ss<<"request["<<this->cpt<<"]: "<<"\n";
    ss<<"->\n"<<atm.toString()<<"\n";
    //========== send request ==============
    this->send(buffer);
  }




  #ifdef PROFILING
  double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
  double WTimeReceiveBegin=  WTimeSendEnd;
  #endif
  //========== get response ==============


  if(this->manager->isMainProcessor()){
    vector<char> data;
    this->receive(data);
    AtmosphereDataResponse atmResponse;
    atmResponse.unSerialize(data);

    ss<<"<-\n"<<atmResponse.toString()<<"\n";
    ss<<"---------------"<<"\n";
    logger->getWriter()<<ss.str();
  }


  #ifdef PROFILING
  double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
  this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
  this->profiler.registerBeginUpdateLevelStamp();
  #endif
  //

  #ifdef PROFILING
  this->profiler.registerEndUpdateLevelStamp();
  #endif

  //
  #ifdef PROFILING
  this->profiler.commitStampRecording();
  if(this->cpt% nbrItrToSaveStatsInFile==0 && this->cpt>0){
    this->profiler.saveToFile();
  }
  #endif
  if( debug && manager->getRank() == 0 ) logger->getWriter() <<  "====================" << "\n";
}
