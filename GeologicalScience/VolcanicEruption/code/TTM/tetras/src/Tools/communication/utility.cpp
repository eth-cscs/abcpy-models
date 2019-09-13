#include "utility.h"

namespace unige_pasc{


/***************************************** util ***********************************************/
MapperUtil::MapperUtil()
{
}

template <typename T>
string  MapperUtil::stringify(T & number)
{
    std::ostringstream o;
    if (!(o << number))
    {
        cout << "error with stringify" << endl;
        exit(1);
    }

    return o.str();
}


template <typename T>
void MapperUtil::unstringify(std::string str, T& value) {
    std::stringstream valueStr(str);
    T tmp = T();
    if (!(valueStr>>tmp)) {
        cout<<(std::string("Cannot read value  ") + str);
        exit(1);
    }
    value = tmp;
}


void MapperUtil::unstringify_pluint(std::string str, pluint& value) {
    std::stringstream valueStr(str);
    pluint tmp;
    if (!(valueStr>>tmp)) {
        cout<<(std::string("Cannot read value  ") + str);
        exit(1);
    }
    value = tmp;
}

/*
#ifdef USE_PARALLEL_MPI
void util::broadCast(std::string& message, int rank, int root, MPI_Comm comm)
{
    int length = (int) message.size();
    // bCast(&length, 1, root, comm);
    MPI_Bcast(static_cast<void*>(&length), 1, MPI_INT, root, comm);
    char* buffer = new char[length+1];
    if (rank==root) {
        std::copy(message.c_str(), message.c_str()+length+1, buffer);
    }
    //bCast(buffer, length+1, root, comm);
    MPI_Bcast(static_cast<void*>(buffer), length+1, MPI_CHAR, root, comm);
    if (rank!=root) {
        message = buffer;
    }
    delete [] buffer;
}
#endif
*/
void* MapperUtil::memcpyEndian(void* dest, const void* src, size_t count) {
#ifdef BIG_ENDIAN_SERIALIZER
    char* dst8 = (char*)dest;
    char* src8 = (char*)src;
    --src8;
    dst8=dst8+count;

    while (count--) {
        *--dst8= *++src8;
    }
    return dest;
#else
    return memcpy(dest,src,count);
#endif
}

/*time_t util::to_time_t(const boost::gregorian::date& date ){
    using namespace boost::posix_time;
    static ptime epoch(boost::gregorian::date(1970, 1, 1));
    time_duration::sec_type secs = (ptime(date,seconds(0)) - epoch).total_seconds();
    return time_t(secs);
}*/

time_t MapperUtil::to_time_t(boost::posix_time::ptime const & t)
{
    using namespace boost::posix_time;
    ptime epoch(boost::gregorian::date(1970,1,1));
    time_duration::sec_type x = (t - epoch).total_seconds();

    return time_t(x);
}

void MapperUtil::serialize(std::vector< BoundaryParticle >  & boundaryParticles, bool isConverged, vector<char> & toSend, int cpt){
    try{
        size_t particleSize= 0;
        if(! boundaryParticles.empty()){
            BoundaryParticle & p= boundaryParticles[0];
            particleSize= sizeof(p.getFamilyId())
                    +3*(sizeof(p.getDisplacement().x_)+sizeof(p.getDisplacement().y_)+ sizeof(p.getDisplacement().z_))
                    +sizeof(p.getTimeStamp());
        }
        int isConvergedSim= (isConverged)?1:0;
        size_t vectParticlesSize = ( particleSize * boundaryParticles.size() )+sizeof(isConvergedSim);
        toSend.resize(vectParticlesSize) ;
        size_t iData=0;
        memcpy( (void *) &toSend[iData], (const void *) &isConvergedSim, sizeof(isConvergedSim)); // "Serialize"
        iData+=sizeof(isConvergedSim);
        if(cpt > 0){
            //stringstream ss;
            //ss<<"["<<cpt<<"]: Plume isconverged="<< isConvergedSim<< " size type"<<sizeof(isConvergedSim)<<"\n";
            //cout<<ss.str();

        }


        // code that could cause exception

        if(! boundaryParticles.empty()){
            for (vector<BoundaryParticle>::iterator it= boundaryParticles.begin(); it != boundaryParticles.end(); it++){
                memcopyIntoVector( *it, toSend, iData);
            }
            //memcpy( (void *) &toSend[iData], (const void *)&boundaryParticles[0], sizeof(BoundaryParticle)*boundaryParticles.size()); // "Serialize"
        }
    }
    catch (const std::exception &exc)
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr <<"#######> serialisation: "<< exc.what()<<endl;
    }

    /* boost::asio::streambuf write_buffer;
    std::ostream streamToSend(&write_buffer);
    boost::archive::text_oarchive oa(streamToSend);
    // add the information about the serialization
    // write class instance to archive
    SerializerParticles boundary(boundaryParticles, isConverged);
    oa << boundary;
    pluint buffSize=   write_buffer.size();
    toSend.resize(buffSize);
    write_buffer.sgetn (&toSend[0],buffSize);*/

}

void MapperUtil::deserialize(std::vector< BoundaryParticle > & newBoundaryParticles, bool & isConverged,  vector<char> & vect, int cpt){
    try{
        int isConvergedSim=0;
        assert(vect.size() >= sizeof(isConvergedSim));
        size_t pBytesSize = vect.size() - sizeof(isConvergedSim);


        BoundaryParticle p;
        size_t particleSize= sizeof(p.getFamilyId())
                +3*(sizeof(p.getDisplacement().x_)+sizeof(p.getDisplacement().y_)+ sizeof(p.getDisplacement().z_))
                +sizeof(p.getTimeStamp());


        assert ( (pBytesSize % particleSize) == 0);
        size_t vectParticlesSize= pBytesSize /particleSize;
        size_t iData=0;
        memcpy((void*)&isConvergedSim, (const void *)&vect[iData], sizeof(isConvergedSim)); // "Serialize"
        iData+=sizeof(isConvergedSim);
        isConverged = (isConvergedSim==0)? false: true;
        if(cpt>0){
            stringstream ss;
            ss<<"["<<cpt<<"]: G isConverged="<< isConverged<<" size type"<<sizeof(isConverged)<<"\n";
            cout<<ss.str();
        }

        if(vectParticlesSize>0){
            //newBoundaryParticles.resize(vectParticlesSize);
            for (size_t i=0; i < vectParticlesSize; i++){
                BoundaryParticle p;
                memcopyIntoParticle( p, vect, iData);
                newBoundaryParticles.push_back(p);
            }
            //memcpy((void*) &newBoundaryParticles[0], (const void*) &vect[iData], sizeof(BoundaryParticle)* vectParticlesSize); // "Serialize"
        }
    }
    catch (const std::exception &exc)
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr <<"#######> deserialisation: "<< exc.what()<<endl;
    }

    /*membuf write_reader(&vect[0],&vect[0]+vect.size());
    std::istream streamReceived(&write_reader);
    boost::archive::text_iarchive ia(streamReceived);
    // read class state from archive
    SerializerParticles boundaryRec(newBoundaryParticles, isConverged);
    ia >> boundaryRec;*/
}

void MapperUtil::memcopyIntoVector( BoundaryParticle & p, vector<char> & buffer, size_t & iData){

    /* int familyId_ = p.getFamilyId();
    Double3 displacement_ = p.getDisplacement();
    Double3 lagrangianSpeed_ = p.getLagrangianSpeed();
    Double3 diffusionSpeed_ = p.getDiffusionSpeed();
    double timeStamp_ = p.getTimeStamp();*/

    memcpy((void*) &buffer[iData], (const void*) &p.familyId_, sizeof(p.familyId_)); // "Serialize"
    iData += sizeof(p.familyId_);

    memcpy((void*) &buffer[iData], (const void*) & (p.displacement_.x_), sizeof(p.displacement_.x_)); // "Serialize"
    iData += sizeof(p.displacement_.x_);
    memcpy((void*) &buffer[iData], (const void*) & (p.displacement_.y_), sizeof(p.displacement_.y_)); // "Serialize"
    iData += sizeof(p.displacement_.y_);
    memcpy((void*) &buffer[iData], (const void*) & (p.displacement_.z_), sizeof(p.displacement_.z_)); // "Serialize"
    iData += sizeof(p.displacement_.z_);


    memcpy((void*) &buffer[iData], (const void*) & (p.lagrangianSpeed_.x_), sizeof(p.lagrangianSpeed_.x_)); // "Serialize"
    iData += sizeof(p.lagrangianSpeed_.x_);
    memcpy((void*) &buffer[iData], (const void*) & (p.lagrangianSpeed_.y_), sizeof(p.lagrangianSpeed_.y_)); // "Serialize"
    iData += sizeof(p.lagrangianSpeed_.y_);
    memcpy((void*) &buffer[iData], (const void*) & (p.lagrangianSpeed_.z_), sizeof(p.lagrangianSpeed_.z_)); // "Serialize"
    iData += sizeof(p.lagrangianSpeed_.z_);

    memcpy((void*) &buffer[iData], (const void*) & (p.diffusionSpeed_.x_), sizeof(p.diffusionSpeed_.x_)); // "Serialize"
    iData += sizeof(p.diffusionSpeed_.x_);
    memcpy((void*) &buffer[iData], (const void*) & (p.diffusionSpeed_.y_), sizeof(p.diffusionSpeed_.y_)); // "Serialize"
    iData += sizeof(p.diffusionSpeed_.y_);
    memcpy((void*) &buffer[iData], (const void*) & (p.diffusionSpeed_.z_), sizeof(p.diffusionSpeed_.z_)); // "Serialize"
    iData += sizeof(p.diffusionSpeed_.z_);

    double timeStamp_ = p.getTimeStamp();
    memcpy((void*) &buffer[iData], (const void*) & (timeStamp_), sizeof(timeStamp_)); // "Serialize"
    iData+=sizeof(timeStamp_);
}


void MapperUtil::memcopyIntoParticle(BoundaryParticle & p, vector<char> & buffer, size_t & iData){

    int familyId_ ;
    Double3 displacement_ ;
    Double3 lagrangianSpeed_ ;
    Double3 diffusionSpeed_ ;
    double timeStamp_;

    double x,y,z;
    memcpy((void*) &familyId_ , (const void*) &buffer[iData], sizeof(familyId_));
    iData += sizeof(familyId_);

    memcpy((void*) &x , (const void*) &buffer[iData], sizeof(x));
    iData +=sizeof(x);
    memcpy((void*) &y , (const void*) &buffer[iData], sizeof(y));
    iData +=sizeof(y);
    memcpy((void*) &z , (const void*) &buffer[iData], sizeof(z)); // "Serialize"
    iData +=sizeof(z);
    displacement_.x_=x;
    displacement_.y_=y;
    displacement_.z_=z;

    memcpy((void*) &x , (const void*) &buffer[iData], sizeof(x));
    iData +=sizeof(x);
    memcpy((void*) &y , (const void*) &buffer[iData], sizeof(y));
    iData +=sizeof(y);
    memcpy((void*) &z , (const void*) &buffer[iData], sizeof(z)); // "Serialize"
    iData +=sizeof(z);
    lagrangianSpeed_.x_=x;
    lagrangianSpeed_.y_=y;
    lagrangianSpeed_.z_=z;

    memcpy((void*) &x , (const void*) &buffer[iData], sizeof(x));
    iData +=sizeof(x);
    memcpy((void*) &y , (const void*) &buffer[iData], sizeof(y));
    iData +=sizeof(y);
    memcpy((void*) &z , (const void*) &buffer[iData], sizeof(z)); // "Serialize"
    iData +=sizeof(z);
    diffusionSpeed_.x_=x;
    diffusionSpeed_.y_=y;
    diffusionSpeed_.z_=z;

    memcpy((void*) &timeStamp_, (const void*) &buffer[iData], sizeof(timeStamp_)); // "Serialize"
    iData+=sizeof(timeStamp_);

    p.setFamilyId(familyId_);
    p.setDisplacement(displacement_);
    p.setLagrangianSpeed(lagrangianSpeed_);
    p.setDiffusionSpeed(diffusionSpeed_);
    p.setTimeStamp(timeStamp_);
}


void MapperUtil::serialize(bool isConverged, vector<char> & toSend){
    memcpy((void*) &toSend[0], (const void*) &isConverged, sizeof(isConverged)); // "Serialize"
}
void MapperUtil::deserialize(bool &isConverged, vector<char> const& toSend){
    memcpy((void*) &isConverged, (const void*) &toSend[0], sizeof(isConverged));
}

}


