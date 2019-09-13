/*
Particles in Advection Field (PIAF)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "../../include/Tools/FileManager/FileManager.hpp"

namespace piaf{

FileManager::FileManager() {
    // TODO Auto-generated constructor stub

}

FileManager::~FileManager() {
    // TODO Auto-generated destructor stub
}


void FileManager::write( std::string filename, std::vector< std::vector<double> > & domainParticlesNumber, Int3 dim ){
    hid_t file_id;
    hsize_t dims[3] = { hsize_t( dim.x_ ), hsize_t( dim.y_ ), hsize_t( dim.z_ ) };
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for(uint fId = 0; fId<domainParticlesNumber.size(); fId++){
        H5LTmake_dataset(file_id,("f"+boost::lexical_cast<std::string>(fId)).c_str(), 3, dims, H5T_NATIVE_DOUBLE, domainParticlesNumber.at(fId).data() );
    }
    H5Fclose(file_id);
}

// doesn't work as expected (data not ordered correctly)

void FileManager::write( std::string filename, AbstractEulerianDomain *dom ){
    // storage in memory and hdf5 not compatible so this output is not correct
    //throw "method not usable";
    auto domainParticlesNumber = dom->getDomainParticlesNumber();
    hid_t file_id;
    hsize_t dims[3] = { hsize_t( dom->getSize().x_ ), hsize_t( dom->getSize().y_ ), hsize_t( dom->getSize().z_ ) };
    //hsize_t dims[3] = { hsize_t( dom->getSize().z_ ), hsize_t( dom->getSize().y_ ), hsize_t( dom->getSize().x_ ) };
    //hsize_t dims[3] = { hsize_t( dom->getSize().y_ ), hsize_t( dom->getSize().x_ ), hsize_t( dom->getSize().z_ ) };
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for(uint fId = 0; fId<domainParticlesNumber.size(); fId++){
        H5LTmake_dataset(file_id,("f"+boost::lexical_cast<std::string>(fId)).c_str(), 3, dims, H5T_NATIVE_DOUBLE, domainParticlesNumber.at(fId).data() );
    }
    double dx = dom->getDx();
    H5LTset_attribute_double( file_id, "/", "dx", &dx, 1);
    H5Fclose(file_id);
}


void FileManager::write( std::string filename, GridTerrain *terrain ){

    std::vector<std::vector<int> > deposition = terrain->getParticleDeposition();

    Int2 size = terrain->getDiscreteSize();
    std::vector< std::shared_ptr< GenericParticleFamily > > families = terrain->getParticleFamilies();
    hid_t file_id;
    hsize_t dims[2] = { hsize_t( size.x_ ), hsize_t( size.y_ ) };
    int *data = new int[size.x_*size.y_];
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    for(auto fam: families){
        GenericParticleFamily *family = &(*fam);
        for(int y=0; y<size.y_; y++){
            for(int x=0; x<size.x_; x++){
                int npart = (deposition[terrain->index2d({x,y})])[family->familyId_];
                //data[terrain->index2d({x,y})] = npart;
                // stored as y,x in order to read x,y from data
                data[y+x*terrain->getYSize()] = npart;
            }
        }
        H5LTmake_dataset(file_id,("f"+boost::lexical_cast<std::string>(fam->familyId_)).c_str(),2,dims,H5T_NATIVE_INT,data);

    }

    H5Fclose(file_id);
    delete[] data;

}

int indexId( const std::vector<int>& ids, int id ){
    for( uint i=0; i<ids.size(); i++ ){
        if( ids[i]==id ) return i;
    }
    return -1;
}

void FileManager::writeVtk( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){
    std::vector< Particle > particles;
    std::string currentFilename;

    for( uint i=0; i<std::get<0>(trackedParticles).size(); i++ ){
        particles = std::get<1>( trackedParticles )[i];
        currentFilename = std::to_string(vtkTrackFileCounter_)+"_"+filename+".vtk";
        std::ofstream vtkfile;
        vtkfile.open(currentFilename);
        vtkfile << "# vtk DataFile Version 3.0" << std::endl;
        vtkfile << "Particle file created by piaf" << std::endl;
        vtkfile << "ASCII" << std::endl;
        vtkfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        vtkfile << "POINTS " << particles.size() << " double" << std::endl;
        for(auto& part : particles){
            vtkfile << part.displacement_.x_ << " " << part.displacement_.y_ << " " << part.displacement_.z_  << std::endl;
        }
        vtkfile << "POINT_DATA " << particles.size() << std::endl;
        vtkfile << "SCALARS FamilyId int" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for(auto& part : particles){
            vtkfile << part.familyId_ << std::endl;
        }
        vtkfile << "SCALARS Scaling double" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for(auto& part : particles){
            vtkfile << part.scaling_ << std::endl;
        }
        vtkfile << "SCALARS ParticleId int" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for(auto& part : particles){
            vtkfile << part.getParticleId() << std::endl;
        }
        vtkfile.close();
        vtkTrackFileCounter_++;
    }
}

void FileManager::writeParticlesCsv( std::string filename,  std::vector< piaf::Particle >& particles ){
    std::string currentFilename;

    currentFilename = std::to_string(csvTrackFileCounter_)+"_"+filename+".csv";
    std::ofstream csvfile;
    //std::cout << "writing file " << currentFilename << std::endl;
    csvfile.open(currentFilename);
    csvfile << "x,y,z,particleId,class,scaling" << std::endl;
    for(auto& part : particles){
         csvfile << part.displacement_.x_ << "," << part.displacement_.y_ << "," << part.displacement_.z_ << "," << part.familyId_ << "," << part.scaling_ << " " << part.getParticleId() << std::endl;
    }
    csvfile.close();
    csvTrackFileCounter_++;
  
}

void FileManager::writeCsv( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){
    std::vector< Particle > particles;
    std::string currentFilename;

    for( uint i=0; i<std::get<0>(trackedParticles).size(); i++ ){
        particles = std::get<1>( trackedParticles )[i];
        currentFilename = std::to_string(csvTrackFileCounter_)+"_"+filename+".csv";
        std::ofstream csvfile;
        //std::cout << "writing file " << currentFilename << std::endl;
        csvfile.open(currentFilename);
        csvfile << "x,y,z,particleId,class,scaling" << std::endl;
        for(auto& part : particles){
            csvfile << part.displacement_.x_ << "," << part.displacement_.y_ << "," << part.displacement_.z_ << "," << part.getParticleId() << "," << part.familyId_ << "," << part.scaling_ << std::endl;
        }
        csvfile.close();
        csvTrackFileCounter_++;
    }
}


void FileManager::writeHdf5( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){

    double time;
    std::vector< Particle > particles;
    hid_t file_id, group_id;
    hsize_t dims[2];

    if(std::get<0>(trackedParticles).size() > 0)
        file_id = H5Fcreate((std::to_string(hdf5TrackFileCounter_)+"_"+filename+".hdf5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    for( uint i=0; i<std::get<0>(trackedParticles).size(); i++ ){
        time = std::get<0>( trackedParticles )[i];
        particles = std::get<1>( trackedParticles )[i];

        std::vector<int> ids;
        std::vector< std::vector < Double3 > > byIdPartCoord;
        for( auto& part : particles ){
            int index = indexId( ids, part.familyId_ );
            if( index==-1 ){
                ids.push_back( part.familyId_ );
                byIdPartCoord.push_back( std::vector< Double3 >() );
                byIdPartCoord[ ids.size() - 1 ].push_back( Double3{ part.displacement_.x_, part.displacement_.y_, part.displacement_.z_ } );
                byIdPartCoord[ ids.size() - 1 ].push_back( Double3{ part.lagrangianSpeed_.x_-part.diffusionSpeed_.x_, part.lagrangianSpeed_.y_-part.diffusionSpeed_.y_, part.lagrangianSpeed_.z_-part.diffusionSpeed_.z_ } );
                byIdPartCoord[ ids.size() - 1 ].push_back( Double3{ part.diffusionSpeed_.x_, part.diffusionSpeed_.y_, part.diffusionSpeed_.z_ } );
            }
            else{
                byIdPartCoord[index].push_back( Double3{ part.displacement_.x_, part.displacement_.y_, part.displacement_.z_ } );
                byIdPartCoord[index].push_back( Double3{ part.lagrangianSpeed_.x_-part.diffusionSpeed_.x_, part.lagrangianSpeed_.y_-part.diffusionSpeed_.y_, part.lagrangianSpeed_.z_-part.diffusionSpeed_.z_ } );
                byIdPartCoord[index].push_back( Double3{ part.diffusionSpeed_.x_, part.diffusionSpeed_.y_, part.diffusionSpeed_.z_ } );
            }
        }

        std::string groupName = std::to_string(time);
        group_id = H5Gcreate2(file_id, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for( uint i=0; i<ids.size(); i++ ){
            //dims[1]=3;
            dims[1]=9;
            dims[0]=byIdPartCoord[i].size()/3;
            H5LTmake_dataset(file_id,(groupName+"/f"+std::to_string(ids[i])).c_str(),2,dims,H5T_NATIVE_DOUBLE,byIdPartCoord[i].data());
        }

        H5Gclose(group_id);
    }

    if(std::get<0>(trackedParticles).size() > 0){
        H5Fclose(file_id);
        hdf5TrackFileCounter_++;
    }
}

void FileManager::writeVelocityField( std::string filename, VelocityField& vc, double t, double dt ){
    auto velocitiesX = vc.getXHdf5Ordering(t, dt);
    auto velocitiesY = vc.getYHdf5Ordering(t, dt);
    auto velocitiesZ = vc.getZHdf5Ordering(t, dt);

    auto file_id_x = H5Fcreate((std::string("x_")+filename).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    auto file_id_y = H5Fcreate((std::string("y_")+filename).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    auto file_id_z = H5Fcreate((std::string("z_")+filename).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    auto group_id_x = H5Gcreate2(file_id_x, "velocities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto group_id_y = H5Gcreate2(file_id_y, "velocities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto group_id_z = H5Gcreate2(file_id_z, "velocities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    uint sx = vc.getDiscreteSize().x_;
    std::cout << "x = " << sx << std::endl;
    uint sy = vc.getDiscreteSize().y_;
    std::cout << "y = " << sy << std::endl;
    uint sz = vc.getDiscreteSize().z_;
    std::cout << "z = " << sz << std::endl;

    hsize_t dimsDataset[] = {sx, sy, sz};

    H5LTmake_dataset(file_id_x,"velocities/x",3,dimsDataset,H5T_NATIVE_DOUBLE,velocitiesX.get());
    H5LTmake_dataset(file_id_y,"velocities/y",3,dimsDataset,H5T_NATIVE_DOUBLE,velocitiesY.get());
    H5LTmake_dataset(file_id_z,"velocities/z",3,dimsDataset,H5T_NATIVE_DOUBLE,velocitiesZ.get());

    H5Gclose(group_id_x);
    H5Gclose(group_id_y);
    H5Gclose(group_id_z);

    H5Fclose(file_id_x);
    H5Fclose(file_id_y);
    H5Fclose(file_id_z);

}

} // namespace piaf
