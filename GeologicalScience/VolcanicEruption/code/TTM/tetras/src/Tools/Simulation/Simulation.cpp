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


#include "Simulation.hpp"

LogLevel GLOBAL_LOG_LEVEL;

double initialTime;

/**************************************** Plume ************************************************/

Plume::Plume(params p, shared_ptr<MPITopology> topology, double U0, double L0 ):
    p_				( p ),
    topology_			( topology ),
    dom_				( NULL ),
    repository_			( NULL ),
    terrain_			( NULL ),
    currentTime_		(0.0),
    isProcess(true),
    frequenceCheckConvergence(100),
    //frequenceCheckConvergence(1),
    previous_particules_deposees(0),
    injectedParticles(0),
    sum_injectedParticles(0),
    sum_outsideParticles(0),
    U0_(U0),
    L0_(L0)
{
    mt_ = std::mt19937(rd_());

    InitSimulation = std::bind(&Plume::InitSimulationReal, this);

    std::cout << "Value of injectAtCrater in Plume constructor FIX THIS PARAMETER, NOT IN USE : " << p_.injectAtCrater_ << std::endl;

    // if(p_.injectAtCrater_){
        InjectParticles = std::bind(&Plume::InjectParticlesAtCrater, this);
    // }
    // else{
        InjectParticles = std::bind(&Plume::InjectParticlesAlongPlume, this);
    // }
}

Plume::Plume(int argc , char ** argv, MPI_Comm globalCommunicator, double U0, double L0):
    argc(argc),
    argv(argv),
    globalCommunicator(globalCommunicator),
    dom_				( NULL ),
    repository_			( NULL ),
    terrain_			( NULL ),
    currentTime_		(0.0),
    isProcess(true),
    frequenceCheckConvergence(100),
    //frequenceCheckConvergence(1),
    previous_particules_deposees(0),
    injectedParticles(0),
    sum_injectedParticles(0),
    sum_outsideParticles(0),
    U0_(U0),
    L0_(L0)
{
    mt_ = std::mt19937(rd_());

    InitSimulation = std::bind(&Plume::InitSimulationReal, this);

    std::cout << "Value of injectAtCrater in Plume constructor FIX THIS PARAMETER, NOT IN USE : " << p_.injectAtCrater_ << std::endl;

    // if(p_.injectAtCrater_){
        InjectParticles = std::bind(&Plume::InjectParticlesAtCrater, this);
    // }
    // else{
    //     InjectParticles = std::bind(&Plume::InjectParticlesAlongPlume, this);
    // }
}



Plume::~Plume() {
    delete dom_;
    delete repository_;
    delete terrain_;
    //this->logger->getWriter().close();
}

MPITopology* Plume::getMpiTopology(){
    return this->topology_.get();
}

void Plume::createVolcanTopology(){
    if(! topology_){
        topology_= TolopogyCreator::getInstance().createTopology(p_, argc, argv, globalCommunicator);
    }
}

void Plume::initLogger(string name)
{
    //------------------------init logger------------------------------
    //assert(argv[0]);
    //string kname(argv[0]);// 1st argument is the kernel name
    string kname= name+".log";
    // plb_ofstream  pcout(kname.c_str());
    this->logger= std::move(std::make_shared<plb_ofstream>(topology_->getRank(), kname.c_str()));
    //this->logger= std::move(std::make_shared<PlumeLogger>(kname, topology_->getRank(), pcout));
    //this->logger->getWriter() <<"hello" << topology_->getSize()<<"\n";
    //---------------------------------------------------------
}
std::vector< std::shared_ptr< SpeedFunctor > > Plume::createWoodsSpeeds( EruptionParameters ep, bool isoDif ){

    auto cp = ep.getColumnProfile();
    auto up = cp->getUCentralProfile();
    auto lp = cp->getLProfile();

    std::vector< std::shared_ptr< SpeedFunctor > > speeds;

    /* creating speeds for Woods model */

    /***********************************************************/
    /* WIND MODEL 1 */
    //std::shared_ptr<WM1> wm1(new WM1(ep.getColumnProfile(), ep.getVentPosition()));
    //wm1->setAtmosphere( ep.getAtmosphere().get() );
    //wm1->setMaximumWindSpeed( ep.getMaxWindSpeed() );
    //wm1->setWindAngle( ep.getWindDirection() );
    //wm1->setSpeedType(EULERIAN_STATIC);
    /***********************************************************/

    /***********************************************************/
    /* WIND MODEL 2 */
    std::shared_ptr<WM2> wm2(new WM2(ep.getColumnProfile(), ep.getVentPosition()));
    wm2->setAtmosphere( ep.getAtmosphere().get() );
    wm2->setMaximumWindSpeed( ep.getMaxWindSpeed() );
    wm2->setWindAngle( ep.getWindDirection() );
    wm2->setSpeedType(EULERIAN_STATIC);
    /***********************************************************/

    /***********************************************************/
    /* DIFFUSION */
    std::shared_ptr<DiffusionWoods> diffusion(new DiffusionWoods( ep.getColumnProfile(), ep.getVentPosition() ) );
    diffusion->dAtm_ = ep.getAtmosphereHorizontalDiffusion();
    diffusion->dCol_ = ep.getColumnHorizontalDiffusion();
    diffusion->setIsDiffusion(true);
    diffusion->setSpeedType(LAGRANGIAN_DYNAMIC);

    std::shared_ptr<DirectionalDiffusionWoods> directionalDiffusion(new DirectionalDiffusionWoods( ep.getColumnProfile(), ep.getVentPosition() ) );
    directionalDiffusion->dVAtm_ = ep.getAtmosphereVerticalDiffusion();
    directionalDiffusion->dHAtm_ = ep.getAtmosphereHorizontalDiffusion();
    directionalDiffusion->dVCol_ = ep.getColumnVerticalDiffusion();
    directionalDiffusion->dHCol_ = ep.getColumnHorizontalDiffusion();
    directionalDiffusion->setIsDiffusion(true);
    directionalDiffusion->setSpeedType(LAGRANGIAN_DYNAMIC);
    /***********************************************************/

    /***********************************************************/
    /* SEDIMENTATION SPEED */
    std::shared_ptr<SedimentationSpeed> sedimentation(new SedimentationSpeed());
    sedimentation->setSpeedType(EULERIAN_STATIC);
    /***********************************************************/

    /***********************************************************/
    /* UMBRELLA CLOUD SPEED */
    std::shared_ptr<UmbrellaCloudSpeed> umbrellaCloud(new UmbrellaCloudSpeed());
    // no wind -> e = 0.0
    umbrellaCloud->setE( ep.getEccentricity() );
    umbrellaCloud->setFocusSource( ep.getFocusSource() );
    umbrellaCloud->setHb( ep.getHb() );
    umbrellaCloud->setHt( ep.getHt() );
    umbrellaCloud->setTropopauseHeight( ep.getAtmosphere()->getTropopauseHeight() );
    umbrellaCloud->setSpeedType( LAGRANGIAN_DYNAMIC );
    umbrellaCloud->setDx( dom_->getDx() );
    //umbrellaCloud->setDt( p_.dt_ );
    /***********************************************************/

    /***********************************************************/
    /* VOLCANIC COLUMN SPEED */
    std::shared_ptr<VolcanicColumnSpeed> volcanicColumnSpeed(new VolcanicColumnSpeed( ep.getColumnProfile(), ep.getVentPosition() ));
    volcanicColumnSpeed->setSpeedType(EULERIAN_STATIC);
    /***********************************************************/

    speeds.push_back( wm2 );
    speeds.push_back( sedimentation );
    speeds.push_back( umbrellaCloud );
    speeds.push_back( volcanicColumnSpeed );

    // diffusion
    if( !p_.nodif_ )
    {
        if( isoDif ) {
            if( topology_->getRank() == 0 ) this->logger->getWriter() << "using isotropic diffusion" << "\n";
            speeds.push_back( diffusion );
        }
        else{
            if( topology_->getRank() == 0 ) this->logger->getWriter() << "using anisotropic diffusion" << "\n";
            speeds.push_back( directionalDiffusion );
        }
    }
    else{
        this->logger->getWriter() << "using no diffusion" << "\n";
    }

    return speeds;

}

std::vector< std::shared_ptr< SpeedFunctor > > Plume::createWimSpeeds( EruptionParameters ep ){

    std::vector< std::shared_ptr< SpeedFunctor > > speeds;

    /***********************************************************/
    /* SEDIMENTATION SPEED */
    std::shared_ptr<SedimentationSpeed> sedimentation(new SedimentationSpeed());
    sedimentation->setSpeedType(EULERIAN_STATIC);
    /***********************************************************/

    std::shared_ptr< WindWimPlumeSpeed > wwp(
                new WindWimPlumeSpeed( ep.getPlume(), ep.getWind(), ep.getVentPosition(),
                                       ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                );
    wwp->setSpeedType(LAGRANGIAN_DYNAMIC);

    std::shared_ptr< DiffusionWim > difwim(
                new DiffusionWim( ep.getPlume(), ep.getVentPosition(),
                                  ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                );
    difwim->setSpeedType(LAGRANGIAN_DYNAMIC);
    difwim->setIsDiffusion(true);

    std::shared_ptr< DiffusionWimUmbrellaCloud > difwimUmbCloud(
                new DiffusionWimUmbrellaCloud( ep.getPlume(), ep.getVentPosition(),
                                  ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion(),
                                  ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_,
                                  ep.getPlume().getZ(0.0)->back() + ep.getVentPosition().z_ )
                );
    difwim->setSpeedType(LAGRANGIAN_DYNAMIC);
    difwim->setIsDiffusion(true);

    std::shared_ptr< DiffusionWimGaussian > difwimgaussian(
                new DiffusionWimGaussian( ep.getPlume(), ep.getVentPosition(),
                                          ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                );
    difwimgaussian->setSpeedType(LAGRANGIAN_DYNAMIC);
    difwimgaussian->setIsDiffusion(true);

    speeds.push_back( sedimentation );
    speeds.push_back( wwp );
    //speeds.push_back( difwim );

    std::shared_ptr< WindWimPlumeSpeedGaussian > wwpg(
                new WindWimPlumeSpeedGaussian( ep.getPlume(), ep.getWind(), ep.getVentPosition(),
                                               ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                );
    wwpg->setSpeedType(LAGRANGIAN_DYNAMIC);
    speeds.push_back( wwpg );

    std::shared_ptr< WindWimPlumeSpeedGaussianAboveHb > wwpgah(
                new WindWimPlumeSpeedGaussianAboveHb( ep.getPlume(), ep.getWind(), ep.getVentPosition(),
                                               ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion(), ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_ )
                );
    wwpgah->setSpeedType(LAGRANGIAN_DYNAMIC);
    speeds.push_back( wwpgah );

    std::shared_ptr< WindWimPlumeSpeedGaussianSedimentation > wwpgs(
                new WindWimPlumeSpeedGaussianSedimentation( ep.getPlume(), ep.getWind(), ep.getVentPosition(),
                                                            ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                );
    wwpgs->setSpeedType(LAGRANGIAN_DYNAMIC);
    speeds.push_back( wwpgs );

    auto plume = ep.getPlume();
    double hb = plume.getZ(0.0)->back() * 0.7;
    double ht = plume.getZ(0.0)->back();
    int indexHb = std::get<0>(plume.indexClosestPoint({ep.getVentPosition().x_, ep.getVentPosition().y_, hb} ,0.0 ,false));

    FakeEruptionParameters fakeEp;
    std::shared_ptr< UmbrellaCloudRevised > umbrellacloudrevised(
        new UmbrellaCloudRevised( ep.getFocusSource(), hb+ep.getVentPosition().z_, ht+ep.getVentPosition().z_, ep.getAtmosphere()->getTropopauseHeight(),fakeEp.getVwind() , fakeEp.getTheta(),
        plume.getU(0.0)->at(indexHb), plume.getR(0.0)->at(indexHb) )
    );
    umbrellacloudrevised->setSpeedType(LAGRANGIAN_DYNAMIC);
    speeds.push_back(umbrellacloudrevised);

    /* umbrella cloud */
    //if(ep.getEccentricity() != -1.0 ){
    for( auto name : ep.getVelocitiesNames() ) {
        if( name.compare("tetras_umbrella_cloud")==0 || name.compare("tetras_umbrella_cloud_gaussian")==0 ) {
            std::shared_ptr<UmbrellaCloudSpeed> umbrellaCloud(new UmbrellaCloudSpeed());
            umbrellaCloud->setE( ep.getEccentricity() );
            umbrellaCloud->setFocusSource( ep.getFocusSource() );
            umbrellaCloud->setHb( ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_);
            umbrellaCloud->setHt( ep.getPlume().getZ(0.0)->back() + ep.getVentPosition().z_ );
            umbrellaCloud->setTropopauseHeight( ep.getAtmosphere()->getTropopauseHeight() );
            umbrellaCloud->setSpeedType( LAGRANGIAN_DYNAMIC );
            umbrellaCloud->setDx( dom_->getDx() );
            speeds.push_back( umbrellaCloud );

            std::shared_ptr<UmbrellaCloudSpeedGaussian> umbrellaCloudGaussian(new UmbrellaCloudSpeedGaussian());
            umbrellaCloudGaussian->setE( ep.getEccentricity() );
            umbrellaCloudGaussian->setFocusSource( ep.getFocusSource() );
            umbrellaCloudGaussian->setHb( ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_);
            umbrellaCloudGaussian->setHt( ep.getPlume().getZ(0.0)->back() + ep.getVentPosition().z_ );
            umbrellaCloudGaussian->setTropopauseHeight( ep.getAtmosphere()->getTropopauseHeight() );
            umbrellaCloudGaussian->setSpeedType( LAGRANGIAN_DYNAMIC );
            umbrellaCloudGaussian->setDx( dom_->getDx() );
            speeds.push_back( umbrellaCloudGaussian );

            std::shared_ptr< WindWimPlumeSpeedGaussianSedimentationThreshold > wwpgst(
                        new WindWimPlumeSpeedGaussianSedimentationThreshold( ep.getPlume(), ep.getWind(), ep.getVentPosition(),
                                                                             ep.getColumnHorizontalDiffusion(), ep.getAtmosphereHorizontalDiffusion() )
                        );
            wwpgst->setHb(ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_);
            wwpgst->setSpeedType(LAGRANGIAN_DYNAMIC);
            speeds.push_back( wwpgst );


            std::shared_ptr<WindWimPlumeSpeedHblimited> windPlumeHbLimited(
                        new WindWimPlumeSpeedHblimited( ep.getPlume(), ep.getWind(), ep.getVentPosition(), ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_ )
                        );
            speeds.push_back(windPlumeHbLimited);

            std::shared_ptr<WindWimPlumeSpeedHbMaxWidth> windPlumeHbMaxWidth(
                        new WindWimPlumeSpeedHbMaxWidth( ep.getPlume(), ep.getWind(), ep.getVentPosition(), ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_ )
                        );
            speeds.push_back(windPlumeHbMaxWidth);
        }
    }
    //}

    if( !p_.nodif_ ) {
        if( topology_->getRank() == 0 ) this->logger->getWriter() << "using isotropic diffusion" << "\n";
        speeds.push_back( difwim );
        speeds.push_back( difwimgaussian );
        speeds.push_back( difwimUmbCloud );
    }
    else{
        this->logger->getWriter() << "using no diffusion" << "\n";
    }

    std::shared_ptr< Backflow > backflow(
                new Backflow( ep.getPlume(),  ep.getVentPosition(), ep.getPlume().getZ(0.0)->back() * 0.7 + ep.getVentPosition().z_ )
                );
    backflow->setSpeedType(LAGRANGIAN_DYNAMIC);

    speeds.push_back( backflow );

    return speeds;
}

void Plume::InitCommon(){

    initialTime = MPI_Wtime();

    totalParticlePerFamily = p_.npart_;

    //std::vector<double> scalingFactors;
    for(Uint i=0; i<eruptionParameters->getParticleFamilies().size(); i++) scalingFactors.push_back( eruptionParameters->getFamilyTotalMass(i) / ( totalParticlePerFamily * eruptionParameters->getParticleFamilies()[i]->mass_ ) );

    if( typeid( *topology_ ) == typeid( MPI3DGridTopology ) ){
        Int3 topo = ((MPI3DGridTopology*)topology_.get())->getTopology();
        if( p_.blockCyclic_ == -1 ){
            //dom_ = new MPI3DBlockCyclicDomain( { eruptionParameters->getVentPosition().x_ - int( p_.domainSizeX_ / 2.0 ) , eruptionParameters->getVentPosition().y_ - int( p_.domainSizeY_ / 2.0 ), 0.0 } , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , (MPI3DGridTopology*)topology_ );
            dom_ = new MPI3DBlockCyclicDomain(  p_.origin_ , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , (MPI3DGridTopology*)topology_.get() );
            Int3 numBlockLoc = dom_->getNumBlockLocal();
            if( topology_->getRank() == 0 ) {
                this->logger->getWriter() << "topology used : MPI3DGrid with block in each dim = 5" << "\n";
                this->logger->getWriter() << "topology = " << topo.x_ << "," << topo.y_ << "," << topo.z_ << "\n";
                this->logger->getWriter() << "num block local = " << numBlockLoc.x_ << "," << numBlockLoc.y_ << "," << numBlockLoc.z_ << "\n";
            }
        }
        else{
            //dom_ = new MPI3DBlockCyclicDomain( { eruptionParameters->getVentPosition().x_ - int( p_.domainSizeX_ / 2.0 ) , eruptionParameters->getVentPosition().y_ - int( p_.domainSizeY_ / 2.0 ), 0.0 } , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , (MPI3DGridTopology*)topology_, p_.blockCyclic_ );
            dom_ = new MPI3DBlockCyclicDomain( p_.origin_ , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , (MPI3DGridTopology*)topology_.get(), p_.blockCyclic_ );
            Int3 numBlockLoc = dom_->getNumBlockLocal();
            if( topology_->getRank() == 0 ){
                this->logger->getWriter() << "topology used : MPI3DGrid with block in each dim = " << p_.blockCyclic_ << "\n";
                this->logger->getWriter() << "topology = " << topo.x_ << "," << topo.y_ << "," << topo.z_ << "\n";
                this->logger->getWriter() << "num block local = " << numBlockLoc.x_ << "," << numBlockLoc.y_ << "," << numBlockLoc.z_ << "\n";
            }
        }
    }
    else if( typeid( *topology_.get() ) == typeid( MPITopology ) ) {
        if( p_.simpleDom_ == -1 ){
            //dom_ = new MPISimpleDomain( { eruptionParameters->getVentPosition().x_ - int( p_.domainSizeX_ / 2.0 ) , eruptionParameters->getVentPosition().y_ - int( p_.domainSizeY_ / 2.0 ), 0.0 } , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , topology_ );
            dom_ = new MPISimpleDomain( p_.origin_ , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , topology_.get() );
            if( topology_->getRank() == 0 ) this->logger->getWriter() << "topology used : simple with dimension splitted : x" << "\n";
        }
        else{
            //dom_ = new MPISimpleDomain( { eruptionParameters->getVentPosition().x_ - int( p_.domainSizeX_ / 2.0 ) , eruptionParameters->getVentPosition().y_ - int( p_.domainSizeY_ / 2.0 ), 0.0 } , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , topology_, p_.simpleDom_ );
            dom_ = new MPISimpleDomain( p_.origin_ , { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ } , p_.dx_ , topology_.get(), p_.simpleDom_ );
            if( topology_->getRank() == 0 && p_.simpleDom_ == 0 ) this->logger->getWriter() << "topology used : simple with dimension splitted : x" << "\n";
            if( topology_->getRank() == 0 && p_.simpleDom_ == 1 ) this->logger->getWriter() << "topology used : simple with dimension splitted : y" << "\n";
            if( topology_->getRank() == 0 && p_.simpleDom_ == 2 ) this->logger->getWriter() << "topology used : simple with dimension splitted : z" << "\n";
        }
    }
    else{
        this->logger->getWriter() << "unknown topology type" << "\n";
        exit(0);
    }

    if( topology_->getRank() == 0 ){
        Int3 numBlockLoc = dom_->getNumBlockLocal();
        Int3 topo = ((MPI3DGridTopology*)topology_.get())->getTopology();
        this->logger->getWriter() << "block size = " << dom_->getBlockSize().x_ << "," << dom_->getBlockSize().y_ << "," << dom_->getBlockSize().z_ << "\n";
        this->logger->getWriter() << "domain size = " << dom_->getGlobalSize().x_ << "," << dom_->getGlobalSize().y_ << "," << dom_->getGlobalSize().z_ << "\n";
        this->logger->getWriter() << "actual domain size = " << numBlockLoc.x_ * topo.x_ * dom_->getBlockSize().x_ << "," << numBlockLoc.y_ * topo.y_ * dom_->getBlockSize().y_ << "," << numBlockLoc.z_ * topo.z_ * dom_->getBlockSize().z_ << "\n";
        double fact = double( numBlockLoc.x_ * topo.x_ * dom_->getBlockSize().x_ * numBlockLoc.y_ * topo.y_ * dom_->getBlockSize().y_ * numBlockLoc.z_ * topo.z_ * dom_->getBlockSize().z_ ) / double( dom_->getGlobalSize().x_ * dom_->getGlobalSize().y_ * dom_->getGlobalSize().z_ );
        this->logger->getWriter() << "actual domain is " << fact << " times larger than expected" << "\n";
    }

    for( auto fam: eruptionParameters->getParticleFamilies() ){
        dom_->addParticleFamily( fam );
    }

    repository_ = new MPIParticleRepository(dom_->getParticleFamilies(),(MPI3DGridTopology*)topology_.get());
    terrain_ = new MPIGridTerrain( { p_.origin_.x_, p_.origin_.y_ }, {p_.domainSizeX_, p_.domainSizeY_},p_.dxT_,dom_->getParticleFamilies(), (MPI3DGridTopology*)topology_.get() );

    terrain_->setActive(true);
    repository_->setActive(true);
    this->logger->getWriter().flush();
}

void Plume::InitSimulationTest(){

    InitCommon();

    std::shared_ptr< piaf::ConstantSpeed > constSpeed( new piaf::ConstantSpeed( { 1.0, 0.0, 0.0 } ) );
    constSpeed->setSpeedType(EULERIAN_STATIC);
    dom_->addSpeed( constSpeed );

}

void Plume::InitSimulationReal(){

    InitCommon();

    if( eruptionParameters->useWoodsModel() ){
        if( topology_->getRank() == 0 ) this->logger->getWriter() << "this simulation will use Woods plume model" << "\n";
    }
    else{
        if( topology_->getRank() == 0 ) this->logger->getWriter() << "this simulation will use Wim plume model" << "\n";
        //eruptionParameters->getPlume().setRadiusLimit(5000.0);
    }

    // create velocities and add them to the domain
    std::vector< std::shared_ptr< SpeedFunctor > > tmpSpeeds;
    std::vector< std::shared_ptr< SpeedFunctor > > speeds;
    if( eruptionParameters->useWoodsModel() ) {
        tmpSpeeds = createWoodsSpeeds( *eruptionParameters, p_.id_ );
    }
    else{
        tmpSpeeds = createWimSpeeds( *eruptionParameters );
    }

    // put all selected velocities in speeds vector
    for( auto name : eruptionParameters->getVelocitiesNames() ) {
        for( auto speed : tmpSpeeds ) {
            if( name.compare( speed->getName() ) == 0 ) speeds.push_back( speed );
        }
    }

    if( topology_->getRank() == 0 ) {
        this->logger->getWriter() << "Velocities used in this simulation : " << "\n";
        for( auto s : speeds  ) this->logger->getWriter() << "  - " << s->getName() << ", isDiffusion() = "<< s->isDiffusion() << "\n";
    }

    // Int3 size = dom_->getSize();
    // std::shared_ptr<GenericParticleFamily> particleFamily = (dom_->getParticleFamilies())[0];
    // VelocityField vf = VelocityField(speeds, dom_->getPosition(), {size.x_*dom_->getDx(), size.y_*dom_->getDx(), size.z_*dom_->getDx()}, dom_->getDx(), particleFamily);
    // fm->MPIFileManager::writeVelocityField(std::string("velocity_field.h5"), vf, 0.0, p_.dt_);
    

    for( std::shared_ptr< SpeedFunctor > s : speeds ) dom_->addSpeed( s );
    this->logger->getWriter().flush();

}

void Plume::F_init()
{
    // create the topology --------------
    this->createVolcanTopology();
    //----------------------------------

    this->logger->getWriter() << "tetras version " + p_.version_ << "\n";

#ifdef PROFILING
    this->profiler.setPrefixFileName("Plume");
    nbrItrToSaveStatsInFile=50;
    this->profiler.setRank(topology_->getRank());
#endif

    this->topology_->init();

    ///simulate
    //    if( topology_->getRank() == 0 ) {
    this->logger->getWriter() << "start simulation with " << topology_->getSize() << " process" << "\n";
    if( p_.ddt_ ) this->logger->getWriter() << "ddt activated, ddt = " << p_.ddtValue_ << "\n";
    else this->logger->getWriter() << "ddt not activated" << "\n";
    //  }

    this->fm.reset(new VolcanoMPIFileManager( topology_.get() ) );
    eruptionParameters.reset(new EruptionParameters (fm->loadEruptionParameters(p_.eruptionFile_, U0_, L0_ )));

    InitSimulation();

    if(this->topology_->getRank() == 0) {
        std::vector<double> diameters, grainSizes, densities, wtpcts;
        for( auto& pf : eruptionParameters->getParticleFamilies() ) {
            //this->logger->getWriter() << "diameter : " << pf->diameter_ << ", density : " << pf->density_ << "\n";
            diameters.push_back(pf->diameter_);
            densities.push_back(pf->density_);
            grainSizes.push_back(pf->grainSize_);
            wtpcts.push_back(pf->originalWeightPercent_);
        }

        this->logger->getWriter() << "phi : [ ";
        for(auto e : grainSizes) this->logger->getWriter() << e << ", ";
        this->logger->getWriter() << "]\n";

        this->logger->getWriter() << "diameter : [ ";
        for(auto e : diameters) this->logger->getWriter() << e << ", ";
        this->logger->getWriter() << "]\n";

        this->logger->getWriter() << "density : [ ";
        for(auto e : densities) this->logger->getWriter() << e << ", ";
        this->logger->getWriter() << "]\n";

        this->logger->getWriter() << "wtpcts : [ ";
        for(auto e : wtpcts) this->logger->getWriter() << e << ", ";
        this->logger->getWriter() << "]\n";

    }

    sim = new MPIExactSimulator(dom_, 0.0 ,p_.dt_,terrain_,repository_, false, topology_->getMpiManager());
    cpt = 0;
    timeStart = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());

    sim->setDdtValue( p_.ddtValue_ );

    if( p_.track_ ) {
        tracker_ = NULL;
        auto pftemp = eruptionParameters->getParticleFamilies();
        std::vector<std::shared_ptr<piaf::GenericParticleFamily> > pf;//( pftemp.begin(), pftemp.end() );
        for(auto partfam : pftemp){
            pf.push_back(std::shared_ptr<piaf::GenericParticleFamily>(partfam.get()));
        }
        if( p_.trackinterval_ > 0.0 ) tracker_ = new MPIParticleTracker( pf, topology_.get(), p_.trackinterval_ );
        else tracker_ = new MPIParticleTracker( pf, topology_.get() );
        sim->setParticleTracker( tracker_ );
    }

    if( p_.seed_ != -1 ) {
        sim->initRand( p_.seed_ );
    }
    if( topology_->getRank() == 0 ){
        this->logger->getWriter() << "running with seed (for process 0, others have seed = seed + rank) = " << sim->getSeed() << "\n";
        this->logger->getWriter() << totalParticlePerFamily*eruptionParameters->getParticleFamilies().size() << " particles will be injected" << "\n";
    }

    totalParticles.resize(eruptionParameters->getParticleFamilies().size(), 0);
    totalParticle = 0;

    //    injectedParticles = InjectParticles();
    //    sum_injectedParticles += injectedParticles;

    if( p_.ddt_ ) sim->enableDdt();
    //if( topology_->getRank() == 0 ) this->logger->getWriter() << "performing simulation without communication overlapping" << "\n";
#ifdef PROFILING
    this->profiler.start();
#endif
    this->logger->getWriter().flush();
}


bool Plume::isConverged(){
    return (!this->isProcess);
}


void Plume::verifyDomainContainsParticles(){
    if( cpt % frequenceCheckConvergence == 0){
        // int totalPart = dom_->countDomainParticles();
        // bool containsPart = dom_->containsParticles();
        // double currentTime = sim->getTime();
        // double endTime = eruptionParameters->getEndOfEruptionEvents();
        if( (! dom_->containsParticles()) && (sim->getTime() > eruptionParameters->getEndOfEruptionEvents()) ) {
            // if(topology_->getMpiManager().getRank() == 0){
            //     std::cout << "===========================================================" << std::endl;
            //     std::cout << "Plume::verifyDomainContainsParticles ... isProcess = FALSE" << std::endl;
            //     std::cout << "dom_->containsParticles() = " << containsPart << std::endl;
            //     std::cout << "sim->getTime() = " << currentTime << std::endl;
            //     std::cout << "eruptionParameters->getEndOfEruptionEvents() = " << endTime << std::endl;
            //     std::cout << "dom_->countDomainParticles() = " << totalPart << std::endl;
            //     std::cout << "cpt = " << cpt << std::endl;
            //     std::cout << "frequenceCheckConvergence = " << frequenceCheckConvergence << std::endl;
            //     std::cout << "===========================================================" << std::endl;
            // }
            this->isProcess=false;
        }  else{
            // if(topology_->getMpiManager().getRank() == 0){
            //     std::cout << "===========================================================" << std::endl;
            //     std::cout << "Plume::verifyDomainContainsParticles ... isProcess = TRUE" << std::endl;
            //     std::cout << "dom_->containsParticles() = " << containsPart << std::endl;
            //     std::cout << "sim->getTime() = " << currentTime << std::endl;
            //     std::cout << "eruptionParameters->getEndOfEruptionEvents() = " << endTime << std::endl;
            //     std::cout << "dom_->countDomainParticles() = " << totalPart << std::endl;
            //     std::cout << "cpt = " << cpt << std::endl;
            //     std::cout << "frequenceCheckConvergence = " << frequenceCheckConvergence << std::endl;
            //     std::cout << "===========================================================" << std::endl;
            // }
            this->isProcess=true;
        }
        // check timeout
        if( MPI_Wtime()-initialTime > p_.timeOut_*60.0 ){
            this->isProcess=false;
            if(topology_->getMpiManager().getRank() == 0){
                std::cout << "End of simulation due to time out ("<< p_.timeOut_ <<" minutes)" << std::endl;
            }
            throw std::string("End of simulation due to time out");
        }
    }
}

int Plume::computeParticleIncrement(int i){
    double mass = 0.0;
    double scaledMass = 0.0;
    int incr = 0;
    if( eruptionParameters->getMer( sim->getTime() ) == -1.0 ) incr = totalParticlePerFamily-totalParticles[i];
    else{
        mass = ( eruptionParameters->getMer( sim->getTime() ) * sim->getDt() ) * (eruptionParameters->getParticleFamilies()[i]->originalWeightPercent_ / 100.0);
        scaledMass = mass / scalingFactors[i];
        incr = scaledMass / eruptionParameters->getParticleFamilies()[i]->mass_;
        if( totalParticles[i]+incr > totalParticlePerFamily ) incr = totalParticlePerFamily - totalParticles[i];
    }
    return incr;
}

int Plume::InjectParticlesTest(){
    int totIncr;
    totIncr = 0;

    // insert particles
    // ...

    return totIncr;
}

int Plume::InjectParticlesAlongPlume(){

    int totIncr = 0;
    int incr = 0;

    //inject randomly linearly particles from vent to ht

    Double3 ventPosition = eruptionParameters->getVentPosition();
    double ventHeight = ventPosition.z_;
    double topHeightAboveVent = eruptionParameters->getPlume().getZ(0.0)->back();

    //  this->logger->getWriter() << "incr = " << incr << ", slices = " << slices << "\n";
    if( totalParticle < totalParticlePerFamily*(int)totalParticles.size() ){
        //std::uniform_real_distribution<double> dist(ventHeight, ventHeight+topHeightAboveVent);
        std::lognormal_distribution<double> dist(3.0,2.0);
        for (Uint i=0; i< totalParticles.size(); i++) {
            incr = computeParticleIncrement(i);
            injectedParticles += incr;
            totalParticles[i] += incr;
            totalParticle += incr;
            totIncr += incr;
            for(int p=0; p<incr; p++){
                // dom_->insertParticles( {ventPosition.x_, ventPosition.y_, dist(mt_)}, 1, i );
                double roll = 100.0 - dist(mt_);
                while(roll>100.0 || roll<0.0) roll = 100.0 - dist(mt_);
                dom_->insertParticles( {ventPosition.x_, ventPosition.y_, ventHeight + (roll/100.0)*topHeightAboveVent}, 1, i, scalingFactors[i] );
            }
        }
    }

    return totIncr;

}


int Plume::InjectParticlesAtCrater(){
    int incr = 0;
    int totIncr = 0;
    //double mass = 0.0;
    //double scaledMass = 0.0;
    injectedParticles = 0;

    // insert particles
    // if( totalParticles[0] < totalParticlePerFamily ) {
    if( totalParticle < totalParticlePerFamily*(int)totalParticles.size() ){
        for (Uint i=0; i< totalParticles.size(); i++) {
            incr = computeParticleIncrement(i);
            injectedParticles += incr;
            totalParticles[i] += incr;
            totalParticle += incr;
            totIncr += incr;
            dom_->insertParticles( eruptionParameters->getVentPosition(), incr,i, scalingFactors[i]);
        }
    }

    return totIncr;
}

void Plume::S(){
#ifdef PROFILING
    this->profiler.createStampRecording();
    this->profiler.registerBeginProcessLevelStamp();
#endif

    injectedParticles = InjectParticles();

    sum_injectedParticles +=injectedParticles;
    this->verifyDomainContainsParticles();

    sim->step();

#ifdef PROFILING
    this->profiler.registerEndProcessLevelStamp();
#endif
}


void Plume::Oi(){

    currentTime_ += sim->getDt();
    cpt++;

    static double lastWrite = 0.0;

    if( (lastWrite == 0 || currentTime_-lastWrite>=p_.trackwriteinterval_) && p_.track_ ) {
        auto removedParticles = tracker_->removeParticles();
        if(p_.trackFileType_ == "hdf5")
            ((MPIFileManager*)(&(*fm)))->writeHdf5(p_.trackFileBaseName_, removedParticles);
        if(p_.trackFileType_ == "vtk")
            ((MPIFileManager*)(&(*fm)))->writeVtk(p_.trackFileBaseName_, removedParticles);
        lastWrite = currentTime_;
    }

    // std::cout << "Plume::Oi() : p_.particlescsv_ = " << p_.particlescsv_ << std::endl; 

     if(p_.particlescsv_ == true){
        // std::cout << "getParticlesinBox" << std::endl; 
        auto particles = dom_->getParticlesInBox(p_.origin_, { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ });
        // std::cout << "writeParticlesCsv" << std::endl;
        ((MPIFileManager*)(&(*fm)))->writeParticlesCsv("track", particles);
        // std::cout << "done" << std::endl; 
    }

}


void Plume::Of(){
#ifdef PROFILING
    this->profiler.end();
    if( topology_->getRank() == 0 ) this->logger->getWriter() << this->profiler.getAveragedStats()<<"\n";
#endif

    std::vector<int> particleCountRepo = repository_->getStoredParticlesCounts();
    std::vector<int> particleCountTera = terrain_->getStoredParticlesCounts();
    std::vector<int> particleCount;
    std::transform(particleCountRepo.begin(), particleCountRepo.end(), particleCountTera.begin(), std::back_inserter(particleCount), std::plus<int>());

    // do the scaling
    //std::vector<double> scalingFactors;
    int totalPart = 0;
    for(Uint i=0; i<totalParticles.size(); i++) {
        //scalingFactors.push_back( eruptionParameters.getFamilyTotalMass(i) / ( totalParticle[i] * eruptionParameters.getParticleFamilies()[i]->mass_ ) );
        totalPart += totalParticles[i];
        //if( topology_->getRank() == 0 )log<LOG_DEBUG>("Required total mas for family %1% : %2% - total mass injected : %3% - factor : %4%") % i % eruptionParameters.getFamilyTotalMass(i) % ( totalParticle[i] * eruptionParameters.getParticleFamilies()[i]->mass_ ) % factors[i];
    }

    if( topology_->getRank() == 0 ) this->logger->getWriter() << "total inserted particles : " << totalPart << "\n" ;

    if( p_.outputFile_.compare("") != 0 ){
        fm->write(p_.outputFile_,terrain_, scalingFactors,*eruptionParameters.get(), totalParticles,  p_.dx_, sim->getDt(), p_.simulatorType_, currentTime_, dom_->getSpeedsName() );
    }

    if( p_.track_ ) {
        auto removedParticles = tracker_->removeParticles();
        //        ((MPIFileManager*)(&(*fm)))->write(p_.trackFileBaseName_+"-"+std::to_string(sim->getTime())+".h5", removedParticles);
        if(p_.trackFileType_ == "hdf5")
            ((MPIFileManager*)(&(*fm)))->writeHdf5(p_.trackFileBaseName_, removedParticles);
        if(p_.trackFileType_ == "vtk")
            ((MPIFileManager*)(&(*fm)))->writeVtk(p_.trackFileBaseName_, removedParticles);
    }

    timeEnd = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());
    timeDif = boost::posix_time::time_duration(timeEnd - timeStart);

    double totalMass = 0.0;
    for(uint family=0; family<particleCount.size(); family++){
        totalMass += particleCount[family] * scalingFactors[family] * eruptionParameters->getParticleFamilies()[family]->mass_;
    }

    if( topology_->getRank() == 0 ){
        this->logger->getWriter() << "Total mass gathered at the end of the simulation : " << totalMass << "\n";
        this->logger->getWriter() << "Number of particles at the end of simulation : ";
        for(auto n : particleCount) this->logger->getWriter() << n << ", ";
        this->logger->getWriter() << "\n";
        this->logger->getWriter() << "Computation time : " << timeDif << " ( " << timeDif.total_seconds() << " [s] )" << " - Number of particle families : " << eruptionParameters->getParticleFamilies().size() << " - Particles per family : " << totalParticlePerFamily << " - Number of process : " << topology_->getSize() << " - dx : " << dom_->getDx() << " - physical elapsed time : " << currentTime_ << "\n";
    }

    delete sim;
    this->logger->getWriter().flush();
    this->logger->getWriter().close();

}



void Plume::serialize(std::vector< BoundaryParticle >  & boundaryParticles, vector<char> & toSend){

    bool isconverged= (!this->isProcess); //this->isConverged();
    MapperUtil::serialize(boundaryParticles, isconverged ,toSend);
}

void Plume::deserialize(std::vector< BoundaryParticle > & newBoundaryParticles, vector<char> & vect){

    bool isRemoteConverged;
    MapperUtil::deserialize(newBoundaryParticles, isRemoteConverged, vect);

}

/**************************************** Continent ************************************************/

Continent::Continent(params p, shared_ptr<MPITopology> topology, double U0, double L0 )
    : Plume( p, topology, U0, L0 )
{
#ifdef PROFILING
    this->profiler.setPrefixFileName("Continent");
#endif
}

Continent::Continent(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)
    : Plume( argc, argv, globalCommunicator, U0, L0 )
{
#ifdef PROFILING
    this->profiler.setPrefixFileName("Continent");
#endif
}


Continent::~Continent() {

}

void Continent::verifyDomainContainsParticles(){
    if( cpt % frequenceCheckConvergence == 0 && cpt>0){
        // int totalPart = dom_->countDomainParticles();
        // bool containsPart = dom_->containsParticles();
        // double currentTime = sim->getTime();
        // double endTime = eruptionParameters->getEndOfEruptionEvents();
        if( (! dom_->containsParticles()) && (sim->getTime() > eruptionParameters->getEndOfEruptionEvents()) ) {
            // if(topology_->getMpiManager().getRank() == 0){
            //     std::cout << "===========================================================" << std::endl;
            //     std::cout << "Continent::verifyDomainContainsParticles ... isProcess = FALSE" << std::endl;
            //     std::cout << "dom_->containsParticles() = " << containsPart << std::endl;
            //     std::cout << "sim->getTime() = " << currentTime << std::endl;
            //     std::cout << "eruptionParameters->getEndOfEruptionEvents() = " << endTime << std::endl;
            //     std::cout << "dom_->countDomainParticles() = " << totalPart << std::endl;
            //     std::cout << "cpt = " << cpt << std::endl;
            //     std::cout << "frequenceCheckConvergence = " << frequenceCheckConvergence << std::endl;
            //     std::cout << "===========================================================" << std::endl;
            // }
            this->isProcess=false;
        }  else{
            // if(topology_->getMpiManager().getRank() == 0){
            //     std::cout << "===========================================================" << std::endl;
            //     std::cout << "Continent::verifyDomainContainsParticles ... isProcess = TRUE" << std::endl;
            //     std::cout << "dom_->containsParticles() = " << containsPart << std::endl;
            //     std::cout << "sim->getTime() = " << currentTime << std::endl;
            //     std::cout << "eruptionParameters->getEndOfEruptionEvents() = " << endTime << std::endl;
            //     std::cout << "dom_->countDomainParticles() = " << totalPart << std::endl;
            //     std::cout << "cpt = " << cpt << std::endl;
            //     std::cout << "frequenceCheckConvergence = " << frequenceCheckConvergence << std::endl;
            //     std::cout << "===========================================================" << std::endl;
            // }
            this->isProcess=true;
        }
    }
}

// To complete by Pierre
void Continent::F_init()
{

    // create the topology --------------
    this->createVolcanTopology();
    //----------------------------------

    this->topology_->init();
#ifdef PROFILING
    this->profiler.start();
#endif
    ///simulate
    if( topology_->getRank() == 0 ) {
        this->logger->getWriter() << "start continent scale simulation with " << topology_->getSize() << " process" << "\n";
        if( p_.ddt_ ) this->logger->getWriter() << "ddt activated, ddt = " << p_.ddtValue_ << "\n";
        else this->logger->getWriter() << "ddt not activated" << "\n";
    }

    //VolcanoMPIFileManager fm;
    this->fm.reset(new VolcanoMPIFileManager( topology_.get() ) );
    eruptionParameters.reset(new EruptionParameters (fm->loadEruptionParameters(p_.eruptionFile_, -1.0, -1.0 )));

    InitSimulation();

    timeStart = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());

    sim = new MPIExactSimulator(dom_,0.0,p_.dt_,terrain_,repository_, false, topology_->getMpiManager());
    cpt = 0;

    sim->setDdtValue( p_.ddtValue_ );

    if( p_.track_ ) {
        tracker_ = NULL;
        auto pftemp = eruptionParameters->getParticleFamilies();
        std::vector<std::shared_ptr<piaf::GenericParticleFamily> > pf;//( pftemp.begin(), pftemp.end() );
        for( auto partfam : pftemp){
            pf.push_back(std::shared_ptr<piaf::GenericParticleFamily>(partfam.get()));
        }
        if( p_.trackinterval_ > 0.0 ) tracker_ = new MPIParticleTracker( pf, topology_.get(), p_.trackinterval_ );
        else tracker_ = new MPIParticleTracker( pf, topology_.get() );
        sim->setParticleTracker( tracker_ );
    }

    if( p_.seed_ != -1 ) {
        sim->initRand( p_.seed_ );
    }
    if( topology_->getRank() == 0 ){
        this->logger->getWriter() << "running with seed (for process 0, others have seed = seed + rank) = " << sim->getSeed() << "\n";
    }

    totalParticles.resize(eruptionParameters->getParticleFamilies().size());

    if( p_.ddt_ ) sim->enableDdt();
    this->logger->getWriter().flush();
}

void Continent::S(){

#ifdef PROFILING
    this->profiler.createStampRecording();
    this->profiler.registerBeginProcessLevelStamp();
#endif

    //std::cout << "Continent::S() : before verifyDomainContainsParticles" << std::endl;
    
    this->verifyDomainContainsParticles();

    //std::cout << "Continent::S() : after verifyDomainContainsParticles" << std::endl;
    
    sim->step();
#ifdef PROFILING
    this->profiler.registerEndProcessLevelStamp();
#endif
}

void Continent::Oi(){

    currentTime_ += sim->getDt();
    cpt++;

    static double lastWrite = 0.0;

    if( (lastWrite == 0 || currentTime_-lastWrite>=p_.trackwriteinterval_) && p_.track_ ) {
        auto removedParticles = tracker_->removeParticles();
        //        ((MPIFileManager*)(&(*fm)))->write(p_.trackFileBaseName_+"-"+std::to_string(currentTime_)+".h5", removedParticles);
        if(p_.trackFileType_ == "hdf5")
            ((MPIFileManager*)(&(*fm)))->writeHdf5(p_.trackFileBaseName_, removedParticles);
        if(p_.trackFileType_ == "vtk")
            ((MPIFileManager*)(&(*fm)))->writeVtk(p_.trackFileBaseName_, removedParticles);
        lastWrite = currentTime_;
    }

    if(p_.particlescsv_ == true){
        // std::cout << "getParticlesinBox" << std::endl; 
        auto particles = dom_->getParticlesInBox(p_.origin_, { p_.domainSizeX_ , p_.domainSizeY_ , p_.domainSizeZ_ });
        // std::cout << "writeParticlesCsv" << std::endl;
        ((MPIFileManager*)(&(*fm)))->writeParticlesCsv("track", particles);
        // std::cout << "done" << std::endl; 
    }

}

void Continent::Of(){

#ifdef PROFILING
    this->profiler.end();
    if( topology_->getRank() == 0 ) this->logger->getWriter() << this->profiler.getAveragedStats()<<"\n";
#endif

    if( topology_->getRank() == 0 ) this->logger->getWriter() << cpt << " iterations completed" << "\n";

    // do the scaling
    int totalPart = 0;
    for(Uint i=0; i<totalParticles.size(); i++) {
        totalPart += totalParticles[i];
    }

    if( topology_->getRank() == 0 ) this->logger->getWriter() << "total inserted particles : " << totalPart << "\n" ;

    if( p_.outputFile_.compare("") != 0 ){
        fm->write(p_.outputFile_,terrain_, scalingFactors,*eruptionParameters.get(), totalParticles,  p_.dx_, sim->getDt(), p_.simulatorType_, currentTime_, dom_->getSpeedsName() );
    }

    if( p_.track_ ) {
        auto removedParticles = tracker_->removeParticles();
        //        ((MPIFileManager*)(&(*fm)))->write(p_.trackFileBaseName_+"-"+std::to_string(sim->getTime())+".h5", removedParticles);
        if(p_.trackFileType_ == "hdf5")
            ((MPIFileManager*)(&(*fm)))->writeHdf5(p_.trackFileBaseName_, removedParticles);
        if(p_.trackFileType_ == "vtk")
            ((MPIFileManager*)(&(*fm)))->writeVtk(p_.trackFileBaseName_, removedParticles);
    }

    timeEnd = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());
    timeDif = boost::posix_time::time_duration(timeEnd - timeStart);

    if( topology_->getRank() == 0 ){
        this->logger->getWriter() << "Computation time : " << timeDif << " ( " << timeDif.total_seconds() << " [s] )" << " - Number of particle families : " << eruptionParameters->getParticleFamilies().size() << " - Particles per family : " << totalParticlePerFamily << " - Number of process : " << topology_->getSize() << " - dx : " << dom_->getDx() << " - physical elapsed time : " << currentTime_ << "\n";
    }

    delete sim;
    this->logger->getWriter().flush();
    this->logger->getWriter().close();

}

/**************************************** LocalPlume ************************************************/

LocalPlume::LocalPlume( params p, shared_ptr<MPITopology> topology, double U0, double L0 )
    : Plume(p,topology, U0, L0)
{
    connector = &LocalVolcanicConnector::getInstance();
    connector->init();
}

LocalPlume::LocalPlume(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)
    : Plume(argc, argv, globalCommunicator, U0, L0)
{
    connector = &LocalVolcanicConnector::getInstance();
    connector->init();
}


LocalPlume::~LocalPlume(){

}

void LocalPlume::U(){

#ifdef PROFILING
    this->profiler.registerBeginGatherLevelStamp();
#endif
    //== remove the particles from the domain ==
    std::vector< BoundaryParticle > newBoundaryParticles = repository_->removeBoundaryParticles();
#ifdef PROFILING
    this->profiler.registerEndGatherLevelStamp();
#endif

    //if(newBoundaryParticles.empty()) return;

    //== display them ==
    //   if( newBoundaryParticles.size() > 0 ){
    //   if( topology_->getRank() == 0 )  this->logger->getWriter() << "Got " << newBoundaryParticles.size() << " boundary particles at iteration " << cpt << ". Simulation time = " << currentTime_ << "\n";
    //   int fid = newBoundaryParticles[0].getFamilyId();
    //   Double3 disp = newBoundaryParticles[0].getDisplacement();
    //   Double3 ls = newBoundaryParticles[0].getLagrangianSpeed();
    //   double ts = newBoundaryParticles[0].getTimeStamp();
    //   Double3 ds = newBoundaryParticles[0].getDiffusionWoods();
    //   if( topology_->getRank() == 0 ) this->logger->getWriter() << "First particle of vector : familyId = " << fid << " - displacement = ("<< disp.x_<<","<<disp.y_<<","<<disp.z_<<")  - lagrangian speed = ("<<ls.x_<<","<<ls.y_<<","<<ls.z_<<")  - diffusion speed = ("<<ds.x_<<","<<ds.y_<<","<<ds.z_<<")  - resultant speed = ("<<ls.x_-ds.x_<<","<<ls.y_-ds.y_<<","<<ls.z_-ds.z_<<") - time stamp = " << ts << "\n";
    // }

#ifdef PROFILING
    double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
#endif
    // exchange through Pipe/File or whatever (NO YET IMPLEMENTED)

    // vector<char> toSend;
    // this->serialize( newBoundaryParticles, toSend);
    // this->connector->send(toSend);
#ifdef PROFILING
    double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
    double WTimeReceiveBegin=  WTimeSendEnd;
#endif

    // vector<char> rec;
    // this->connector->receive(rec);

#ifdef PROFILING
    double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
    this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
    this->profiler.registerBeginUpdateLevelStamp();
    this->profiler.registerEndUpdateLevelStamp();
#endif
    //int totalParticlesInDomain= this->dom_->countDomainParticles();// per iteration
    int totalParticlesInTerrain=this->terrain_->countTerrainParticles();// total since begin of simultion

    int particules_sorties=newBoundaryParticles.size();// per iteration
    this->sum_outsideParticles +=particules_sorties;
    //int particules_deposees= totalParticlesInTerrain - previous_particules_deposees;// // per iteration
    previous_particules_deposees =totalParticlesInTerrain;
    // to verify
    //totalParticlesInDomain= Sum(injectedParticles) - (totalParticlesInTerrain+sum(particules_sorties) )
    //if(topology_->getRank() == 0 ) cout <<"["<< cpt<<"]: nbrTotalOfDomainParticles("<<totalParticlesInDomain<<") = sum_injectedParticles("<<sum_injectedParticles<<") - [ nbrParticlesInTerrain("<<totalParticlesInTerrain<<") + sum_outsideParticles("<<sum_outsideParticles<<") ]"<<"\n";

#ifdef PROFILING
    //registerParticles (int injected, int deposed, int outside){
    this->profiler.registerParticles (injectedParticles, particules_deposees, particules_sorties,
                                      totalParticlesInDomain, this->sum_injectedParticles, totalParticlesInTerrain,  this->sum_outsideParticles);

    this->profiler.commitStampRecording();
    if(this->cpt%nbrItrToSaveStatsInFile==0 && this->cpt>0){
        this->profiler.saveToFile();
    }
#endif

}

// inline bool fileExists(const std::string& name) {
//     ifstream f(name.c_str());
//     return f.good();
// }

void LocalPlume::mainLoop(){

    this->createVolcanTopology();

    // init the logger with name tetras-date-suffixe
    // check if a file with the same name exists until
    // finding one that doesnt exist (by testing different suffix)
    // there is a risk of race condition with the file system
    auto now = std::chrono::system_clock::now();
    time_t in_time_t = std::chrono::system_clock::to_time_t(now);
    std::string prefix("tetras");
    char buff[21];
    strftime(buff, 21, "%d-%m-%Y_%Hh%Mm%Ss", localtime(&in_time_t));
    std::string date(buff);
    int s = 0;
    std::string suffix(to_string(s));
    auto filename = prefix+"_"+date+"_"+suffix;
    while(fileExists(filename)){
        s++;
        std::string suffix(to_string(s));
        filename = prefix+"_"+date+"_"+suffix;
    }
    if(topology_->getRank() == 0) std::cout << "logfile " << filename << ".log" << std::endl;

    this->initLogger(filename);

    this->simulate();

}

/**************************************** LocalContinent ************************************************/

LocalContinent::LocalContinent(params p, shared_ptr<MPITopology> topology, double U0, double L0 )
    : Continent(p,topology, U0, L0)
{
    connector = &LocalVolcanicConnector::getInstance();
    connector->init();
}

LocalContinent::LocalContinent(int argc, char ** argv, MPI_Comm globalCommunicator, double U0, double L0)
    : Continent(argc, argv, globalCommunicator, U0, L0)
{
    connector = &LocalVolcanicConnector::getInstance();
    connector->init();
}


LocalContinent::~LocalContinent(){}

void LocalContinent::mainLoop(){
    this->createVolcanTopology();
    this->initLogger("Continent");
    this->simulate();
}

void LocalContinent::U(){

    // TO DO BY Pierre, for example:
    //      1. receive/send from/to somewhere
    //      2. update Domain
    //      4. send/receive to/from/ somewhere
#ifdef PROFILING
    this->profiler.registerBeginGatherLevelStamp();
#endif
    std::vector< BoundaryParticle > newBoundaryParticles = repository_->removeBoundaryParticles();
#ifdef PROFILING
    this->profiler.registerEndGatherLevelStamp();
#endif
    
#ifdef PROFILING
    double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
    double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
    double WTimeReceiveBegin=  WTimeSendEnd;
    double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
    this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
    this->profiler.registerBeginUpdateLevelStamp();
    this->profiler.registerEndUpdateLevelStamp();
#endif


    //int totalParticlesInDomain= this->dom_->countDomainParticles();// per iteration
    int totalParticlesInTerrain=this->terrain_->countTerrainParticles();// total since begin of simultion

    int particules_sorties=newBoundaryParticles.size();// per iteration
    this->sum_outsideParticles +=particules_sorties;
    //int particules_deposees= totalParticlesInTerrain - previous_particules_deposees;// // per iteration
    previous_particules_deposees =totalParticlesInTerrain;
    // to verify
    //totalParticlesInDomain= Sum(injectedParticles) - (totalParticlesInTerrain+sum(particules_sorties) )
    // if(topology_->getRank() == 0 ) cout <<"["<< cpt<<"]: nbrTotalOfDomainParticles("<<totalParticlesInDomain<<") = sum_injectedParticles("<<sum_injectedParticles<<") - [ nbrParticlesInTerrain("<<totalParticlesInTerrain<<") + sum_outsideParticles("<<sum_outsideParticles<<") ]"<<"\n";

#ifdef PROFILING
    //registerParticles (int injected, int deposed, int outside){
    this->profiler.registerParticles (injectedParticles, particules_deposees, particules_sorties,
                                      totalParticlesInDomain, this->sum_injectedParticles, totalParticlesInTerrain,  this->sum_outsideParticles);
    this->profiler.commitStampRecording();
    if(this->cpt%nbrItrToSaveStatsInFile==0 && this->cpt>0){
        cout<<"cpt="<<cpt<<"---";
        this->profiler.saveToFile();
    }
#endif

}

/*********************************** LocalConnector *****************************************************/

LocalVolcanicConnector::LocalVolcanicConnector(): isInitialized (false){}

LocalVolcanicConnector::~LocalVolcanicConnector(){}

void LocalVolcanicConnector::init (){
}

void LocalVolcanicConnector::send(vector<char> const & buffer){

    cout<<" sending buffer: "<<buffer.size()<<"\n";
}

void  LocalVolcanicConnector::receive(vector<char> & buffer){
    cout<<" receiving buffer: "<< buffer.size()<<"\n";

}

int LocalVolcanicConnector::getSubmodel_ID(){return 0;}
void LocalVolcanicConnector::end(){}

/**************************************** DistributedSimulation ************************************************/
DistributedPlume::DistributedPlume(params p, shared_ptr<MPITopology> topology, double U0, double L0)
    : Plume (p, topology, U0, L0){
    //if( topology->getRank() == 0 ) cout <<" This is the distributed Plume " <<"\n";

    this->isRemoteConverged=false;

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.type=MapperType::SUBMODEL;
    this->info.dt=p.dt_;
    this->info.dx=p.dx_;
    //this->info.sectionId= this->getSubmodel_ID();
#ifdef PROFILING
    this->profiler.setPrefixFileName("distributedPlume");
#endif
}


DistributedPlume::DistributedPlume(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : Plume (argc, argv, globalCommunicator, U0, L0)
{
    //if( topology->getRank() == 0 ) cout <<" This is the distributed Plume " <<"\n";

    this->isRemoteConverged=false;

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.type=MapperType::SUBMODEL;

    //this->info.sectionId= this->getSubmodel_ID();
#ifdef PROFILING
    this->profiler.setPrefixFileName("distributedPlume");
#endif
}




DistributedPlume::~DistributedPlume(){}

void DistributedPlume::synchronizeCommunication(){

    this->info.dt=p_.dt_;
    this->info.dx=p_.dx_;

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);

    if( topology_->getRank() == 0 ) { //requires only on

        vector<char> toSend;
        this->serializeHeader(toSend, this->info);
        this->send(toSend);//HERE: send toSend vector through MUSCLE

        vector<char> toReceive;
        this->receive(toReceive);// <------ Attentin: data only exists on rank 0 -> you need to broadcast
        INFO_SMART_MAPPER info_tmp;
        this->unserializeHeader(toReceive, info_tmp);
        //if( topology_->getRank() == 0 ) this->logger->getWriter() <<"[sync] -> synchronized with the  remote SubmodelID  "<< info_tmp.sectionId <<"\n";
    }
}

void DistributedPlume::F_init(){
    // call parent F_init
    Plume::F_init();
    this->info.sectionId= this->getSubmodel_ID();
    this->synchronizeCommunication();
}

bool DistributedPlume::isConverged(){
    return (Plume::isConverged() && this->isRemoteConverged);
}
void DistributedPlume::U(){

    try{
        //bool debug=true;
        bool debug=false;
        //== remove the particles from the domain ==

#ifdef PROFILING
        this->profiler.registerBeginGatherLevelStamp();
#endif
        std::vector< BoundaryParticle > boundaryParticles = repository_->removeBoundaryParticles();
#ifdef PROFILING
        this->profiler.registerEndGatherLevelStamp();
#endif
        //if(boundaryParticles.empty()) return;
        if( debug && topology_->getRank() == 0 ) this->logger->getWriter() <<  "== Update Boundary -> cpt:"<< this->cpt<< " ==" << "\n";
#ifdef PROFILING
        double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
#endif
        this->logger->getWriter().flush();
        std::vector< BoundaryParticle > newBoundaryParticles;
        if( topology_->getRank() == 0 ) { // serialize/deserialize and send/receive only on rank 0
            //== serialize particles ==
            vector<char> toSend;
            this->serialize( boundaryParticles, toSend);
            //====== send ===========
            this->sendConvergence();

            this->send(toSend);//HERE: send toSend vector through MUSCLE
            if( debug && topology_->getRank() == 0 ) {
                this->logger->getWriter() <<  "local convergence of DistributedPlume : " << !this->isProcess <<"\n";
                this->logger->getWriter() <<  "sending vector<boundaryParticles> of size : " <<boundaryParticles.size() <<"\n";
            }

        }
        this->logger->getWriter().flush();
#ifdef PROFILING
        double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
        double WTimeReceiveBegin=  WTimeSendEnd;
#endif
        if( topology_->getRank() == 0 ) {

            //== deserialize particles ==
            vector<char> toReceive;
            // HERE: receive from MUSCLE
            this->receive(toReceive); // <------ Attention: data only exists on rank 0 -> you need to broadcast
            //toReceive.swap(toSend);
            this->deserialize(newBoundaryParticles,toReceive);
            if( debug && topology_->getRank() == 0 )this->logger->getWriter() <<  "receiving vector<BoundaryParticle> of size: " <<newBoundaryParticles.size() <<"\n";
            vector<char>().swap(toReceive);
        }
        this->logger->getWriter().flush();
#ifdef PROFILING
        double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
        this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
        this->profiler.registerBeginUpdateLevelStamp();
#endif
        MpiManager & mpiManager= this->topology_->getMpiManager();

        // broadcast newBoundaryParticles <- MPI operations
        int localSize = newBoundaryParticles.size()*sizeof(BoundaryParticle);
        mpiManager.bCast<int>(&localSize, 1,   mpiManager.bossId());
        if( topology_->getRank() > 0 ) newBoundaryParticles.resize(localSize/sizeof(BoundaryParticle));
        mpiManager.bCast<char>((char *) &newBoundaryParticles[0], localSize, mpiManager.bossId());
        // braodcast the isRemoteConverged boolean received from mapper: only rank zero got it
        mpiManager.bCast<char>((char*)&isRemoteConverged, sizeof(isRemoteConverged), mpiManager.bossId());

        if( newBoundaryParticles.size() > 0 ){
            if( debug && topology_->getRank() == 0 )  this->logger->getWriter() << "Got " << newBoundaryParticles.size() << " boundary particles at iteration " << cpt << ". Simulation time = " << currentTime_ << "\n";
            int fid = newBoundaryParticles[0].getFamilyId();
            Double3 disp = newBoundaryParticles[0].getDisplacement();
            Double3 ls = newBoundaryParticles[0].getLagrangianSpeed();
            double ts = newBoundaryParticles[0].getTimeStamp();
            Double3 ds = newBoundaryParticles[0].getDiffusionSpeed();
            if(debug &&  topology_->getRank() == 0) this->logger->getWriter() << "First particle of vector : familyId = " << fid << " - displacement = ("<< disp.x_<<","<<disp.y_<<","<<disp.z_<<")  - lagrangian speed = ("<<ls.x_<<","<<ls.y_<<","<<ls.z_<<")  - diffusion speed = ("<<ds.x_<<","<<ds.y_<<","<<ds.z_<<")  - resultant speed = ("<<ls.x_-ds.x_<<","<<ls.y_-ds.y_<<","<<ls.z_-ds.z_<<") - time stamp = " << ts << "\n";
        }

#ifdef PROFILING
        this->profiler.registerEndUpdateLevelStamp();
#endif

        //int totalParticlesInDomain= this->dom_->countDomainParticles();// per iteration
        int totalParticlesInTerrain=this->terrain_->countTerrainParticles();// total since begin of simultion

        int particules_sorties=boundaryParticles.size();// per iteration
        this->sum_outsideParticles +=particules_sorties;
        //int particules_deposees= totalParticlesInTerrain - previous_particules_deposees;// // per iteration
        previous_particules_deposees =totalParticlesInTerrain;
        // to verify
        //totalParticlesInDomain= Sum(injectedParticles) - (totalParticlesInTerrain+sum(particules_sorties) )
        //if(topology_->getRank() == 0 ) cout <<"["<< cpt<<"]: nbrTotalOfDomainParticles("<<totalParticlesInDomain<<") = sum_injectedParticles("<<sum_injectedParticles<<") - [ nbrParticlesInTerrain("<<totalParticlesInTerrain<<") + sum_outsideParticles("<<sum_outsideParticles<<") ]"<<"\n";

#ifdef PROFILING
        //registerParticles (int injected, int deposed, int outside){
        this->profiler.registerParticles (injectedParticles, particules_deposees, particules_sorties,
                                          totalParticlesInDomain, this->sum_injectedParticles, totalParticlesInTerrain,  this->sum_outsideParticles);
        //void registerParticles (int injected, int deposed, int outside,
        //                      int totalParticlesInDomain, int sum_injectedParticles, int totalParticlesInTerrain, int sum_outsideParticles){

        this->profiler.commitStampRecording();
        if(this->cpt% nbrItrToSaveStatsInFile==0 && this->cpt>0){
            this->profiler.saveToFile();
        }
#endif
        if( debug && topology_->getRank() == 0 ) this->logger->getWriter() <<  "====================" << "\n";
        this->logger->getWriter().flush();
    }catch (const std::exception &exc)
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr <<"#######> plume sendung: "<< exc.what()<<endl;
    }
}

void DistributedPlume::Of(){
    // call parent
    Plume::Of();

    this->end();
}

void DistributedPlume::serialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & toSend){
    bool isconverged= (!this->isProcess) ;//Plume::isConverged();// pay attention here we use the isConverged of Plume
    // stringstream ss;
    // ss<<"["<<this->cpt<<"]: Plume isconverged="<< isconverged<< " isRemoteConverged="<< this->isRemoteConverged<<"\n";
    // cout<<ss.str();
    // this->logger->getWriter() << ss.str();
    //  this->logger->getWriter().flush();
    MapperUtil::serialize(boundaryParticles, isconverged ,toSend, this->cpt);
}


void DistributedPlume::deserialize(std::vector< BoundaryParticle > & newBoundaryParticles, vector<char> & vect){
    MapperUtil::deserialize(newBoundaryParticles, this->isRemoteConverged, vect);
}


//*************************************************************************************************************/
DistributedAtmospheredPlume::DistributedAtmospheredPlume(params p, shared_ptr<MPITopology> topology, double U0, double L0)
    : DistributedPlume (p, topology, U0, L0){
    lastTime = -1.0;
}

DistributedAtmospheredPlume::DistributedAtmospheredPlume(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedPlume (argc, argv, globalCommunicator, U0, L0){
    lastTime = -1.0;
}

DistributedAtmospheredPlume::~DistributedAtmospheredPlume(){}


void DistributedAtmospheredPlume::synchronizeCommunication(){

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.dt=p_.dt_;
    this->info.dx=p_.dx_;

    if( topology_->getRank() == 0 ) { //requires only on

        vector<char> toSend;
        this->serializeHeader(toSend, this->info);
        this->send(toSend);//HERE: send toSend vector through MUSCLE

        vector<char> toReceive;
        this->receive(toReceive);// <------ Attentin: data only exists on rank 0 -> you need to broadcast
        INFO_SMART_MAPPER info_tmp, info_atm;
        this->unserializeHeader(toReceive, info_tmp);
        //if( topology_->getRank() == 0 ) this->logger->getWriter() <<"[sync] -> synchronized with the  remote SubmodelID  "<< info_tmp.sectionId <<"\n";
        // sync with Atmosphere
        this->sendToAtmosphere(toSend);//HERE: send toSend vector through MUSCLE
        vector<char> toReceiveAtm;
        this->receiveFromAtmosphere(toReceiveAtm);// <------ Attentin: data only exists on rank 0 -> you need to broadcast
        this->unserializeHeader(toReceiveAtm, info_atm);
    }
}

AtmosphereDataRequest DistributedAtmospheredPlume::createAtmosphereDataRequest(){
    Coord3D craterPosition(0.0,0.0,0.0);
    Coord3D domainSize (200000.0, 200000.0, 200000.0);
    Coord3D domainPosition(-100000.0, -100000.0, -100000.0);
    CoordGeo craterGeoPosition(45.989286,11.531441);
    string date="20000101";
    double time= 0.0;
    AtmosphereDataRequest request(craterPosition, domainSize, domainPosition, craterGeoPosition, date,  time);
    return request;
}

void DistributedAtmospheredPlume::sendAtmosphereDataRequest(){

    if(this->topology_->getRank() == 0){

        stringstream ss;

        AtmosphereDataRequest request=createAtmosphereDataRequest();
        vector<char> buffer;

        vector<char> tmp;
        bool isconv= isConverged();
        buffer.resize(sizeof(isconv));
        memcpy(&buffer[0], &isconv, sizeof(isconv) );
        request.serialize(tmp);
        buffer.insert(buffer.end(), tmp.begin(), tmp.end());

        ss<<"---------------"<<"\n";
        ss<<"request["<<this->cpt<<"]: "<<"\n";
        ss<<"->\n"<<request.toString()<<"\n";
        cout<<ss.str();
        //========== send request ==============
        this->sendToAtmosphere(buffer);
    }

}

void DistributedAtmospheredPlume::receiveAtmosphereDataResponse(AtmosphereDataResponse &atmResponse){

    if(this->topology_->getRank() == 0){
        stringstream ss;
        vector<char> data;
        this->receiveFromAtmosphere(data);
        atmResponse.unSerialize(data);

        ss<<"<-\n"<<atmResponse.toString()<<"\n";
        ss<<"---------------"<<"\n";
        cout<<ss.str();
    }

}
void DistributedAtmospheredPlume::U(){
    // do treatement
    DistributedPlume::U();
    // get atmosphere data every hour of simulation
    if( lastTime == -1.0 || (int(sim->getTime()) % 3600) < int(lastTime)  ){
        sendAtmosphereDataRequest();
        AtmosphereDataResponse atmResponse;
        receiveAtmosphereDataResponse(atmResponse);
    }
    lastTime = sim->getTime();
}


/**************************************** DistributedContinent ************************************************/


DistributedContinent::DistributedContinent(params p, shared_ptr<MPITopology> topology, double U0, double L0)
    : Continent (p, topology, U0, L0)
{
    this->isRemoteConverged=false;
    // if( topology->getRank() == 0 ) cout <<" This is the distributed continent " <<"\n";


    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.type=MapperType::SUBMODEL;
    this->info.dt=p.dt_;
    this->info.dx=p.dx_;
    //this->info.sectionId=getSubmodel_ID();

#ifdef PROFILING
    this->profiler.setPrefixFileName("distributedContinent");
#endif
}
DistributedContinent::DistributedContinent(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : Continent (argc, argv, globalCommunicator, U0, L0)
{
    this->isRemoteConverged=false;
    // if( topology->getRank() == 0 ) cout <<" This is the distributed continent " <<"\n";


    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.type=MapperType::SUBMODEL;

    //this->info.sectionId=getSubmodel_ID();

#ifdef PROFILING
    this->profiler.setPrefixFileName("distributedContinent");
#endif
}


DistributedContinent::~DistributedContinent(){}


void DistributedContinent::synchronizeCommunication(){

    this->info.dt=p_.dt_;
    this->info.dx=p_.dx_;
    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);

    if( topology_->getRank() == 0 ) { //requires only on

        vector<char> toSend;
        this->serializeHeader(toSend, this->info);
        this->send(toSend);//HERE: send toSend vector through MUSCLE
        vector<char> toReceive;
        this->receive(toReceive);// <------ Attention: data only exists on rank 0 -> you need to broadcast
        INFO_SMART_MAPPER info_tmp;
        this->unserializeHeader(toReceive, info_tmp);
        //if( topology_->getRank() == 0 ) this->logger->getWriter() <<"[sync] -> synchronized with the  remote SubmodelID  "<< info_tmp.sectionId <<"\n";
    }
}

void DistributedContinent::F_init(){
    // call parent F_init
    Continent::F_init();
    this->info.sectionId= this->getSubmodel_ID();
    this->synchronizeCommunication();
}

void DistributedContinent::Of(){
    // call parent F_init
    Continent::Of();
    this->end();
}


void DistributedContinent::U(){
    //bool debug=true;
    bool debug=false;
    try{
#ifdef PROFILING
        this->profiler.registerBeginGatherLevelStamp();
#endif
        //== remove the particles from the domain ==
        std::vector< BoundaryParticle > boundaryParticles = repository_->removeBoundaryParticles();
#ifdef PROFILING
        this->profiler.registerEndGatherLevelStamp();
#endif

        //if(boundaryParticles.empty()) return;
        if( debug && topology_->getRank() == 0 ) this->logger->getWriter() <<  "== Update Boundary-> cpt:"<< this->cpt <<" ==" << "\n";

#ifdef PROFILING
        double WTimeSendBegin=unige_pasc::getTimeStamp<double>();
#endif
        std::vector< BoundaryParticle > newBoundaryParticles;
        if( topology_->getRank() == 0 ) { // serialize/deserialize and send/receive only on rank 0

            //== serialize particles ==
            vector<char> toSend;
            //if (debug) if( topology_->getRank() == 0 ) this->logger->getWriter() <<  "start serialisation";
            this->serialize( boundaryParticles, toSend);
            // if (debug) if( topology_->getRank() == 0 ) this->logger->getWriter() <<  " -> [OK]" << "\n";
            if (debug) if( topology_->getRank() == 0 ){
                this->logger->getWriter() <<  "local convergence : " << !this->isProcess <<"\n";
                this->logger->getWriter() <<  "sending vector<boundaryParticles> of size : " <<boundaryParticles.size() <<"\n";
            }
            //== send =0
            this->sendConvergence();
            this->send(toSend);//HERE: send toSend vector through MUSCLE

        }

#ifdef PROFILING
        double WTimeSendEnd=unige_pasc::getTimeStamp<double>();
        double WTimeReceiveBegin=  WTimeSendEnd;
#endif
        //if (debug) if( topology_->getRank() == 0 ) this->logger->getWriter() <<  "sending of char size " <<toSend.size() <<"\n";
        if( topology_->getRank() == 0 ) {

            //== deserialize particles ==
            vector<char> toReceive;
            // HERE: receive from MUSCLE
            this->receive(toReceive); // <------ Attentoin: data only exists on rank 0 -> you need to broadcast
            //toReceive.swap(toSend);
            // if (debug) if( topology_->getRank() == 0 ) this->logger->getWriter() <<  "receiving size char: "<< toReceive.size()<<"\n";
            this->deserialize(newBoundaryParticles,toReceive);
            if( debug && topology_->getRank() == 0 ) this->logger->getWriter() <<  "receiving vector<BoundaryParticle> of size : "<< newBoundaryParticles.size()<<"\n";

        }

#ifdef PROFILING
        double WTimeReceiveEnd=unige_pasc::getTimeStamp<double>();
        this->profiler.registerSendReceiveLevelStamp( WTimeSendBegin, WTimeSendEnd, WTimeReceiveBegin, WTimeReceiveEnd);
        this->profiler.registerBeginUpdateLevelStamp();
#endif
        MpiManager & mpiManager= this->topology_->getMpiManager();

        // broadcast newBoundaryParticles <- MPI operations
        int localSize = newBoundaryParticles.size()*sizeof(BoundaryParticle);
        mpiManager.bCast<int>(&localSize, 1, mpiManager.bossId());
        if( topology_->getRank() > 0 ) newBoundaryParticles.resize(localSize/sizeof(BoundaryParticle));
        mpiManager.bCast<char>((char*)&newBoundaryParticles[0], localSize, mpiManager.bossId());
        // braodcast the isRemoteConverged boolean received from mapper: only rank zero got it
        mpiManager.bCast<char>((char*)&isRemoteConverged, sizeof(isRemoteConverged), mpiManager.bossId());

        if( newBoundaryParticles.size() > 0 ){
            if( debug && topology_->getRank() == 0 )  this->logger->getWriter() << "Got " << newBoundaryParticles.size() << " boundary particles at iteration " << cpt << ". Simulation time = " << currentTime_ << "\n";
            int fid = newBoundaryParticles[0].getFamilyId();
            Double3 disp = newBoundaryParticles[0].getDisplacement();
            Double3 ls = newBoundaryParticles[0].getLagrangianSpeed();
            double ts = newBoundaryParticles[0].getTimeStamp();
            Double3 ds = newBoundaryParticles[0].getDiffusionSpeed();
            if( debug && topology_->getRank() == 0 ) this->logger->getWriter() << "First particle of vector : familyId = " << fid << " - displacement = ("<< disp.x_<<","<<disp.y_<<","<<disp.z_<<")  - lagrangian speed = ("<<ls.x_<<","<<ls.y_<<","<<ls.z_<<")  - diffusion speed = ("<<ds.x_<<","<<ds.y_<<","<<ds.z_<<")  - resultant speed = ("<<ls.x_-ds.x_<<","<<ls.y_-ds.y_<<","<<ls.z_-ds.z_<<") - time stamp = " << ts << "\n";
            for( auto part : newBoundaryParticles ) dom_->putParticle( part );
            this->sum_injectedParticles+=newBoundaryParticles.size();
        }
#ifdef PROFILING
        this->profiler.registerEndUpdateLevelStamp();
#endif

        //int totalParticlesInDomain= this->dom_->countDomainParticles();// per iteration
        int totalParticlesInTerrain=this->terrain_->countTerrainParticles();// total since begin of simultion

        int particules_sorties=boundaryParticles.size();// per iteration
        this->sum_outsideParticles +=particules_sorties;
        //int particules_deposees= totalParticlesInTerrain - previous_particules_deposees;// // per iteration
        previous_particules_deposees =totalParticlesInTerrain;
        // to verify
        //totalParticlesInDomain= Sum(injectedParticles) - (totalParticlesInTerrain+sum(particules_sorties) )
        // if(topology_->getRank() == 0 ) cout <<"["<< cpt<<"]: nbrTotalOfDomainParticles("<<totalParticlesInDomain<<") = sum_injectedParticles("<<sum_injectedParticles<<") - [ nbrParticlesInTerrain("<<totalParticlesInTerrain<<") + sum_outsideParticles("<<sum_outsideParticles<<") ]"<<"\n";


#ifdef PROFILING
        //registerParticles (int injected, int deposed, int outside){
        this->profiler.registerParticles (injectedParticles, particules_deposees, particules_sorties,
                                          totalParticlesInDomain, this->sum_injectedParticles, totalParticlesInTerrain,  this->sum_outsideParticles);
        this->profiler.commitStampRecording();

        if(this->cpt%nbrItrToSaveStatsInFile==0 && this->cpt>0){
            // cout<<"cpt="<<cpt<<"\n";
            this->profiler.saveToFile();
        }
#endif
        if( debug && topology_->getRank() == 0 ) this->logger->getWriter() <<  "====================" << "\n";
        this->logger->getWriter().flush();
    }catch (const std::exception &exc)
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr <<"#######> contient  "<< exc.what()<<endl;
    }
}

bool DistributedContinent::isConverged(){
    return (Continent::isConverged() && this->isRemoteConverged);
}


void DistributedContinent::serialize(std::vector< BoundaryParticle > & boundaryParticles, vector<char> & toSend){
    bool isconverged= (!this->isProcess) ;//Continent::isConverged();// pay attention here we use the isConverged of Continent
    MapperUtil::serialize(boundaryParticles, isconverged ,toSend);
}


void DistributedContinent::deserialize(std::vector< BoundaryParticle > & newBoundaryParticles, vector<char> & vect){
    MapperUtil::deserialize(newBoundaryParticles, this->isRemoteConverged, vect);
}



//********************************************************* DistributedAtmospheredContinent ****************************************************/

DistributedAtmospheredContinent::DistributedAtmospheredContinent(  params p, shared_ptr<MPITopology> topology, double U0, double L0):
    DistributedContinent (p, topology, U0, L0){
    lastTime = -1.0;
}
DistributedAtmospheredContinent::DistributedAtmospheredContinent(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedContinent (argc, argv, globalCommunicator, U0, L0){
    lastTime = -1.0;
}

DistributedAtmospheredContinent::~DistributedAtmospheredContinent(){}


void DistributedAtmospheredContinent::synchronizeCommunication(){

    boost::posix_time::ptime tm(boost::posix_time::microsec_clock::local_time());
    this->info.timeStamp=MapperUtil::to_time_t(tm);
    this->info.dt=p_.dt_;
    this->info.dx=p_.dx_;

    if( topology_->getRank() == 0 ) { //requires only on

        vector<char> toSend;
        this->serializeHeader(toSend, this->info);
        this->send(toSend);//HERE: send toSend vector through MUSCLE

        vector<char> toReceive;
        this->receive(toReceive);// <------ Attentin: data only exists on rank 0 -> you need to broadcast
        INFO_SMART_MAPPER info_tmp, info_atm;
        this->unserializeHeader(toReceive, info_tmp);
        //if( topology_->getRank() == 0 ) this->logger->getWriter() <<"[sync] -> synchronized with the  remote SubmodelID  "<< info_tmp.sectionId <<"\n";
        // sync with Atmosphere
        this->sendToAtmosphere(toSend);//HERE: send toSend vector through MUSCLE
        vector<char> toReceiveAtm;
        this->receiveFromAtmosphere(toReceiveAtm);// <------ Attentin: data only exists on rank 0 -> you need to broadcast
        this->unserializeHeader(toReceiveAtm, info_atm);
    }
}

AtmosphereDataRequest DistributedAtmospheredContinent::createAtmosphereDataRequest(){
    Coord3D craterPosition(0.0,0.0,0.0);
    Coord3D domainSize (200000.0, 200000.0, 200000.0);
    Coord3D domainPosition(-100000.0, -100000.0, -100000.0);
    CoordGeo craterGeoPosition(45.989286,11.531441);
    string date="20000101";
    double time= 0.0;
    AtmosphereDataRequest request(craterPosition, domainSize, domainPosition, craterGeoPosition, date,  time);
    return request;
}

void DistributedAtmospheredContinent::sendAtmosphereDataRequest(){

    if(this->topology_->getRank() == 0){

        stringstream ss;

        AtmosphereDataRequest request=createAtmosphereDataRequest();
        vector<char> buffer;

        vector<char> tmp;
        bool isconv= isConverged();
        buffer.resize(sizeof(isconv));
        memcpy(&buffer[0], &isconv, sizeof(isconv) );
        request.serialize(tmp);
        buffer.insert(buffer.end(), tmp.begin(), tmp.end());

        ss<<"---------------"<<"\n";
        ss<<"request["<<this->cpt<<"]: "<<"\n";
        ss<<"->\n"<<request.toString()<<"\n";
        cout<<ss.str();
        //========== send request ==============
        this->sendToAtmosphere(buffer);
    }

}

void DistributedAtmospheredContinent::receiveAtmosphereDataResponse(AtmosphereDataResponse &atmResponse){

    if(this->topology_->getRank() == 0){
        stringstream ss;
        vector<char> data;
        this->receiveFromAtmosphere(data);
        atmResponse.unSerialize(data);

        ss<<"<-\n"<<atmResponse.toString()<<"\n";
        ss<<"---------------"<<"\n";
        cout<<ss.str();
    }

}
void DistributedAtmospheredContinent::U(){
    // do treatement
    DistributedContinent::U();
    // get atmosphere data every hour of simulation
    if( lastTime == -1.0 || (int(sim->getTime()) % 3600) < int(lastTime)  ){
        sendAtmosphereDataRequest();
        AtmosphereDataResponse atmResponse;
        receiveAtmosphereDataResponse(atmResponse);
    }
    lastTime = sim->getTime();
}

/*****************************************************LightDistributedPlume ***********************************************************/

LightDistributedPlume::LightDistributedPlume(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedPlume(argc, argv, globalCommunicator, U0, L0), AsynchronousRemoteMpiKernel(){

}

LightDistributedPlume::~LightDistributedPlume(){}

int LightDistributedPlume::getSubmodel_ID(){
    return 0;
}

void LightDistributedPlume::send(vector<char> const& toSend){
    //f_out->send((char*)toSend.data(), toSend.size());
    this->iSend(f_out, (char*)toSend.data(), toSend.size());
}

void LightDistributedPlume::receive(vector<char> & data){
    int count;
    //char * buffer= this->f_in->receive(count);
    char * buffer= this->iRecvAndWaitAll(f_in, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}


void LightDistributedPlume::sendConvergence(){
    bool islocalConverged= !this->isProcess;
    f_out->send((char*) &islocalConverged, sizeof(islocalConverged));
}

bool LightDistributedPlume::receiveConvergence(){
    int count=0;
    char * buffer= this->f_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete[] buffer;
    MapperUtil::deserialize(this->isRemoteConverged, vect);
    vector<char>().swap(vect);
    return (count)? true:false;
}


void LightDistributedPlume::end(){}


void LightDistributedPlume::mainLoop(){
    this->createVolcanTopology();
    this->topology_->setMpiManager(this->mpiManager);
    this->initLogger(this->getName());
    f_out=this->getConduitEntrance<char>("f_out");
    f_in=this->getConduitExit<char>("f_in");

    std::cout << "Hey from LightDistributedPlume::mainLoop(), my rank is " << this->getMpiManager()->getRank() << std::endl;

    this->simulate();
}

void LightDistributedPlume::createVolcanTopology(){
    this->argc=this->getKernelArgc();
    this->argv=this->getKernelArgv();
    //cout<<" ++++++++> "<<this->getKernelArgv()[0] <<" "<<this->getKernelArgv()[1]<<"\n";
    Plume::createVolcanTopology();
}



/*****************************************************LightDistributedAtmospheredPlume ***********************************************************/

LightDistributedAtmospheredPlume::LightDistributedAtmospheredPlume(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedAtmospheredPlume(argc, argv, globalCommunicator, U0, L0), AsynchronousRemoteMpiKernel(){
}

LightDistributedAtmospheredPlume::~LightDistributedAtmospheredPlume(){}

int LightDistributedAtmospheredPlume::getSubmodel_ID(){
    return 0;
}

void LightDistributedAtmospheredPlume::send(vector<char> const& toSend){
    this->iSend(f_out, (char*)toSend.data(), toSend.size());
}

void LightDistributedAtmospheredPlume::receive(vector<char> & data){
    int count;
    char * buffer= this->iRecvAndWaitAll(f_in, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}


void LightDistributedAtmospheredPlume::sendConvergence(){
    bool islocalConverged= !this->isProcess;
    f_out->send((char*) &islocalConverged, sizeof(islocalConverged));
}

bool LightDistributedAtmospheredPlume::receiveConvergence(){
    int count=0;
    char * buffer= this->f_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete[] buffer;
    MapperUtil::deserialize(this->isRemoteConverged, vect);
    vector<char>().swap(vect);
    return (count)? true:false;
}


void LightDistributedAtmospheredPlume::sendToAtmosphere(vector<char> const& toSend){
    this->iSend(f_out_atm, (char*)toSend.data(), toSend.size());
}

void LightDistributedAtmospheredPlume::receiveFromAtmosphere(vector<char> & data){
    int count;
    char * buffer= this->iRecvAndWaitAll(f_in_atm, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

void LightDistributedAtmospheredPlume::end(){}


void LightDistributedAtmospheredPlume::mainLoop(){
    this->createVolcanTopology();
    this->topology_->setMpiManager(this->mpiManager);
    this->initLogger(this->getName());
    f_out = this->getConduitEntrance<char>("f_out");
    f_in = this->getConduitExit<char>("f_in");
    f_out_atm = this->getConduitEntrance<char>("f_out_atm");
    f_in_atm = this->getConduitExit<char>("f_in_atm");
    this->simulate();
}

void LightDistributedAtmospheredPlume::createVolcanTopology(){
    this->argc=this->getKernelArgc();
    this->argv=this->getKernelArgv();
    Plume::createVolcanTopology();
}


/*****************************************************LightDistributedAtmospheredContinent ***********************************************************/

LightDistributedAtmospheredContinent::LightDistributedAtmospheredContinent(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedAtmospheredContinent( argc, argv, globalCommunicator, U0, L0), AsynchronousRemoteMpiKernel(){
}

LightDistributedAtmospheredContinent::~LightDistributedAtmospheredContinent(){}

int LightDistributedAtmospheredContinent::getSubmodel_ID(){
    return 0;
}

void LightDistributedAtmospheredContinent::send(vector<char> const& toSend){
    this->iSend(f_out, (char*)toSend.data(), toSend.size());
}

void LightDistributedAtmospheredContinent::receive(vector<char> & data){
    int count;
    char * buffer= this->iRecvAndWaitAll(f_in, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

void LightDistributedAtmospheredContinent::sendConvergence(){
    bool islocalConverged= !this->isProcess;
    f_out->send((char*) &islocalConverged, sizeof(islocalConverged));
}

bool LightDistributedAtmospheredContinent::receiveConvergence(){
    int count=0;
    char * buffer= this->f_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete[] buffer;
    MapperUtil::deserialize(this->isRemoteConverged, vect);
    vector<char>().swap(vect);
    return (count)? true:false;
}

void LightDistributedAtmospheredContinent::sendToAtmosphere(vector<char> const& toSend){
    this->iSend(f_out_atm, (char*)toSend.data(), toSend.size());
}

void LightDistributedAtmospheredContinent::receiveFromAtmosphere(vector<char> & data){
    int count;
    char * buffer= this->iRecvAndWaitAll(f_in_atm, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

void LightDistributedAtmospheredContinent::end(){}


void LightDistributedAtmospheredContinent::mainLoop(){
    this->createVolcanTopology();
    this->topology_->setMpiManager(this->mpiManager);
    this->initLogger(this->getName());
    f_out = this->getConduitEntrance<char>("f_out");
    f_in = this->getConduitExit<char>("f_in");
    f_out_atm = this->getConduitEntrance<char>("f_out_atm");
    f_in_atm = this->getConduitExit<char>("f_in_atm");
    this->simulate();
}


//---------------------------------------------------------------------------------------------

LightDistributedContinent::LightDistributedContinent(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0)
    : DistributedContinent(argc, argv, globalCommunicator, U0, L0), AsynchronousRemoteMpiKernel() {
}


LightDistributedContinent::~LightDistributedContinent(){}
int LightDistributedContinent::getSubmodel_ID(){
    return 2;
}

void LightDistributedContinent::send(vector<char> const& toSend){
    //f_out->send((char*)toSend.data(), toSend.size());
    this->iSend(f_out, (char*)toSend.data(), toSend.size());
}

void LightDistributedContinent::receive(vector<char> & data){


    int count;
    //char * buffer= this->f_in->receive(count);
    char * buffer= this->iRecvAndWaitAll(f_in, count);
    std::vector<char> vect(buffer, buffer + count);
    vect.swap(data);
    delete[] buffer;
}

void LightDistributedContinent::sendConvergence(){
    bool islocalConverged= !this->isProcess;
    f_out->send((char*) &islocalConverged, sizeof(islocalConverged));
}

bool LightDistributedContinent::receiveConvergence(){
    int count=0;
    char * buffer= this->f_in->receive(count);
    std::vector<char> vect(buffer, buffer + count);
    delete[] buffer;
    MapperUtil::deserialize(this->isRemoteConverged, vect);
    vector<char>().swap(vect);
    return (count)? true:false;
}

void LightDistributedContinent::end(){

}

void LightDistributedContinent::mainLoop(){
    this->createVolcanTopology();
    this->topology_->setMpiManager(this->mpiManager);
    this->initLogger(this->getName());
    f_out=this->getConduitEntrance<char>("f_out");
    f_in=this->getConduitExit<char>("f_in");
    this->simulate();
}

void LightDistributedContinent::createVolcanTopology()
{
    this->argc=this->getKernelArgc();
    this->argv=this->getKernelArgv();
    Continent::createVolcanTopology();
}

//********************************************* TolopogyCreator ******************************************************

TolopogyCreator::TolopogyCreator(){}
shared_ptr<MPITopology> TolopogyCreator::createTopology(params &p, int argc, char **argv, MPI_Comm globalCommunicator){
    GLOBAL_LOG_LEVEL = LOG_DEBUG;

    p.version_ = VERSION;
    ArgParser & parser= getArgParser();

    try{
        parser.parseArgs(argc, argv, p);
    }catch(exception& e){
        cout << "Standard exception: " << e.what() << "\n";
    }

    shared_ptr<MPITopology> topology;
    shared_ptr<MpiManager> mpiManager(std::make_shared<MpiManager>()) ;//getManagerMPI();
    mpiManager->init(& argc, &argv, globalCommunicator);

    if( p.simpleDom_ != -1 ){
        topology = std::move(std::make_shared<MPITopology>( argc, argv, mpiManager ));

    }else{
        topology= std::move(std::make_shared<MPI3DGridTopology>( argc, argv, mpiManager ));
        //  topology = new MPI3DGridTopology( argc, argv, mpiManager );
    }

    if( topology->getRank() == 0 ) cout << "tetras version " + p.version_ << "\n";
    return topology;
}


/************************************** FactoryLocalPlume *********************************************/
FactoryLocalPlume::FactoryLocalPlume(){}
FactoryLocalPlume::~FactoryLocalPlume(){}
/*shared_ptr<MPITopology> FactoryLocalPlume::getTopology(params &p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/
LocalPlume * FactoryLocalPlume::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (LocalPlume *) new LocalPlume(argc, argv, globalCommunicator, U0, L0);
}

/************************************** FactoryLocalContinent *********************************************/
FactoryLocalContinent::FactoryLocalContinent(){}
FactoryLocalContinent::~FactoryLocalContinent(){}
/*shared_ptr<MPITopology> FactoryLocalContinent::getTopology(params &p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/
LocalContinent * FactoryLocalContinent::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (LocalContinent *) new LocalContinent( argc, argv, globalCommunicator, U0, L0);
}

//********************************************* FactoryDistributedPlume ******************************************************
FactoryDistributedPlume::FactoryDistributedPlume(){}
FactoryDistributedPlume::~FactoryDistributedPlume(){}

/*shared_ptr<MPITopology> FactoryDistributedPlume::getTopology(params & p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/
////
FactoryDistributedAtmospheredPlume::FactoryDistributedAtmospheredPlume(){}
FactoryDistributedAtmospheredPlume::~FactoryDistributedAtmospheredPlume(){}
/*shared_ptr<MPITopology> FactoryDistributedAtmospheredPlume::getTopology(params & p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/

//********************************************* FactoryLightDistributedPlume ******************************************************
FactoryLightDistributedPlume::FactoryLightDistributedPlume():FactoryDistributedPlume(){}
FactoryLightDistributedPlume::~FactoryLightDistributedPlume(){}

DistributedPlume * FactoryLightDistributedPlume::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (DistributedPlume *) new LightDistributedPlume( argc,  argv, globalCommunicator, U0, L0);
}

//********************************************* FactoryLightDistributedAmospheredPlume ******************************************************
FactoryLightDistributedAtmospheredPlume::FactoryLightDistributedAtmospheredPlume():FactoryDistributedAtmospheredPlume(){}
FactoryLightDistributedAtmospheredPlume::~FactoryLightDistributedAtmospheredPlume(){}

DistributedAtmospheredPlume *FactoryLightDistributedAtmospheredPlume::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (DistributedAtmospheredPlume *) new LightDistributedAtmospheredPlume( argc,  argv, globalCommunicator, U0, L0);
}
//********************************************* FactoryDistributedContinent ******************************************************
FactoryDistributedContinent::FactoryDistributedContinent(){}
FactoryDistributedContinent::~FactoryDistributedContinent(){}

/*shared_ptr<MPITopology> FactoryDistributedContinent::getTopology(params & p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/
////
FactoryDistributedAtmospheredContinent::FactoryDistributedAtmospheredContinent(){}
FactoryDistributedAtmospheredContinent::~FactoryDistributedAtmospheredContinent(){}
/*shared_ptr<MPITopology> FactoryDistributedAtmospheredContinent::getTopology(params & p, int argc, char **argv){
    return TolopogyCreator::getInstance().createTopology(p, argc, argv);
}*/

//********************************************* FactoryLightDistributedContinent ******************************************************
FactoryLightDistributedContinent::FactoryLightDistributedContinent():FactoryDistributedContinent(){}
FactoryLightDistributedContinent::~FactoryLightDistributedContinent(){}
DistributedContinent * FactoryLightDistributedContinent::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (DistributedContinent *) new LightDistributedContinent(argc,  argv, globalCommunicator, U0, L0);
}
//********************************************* FactoryLightDistributedAmospheredContinent ******************************************************
FactoryLightDistributedAtmospheredContinent::FactoryLightDistributedAtmospheredContinent():FactoryDistributedAtmospheredContinent(){}
FactoryLightDistributedAtmospheredContinent::~FactoryLightDistributedAtmospheredContinent(){}

DistributedAtmospheredContinent *FactoryLightDistributedAtmospheredContinent::newInstance(int argc, char **argv, MPI_Comm globalCommunicator, double U0, double L0){
    return  (DistributedAtmospheredContinent *) new LightDistributedAtmospheredContinent(argc,  argv, globalCommunicator, U0, L0);
}
