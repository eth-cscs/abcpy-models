/*
TEphra TRAsport Simulator (TETRAS)
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


#ifndef VOLCANICTEPHRASPEED_H_
#define VOLCANICTEPHRASPEED_H_

#include "VolcanicTephraParticleFamilies.hpp"
#include "Tools/Atmosphere.hpp"
#include "Tools/WindProfile.hpp"
#include "Tools/WimPlume.hpp"
#include "Tools/VolcanicColumnProfile.hpp"
#include <boost/math/constants/constants.hpp>

#include "piaf.hpp"

using namespace piaf;

static const double pivalue = boost::math::constants::pi<double>();

/* sedimentation from "volcanic ash" book eq. 25 for CD */

// for spherical particles
inline double FS(){
    return 1.0;
}

// for spherical particles
inline double FN(){
    return 1.0;
}

inline double kS(double FS){
    return 0.5*(pow(FS, 1.0/3.0) + pow(FS, -1.0/3.0));
}

inline double kN(double FN){
    return pow(10.0, 0.45*pow(-log10(FN), 0.99));
}

inline double CD(double re){
    return (24.0*kS(FS()) / re) * (1.0+0.125*pow( re*kN(FN()) / kS(FS()) , 2.0/3.0 )) + 0.46*kN(FN()) / (1.0 + 5330.0 / (re* kN(FN()) / kS(FS())));
}


inline double sedimentation(double altitude, double diameter, double density){
    double oldRe = 0.1;
    double re = 1.0;
    double rhoP = density;
    double d = diameter;
    double v = 0.0;

    /*!
    * \todo ca pourrait etre repris de l'atmosphere
    */
    // const double rhoA = 1.25*exp( - (altitude/1000.0)/8.2 );
    // viscosity of air
    // const double mu = 1.75*pow(10.0, -5.0);

    const double rhoA = 1.177;
    // viscosity of air
    //const double mu = 1.75*pow(10.0, -5.0);
    const double mu = 1.846e-5;

    while(relativeDif(oldRe,re)>0.01){
        v = sqrt( (4.0*9.81*d*(rhoP-rhoA)) / (3.0*rhoA*CD(re)) );
        oldRe = re;
        re = rhoA * abs(v) * d / mu;
    }

    return v;
}


/*=====================================================*/

// compute speed from equation A1 (Bonadonna 2003)
inline double VA1(double re, double rhoP, double rhoA, double d){
    double cd = 0.0;
    if(re<6.0) cd = 24.0/re;
    else if((re>=6.0) && (re<500.0)) cd = 10.0/(pow(re,0.5));
    else if(re>=500.0) cd = 0.43;

    double v = sqrt((4.0*9.81*d*(rhoP - rhoA)) / (3.0*rhoA*cd));

    return v;
}

// compute speed from equation A4 (Bonadonna 2003)
inline double VA4(double re, double rhoP, double rhoA, double d, double mu){
    if(re<6.0) return (9.81*d*d*(rhoP-rhoA))/(18.0*mu);
    else if((re>=6.0) && (re<500.0)) return d*pow( ( ( 4.0*9.81*9.81*pow( ( rhoP-rhoA ),2.0) ) / ( 225.0*rhoA*mu )), ( 1.0/3.0 ) );
    else return pow(((3.1*9.81*d*(rhoP-rhoA))/(rhoA)),1.0/2.0);
}

inline double RE(double v, double rhoA, double d, double mu){
    return (d*rhoA*v)/mu;
}

inline double sedimentationOld(double altitude, double diameter, double density){
    // from kae's code
    /*!
    * \todo ca pourrait etre repris de l'atmosphere
    */
    const double rhoA = 1.25*exp( - (altitude/1000.0)/8.2 );
    // viscosity of air
    const double mu = 1.75*pow(10.0, -5.0);

    double v = 0.0;
    double d;
    double rhoP;
    double re;
    double oldRe;

    d = diameter;
    rhoP = density;

    re = 1.0;
    oldRe = 0.1;

    /* \todo instabilite numerique possible, voir quelle est la facon la plus efficace de traiter ca (relative diff ou memoriser avant derniere valeur) */
    //v = VA1(re,rhoP,rhoA,d);
    while(relativeDif(re,oldRe)>0.01){
        v = VA1(re,rhoP,rhoA,d);
        oldRe = re;
        re = RE(v, rhoA, d, mu);
        //v = VA1(re,rhoP,rhoA,d);
    }

    return v;
}

class Backflow: public SpeedFunctor{
public:

    Backflow( WimPlume plume, Double3 ventPosition, double hb): SpeedFunctor("tetras_backflow"),
        plume_( plume ), ventPosition_( ventPosition ), hb_( hb ){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:

    WimPlume plume_;
    Double3 ventPosition_;
    double hb_;

};
inline Double3 Backflow::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double threshold = 500.0;
    //double v = 50.0;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( (! std::get<1>(indexClosestPoint)) && pos.z_ < hb_ ) {
        Double3 closestPoint = (*plume_.getCenterLinePoints( t ))[std::get<0>(indexClosestPoint)];
        // we are outside the plume, compute distance to the border
        double distBorder = dist3D( closestPoint , pos ) - (*plume_.getR( t ))[std::get<0>(indexClosestPoint)];
        if( distBorder <= threshold ) {
            Double3 backflowV = (closestPoint - pos) / dist3D( closestPoint, pos );
            //backflowV = backflowV * v;
            backflowV = backflowV * distBorder;
            return backflowV;
        }
    }
    return {0.0, 0.0, 0.0};
}

// For Wim's plume model, wind and plume are packed in the same structure, and sedimentation for this version.
// Vertical plume velocity is no more applied after plume corner in this version
class WindWimPlumeSpeedGaussianSedimentationThreshold: public SpeedFunctor{
public:

    WindWimPlumeSpeedGaussianSedimentationThreshold( WimPlume plume, WindProfile wind, Double3 ventPosition, double dPlume, double dAtm): SpeedFunctor("tetras_wind_and_wim_plume_gaussian_sedimentation_triggered_threshold"),
        plume_( plume ), wind_( wind ),
        ventPosition_( ventPosition ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    void setHb(double hb){ hb_ = hb; }

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;
    double hb_;


};
inline Double3 WindWimPlumeSpeedGaussianSedimentationThreshold::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    std::tuple<int, bool> indexClosestPointToCorner = plume_.indexClosestPoint( {ventPosition_.x_, ventPosition_.y_, hb_}, t, false );
    double thresholdRadius = plume_.getR(t)->at(get<0>(indexClosestPointToCorner));
    // double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // compute sedimentation speed
    std::shared_ptr<TephraParticleFamily> *fam;
    fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
    double sed = sedimentation(pos.z_, (*fam)->diameter_, (*fam)->density_);
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    // to be more precise we should test on the norm of the plume speed and we should be able differentiate wind speed and plume speed inside the plume
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
        double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
        double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        if(r > thresholdRadius) plumeSpeed = {0.0, 0.0, 0.0};
        else plumeSpeed = 2.0 * plumeSpeed * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
    }
    //    else{
    std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
    auto windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    //    }
    if(sed > plumeSpeed.z_){
        windSpeed.z_ -= sed;
        return windSpeed;
    }
    return plumeSpeed;
}

// For Wim's plume model, wind and plume are packed in the same structure, and sedimentation for this version
// either we apply vertical plume speed, either sedimentation speed
class WindWimPlumeSpeedGaussianSedimentation: public SpeedFunctor{
public:

    WindWimPlumeSpeedGaussianSedimentation( WimPlume plume, WindProfile wind, Double3 ventPosition, double dPlume, double dAtm): SpeedFunctor("tetras_wind_and_wim_plume_gaussian_sedimentation_triggered"),
        plume_( plume ), wind_( wind ),
        ventPosition_( ventPosition ){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;


};
inline Double3 WindWimPlumeSpeedGaussianSedimentation::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    // double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // compute sedimentation speed
    std::shared_ptr<TephraParticleFamily> *fam;
    fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
    double sed = sedimentation(pos.z_, (*fam)->diameter_, (*fam)->density_);
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    // to be more precise we should test on the norm of the plume speed and we should be able differentiate wind speed and plume speed inside the plume
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
        double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
        double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        plumeSpeed = 2.0 * plumeSpeed * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
    }
    //    else{
    std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
    auto windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    //    }
    if(sed > plumeSpeed.z_){
        windSpeed.z_ -= sed;
        return windSpeed;
    }
    return plumeSpeed;
}

// For Wim's plume model, wind and plume are packed in the same structure
class WindWimPlumeSpeedGaussianAboveHb: public SpeedFunctor{
public:

    WindWimPlumeSpeedGaussianAboveHb( WimPlume plume, WindProfile wind, Double3 ventPosition, double dPlume, double dAtm, double hb): SpeedFunctor("tetras_wind_and_wim_plume_gaussian_above_hb"),
        plume_( plume ), wind_( wind ), hb_( hb ),
        ventPosition_( ventPosition ){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:

    WimPlume plume_;
    WindProfile wind_;
    double hb_;
    Double3 ventPosition_;


};
inline Double3 WindWimPlumeSpeedGaussianAboveHb::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    Double3 windSpeed = {0.0, 0.0, 0.0};
    //	double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        if(pos.z_ > hb_){
            plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
            double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
            double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
            plumeSpeed = 2.0 * plumeSpeed * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
        } else {
            plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
        }
    }
    else{
        std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
        windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    }
    return plumeSpeed + windSpeed;
}

// For Wim's plume model, wind and plume are packed in the same structure
class WindWimPlumeSpeedGaussian: public SpeedFunctor{
public:

    WindWimPlumeSpeedGaussian( WimPlume plume, WindProfile wind, Double3 ventPosition, double dPlume, double dAtm): SpeedFunctor("tetras_wind_and_wim_plume_gaussian"),
        plume_( plume ), wind_( wind ),
        ventPosition_( ventPosition ){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;


};
inline Double3 WindWimPlumeSpeedGaussian::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    Double3 windSpeed = {0.0, 0.0, 0.0};
    //	double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
        double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
        double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        plumeSpeed = 2.0 * plumeSpeed * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
    }
    else{
        std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
        windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    }
    return plumeSpeed + windSpeed;
}

// For Wim's plume model, wind and plume are packed in the same structure
class WindWimPlumeSpeed: public SpeedFunctor{
public:

    WindWimPlumeSpeed( WimPlume plume, WindProfile wind, Double3 ventPosition, double dPlume, double dAtm): SpeedFunctor("tetras_wind_and_wim_plume"),
        plume_( plume ), wind_( wind ),
        ventPosition_( ventPosition ){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;


};
inline Double3 WindWimPlumeSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    Double3 windSpeed = {0.0, 0.0, 0.0};
    //	double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
    }
    else{
        std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
        windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    }
    return plumeSpeed + windSpeed;
}




class WindWimPlumeSpeedHblimited: public SpeedFunctor{
public:

    WindWimPlumeSpeedHblimited( WimPlume plume, WindProfile wind, Double3 ventPosition, double hb): SpeedFunctor("tetras_wind_and_wim_plume_hblimited"),
        plume_( plume ), wind_( wind ), ventPosition_( ventPosition ), hb_( hb ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
    //void setHb(double hb){ hb_ = hb; }

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;
    double hb_;

};
inline Double3 WindWimPlumeSpeedHblimited::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    Double3 windSpeed = {0.0, 0.0, 0.0};
    //	double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( get<1>(indexClosestPoint) ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
    }
    else{
        std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
        windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    }
    if(pos.z_ > hb_) plumeSpeed = {0.0, 0.0, 0.0};
    return plumeSpeed + windSpeed;
}

class WindWimPlumeSpeedHbMaxWidth: public SpeedFunctor{
public:

    WindWimPlumeSpeedHbMaxWidth( WimPlume plume, WindProfile wind, Double3 ventPosition, double hb): SpeedFunctor("tetras_wind_and_wim_plume_hbmaxwidth"),
        plume_( plume ), wind_( wind ), ventPosition_( ventPosition ), hb_( hb ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
    //void setHb(double hb){ hb_ = hb; }

private:

    WimPlume plume_;
    WindProfile wind_;
    Double3 ventPosition_;
    double hb_;

};
inline Double3 WindWimPlumeSpeedHbMaxWidth::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    Double3 plumeSpeed = {0.0, 0.0, 0.0};
    Double3 windSpeed = {0.0, 0.0, 0.0};
    //	double radius = (plume_.getX( t ))->back();
    std::tuple<int, bool> indexClosestPoint;
    auto posTest = pos;
    if( pos.z_ > hb_ ) posTest.z_ = hb_;
    bool inPlume = get<1>(plume_.indexClosestPoint( posTest, t, false ));
    // precise plume test
    indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( inPlume ){
        // set plume speed
        plumeSpeed  = plume_.getUVect(t)->at(get<0>(indexClosestPoint));
    }
    else{
        std::vector< piaf::Double3 >* windLinearSampling = wind_.getLinearSampling( t, 1000.0 );
        windSpeed = windLinearSampling->at( pos.z_ / 1000.0 );
    }
    return plumeSpeed + windSpeed;
}





inline Double3 anisotropicDiffusion( double DH, double DV, double dt ){
    double urV = sqrt( ( 2.0 * DV ) / dt );
    double urH = sqrt( ( 4.0 * DH ) / dt );
    double alpha = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
    urV = urV * sin(alpha) + urV * cos(alpha);
    alpha = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
    return{ urH * cos( alpha ), urH * sin( alpha ), urV };
}

inline Double3 isotropicDiffusion( double D, double dt ){
    double ur = sqrt( ( 6.0 * D ) / dt );
    double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * boost::math::constants::pi<double>();
    double phi = acos( ( double( rand() ) / double( RAND_MAX ) ) * 2.0 - 1.0 );
    return { ur * cos( theta ) * sin( phi ), ur * sin( theta ) * sin( phi ), ur * cos( phi ) };
}

class DiffusionWimUmbrellaCloud: public SpeedFunctor{
public:

    DiffusionWimUmbrellaCloud( WimPlume plume, Double3 ventPosition, double dPlume, double dAtm, double hb, double ht):SpeedFunctor("tetras_diffusion_wim_umbrella_cloud"),
        dPlume_( dPlume ), dAtm_( dAtm ), hb_( hb ), ht_( ht ), plume_( plume ), ventPosition_( ventPosition ) {}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    double dPlume_;
    double dAtm_;
    double hb_;
    double ht_;

private:

    WimPlume plume_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWimUmbrellaCloud::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double D = dAtm_;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );;
    if( std::get<1>(indexClosestPoint) || (pos.z_ > hb_ && pos.z_ < ht_) ) D = dPlume_;
    return isotropicDiffusion( D, dt );
}

class DiffusionWim: public SpeedFunctor{
public:

    DiffusionWim( WimPlume plume, Double3 ventPosition, double dPlume, double dAtm):SpeedFunctor("tetras_diffusion_wim"),
        dPlume_( dPlume ), dAtm_( dAtm ), plume_( plume ), ventPosition_( ventPosition ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    double dPlume_;
    double dAtm_;

private:

    WimPlume plume_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWim::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double D = dAtm_;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );
    if( std::get<1>(indexClosestPoint) ) D = dPlume_;
    return isotropicDiffusion( D, dt );
}

class DiffusionWimGaussian: public SpeedFunctor{
public:

    DiffusionWimGaussian( WimPlume plume, Double3 ventPosition, double dPlume, double dAtm):SpeedFunctor("tetras_diffusion_wim_gaussian"),
        dPlume_( dPlume ), dAtm_( dAtm ), plume_( plume ), ventPosition_( ventPosition ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    double dPlume_;
    double dAtm_;

private:

    WimPlume plume_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWimGaussian::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double D = dAtm_;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );;
    if( std::get<1>(indexClosestPoint) ){
        //double b = plume_.getR(t)->at(get<0>(indexClosestPoint)) / 2.0;
        double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
        double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        //D = dPlume_ * pow(boost::math::constants::e<double>(),-((r*r)/(b*b)));
        // diffusion is the gaussian diffusion of the plume + the constant diffusion of atmosphere
        D = dAtm_ + 2.0 * dPlume_ * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
    }
    return isotropicDiffusion( D, dt );
}

class DiffusionWimGaussianTriggered: public SpeedFunctor{
public:

    DiffusionWimGaussianTriggered( WimPlume plume, Double3 ventPosition, double dPlume, double dAtm):SpeedFunctor("tetras_diffusion_wim_gaussian_triggered"),
        dPlume_( dPlume ), dAtm_( dAtm ), plume_( plume ), ventPosition_( ventPosition ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    double dPlume_;
    double dAtm_;

private:

    WimPlume plume_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWimGaussianTriggered::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double D = dAtm_;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );;
    if( std::get<1>(indexClosestPoint) ){
        std::shared_ptr<TephraParticleFamily> *fam;
        fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
        // double b = plume_.getR(t)->at(get<0>(indexClosestPoint)) / 2.0;
        double b2 = pow(plume_.getR(t)->at(get<0>(indexClosestPoint)), 2.0) / 2.0;
        double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        // D = dPlume_ * pow(boost::math::constants::e<double>(),-((r*r)/(b*b)));
        D = dPlume_ * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
        if( sedimentation(pos.z_, (*fam)->diameter_, (*fam)->density_) >  plume_.getU(t)->at(get<0>(indexClosestPoint)) ) D = 0.0;

    }
    return isotropicDiffusion( D, dt );
}

class DiffusionWimTriggered: public SpeedFunctor{
public:

    DiffusionWimTriggered( WimPlume plume, Double3 ventPosition, double dPlume, double dAtm):SpeedFunctor("tetras_diffusion_wim_triggered"),
        dPlume_( dPlume ), dAtm_( dAtm ), plume_( plume ), ventPosition_( ventPosition ){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    double dPlume_;
    double dAtm_;

private:

    WimPlume plume_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWimTriggered::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double D = dAtm_;
    std::tuple<int, bool> indexClosestPoint = plume_.indexClosestPoint( pos, t, false );;
    if( std::get<1>(indexClosestPoint) ){
        std::shared_ptr<TephraParticleFamily> *fam;
        fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
        //double b = plume_.getR(t)->at(get<0>(indexClosestPoint)) / 2.0;
        //double r = dist3D(pos,plume_.getCenterLinePoints(t)->at(get<0>(indexClosestPoint)));
        D = dPlume_;
        if( sedimentation(pos.z_, (*fam)->diameter_, (*fam)->density_) >  plume_.getU(t)->at(get<0>(indexClosestPoint)) ) D = 0.0;

    }
    return isotropicDiffusion( D, dt );
}

/*==========================================================================================================*/
/*!
* \brief speed used to define diffusion
*/
class DiffusionWoods: public SpeedFunctor{
public:

    DiffusionWoods(std::shared_ptr<VolcanicColumnProfile> columnProfile, Double3 ventPosition): SpeedFunctor("tetras_diffusion_woods"),
        columnProfile_(columnProfile), ventPosition_(ventPosition){}
    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
    double dAtm_;
    double dCol_;

private:

    std::shared_ptr<VolcanicColumnProfile> columnProfile_;
    Double3 ventPosition_;

};
inline Double3 DiffusionWoods::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double r = sqrt( pow( ( pos.x_ - ventPosition_.x_ ), 2.0 ) + pow( ( pos.y_ - ventPosition_.y_ ), 2.0 ) );
    double b = ( columnProfile_->getL( pos.z_ ) ) / 2.0;
    double D = dCol_;
    if( r > 2.0 * b ) D = dAtm_;
    return isotropicDiffusion( D, dt );
}
/*==========================================================================================================*/

class DirectionalDiffusionWoods: public SpeedFunctor{
public:
    DirectionalDiffusionWoods(std::shared_ptr<VolcanicColumnProfile> columnProfile, Double3 ventPosition):SpeedFunctor("tetras_anisotropic_diffusion_woods"),
        columnProfile_(columnProfile), ventPosition_(ventPosition){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
    // diffusion coefficient
    double dVAtm_;
    double dHAtm_;
    double dVCol_;
    double dHCol_;

private:
    std::shared_ptr<VolcanicColumnProfile> columnProfile_;
    Double3 ventPosition_;
};
inline Double3 DirectionalDiffusionWoods::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){

    Double3 speedAtm;
    Double3 speedCol;

    double urVAtm = sqrt( ( 2.0 * dVAtm_ ) / dt );
    double urHAtm = sqrt( ( 4.0 * dHAtm_ ) / dt );
    double urVCol = sqrt( ( 2.0 * dVCol_ ) / dt );
    double urHCol = sqrt( ( 4.0 * dHCol_ ) / dt );

    double r;
    double b;

    r = sqrt( pow( ( pos.x_ - ventPosition_.x_ ), 2.0 ) + pow( ( pos.y_ - ventPosition_.y_ ), 2.0 ) );
    b = ( columnProfile_->getL( pos.z_ ) ) / 2.0;

    double alpha = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
    double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
    urVAtm = urVAtm * sin(alpha) + urVAtm * cos(alpha);
    speedAtm = { urHAtm * cos( theta ), urHAtm * sin( theta ), urVAtm };

    if( r > 2.0 * b ) speedCol = { 0.0, 0.0, 0.0 };
    else{
        alpha = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
        theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
        urVCol = urVCol * sin(alpha) + urVCol * cos(alpha);
        speedCol = { urHCol * cos( theta ), urHCol * sin( theta ), urVCol };
    }

    return { speedAtm.x_ + speedCol.x_, speedAtm.y_ + speedCol.y_, speedAtm.z_ + speedCol.z_ };

}

/*==========================================================================================================*/
/*!
* \brief wind model 1 from Carey and Sparks (1986)
*/
class WM1: public SpeedFunctor{
public:
    WM1(std::shared_ptr<VolcanicColumnProfile> columnProfile,Double3 ventPosition):SpeedFunctor("tetras_wm1"),
        columnProfile_(columnProfile),ventPosition_(ventPosition),tropopauseHeight_(-1.0),maximumWindSpeed_(-1.0),windAngle_(0.0){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    void setAtmosphere(Atmosphere* atmosphere);
    void setMaximumWindSpeed(double speed);
    void setWindAngle(double angle);

private:
    std::shared_ptr<VolcanicColumnProfile> columnProfile_;
    Double3 ventPosition_;
    double tropopauseHeight_;
    double maximumWindSpeed_;
    // geographic oritentation (from north, rotation anti-trigo)
    double windAngle_;
    Atmosphere* atmosphere_;
};
inline void WM1::setAtmosphere(Atmosphere* atmosphere){
    atmosphere_ = atmosphere;
    tropopauseHeight_ = atmosphere->getTropopauseHeight();
}
inline void WM1::setMaximumWindSpeed(double speed){maximumWindSpeed_ = speed;}
inline void WM1::setWindAngle(double angle){windAngle_ = angle;}
inline Double3 WM1::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    // defining the base vector
    Double3 speed;
    speed = {0.0,0.0,0.0};

    double r;
    double b;
    r = sqrt( pow( ( pos.x_ - ventPosition_.x_ ), 2.0 ) + pow( ( pos.y_ - ventPosition_.y_ ), 2.0 ) );
    b = ( columnProfile_->getL( pos.z_ ) ) / 2.0;

    if( r > 2.0 * b ){
        if(tropopauseHeight_ != -1 || maximumWindSpeed_ != -1){
            if((pos.z_/tropopauseHeight_)<=1.0) speed = {0.0,(pos.z_/tropopauseHeight_)*maximumWindSpeed_,0.0};
            else speed ={0.0,0.75*maximumWindSpeed_,0.0};
        }
        speed = {cos((360.0-windAngle_)*(pivalue/180.0))*speed.x_-sin((360.0-windAngle_)*(pivalue/180.0))*speed.y_,sin((360.0-windAngle_)*(pivalue/180.0))*speed.x_+cos((360.0-windAngle_)*(pivalue/180.0))*speed.y_,speed.z_};
    }
    return speed;
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/*!
* \brief wind model 2 from : Holasek and Self, 1995 ; Carey and Sigurdsson, 1982; Holasek and Self, 1995;
* Woods et al., 1995; Bonadonna et al., 2002a
*/
class WM2: public SpeedFunctor{
public:
    WM2(std::shared_ptr<VolcanicColumnProfile> columnProfile,Double3 ventPosition):SpeedFunctor("tetras_wm2"),
        columnProfile_(columnProfile),ventPosition_(ventPosition),tropopauseHeight_(-1.0),topTropopauseHeight_(-1.0),maximumWindSpeed_(-1.0),windAngle_(0.0){}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    void setAtmosphere(Atmosphere* atmosphere);
    void setMaximumWindSpeed(double speed);
    void setWindAngle( double angle );

private:
    std::shared_ptr<VolcanicColumnProfile> columnProfile_;
    Double3 ventPosition_;
    double tropopauseHeight_;
    double topTropopauseHeight_;
    double maximumWindSpeed_;
    // geographic oritentation (from north, rotation anti-trigo)
    double windAngle_;
    Atmosphere* atmosphere_;
};
inline void WM2::setAtmosphere(Atmosphere* atmosphere){
    atmosphere_ = atmosphere;
    tropopauseHeight_ = atmosphere->getTropopauseHeight();
    topTropopauseHeight_ = atmosphere->getTopSteadyHeight();
}
inline void WM2::setMaximumWindSpeed(double speed){maximumWindSpeed_ = speed;}
inline void WM2::setWindAngle( double angle ){windAngle_ = angle;}
inline Double3 WM2::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    // defining the base vector
    Double3 speed;
    speed = {0.0,0.0,0.0};

    double r;
    double b;
    r = sqrt( pow( ( pos.x_ - ventPosition_.x_ ), 2.0 ) + pow( ( pos.y_ - ventPosition_.y_ ), 2.0 ) );
    b = ( columnProfile_->getL( pos.z_ ) ) / 2.0;

    if( r > 2.0 * b ){
        if(tropopauseHeight_ != -1 || maximumWindSpeed_ != -1){
            if((pos.z_/tropopauseHeight_)<=1.0) speed = {0.0,(pos.z_/tropopauseHeight_)*maximumWindSpeed_,0.0};
            else{
                speed = {0.0,(0.9*(1-((pos.z_ - tropopauseHeight_)/(topTropopauseHeight_ - tropopauseHeight_))) + 0.1)*maximumWindSpeed_,0.0};
                if(pos.z_ > topTropopauseHeight_) speed = {0.0, 0.1*maximumWindSpeed_, 0.0};
            }
        }
        speed = { cos( (360.0-windAngle_)*(pivalue/180.0) ) * speed.x_ - sin( (360.0-windAngle_)*(pivalue/180.0) ) * speed.y_, sin( (360.0-windAngle_)*(pivalue/180.0) ) * speed.x_ + cos( (360.0-windAngle_)*(pivalue/180.0) ) * speed.y_, speed.z_ };
    }
    return speed;
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/*!
* \brief sedimentation speed , bagheri
*/
class SedimentationSpeed: public SpeedFunctor{
public:
    SedimentationSpeed():SpeedFunctor("tetras_sedimentation"){};
    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
};



inline Double3 SedimentationSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    std::shared_ptr<TephraParticleFamily> *fam;
    fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
    return {0.0, 0.0, -sedimentation(pos.z_, (*fam)->diameter_, (*fam)->density_) };
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/*!
* \brief sedimentation speed , Bonadonna and Phillips, 2003
*/
class SedimentationSpeedOld: public SpeedFunctor{
public:
    SedimentationSpeedOld():SpeedFunctor("tetras_sedimentation_old"){};
    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);
};



inline Double3 SedimentationSpeedOld::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    std::shared_ptr<TephraParticleFamily> *fam;
    fam = reinterpret_cast<std::shared_ptr<TephraParticleFamily>* >(&family);
    return {0.0, 0.0, -sedimentationOld(pos.z_, (*fam)->diameter_, (*fam)->density_) };
}
/*==========================================================================================================*/

class UmbrellaCloudRevised: public SpeedFunctor{
public:
    UmbrellaCloudRevised(Double2 focusSource, double hb, double ht, double tropopauseHeight, double Vwind, double theta, double Ug, double R): SpeedFunctor("tetras_umbrella_cloud_revised"),
    focusSource_(focusSource), hb_(hb), ht_(ht), tropopauseHeight_(tropopauseHeight), Vwind_(Vwind), theta_(theta), Ug_(Ug), R_(R) {}

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:
    Double2 focusSource_;
    double hb_;
    double ht_;
    double tropopauseHeight_;
    double Vwind_;
    double theta_;
    double Ug_;
    double R_;

};

inline Double3 UmbrellaCloudRevised::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){

    Double3 velocity = {0.0, 0.0, 0.0};
    double lambda = 1.0;
    if(pos.z_<=ht_ && pos.z_ >=hb_ && tropopauseHeight_!=-1.0){
        //std::cout << pos.z_ << ", ";
        double N = 0.02;
        if(pos.z_ < tropopauseHeight_) N = 0.01;
        Double3 rVector = {pos.x_-focusSource_.x_, pos.y_-focusSource_.y_};
        double rNorm = sqrt(pow(rVector.x_,2.0) + pow(rVector.y_,2.0));
        Double3 rUnitVector = rVector / rNorm;

        // avoid divide by 0 on r=0, function linear near 0
        double velocityNorm;

        if(rNorm < 0.1) velocityNorm = rNorm;
        else{
            double A = lambda * N * Vwind_ * cos(theta_) * rNorm;
            double B = Ug_ * lambda * N * pow(R_, 2.0);
            double b2 = 2.0*pow(R_, 2.0);
            double C = -(pow(rNorm, 2.0) / b2);
            double D = exp(C);
            double E = 1.0 - D;
            double F = A/2.0 + (B/4.0) * (E/rNorm);
            velocityNorm = pow(F, 1.0/2.0); 
        }
        velocity = velocityNorm * rUnitVector;

    } 
    return velocity;
}

/*==========================================================================================================*/
/*!
* \brief umbrella cloud velocity , Bonadonna and Phillips, 2003
*/
class UmbrellaCloudSpeed: public SpeedFunctor{
public:
    UmbrellaCloudSpeed(): SpeedFunctor("tetras_umbrella_cloud"),
        focusSource_( { 0.0, 0.0 } ), hb_( 0.0 ), ht_( -1.0 ), tropopauseHeight_( -1.0 ), e_( 0.5 ){
        epsilon_ = 2.0 * pivalue * sqrt( ( 2.0 - pow( e_, 2.0 ) ) / ( 2.0 * pow( ( 1.0 + e_ ) , 2.0 ) ) );
    };

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    // mandatory set functions
    void setFocusSource(Double2 source);
    void setHb(double hb);
    void setHt(double ht);
    void setTropopauseHeight(double height);

    void setDx( double dx );
    void setDt( double dt );

    // optional set functions
    void setE(double e);

private:
    Double2 focusSource_;
    double hb_;
    double ht_;
    double tropopauseHeight_;
    double e_;

    double epsilon_;
    double q_;

    double dx_;
};
inline void UmbrellaCloudSpeed::setFocusSource(Double2 source){
    focusSource_ = source;
}
inline void UmbrellaCloudSpeed::setHb(double hb){
    hb_ = hb;
}
// bursik et al 1992 (eq. 6)
inline void UmbrellaCloudSpeed::setHt( double ht ){
    ht_ = ht;
    q_ = pow( ( ( ht_ / 1000.0 ) / 0.287 ) , 5.2632 );
}

inline void UmbrellaCloudSpeed::setTropopauseHeight(double height){
    assert(height <= 20000.0);
    tropopauseHeight_ = height;
}
inline void UmbrellaCloudSpeed::setE( double e ){
    assert( e >= 0.0 && e <= 1.0 );
    e_ = e;
    epsilon_ = 2.0 * pivalue * sqrt( ( 2.0 - pow( e_, 2.0 ) ) / ( 2.0 * pow( ( 1.0 + e_ ) , 2.0 ) ) );
}

inline void UmbrellaCloudSpeed::setDx( double dx ) { dx_ = dx; }

inline Double3 UmbrellaCloudSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    const double lambda = 0.15;
    Double3 speed;
    double speedNorm;
    double N;
    Double2 rVector;
    Double2 rUnitVector;
    double rNorm;

    if(pos.z_<=ht_ && pos.z_ >=hb_ && tropopauseHeight_!=-1.0){
        if(pos.z_ < tropopauseHeight_) N = 0.01;
        else N = 0.02;
        rVector = {pos.x_-focusSource_.x_, pos.y_-focusSource_.y_};
        rNorm = sqrt(pow(rVector.x_,2.0) + pow(rVector.y_,2.0));
        // on limite la vitesse des particules trop pres du centre
        if( rNorm > 10.0 ){
            // compute the speed
            speedNorm = sqrt( ( lambda * N * q_ ) / ( epsilon_ * rNorm ) );
            rUnitVector = rVector / rNorm;
            speed.x_ = rUnitVector.x_ * speedNorm;
            speed.y_ = rUnitVector.y_ * speedNorm;
            // comment for 3d simulation !!!
            //speed.y_ = 0.0;
            speed.z_ = 0.0;
        }
        // if rNorm == 0.0, take the maximum acceptable speed (ie. dx) with a random orientation
        // we take speedNorm<dx because of the speed assertion
        else {
            speed = { ( dx_ / dt ) / 2.0, 0.0, 0.0 };
            double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
            // rotation (around y)
            //speed ={cos(theta)*speed.x_+sin(theta)*speed.z_, speed.y_, -sin(theta)*speed.x_+cos(theta)*speed.z_};
            // rotation around z
            speed = { cos( theta ) * speed.x_ - sin( theta ) * speed.y_, sin( theta ) * speed.x_ + cos( theta ) * speed.y_, speed.z_ };
            //speed = {0.0, 0.0, 0.0};
            //speed = {10.0, 0.0, 0.0};
        }
    }
    else{
        // not in the umbrella cloud region
        speed = { 0.0 , 0.0 , 0.0 };
    }
    return speed;
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/*!
* \brief umbrella cloud velocity , Bonadonna and Phillips, 2003
* with a gaussian profile
*/
class UmbrellaCloudSpeedGaussian: public SpeedFunctor{
public:
    UmbrellaCloudSpeedGaussian(): SpeedFunctor("tetras_umbrella_cloud_gaussian"),
        focusSource_( { 0.0, 0.0 } ), hb_( 0.0 ), ht_( -1.0 ), tropopauseHeight_( -1.0 ), e_( 0.5 ){
        epsilon_ = 2.0 * pivalue * sqrt( ( 2.0 - pow( e_, 2.0 ) ) / ( 2.0 * pow( ( 1.0 + e_ ) , 2.0 ) ) );
    };

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

    // mandatory set functions
    void setFocusSource(Double2 source);
    void setHb(double hb);
    void setHt(double ht);
    void setTropopauseHeight(double height);

    void setDx( double dx );
    void setDt( double dt );

    // optional set functions
    void setE(double e);

private:
    Double2 focusSource_;
    double hb_;
    double ht_;
    double tropopauseHeight_;
    double e_;

    double epsilon_;
    double q_;

    double dx_;
};
inline void UmbrellaCloudSpeedGaussian::setFocusSource(Double2 source){
    focusSource_ = source;
}
inline void UmbrellaCloudSpeedGaussian::setHb(double hb){
    hb_ = hb;
}
// bursik et al 1992 (eq. 6)
inline void UmbrellaCloudSpeedGaussian::setHt( double ht ){
    ht_ = ht;
    q_ = pow( ( ( ht_ / 1000.0 ) / 0.287 ) , 5.2632 );
}

inline void UmbrellaCloudSpeedGaussian::setTropopauseHeight(double height){
    assert(height <= 20000.0);
    tropopauseHeight_ = height;
}
inline void UmbrellaCloudSpeedGaussian::setE( double e ){
    assert( e >= 0.0 && e <= 1.0 );
    e_ = e;
    epsilon_ = 2.0 * pivalue * sqrt( ( 2.0 - pow( e_, 2.0 ) ) / ( 2.0 * pow( ( 1.0 + e_ ) , 2.0 ) ) );
}

inline void UmbrellaCloudSpeedGaussian::setDx( double dx ) { dx_ = dx; }

inline Double3 UmbrellaCloudSpeedGaussian::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    const double lambda = 0.15;
    Double3 speed;
    double speedNorm;
    double N;
    Double2 rVector;
    Double2 rUnitVector;
    double rNorm;

    if(pos.z_<=ht_ && pos.z_ >=hb_ && tropopauseHeight_!=-1.0){
        // compute N
        if(pos.z_ < tropopauseHeight_) N = 0.01;
        else N = 0.02;
        // compute r
        rVector = {pos.x_-focusSource_.x_, pos.y_-focusSource_.y_};
        rNorm = sqrt(pow(rVector.x_,2.0) + pow(rVector.y_,2.0));
        // on limite la vitesse des particules trop pres du centre
        if( rNorm > 10.0 ){
            // compute the speed
            speedNorm = sqrt( ( lambda * N * q_ ) / ( epsilon_ * rNorm ) );
            rUnitVector = rVector / rNorm;
            speed.x_ = rUnitVector.x_ * speedNorm;
            speed.y_ = rUnitVector.y_ * speedNorm;
            // comment for 3d simulation !!!
            //speed.y_ = 0.0;
            speed.z_ = 0.0;
        }
        // if rNorm == 0.0, take the maximum acceptable speed (ie. dx) with a random orientation
        // we take speedNorm<dx because of the speed assertion
        else {
            speed = { ( dx_ / dt ) / 2.0, 0.0, 0.0 };
            double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * pivalue;
            // rotation (around y)
            //speed ={cos(theta)*speed.x_+sin(theta)*speed.z_, speed.y_, -sin(theta)*speed.x_+cos(theta)*speed.z_};
            // rotation around z
            speed = { cos( theta ) * speed.x_ - sin( theta ) * speed.y_, sin( theta ) * speed.x_ + cos( theta ) * speed.y_, speed.z_ };
            //speed = {0.0, 0.0, 0.0};
            //speed = {10.0, 0.0, 0.0};
        }
    }
    else{
        // not in the umbrella cloud region
        speed = { 0.0 , 0.0 , 0.0 };
    }
    // just multiplied by a gaussian
    double r = fabs(pos.z_ - (hb_ + (ht_ - hb_) / 2.0));
    double b2 = pow((ht_ - hb_) / 2.0, 2.0) / 2.0;
    speed = speed * pow(boost::math::constants::e<double>(),-((r*r)/(b2)));
    return speed;
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/*!
* \brief volcanic column speed
*/
class VolcanicColumnSpeed: public SpeedFunctor{
public:
    VolcanicColumnSpeed(std::shared_ptr<VolcanicColumnProfile> columnProfile, Double3 ventPosition): SpeedFunctor("tetras_woods_plume"),
        columnProfile_(columnProfile), ventPosition_(ventPosition){};

    Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

private:
    std::shared_ptr<VolcanicColumnProfile> columnProfile_;
    Double3 ventPosition_;
};
inline Double3 VolcanicColumnSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
    double zSpeed, U, r, b;
    U = columnProfile_->getU(pos.z_);
    r = sqrt( pow( ( pos.x_ - ventPosition_.x_ ), 2.0 ) + pow( ( pos.y_ - ventPosition_.y_ ), 2.0 ) );
    b = ( columnProfile_->getL( pos.z_ ) ) / 2.0;
    zSpeed = U * exp( -( pow( ( r / b ) , 2.0 ) ) );
    //if(r>2*b) zSpeed = 0.0;
    // could happen for particles pos.z>ht and pos.x, pos.y == ventPosition
    if(std::isnan(zSpeed)){
        zSpeed = 0.0;
    }
    return { 0.0, 0.0 , zSpeed };
}
/*==========================================================================================================*/

#endif
