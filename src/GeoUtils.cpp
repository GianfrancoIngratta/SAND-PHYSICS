#include <iostream>

#include "GeoUtils.h"
#include "TString.h"

bool GeoUtils::IsInSANDInnerVol(std::string units, double x, double y, double z){
    // conversion to mm
    if(units == "m"){
        x = x*1e3; y = y*1e3; z = z*1e3;
    }else if(units == "cm"){
        x = x*10.; y = y*10.; z = z*10.;
    }else{}
    bool pass_x = (fabs(x - GeoUtils::DRIFT::SAND_CENTER_X) <= GeoUtils::SAND_INNER_VOL_X_LENGTH/2.);
    bool pass_zy = ((y - GeoUtils::DRIFT::SAND_CENTER_Y)*(y - GeoUtils::DRIFT::SAND_CENTER_Y) + 
                    (z - GeoUtils::DRIFT::SAND_CENTER_Z)*(z - GeoUtils::DRIFT::SAND_CENTER_Z) < 
                    (GeoUtils::SAND_INNER_VOL_DIAMETER/2.)*(GeoUtils::SAND_INNER_VOL_DIAMETER/2.));
    return pass_x * pass_zy;
}

bool GeoUtils::DRIFT::IsInSMODFiducialVol(int smod, double x, double y, double z){
    // !!! ALL LENGTHS REQUIRED IN MM !!!
    bool pass_x = (fabs(x - GeoUtils::DRIFT::SAND_CENTER_X) <= GeoUtils::DRIFT::SAND_TRACKER_X_LENGTH/2. - GeoUtils::DRIFT::FIDUCIAL_CUT);
    bool pass_zy = (fabs(y - GeoUtils::DRIFT::SAND_CENTER_Y) <= GeoUtils::DRIFT::SUPERMOD_Y_HEIGHT[smod]/2. - GeoUtils::DRIFT::FIDUCIAL_CUT);
    return pass_x * pass_zy;
}

bool GeoUtils::DRIFT::IsInFiducialVolume(std::string volName, std::string units, double x, double y, double z){
    // conversion to mm
    if(units == "m"){
        x = x*1e3; y = y*1e3; z = z*1e3;
    }else if(units == "cm"){
        x = x*10.; y = y*10.; z = z*10.;
    }else{}
    TString volName_ = volName;
    if(volName_.Contains("Frame")){
        return false;
    }else if(volName_.Contains("_A")){ // one of 2 supermod A
        return GeoUtils::DRIFT::IsInSMODFiducialVol(0, x, y, z);
    }else if(volName_.Contains("_B")){ // one of 2 supermod B
        return GeoUtils::DRIFT::IsInSMODFiducialVol(1, x, y, z);
    }else if(volName_.Contains("_C")){ // one of 2 supermod C
        return GeoUtils::DRIFT::IsInSMODFiducialVol(2, x, y, z);
    }else if(volName_.Contains("_X0")){ // dwstream supermod X0
        return GeoUtils::DRIFT::IsInSMODFiducialVol(3, x, y, z);
    }else if(volName_.Contains("_X1")){ // dwstream supermod X1
        return GeoUtils::DRIFT::IsInSMODFiducialVol(4, x, y, z);
    }else{
        return false;
    }
}