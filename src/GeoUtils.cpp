#include <iostream>
#include <TObjString.h>

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

std::string GeoUtils::InteractionVolume_short(const std::string& detailed_name) {
    if (detailed_name.find("C3H6Target") != std::string::npos) {
        return "C3H6_Target";
    } else if (detailed_name.find("CTarget") != std::string::npos) {
        return "C_Target";
    } else if (detailed_name.find("Mylar") != std::string::npos) {
        return "Mylar";
    } else if (detailed_name.find("Yoke") != std::string::npos) {
        return "Yoke";
    } else if (detailed_name.find("Solenoid") != std::string::npos) {
        return "Yoke";
    } else if (detailed_name.find("ECAL") != std::string::npos) {
        return "ECAL";
    } else if (detailed_name.find("GRAIN") != std::string::npos) {
        if (detailed_name.find("LAr") != std::string::npos) {
            return "GRAIN_LAr";
        } else {
            return "GRAIN vessels";
        }
    } else if (detailed_name.find("Frame") != std::string::npos) {
        return "Supermodules_Frames";
    } else if (detailed_name.find("Drift") != std::string::npos) {
        return "Drift_gas";
    } else if (detailed_name.find("volSAND") != std::string::npos) {
        return "Supermodules_Frames";
    } else {
        return "Other";
    }
}

// ECAL functions ____________________

bool GeoUtils::ECAL::is_ecal_barrel(const TString& volume_name){
  // something like: volECALActiveSlab_21_PV_0
  return volume_name.Contains("volECAL") == true &&
         volume_name.Contains("Active") == true &&
         volume_name.Contains("end") == false;
}


double GeoUtils::ECAL::TfromTDC(const double tdc1, const double tdc2, const double length){
    // vlfb 5.85 [ns/mm]
    return 0.5 * (tdc1 + tdc2 - GeoUtils::ECAL::vlfb * length);
}

double GeoUtils::ECAL::XfromTDC(const double tdc1, const double tdc2){
    // vlfb 5.85 [ns/mm]
    return 0.5 * (tdc1 - tdc2) / GeoUtils::ECAL::vlfb; // mm
}

double GeoUtils::ECAL::EfromADCsingle(double adc, double f){
  return adc / (f * GeoUtils::ECAL::attpassratio * GeoUtils::ECAL::pe2ADC *
                GeoUtils::ECAL::e2pe);
}

double GeoUtils::ECAL::AttenuationFactor(double d, int planeID){
  /*
       dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
         d    distance from photocatode - 2 cells/cell; d1 and d2
        atl1  50. cm
        atl2  430 cm planes 1-2    innermost plane is 1
              380 cm plane 3
              330 cm planes 4-5
         p1   0.35
  */
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    case 1:
      atl2 = GeoUtils::ECAL::atl2_01;
      break;

    case 2:
      atl2 = GeoUtils::ECAL::atl2_2;
      break;

    case 3:
    case 4:
      atl2 = GeoUtils::ECAL::atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
    
    }
    return GeoUtils::ECAL::p1 * TMath::Exp(-d / GeoUtils::ECAL::atl1) + (1. - GeoUtils::ECAL::p1) * TMath::Exp(-d / atl2);
}

double GeoUtils::ECAL::EfromADC(double adc1, double adc2, double d1,
                                       double d2, int planeID){
  double f1 = GeoUtils::ECAL::AttenuationFactor(d1, planeID);
  double f2 = GeoUtils::ECAL::AttenuationFactor(d2, planeID);

  double const attpassratio = 0.187; //new
  double e = 0.5 * (adc1 / f1 + adc2 / f2) / (attpassratio * GeoUtils::ECAL::pe2ADC * GeoUtils::ECAL::e2pe); 
  return e;
}

bool GeoUtils::ECAL::isCellComplete(const dg_cell& cell, int& side){
    if(cell.ps1.size()==0){
        side = 1;
        return false;
    }else if(cell.ps2.size()==0){
        side = 2;
        return false;
    }else{
        side = 0;
        return true;
    }
}

TVector3 GeoUtils::MasterToSAND(const TVector3& point){
    return point - GeoUtils::SAND_CENTER;
}

int GeoUtils::ECAL::GetModuleIdFromPoint(std::string units, TVector3 point){
    /*
        Given a point in the space, find the ECAL module
        that contains the point
    */
    // convert to mm
    if(units == "m"){
        point *= 1e3;
    }else if(units == "cm"){
        point *= 10.;
    }else{}
    // first get coordinates wrt SAND center
    TVector3 Point2SAND = GeoUtils::MasterToSAND(point);
    // Check if point is in one endcap
    if(Point2SAND.X() <= - GeoUtils::SAND_INNER_VOL_X_LENGTH / 2.0){
        return 40;
    } else if (Point2SAND.X() >= GeoUtils::SAND_INNER_VOL_X_LENGTH / 2.0){
        return 30;
    } else {
        // find the coordinates wrt an axis that starts from module 0
        double point_z = Point2SAND.Z() * cos(GeoUtils::ECAL::Module_0_starting_angle) + Point2SAND.Y() * sin(GeoUtils::ECAL::Module_0_starting_angle);
        double point_y = - Point2SAND.Z() * sin(GeoUtils::ECAL::Module_0_starting_angle) + Point2SAND.Y() * cos(GeoUtils::ECAL::Module_0_starting_angle);
        // get angle point with z axis from [0,2pi]
        double angle = (point_y >= 0.) ? 
                            TMath::ATan2(point_y, point_z) : TMath::ATan2(point_y, point_z) + 2 * TMath::Pi();
        
        int module = angle / GeoUtils::ECAL::Module_angle;
        return module;
    }
}

// ECAL functions ____________________

bool GeoUtils::DRIFT::IsInSMODFiducialVol(int smod, double x, double y){
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
        return GeoUtils::DRIFT::IsInSMODFiducialVol(0, x, y);
    }else if(volName_.Contains("_B")){ // one of 2 supermod B
        return GeoUtils::DRIFT::IsInSMODFiducialVol(1, x, y);
    }else if(volName_.Contains("_C")){ // one of 2 supermod C
        return GeoUtils::DRIFT::IsInSMODFiducialVol(2, x, y);
    }else if(volName_.Contains("_X0")){ // dwstream supermod X0
        return GeoUtils::DRIFT::IsInSMODFiducialVol(3, x, y);
    }else if(volName_.Contains("_X1")){ // dwstream supermod X1
        return GeoUtils::DRIFT::IsInSMODFiducialVol(4, x, y);
    }else{
        return false;
    }
}