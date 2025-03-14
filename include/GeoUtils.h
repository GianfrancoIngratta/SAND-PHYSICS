#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#include "TGeoManager.h"
#include "TLorentzVector.h"
#include "struct.h"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

namespace GeoUtils{//GeoUtils

// SAND inner volume
const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

const TVector3 SAND_CENTER = {SAND_CENTER_X, SAND_CENTER_Y, SAND_CENTER_Z};

const double SAND_INNER_VOL_X_LENGTH = 3380.0;

const double SAND_INNER_VOL_DIAMETER = 4000.0;

extern TGeoManager* geo;

bool IsInSANDInnerVol(std::string units, double x, double y, double z);

TVector3 MasterToSAND(const TVector3& point);

double IntersectWithInfCylinder(const TVector3& vtx, 
                                          const TVector3& direction, 
                                          double diameter);

double IntersectWithTube(std::string units, 
                         const TVector3& vertex, 
                         const TVector3& direction,
                         double tube_diameter,
                         double tube_length);

std::string InteractionVolume_short(const std::string& detailed_name);

ROOT::VecOps::RVec<std::string> InteractionVolume_short_v(const ROOT::VecOps::RVec<std::string>& detailed_name);

namespace ECAL{ // ECAL

const double Module_angle = TMath::Pi() / 12.;

const double Module_0_starting_angle = (TMath::Pi() - GeoUtils::ECAL::Module_angle) * 0.5;

const double vlfb = 5.85 * 1e-3; // ns/mm velocity signal in fiber, 5.85 ns/m

double const attpassratio = 0.187;

// photoelectron/counts = 0.25
const double pe2ADC = 1/.25;

// Average number of photoelectrons = 25*Ea(MeV)
// corrected to 18.5 to have mean number of pe of 40
// for mip crossing in the middle of barrel module
const double e2pe = 18.5;

const double thickness = 230; // mm

// get fiber attenuation factor.
// It depends on distance from pmt (d)
// and planeID (planes have different fibers)
double AttenuationFactor(double d, int planeID);
const double p1 = 0.35;
const double atl1 = 500.;
const double atl2_01 = 4300.0;
const double atl2_2 = 3800.0;
const double atl2_34 = 3300.0;

double XfromTDC(const double tdc1, const double tdc2);

double TfromTDC(const double tdc1, const double tdc2, const double length);

double EfromADC(double adc1, double adc2, double d1,
                                       double d2, int planeID);

double EfromADCsingle(double adc, double f);

int GetModuleIdFromPoint(std::string units, TVector3 point);

bool isCellComplete(const dg_cell& cell, int& side);

bool is_ecal_barrel(const TString& volume_name);

} // ECAL

namespace STT{//STT

}//STT

namespace DRIFT{//DRIFT

// does not include frames
const double SAND_TRACKER_X_LENGTH = 3220.0;

// mm (approx), this includes 9 C3H6 mod + 1 C mod + clearances
// const double SAND_SUPEMOD_Z_THICK = 364.6;

// const double SAND_TRACKER_Z_START = 22952.14; // mm

const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

// !!! FIDUCIAL CUT

// distance in mm from SUPERMOD frames
const double FIDUCIAL_CUT = 300.; // mm

// height in mm
const double SUPERMOD_Y_HEIGHT[5] = {3755.16996258, // A1, A2
                                  3550.97659024, // B1, B2
                                  3129.97277395, // C1, C2
                                  1262.30218417,  // X1    
                                  2466.05424949 // X0    

};

bool IsInFiducialVolume(std::string volName, std::string units, double x, double y, double z);

bool IsInSMODFiducialVol(int smod, double x, double y);

}//DRIFT

}//GeoUtils

#endif