#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#include "TGeoManager.h"
#include "TLorentzVector.h"

namespace GeoUtils{//GeoUtils

// SAND inner volume
const double SAND_INNER_VOL_X_LENGTH = 3380.0;
const double SAND_INNER_VOL_DIAMETER = 4000.0;

extern TGeoManager* geo;

bool IsInSANDInnerVol(std::string units, double x, double y, double z);

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
const double FIDUCIAL_CUT = 50.;

// height in mm
const double SUPERMOD_Y_HEIGHT[5] = {3755.16996258, // A1, A2
                                  3550.97659024, // B1, B2
                                  3129.97277395, // C1, C2
                                  1262.30218417,  // X1    
                                  2466.05424949 // X0    

};

bool IsInFiducialVolume(std::string volName, std::string units, double x, double y, double z);

bool IsInSMODFiducialVol(int smod, double x, double y, double z);

}//DRIFT

}//GeoUtils

#endif