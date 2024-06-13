#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#include "TGeoManager.h"
#include "TLorentzVector.h"

namespace GeoUtils{//GeoUtils

extern TGeoManager* geo;

namespace STT{//STT

}//STT

namespace DRIFT{//DRIFT

// !! HARDCODED (I know you are judging me, stop it)
// does not include frames
const double SAND_TRACKER_X_LENGTH = 3220.0;

// mm (approx), this includes 9 C3H6 mod + 1 C mod + clearances
const double SAND_SUPEMOD_Z_THICK = 364.6;

const double SAND_TRACKER_Z_START = 22952.14; // mm

const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

const std::vector<double> SUPERMOD_LENGTHS = {3129.97277395, // A1
                                              3550.97659024, // B1
                                              3755.16996258, // C1
                                              3755.16996258, // C2
                                              3550.97659024, // B2
                                              3129.97277395, // A2
                                              2466.05424949, // D
                                              1262.30218417};// F

bool IsInFiducialVolume(double x, double y, double z);

}//DRIFT

}//GeoUtils

#endif