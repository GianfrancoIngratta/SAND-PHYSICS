#include "GeoUtils.h"

bool GeoUtils::DRIFT::IsInFiducialVolume(double x, double y, double z){
    /*
        select events in the sand fiducial volume at 10 cm far from the frames:
        - 5 cm along x from the frame edge;
    */
    // convert m to mm
    x = x*1000;
    y = y*1000;
    z = z*1000;
    bool pass_x = (fabs(x - GeoUtils::DRIFT::SAND_CENTER_X) <= GeoUtils::DRIFT::SAND_TRACKER_X_LENGTH/2. - 50);
    // get supermod number
    int smod = (z - GeoUtils::DRIFT::SAND_TRACKER_Z_START)/GeoUtils::DRIFT::SAND_SUPEMOD_Z_THICK;
    // event should be 5 cm far from supermod frame
    bool pass_zy = (fabs(y - GeoUtils::DRIFT::SAND_CENTER_Y) <= GeoUtils::DRIFT::SUPERMOD_LENGTHS[smod]/2. - 50);
    return pass_x * pass_zy;
}