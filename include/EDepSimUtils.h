#ifndef EDEP_UTILS_H
#define EDEP_UTILS_H

#include "TG4HitSegment.h"
#include "TLorentzVector.h"

namespace EDepUtils{


struct track_hits
{   
    int track_id = -999;

    int pdg = -999;    

    // list of all energy depositions for each single hit
    std::vector<double>  hit_edep = {};
    
    // time and position associated to each hit
    std::vector<TLorentzVector> hit_LorentzVector = {};
};

}

#endif
