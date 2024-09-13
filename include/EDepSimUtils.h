#ifndef EDEP_UTILS_H
#define EDEP_UTILS_H

#include "TG4HitSegment.h"
#include "TLorentzVector.h"

namespace EDepUtils{

struct track_hits
{   
    int track_id = -999;

    int pdg = -999;

    // list of all hit hindex associated with the track
    std::vector<int> h_indices = {};

    // list of all energy depositions for each single hit
    std::vector<double>  hit_edep = {};
    
    // time and position associated to each hit
    std::vector<TLorentzVector> hit_LorentzVector = {};

    std::vector<std::string> hit_volName = {};
};

}

#endif
