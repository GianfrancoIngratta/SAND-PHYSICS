#ifndef EDEP_UTILS_H
#define EDEP_UTILS_H

#include "TG4HitSegment.h"
#include "TLorentzVector.h"

namespace EDepUtils{

struct track_hits
{   
    int vertex_number = -999;

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

namespace SPILL{
    const int nof_batches_in_spill = 6;
    const int nof_bunches_in_batch = 84;
    const double bunch_duration = 1.;
    const double bunch_separation = 19.;
    
    const double batch_duration = bunch_separation * nof_bunches_in_batch; // 1.6 ms
    const double spill_duration = nof_batches_in_spill * batch_duration; // 9.6 ms
}

#endif
