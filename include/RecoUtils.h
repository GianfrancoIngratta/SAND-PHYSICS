#ifndef SANDRECO_UTILS_H
#define SANDRECO_UTILS_H

#include <string>

#include "TString.h"
#include "struct.h"
#include "TLorentzVector.h"

namespace RECO{

const double DEFAULT_DOUBLE = 0.;

const double DEFAULT_NEGATIVE = -9999999.;

const TLorentzVector DEFAULT_TLV = {DEFAULT_DOUBLE,
                                    DEFAULT_DOUBLE,
                                    DEFAULT_DOUBLE,
                                    DEFAULT_DOUBLE};

namespace SANDRECO{
}

namespace FASTRECO{
}

struct RecoPerformance
{
    std::vector<double> muon_momentum_residuals;
    
    double muon_mometum_resolution  = DEFAULT_NEGATIVE;

    double muon_reco_efficiency     = DEFAULT_NEGATIVE;

    double proton_reco_efficiency   = DEFAULT_NEGATIVE;

    double neutron_reco_efficiency  = DEFAULT_NEGATIVE;
};

}


#endif