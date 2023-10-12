#ifndef GENIE_UTILS_H
#define GENIE_UTILS_H

#include <string>

#include "TString.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

namespace GenieUtils{
std::string PDG2Name(int pdg);
}

#endif