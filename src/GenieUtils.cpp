#include "GenieUtils.h"

std::string GenieUtils::PDG2Name(int pdg){
    std::string name = genie::PDGLibrary::Instance()->Find(pdg)->GetName();
    return name;
}
