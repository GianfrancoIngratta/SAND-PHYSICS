#ifndef GENIE_UTILS_H
#define GENIE_UTILS_H

#include <string>
#include <mutex>

#include "TString.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

namespace GenieUtils{

static auto database = genie::PDGLibrary::Instance();

std::string PDG2Name(int pdg);

struct event_topology
{
    int NofProtons = 0;

    int NofNeutrons = 0;
    // pi+- pi0
    int NofPiP = 0;

    int NofPiM = 0;

    int NofPi0 = 0;

    // e+- gamma
    int NofEM = 0;

    // mesons : K0, KL, K+- ...
    int NofKaons = 0;

    // baryonResonance : Delta, N,
    /*
        see https://github.com/GENIE-MC/Generator/blob/
        d43b421d7f80fa6ef23f5da99ac425fe527cf324/src/Framework
        /ParticleData/BaryonResUtils.cxx#L299
    */
    int NofResonance = 0;

    // other;
    // C12, 016, Recoiled Nuclei ...
    int NofOther = 0;

    // long int unique_code(){
    int unique_code(){
    /*
        encode event topology with a long int
    */                                  
        // long int code = 0;
        int code = 0;
        // long int factor = 1;
        int factor = 1;
        for(auto n : {NofProtons, // factor : 1
                      NofNeutrons, // factor : 10
                      NofPiP, // factor : 10^2
                      NofPiM, // // factor : 10^3
                      NofPi0,
                      NofKaons,
                      NofResonance,
                      NofOther}){
        code += (n > 9) ? 9 * factor : n * factor;
        factor *= 10;
        }
        return code;
    }

std::string Name(){
    long int c = unique_code();
    std::string n;
    switch(c){
        case 0:
            n = "no particles";
            break;
        case 1:
            n = "p";
            break;
        case 10:
            n = "n";
            break;
        case 11:
            n = "p + n";
            break;
        case 101:
            n = "p + pi+";
            break;
        case 1001:
            n = "p + pi-";
            break;
        case 10001:
            n = "p + pi0";
            break;
        case 110:
            n = "n + pi+";
            break;
        case 1010:
            n = "n + pi-";
            break;
        case 10010:
            n = "n + pi0";
            break;
        case 111:
            n = "p + n + pi+";
            break;
        case 1011:
            n = "p + n + pi-";
            break;
        case 10011:
            n = "p + n + pi0";
            break;
        default:
            n = "Other";
            break;
    }
    return n;
}

};

void UpdateTopology(event_topology& t, int pdg);

}

#endif