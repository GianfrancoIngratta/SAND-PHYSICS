#include "RDataFrameUtils.h"

#include "TChain.h"
#include "TMath.h"

TGeoManager* geo = nullptr;

//RDFUtils______________________________________________________________________

ROOT::RDataFrame RDFUtils::InitDF(TString production, 
                                  const char* tree_name, 
                                  unsigned int file_index_start, 
                                  unsigned int file_index_stop){
    
    TChain* chain = new TChain(tree_name, tree_name);
    TChain* chain_friend = new TChain("gtree", "gtree");

    if(production.Contains("*")){
        // chain a given number of file of a production
        for(unsigned int i = file_index_start; i < file_index_stop; i++){
            auto file = TString::Format(production.ReplaceAll("*","%d"), i, i);
            if (!gSystem->AccessPathName(file.Data())) chain->Add(file);
            if (production.Contains("gtrac")){ // processing genie output
                auto file_friend = TString::Format(file.ReplaceAll(".gtrac.root",".ghep.overlayed.root"));
                if (!gSystem->AccessPathName(file_friend.Data())) chain_friend->Add(file_friend);
            }
        }
    }else{
        chain->Add(production.Data());
        if (production.Contains("gtrac")){
            auto file_friend = TString::Format(production.ReplaceAll(".gtrac.root",".ghep.overlayed.root"));
            LOG("ii", TString::Format("file friend %s",file_friend.Data()));
            if (!gSystem->AccessPathName(file_friend.Data())){
                chain_friend->Add(file_friend);
                LOG("ii", "file_friend Added to chain_friend");
                chain->AddFriend(chain_friend, "gtree");
                LOG("ii", "chain_friend ghep Added to chai gtrac");
            }
        }
    }
    std::cout << "tree name : " << tree_name << " number of events : " << chain ->GetEntries() << "\n";
    
    ROOT::RDataFrame* df = new ROOT::RDataFrame(*chain);

    if(!df){
        std::cout<<"DF NOT initialized in file "<<__FILE__<<" line "<<__LINE__<<"\n";
        throw "";
        }

    return *df;
}

void RDFUtils::PrintColumns(ROOT::RDataFrame& df){

    auto colNames = df.GetColumnNames();
    
    for (auto &&colName : colNames){
        std::cout << "colName : " <<colName<<"\n";
        // if((colName == "gmcrec")) 
        //     std::cout << "colName : " <<colName<<", colType : " <<df.GetColumnType(colName)<<"\n";
        // }
}
}

template<int coord>
ROOT::VecOps::RVec<double> RDFUtils::GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL){

    ROOT::VecOps::RVec<double> v;

    for(auto& TL : vTL) v.push_back(TL[coord]);

    return v;
}

double RDFUtils::GetColumnSum(const ROOT::VecOps::RVec<double>& v){
    double sum = 0.;
    for(auto& value : v) sum+=value;
    return sum;
}

ROOT::VecOps::RVec<double> RDFUtils::VectorDifference(const ROOT::VecOps::RVec<double>& v1,
                                                      const ROOT::VecOps::RVec<double>& v2){
    ROOT::VecOps::RVec<double> diff;
    for(auto i=0u; i<v1.size(); i++) diff.push_back(v1[i]-v2[i]);
    return diff;                                                        
}

TVector3 RDFUtils::TLVectorCrossProduct(const TLorentzVector& v1,
                                        const TLorentzVector& v2){
    TVector3 v = v1.Vect().Cross(v2.Vect());
    return v*(1./v.Mag());                                 
}

ROOT::VecOps::RVec<double> RDFUtils::VectorSubtractConst(const ROOT::VecOps::RVec<double>& v1, double c){
    ROOT::VecOps::RVec<double> diff;
    for(auto i=0u; i<v1.size(); i++) diff.push_back(v1[i]-c);
    return diff;                                                        
}

TLorentzVector RDFUtils::VectorFilterByHighest(const ROOT::VecOps::RVec<double>& filter,
                                               const ROOT::VecOps::RVec<TLorentzVector>& v){
    // return value of v at index position corresponding to index of filter's highest
    if(v.size()==0){
        return {0.,0.,0.,0.};
    }else{                                      
        auto pointer2max = std::max_element(filter.begin(),filter.end());
        auto index2max = std::find(filter.begin(),filter.end(),*pointer2max);
        return v[index2max-filter.begin()];
    }
}

// ROOT::VecOps::RVec<double> RDFUtils::GetResolution(const ROOT::VecOps::RVec<double>& reco,
//                                                    const ROOT::VecOps::RVec<double>& true){

//     ROOT::VecOps::RVec<double> resolutions;
    
//     for (auto i = 0u; i < reco.size(); i++){
//         if(reco[i]==RDFUtils::DEFAULT_DOUBLE){
//             resolutions.push_back(RDFUtils::DEFAULT_NEGATIVE);
//         }else{
//             resolutions.push_back(1. - true[i] / reco[i]);
//         }
//     }
    
//     return resolutions;
// }

//RDFUtils::GENIE_________________________________________________________________

std::string RDFUtils::GENIE::InteractionTarget(const ROOT::VecOps::RVec<int>& pdg){
    // return target with which neutrino interact
    return GenieUtils::PDG2Name(pdg[1]);
}

int RDFUtils::GENIE::NeutrinoFlavor(const ROOT::VecOps::RVec<int>& pdg){
    // return interacting neutrino flavor
    return pdg[0];
}

int RDFUtils::GENIE::NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg){
    return pdg.size();
}

std::string RDFUtils::GENIE::EventType(TObjString& s){
    
    if(s.GetString().Contains("DIS")){
        return "DIS";
    }else if(s.GetString().Contains("RES")){
        return "RES";
    }else if(s.GetString().Contains("QES")){
        return "QES";
    }else if(s.GetString().Contains("COH")){
        return "COH";
    }else if(s.GetString().Contains("MEC")){
        return "MEC";
    }else if(s.GetString().Contains("Unknown")){
        return "Unknown";
    }
    else{
        return "Other";
    }
}

bool RDFUtils::GENIE::IsCC(TObjString& s){
    if(s.GetString().Contains("[NC]")){
        return false;
    }else{
        return true;
    }
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                                           const ROOT::VecOps::RVec<int>& status,
                                                                           const ROOT::VecOps::RVec<double>& P4,
                                                                           const ROOT::VecOps::RVec<double>& X4,
                                                                           const ROOT::VecOps::RVec<int>& FirstMother,
                                                                           const ROOT::VecOps::RVec<int>& LastMother,
                                                                           const ROOT::VecOps::RVec<int>& FirstDaughter,
                                                                           const ROOT::VecOps::RVec<int>& LastDaughter){
    /*
        Construct an istance genie::GHepParticle for each particle in the genie output.
        The class genie::GHepParticle is more helpfull for the analysis
        class reference: https://internal.dunescience.org/doxygen/GHepParticle_8h_source.html#l00177
    */                                                                            
    ROOT::VecOps::RVec<genie::GHepParticle> particles;

    for(auto i=0u; i<pdg.size(); i++){
        
        TLorentzVector momentum = {P4[4*i],P4[4*i+1],P4[4*i+2],P4[4*i+3]};
        TLorentzVector position = {X4[4*i],X4[4*i+1],X4[4*i+2],X4[4*i+3]};

        auto flag_status = static_cast<genie::GHepStatus_t>(status[i]);

        genie::GHepParticle particle(pdg[i], 
                                    flag_status,
                                    FirstMother[i],
                                    LastMother[i],
                                    FirstDaughter[i],
                                    LastDaughter[i],
                                    momentum,
                                    position);
        
        particles.push_back(particle);

        // std::cout<< "particle : status " << particle.Status()
        //         << " mother : " << particle.FirstMother()
        //         << " dougther : " << particle.FirstDaughter() << "\n";
    }

    return particles;                                                                            
}

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    /*
        return a vector of genie::GHepParticle with the same STATUS.
        Reference : https://github.com/GENIE-MC/Generator/blob/master/src/Framework/GHEP/GHepStatus.h
        typedef enum EGHepStatus {
            kIStUndefined                  = -1,
            kIStInitialState               =  0,  generator-level initial state 
            kIStStableFinalState           =  1,  generator-level final state: particles to be tracked by detector-level MC 
            kIStIntermediateState          =  2,
            kIStDecayedState               =  3,
            kIStCorrelatedNucleon          = 10,
            kIStNucleonTarget              = 11,
            kIStDISPreFragmHadronicState   = 12,
            kIStPreDecayResonantState      = 13,
            kIStHadronInTheNucleus         = 14,  hadrons inside the nucleus: marked for hadron transport modules to act on 
            kIStFinalStateNuclearRemnant   = 15,  low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle 
            kIStNucleonClusterTarget       = 16   for composite nucleons before phase space decay
        }
    */
    ROOT::VecOps::RVec<genie::GHepParticle> selected;

    for(auto& p : particles){
        if(p.Status() == STATUS) selected.push_back(p);
    }

    return selected;
}

GenieUtils::event_topology RDFUtils::GENIE::GetInteractionTopology(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    
    GenieUtils::event_topology topology;

    if(particles.size()==0){
        // std::cout<<"topology unique code: "<<topology.unique_code()<<" \n";
        return topology;
    }else{
        for(auto const& particle : particles){
            GenieUtils::UpdateTopology(topology, particle.Pdg());
        }
        // std::cout<<"topology unique code: "<<topology.unique_code()<<" \n";
    }
    return topology;                  
}


template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetParticlesWithPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    // return GHepParticle whose pdg is PDG
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        if(p.Pdg()==PDG) filtered_particles.push_back(p);
    }
    
    return filtered_particles;
}

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::ExcludePDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    // return input vector without particles whose PDG is fabs(14)
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        if(abs(p.Pdg())!=PDG) filtered_particles.push_back(p);
    }
    
    return filtered_particles;
}

ROOT::VecOps::RVec<std::string> RDFUtils::GENIE::PDG2Name(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    // convert particles pdg into names
    ROOT::VecOps::RVec<std::string> names;
    for(const auto& p : particles){
        names.push_back(GenieUtils::PDG2Name(p.Pdg()));
    }
    return names;
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::FinalStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    /*
    Extract the final state hadronic system originating from the hadronic vertex
    (after the intranuclear rescattering step).
    
    Reference: https://github.com/GENIE-MC/Generator/blob/975179a5f2595790c5bb722a2515d18bda67af6f/src/Apps/gNtpConv.cxx#L775

    EXCLUDE from final state:
        - particles with 0 momentum
        - pseudo particles
        - primary lepton
        - e+ e- gamma coming from neutrino vertex, nuclear de-excitation
    */
    ROOT::VecOps::RVec<genie::GHepParticle> final_had_syst;
    for (const auto& particle : particles){
        if(particle.P4()->P() == 0) continue; // exclude particles with fs 0 momentum
        if(particle.Status() == genie::kIStStableFinalState){
            int pdg = particle.Pdg();
            int mother = particle.FirstMother();
            // exclude e+ e- gamma with mother = -1 (coming from neutrino vertex, nuclear de-excitation)
            if((pdg == genie::kPdgGamma) | (pdg == genie::kPdgElectron) | (pdg == genie::kPdgPositron)){
                if(mother != -1) final_had_syst.push_back(particle);
            }else if(genie::pdg::IsPseudoParticle(pdg)){ // exclude possible pseudo particles
                continue;
            }else if(mother == 0){ // exclude primary lepton
                continue;
            }else{
                final_had_syst.push_back(particle);
            }
        }
    }
    return final_had_syst;                                                                                   
}

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::PrimaryStateHadronicSystem_status(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    /*
    Called inside: RDFUtils::GENIE::PrimaryStateHadronicSystem(particles, event_type)
    WARNING : 
        Found exception when calling this function with RES;[res9] events where the
        reuested state should be genie::kIStDecayedState and in the final state the 
        unstable status found is 13 genie::kIStPREDecayedState which is not handled 
        in the reference. The particle with state 13 is found having PDG 0 so clearly
        a class of events that has some problems and we wanna just not consider.
    */
    // std::cout<< "request partile with status : "<< STATUS << "\n";
    ROOT::VecOps::RVec<genie::GHepParticle> prim_had_syst;
    int particle_index = -1;
    int ist_store = -10;
    for(const auto& particle : particles){
        particle_index ++;
        int ist_comp  = particle.Status();
        if(ist_comp == STATUS){
            ist_store = particle_index;    //store this mother
            continue;
        }
        if(particle.FirstMother() == ist_store) prim_had_syst.push_back(particle);
    }
    // if(prim_had_syst.size() == 0){
    //     std::cout << "primary hadronic system with status: "<<STATUS<<" has no particles \n";
    //     // throw "";
    // }
    return prim_had_syst;
}   

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::PrimaryStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles,
                                                                                    std::string event_type){
    /*
    Extract info on the primary hadronic system (before any intranuclear rescattering)
    looking for particles with status_code == kIStHadronInTheNucleus
    An exception is the coherent production and scattering off free nucleon targets
    (no intranuclear rescattering) in which case primary hadronic system is set to be
    'identical' with the final  state hadronic system.
    For event_type refer to RDFUtils::GENIE::EventType
    Reference: https://github.com/GENIE-MC/Generator/blob/975179a5f2595790c5bb722a2515d18bda67af6f/src/Apps/gNtpConv.cxx#L775
    */
    
    ROOT::VecOps::RVec<genie::GHepParticle> prim_had_syst;
    
    genie::GHepParticle target = particles[1];

    // if coherent or free nucleon target set primary states equal to final states
    if(!genie::pdg::IsIon(target.Pdg()) || (event_type == "COH")){
        return RDFUtils::GENIE::FinalStateHadronicSystem(particles);
    }else{
        if(event_type == "RES"){
            prim_had_syst = RDFUtils::GENIE::PrimaryStateHadronicSystem_status<genie::kIStDecayedState>(particles);
        }else if(event_type == "DIS"){
            prim_had_syst = RDFUtils::GENIE::PrimaryStateHadronicSystem_status<genie::kIStDISPreFragmHadronicState>(particles);
        }else if(event_type == "QES"){
            prim_had_syst = RDFUtils::GENIE::PrimaryStateHadronicSystem_status<genie::kIStNucleonTarget>(particles);
        }else if(event_type == "MEC"){
            prim_had_syst = RDFUtils::GENIE::PrimaryStateHadronicSystem_status<genie::kIStDecayedState>(particles);
        }else{
        }
    }
    // also include gammas from nuclear de-excitations (appearing in the daughter list of the 
	// hit nucleus, earlier than the primary hadronic system extracted above)
    for(int i = target.FirstDaughter(); i <= target.LastDaughter(); i++) {
        if(i<0) continue;
        if(particles[i].Status()==genie::kIStStableFinalState) { 
            prim_had_syst.push_back(particles[i]); 
        }
	}
    // if(prim_had_syst.size() == 0){
    //     std::cout << "primary hadronic system has no particles \n";
    //     throw "";
    // }
    return prim_had_syst;
    }


ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetExhoticMesons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        int pdg = p.Pdg();
        if((pdg>100 || pdg<600)&&(pdg!=genie::kPdgPiP && pdg!=genie::kPdgPi0)) filtered_particles.push_back(p);
    }

    return filtered_particles;
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetExhoticHadrons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        int pdg = p.Pdg();
        if((pdg>1000 || pdg<50000)&&(pdg!=genie::kPdgProton && pdg!=genie::kPdgNeutron)) filtered_particles.push_back(p);
    }

    return filtered_particles;
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetRecoiledNuclei(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        int pdg = p.Pdg();
        if((pdg>=genie::kPdgTgtFreeP)) filtered_particles.push_back(p);
    }

    return filtered_particles;
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetFinalHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        int pdg = abs(p.Pdg());
        if((pdg!=genie::kPdgMuon || pdg!=genie::kPdgElectron)) filtered_particles.push_back(p);
    }

    return filtered_particles;
}

ROOT::VecOps::RVec<int> RDFUtils::GENIE::GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    ROOT::VecOps::RVec<int> pdgs;
    for(auto& p : particles) pdgs.push_back(p.Pdg());
    return pdgs;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::GENIE::GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    // return momentum (LorentzVector) of each of inout GHepParticle
    ROOT::VecOps::RVec<TLorentzVector> momenta;
    for(auto& p : particles){
        TLorentzVector momentum = *p.GetP4();
        momenta.push_back(momentum);
        }
    return momenta;
}

double RDFUtils::GENIE::GetInitialNucleonMomentum(const TLorentzVector& FinalStateMomentum,
                                                  const TVector3& FinalStateLongitudinalMomentum,
                                                  const TVector3& FinalStateDeltaPT,
                                                  int InteractionTarget){
/*total initial nucleon momentum as calculated in arXiv:2102.03346v1 (T2K collaboration)*/

    double target_nucleus_mass = GenieUtils::database->Find(InteractionTarget)->Mass();
    double proton_mass = GenieUtils::database->Find(genie::kPdgProton)->Mass();
    double proton_mean_excitation_energy = 26.1 * 1E-3; // 26.1 MeV for Carbon
    double residual_nucleus_mass = target_nucleus_mass - proton_mass + proton_mean_excitation_energy;
    double final_state_energy = FinalStateMomentum.T();
    double final_state_long_mom = FinalStateLongitudinalMomentum.Mag();
    double FinalStateDeltaPT2 = FinalStateDeltaPT.Mag() * FinalStateDeltaPT.Mag(); 
    double k = target_nucleus_mass + final_state_long_mom - final_state_energy;
    double nucleon_long_mom = 0.5*(k - (FinalStateDeltaPT2 + residual_nucleus_mass * residual_nucleus_mass) / k);

    return sqrt(FinalStateDeltaPT2 + nucleon_long_mom*nucleon_long_mom);
}

TLorentzVector RDFUtils::GENIE::SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL){
    // return the sum of N TLorentzVectors
    TLorentzVector sum = {0.,0.,0.,0.};
    for(auto& TL : VTL) sum += TL;
    return sum;
}

ROOT::RDF::RNode RDFUtils::AddConstantsToDF(ROOT::RDataFrame& df){
    return df.Define("SomeUsefullConstant","1.");
}

// bool RDFUtils::GENIE::IsUnphysical(genie::NtpMCEventRecord& m){
//     return m.event->IsUnphysical();
// }

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDF::RNode& df){
    return df
             // ghep columns
            //  .Define("IsUnphysical", RDFUtils::GENIE::IsUnphysical, {"gmcrec"})
            // exclude some unphysical events under name unknown
             .Filter([](TObjString& ev){return !(ev.GetString().Contains("Unknown"));}, {"EvtCode"})
             .Define("isCCEvent",                   RDFUtils::GENIE::IsCC, {"EvtCode"})
             .Define("Interaction_vtxX","EvtVtx[0]")
             .Define("Interaction_vtxY","EvtVtx[1]")
             .Define("Interaction_vtxZ","EvtVtx[2]")
             .Define("Interaction_vtxT","EvtVtx[3]")
             .Define("EventType",                   RDFUtils::GENIE::EventType, {"EvtCode"})
             .Define("NeutrinoFlavor",              RDFUtils::GENIE::NeutrinoFlavor, {"StdHepPdg"})
             .Define("InteractionTarget",           RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             .Define("InteractionTargetPDG",        "StdHepPdg[1]")
             .Define("InteractionTargetFromGEO",    RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
             .Define("InteractionVolume",           RDFUtils::GEO::GetVolumeFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
             // all particles: hadrons laptons and nuclei
             .Define("Particles",                   RDFUtils::GENIE::AllGenieParticles, {"StdHepPdg",
                                                                                         "StdHepStatus",
                                                                                         "StdHepP4", // momentum
                                                                                         "StdHepX4", // position
                                                                                         "StdHepFm", // first mother
                                                                                         "StdHepLm", // last mother
                                                                                         "StdHepFd", // first Daughter
                                                                                         "StdHepLd", // last Daughter
                                                                                         })
             /*
                initial state ___________________________________________________________________________________________________
                generator-level initial state
             */
             .Define("InitialStateParticles",       RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStInitialState>, {"Particles"})
             .Define("InitialStateParticlesPDG",    RDFUtils::GENIE::GetPDG, {"InitialStateParticles"})
             .Define("InitialStateParticlesNames",  RDFUtils::GENIE::PDG2Name, {"InitialStateParticles"})
             .Define("InitialStateParticlesP4",     RDFUtils::GENIE::GetMomentum, {"InitialStateParticles"})
             .Define("InitialStateP4",              RDFUtils::GENIE::SumLorentzVectors, {"InitialStateParticlesP4"})
             /* incoming neutrinos */
             .Define("IncomingNuMu",                RDFUtils::GENIE::GetParticlesWithPDG<14>, {"InitialStateParticles"})
             .Define("IncomingNuMu_P4",             RDFUtils::GENIE::GetMomentum, {"IncomingNuMu"})
             .Define("IncomingAntiNuMu",            RDFUtils::GENIE::GetParticlesWithPDG<-14>, {"InitialStateParticles"})
             .Define("IncomingAntiNuMu_P4",         RDFUtils::GENIE::GetMomentum, {"IncomingAntiNuMu"})
             /* nucleon target (free proton - i.e. interaction on H - should NOT be include in the category by GENIE) */
             .Define("NucleonTarget",               RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStNucleonTarget>, {"Particles"})
             .Define("NucleonTargetName",           RDFUtils::GENIE::PDG2Name, {"NucleonTarget"})
             .Define("NucleonTargetP4",             RDFUtils::GENIE::GetMomentum, {"NucleonTarget"})
             /*
                primary hadronic system _________________________________________________________________________________________
                hadrons generated in the nucleus (hadrons produced in the primary interaction of neutrinos) before any FSI
             */
             .Define("HadronsInTheNucleus",       RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStHadronInTheNucleus>, {"Particles"})
             .Define("HadronsInTheNucleusPDG",    RDFUtils::GENIE::GetPDG, {"HadronsInTheNucleus"})
             .Define("HadronsInTheNucleusNames",  RDFUtils::GENIE::PDG2Name, {"HadronsInTheNucleus"})
             .Define("HadronsInTheNucleusP4",     RDFUtils::GENIE::GetMomentum, {"HadronsInTheNucleus"})
             /*
                stable final state particles ____________________________________________________________________________________
                generator-level final state: particles to be tracked by detector-level MC
             */
             .Define("StableFinalStateParticles",   RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStStableFinalState>, {"Particles"})
             .Define("StableFinalStateParticlesPDG",RDFUtils::GENIE::GetPDG, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesName",RDFUtils::GENIE::PDG2Name, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesP4", RDFUtils::GENIE::GetMomentum, {"StableFinalStateParticles"})
             .Define("StableFinalStateP4",          RDFUtils::GENIE::SumLorentzVectors, {"StableFinalStateParticlesP4"})
             .Define("StableFinalStateTopology",    RDFUtils::GENIE::GetInteractionTopology, {"StableFinalStateParticles"})
             /*
                final state nuclear remnant _____________________________________________________________________________________
                low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle
             */
             .Define("FinalStateNuclearRemnant",    RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStFinalStateNuclearRemnant>, {"Particles"})
             .Define("FinalStateNuclearRemnantPDG", RDFUtils::GENIE::GetPDG, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearRemnantP",   RDFUtils::GENIE::GetMomentum, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearP4",         RDFUtils::GENIE::SumLorentzVectors, {"FinalStateNuclearRemnantP"})
             /*
                final state hadronic system
             */
             .Define("FinalStateHadronicSystem",    RDFUtils::GENIE::FinalStateHadronicSystem, {"Particles"})
             .Define("FinalStateHadronicSystemPDG", RDFUtils::GENIE::GetPDG, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemName",RDFUtils::GENIE::PDG2Name, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemP",   RDFUtils::GENIE::GetMomentum, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemTopology", RDFUtils::GENIE::GetInteractionTopology, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemTopology_name", [](GenieUtils::event_topology t){return t.Name();}, {"FinalStateHadronicSystemTopology"})
             /*
                primary state hadronic system
             */
             .Define("PrimaryStateHadronicSystem",  RDFUtils::GENIE::PrimaryStateHadronicSystem, {"Particles", "EventType"})
             .Define("PrimaryStateHadronicSystemPDG", RDFUtils::GENIE::GetPDG, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemName",RDFUtils::GENIE::PDG2Name, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemP",   RDFUtils::GENIE::GetMomentum, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemTopology", RDFUtils::GENIE::GetInteractionTopology, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemTopology_name", [](GenieUtils::event_topology t){return t.Name();}, {"PrimaryStateHadronicSystemTopology"})
             /*
                EVENT TOPOLOGY __________________________________________________________________________________________________
                organize final states
             */
             /*mu*/
             .Define("FinalStateMuons",             RDFUtils::GENIE::GetParticlesWithPDG<13>, {"StableFinalStateParticles"})
             .Define("FinalStateMuonsP",            RDFUtils::GENIE::GetMomentum, {"FinalStateMuons"})
             .Define("FinalStateMuonsP4",           RDFUtils::GENIE::SumLorentzVectors, {"FinalStateMuonsP"}) // just to obtain a TLVector instead of vector<TLVector>
             /*antimu*/
             .Define("FinalStateAntiMuons",         RDFUtils::GENIE::GetParticlesWithPDG<-13>, {"StableFinalStateParticles"})
             .Define("FinalStateAntiMuonsP",        RDFUtils::GENIE::GetMomentum, {"FinalStateAntiMuons"})
             .Define("FinalStateAntiMuonsP4",       RDFUtils::GENIE::SumLorentzVectors, {"FinalStateAntiMuonsP"}) // just to obtain a TLVector instead of vector<TLVector>
             /*protons*/
             .Define("FinalStateProtons",           RDFUtils::GENIE::GetParticlesWithPDG<2212>, {"StableFinalStateParticles"})
             .Define("FinalStateProtonsP",          RDFUtils::GENIE::GetMomentum, {"FinalStateProtons"})
             /*neutrons*/
             .Define("FinalStateNeutrons",          RDFUtils::GENIE::GetParticlesWithPDG<2212>, {"StableFinalStateParticles"})
             .Define("FinalStateNeutronsP",         RDFUtils::GENIE::GetMomentum, {"FinalStateNeutrons"})
             /*charged pions*/
             .Define("FinalStateChargedPi",         RDFUtils::GENIE::GetParticlesWithPDG<211>, {"StableFinalStateParticles"})
             .Define("FinalStateChargedPiP",        RDFUtils::GENIE::GetMomentum, {"FinalStateChargedPi"})
             /*neutral pions*/
             .Define("FinalStateNeutralPi",         RDFUtils::GENIE::GetParticlesWithPDG<111>, {"StableFinalStateParticles"})
             .Define("FinalStateNeutralPiP",        RDFUtils::GENIE::GetMomentum, {"FinalStateNeutralPi"})
             /*electron positons*/
             .Define("FinalStateElectronPositrons", RDFUtils::GENIE::GetParticlesWithPDG<11>, {"StableFinalStateParticles"})
             .Define("FinalStateElectronPositronsP",RDFUtils::GENIE::GetMomentum, {"StableFinalStateParticles"})
             /*photons*/
             .Define("FinalStatePhotons",           RDFUtils::GENIE::GetParticlesWithPDG<22>, {"StableFinalStateParticles"})
             .Define("FinalStatePhotonsP",          RDFUtils::GENIE::GetMomentum, {"FinalStatePhotons"})
             /*exhotic mesons : K0,KL,KS,K+-...*/ 
             .Define("FinalStateExhoticMesons",     RDFUtils::GENIE::GetExhoticMesons, {"StableFinalStateParticles"})
             .Define("FinalStateExhoticMesonsP",    RDFUtils::GENIE::GetMomentum, {"FinalStateExhoticMesons"})
             /*exhotic hadrons : sigma, lambda, ...*/ 
             .Define("FinalStateExhoticHadrons",    RDFUtils::GENIE::GetExhoticHadrons, {"StableFinalStateParticles"})
             .Define("FinalStateExhoticHadronsP",   RDFUtils::GENIE::GetMomentum, {"FinalStateExhoticHadrons"})
             /*recoiled nuclei*/ 
             .Define("FinalStateNuclei",            RDFUtils::GENIE::GetRecoiledNuclei, {"StableFinalStateParticles"})
             .Define("FinalStateNucleiP",           RDFUtils::GENIE::GetMomentum, {"FinalStateNuclei"})
             /* hadronic system : sum of stable final state expet muons, electromagnetic stuff, and nuclear remnant*/
             .Define("FinalHadronicSystem",         RDFUtils::GENIE::GetFinalHadronicSystem, {"StableFinalStateParticles"})
             .Define("FinalHadronicSystemP",        RDFUtils::GENIE::GetMomentum, {"FinalHadronicSystem"})
             .Define("FinalHadronicSystemP4",       RDFUtils::GENIE::SumLorentzVectors, {"FinalHadronicSystemP"}) 
             ;
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(ROOT::RDF::RNode& df){
    /* These variables allows the separation of interacions occurring on free proton (H)
    from those occuring on Carbon. Variables are explained in arXiv:2102.03346v1 and arXiv:1809.08752v2[Petti] */
    //____________________________________________________________________________
    return df
            // FinalHadronicSystemP4_TT : projection of FinalHadronicSystemP4 onto DoubleTransverseAxis perpendicular to neutrino direction
             .Define("DoubleTransverseAxis",       RDFUtils::TLVectorCrossProduct, {"InitialStateNuMu_P4","FinalStateMuonsP4"})
             .Define("FinalHadronicSystemP4_TT",   "DoubleTransverseAxis.Dot(FinalHadronicSystemP4.Vect())")
            // Initial Nucleon Momentum "P_N"
             .Define("BeamDirection",              "InitialStateNuMu_P4.Vect() * (1./InitialStateNuMu_P4.Vect().Mag())")
             .Define("BeamDirectionX", "BeamDirection.X()") // FOR CHECK
             .Define("BeamDirectionY", "BeamDirection.Y()") // FOR CHECK
             .Define("BeamDirectionZ", "BeamDirection.Z()") // FOR CHECK
             .Define("FinalHadronicSystemPL",      "BeamDirection * (FinalHadronicSystemP4.Vect().Dot(BeamDirection))") // component parallel to beam direction
             .Define("FinalStateMuonsPL",          "BeamDirection * (FinalStateMuonsP4.Vect().Dot(BeamDirection))") // component parallel to beam direction
             .Define("FinalHadronicSystemPT",      "FinalHadronicSystemP4.Vect() - FinalHadronicSystemPL")
             .Define("FinalStateMuonsPT",          "FinalStateMuonsP4.Vect() - FinalStateMuonsPL")
             .Define("FinalStatePL",               "FinalHadronicSystemPL + FinalStateMuonsPL")
             .Define("FinalStateDeltaPT",          "FinalStateMuonsPT + FinalHadronicSystemPT")
             .Define("FinalStateP4",               "FinalHadronicSystemP4 + FinalStateMuonsP4")
            //  .Define("InitialNucleonMomentum",     RDFUtils::GENIE::GetInitialNucleonMomentum, {"FinalStateP4", "FinalStatePL","FinalStateDeltaPT","InteractionTargetPDG"})
             /*angle of emission*/
             .Define("FinalStateMuonEmissionAngle",             "FinalStateMuonsP4.Vect().Angle(BeamDirection)")
             .Define("FinalStateHadronicSystemEmissionAngle",   "FinalHadronicSystemP4.Vect().Angle(BeamDirection)")
             .Define("FinalStateHadronicSystemVSMuonsAngle",    "FinalHadronicSystemP4.Vect().Angle(FinalStateMuonsP4.Vect())")
              /*Transverse boosting angle*/
             .Define("TransverseBoostingAngle",    "TMath::ACos((((-1.)*FinalStateMuonsPT).Dot(FinalStateDeltaPT))/FinalStateMuonsPT.Mag()/FinalStateDeltaPT.Mag())")
             // Transverse boosting angle "\deltaPt"
             .Define("MissingTransverseMomentum",  "FinalStateDeltaPT")
             .Define("Asimmetry_RmH",              "(MissingTransverseMomentum.Mag()-FinalHadronicSystemPT.Mag())/(MissingTransverseMomentum.Mag()+FinalHadronicSystemPT.Mag())")
             ;
}

//RDFUtils::GEO______________________________________________________________________

std::string RDFUtils::GEO::GetVolumeFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }
   
    return geo->FindNode(x*100., y*100., z*100.)->GetVolume()->GetName();
}

std::string RDFUtils::GEO::GetMaterialFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }

    auto navigator = geo->GetCurrentNavigator();

    if(!navigator) navigator = geo->AddNavigator();
   
    return geo->FindNode(x*100., y*100., z*100.)->GetVolume()->GetMaterial()->GetName();
}

//RDFUtils::EDEPSIM_______________________________________________________________

int RDFUtils::EDEPSIM::NofEventsPerSpill(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices){
    return vertices.size();
}

template<int coordinate>
ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::PrimariesVertex(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices)
{
    ROOT::VecOps::RVec<double> vertices_positions;
    for(auto& vertex: vertices)
    {
        auto vertex_position = vertex.GetPosition();
        vertices_positions.emplace_back(vertex_position[coordinate]);
    }
    return vertices_positions;
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::PrimariesVertexTimeDiff(ROOT::VecOps::RVec<double>& primaries_time)
{
    ROOT::VecOps::RVec<double> time_diff;

    for(auto i=1u; i<primaries_time.size(); i++)
    {
        time_diff.emplace_back(primaries_time[i]-primaries_time[i-1]);
    }
    return time_diff;
}

ROOT::RDF::RNode RDFUtils::EDEPSIM::AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df){
    return df
             .Define("NofEventsPerSpill", RDFUtils::EDEPSIM::NofEventsPerSpill, {"Primaries"})
             .Define("PrimariesVertexX", RDFUtils::EDEPSIM::PrimariesVertex<0>,{"Primaries"})
             .Define("PrimariesVertexY", RDFUtils::EDEPSIM::PrimariesVertex<1>,{"Primaries"})
             .Define("PrimariesVertexZ", RDFUtils::EDEPSIM::PrimariesVertex<2>,{"Primaries"})
             .Define("PrimariesVertexT", RDFUtils::EDEPSIM::PrimariesVertex<3>,{"Primaries"})
             .Define("PrimariesVertexTimeDiff",RDFUtils::EDEPSIM::PrimariesVertexTimeDiff,{"PrimariesVertexT"})
             ;
}