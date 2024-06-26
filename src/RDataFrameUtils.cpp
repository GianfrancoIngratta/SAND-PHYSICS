#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "RecoUtils.h"
#include "GeoUtils.h"

#include "TChain.h"
#include "TMath.h"

TGeoManager* geo = nullptr;

//RDFUtils______________________________________________________________________

ROOT::RDataFrame RDFUtils::InitDF(TString production, 
                                  const char* tree_name, 
                                  unsigned int file_index_start, 
                                  unsigned int file_index_stop){
    
    TChain* chain = new TChain(tree_name, tree_name);
    if(production.Contains("*")){
        // chain a given number of file of a production
        for(unsigned int i = file_index_start; i < file_index_stop; i++){
            auto file = TString::Format(production.ReplaceAll("*","%d"), i, i);
            if (!gSystem->AccessPathName(file.Data())) chain->Add(file);
        }
    }else{
        chain->Add(production.Data());
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
        // std::cout << "colName : " <<colName<<", colType : " <<df.GetColumnType(colName)<<"\n";
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
    }else{
        return "Other";
    }
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                                           const ROOT::VecOps::RVec<int>& status,
                                                                           const ROOT::VecOps::RVec<double>& P4){
    /*
        Construct an istance genie::GHepParticle for each particle in the genie output.
        The class genie::GHepParticle is more helpfull for the analysis
    */                                                                            
    ROOT::VecOps::RVec<genie::GHepParticle> particles;

    for(auto i=0u; i<pdg.size(); i++){
        
        TLorentzVector momentum = {P4[4*i],P4[4*i+1],P4[4*i+2],P4[4*i+3]};

        auto flag = static_cast<genie::GHepStatus_t>(status[i]);
        
        auto particle = GenieUtils::GenieParticle(pdg[i], flag, momentum);
        
        particles.push_back(particle);
    }

    return particles;                                                                            
}

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> selected;

    for(auto& p : particles){
        if(p.Status() == STATUS) selected.push_back(p);
    }

    return selected;
}

GenieUtils::event_topology RDFUtils::GENIE::GetFinalStateTopology(const ROOT::VecOps::RVec<int>& pdgs)
{
    
    GenieUtils::event_topology t;

    for(auto& p : pdgs)
    {
        auto pdg = abs(p);
        if(pdg==genie::kPdgMuon){
                t.NofMuons++;
        }else if(pdg==genie::kPdgElectron || pdg==genie::kPdgGamma){
                t.NofElPosGamma++;
        }else if(pdg==genie::kPdgProton){
                t.NofProtons++;
        }else if(pdg==genie::kPdgNeutron){
                t.NofNeutrons++;
        }else if(pdg==genie::kPdgPiP || pdg==genie::kPdgPi0){
                t.NofPions++;
        }else if((pdg>100 || pdg<600)&&(pdg!=genie::kPdgPiP && pdg!=genie::kPdgPi0)){
                t.NofExhotic++;
        }else if((pdg>1000 || pdg<50000)&&(pdg!=genie::kPdgProton && pdg!=genie::kPdgNeutron)){
                t.NofExhotic++;
        }else if((pdg>=genie::kPdgTgtFreeP)){
                t.NofRecoiledNuclei++;
        }else{
                continue;
        }
    }

    return t;                  
}

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetParticlesWithPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    // return GHepParticle whose pdg is PDG
    
    ROOT::VecOps::RVec<genie::GHepParticle> filtered_particles;

    for(auto& p : particles){
        if(abs(p.Pdg())==PDG) filtered_particles.push_back(p);
    }
    
    return filtered_particles;
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

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDF::RNode& df){
    return df.Define("Interaction_vtxX","EvtVtx[0]")
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
             .Define("Particles",                   RDFUtils::GENIE::AllGenieParticles, {"StdHepPdg","StdHepStatus","StdHepP4"})
             /*
                initial state ___________________________________________________________________________________________________
                generator-level initial state
             */
             .Define("InitialStateParticles",       RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStInitialState>, {"Particles"})
             .Define("InitialStateParticlesPDG",    RDFUtils::GENIE::GetPDG, {"InitialStateParticles"})
             .Define("InitialStateParticlesP",      RDFUtils::GENIE::GetMomentum, {"InitialStateParticles"})
             .Define("InitialStateParticlesPx",     RDFUtils::GetComponent<0>, {"InitialStateParticlesP"})
             .Define("InitialStateParticlesPy",     RDFUtils::GetComponent<1>, {"InitialStateParticlesP"})
             .Define("InitialStateParticlesPz",     RDFUtils::GetComponent<2>, {"InitialStateParticlesP"})
             .Define("InitialStateParticlesE",      RDFUtils::GetComponent<3>, {"InitialStateParticlesP"})
             .Define("InitialStateP4",              RDFUtils::GENIE::SumLorentzVectors, {"InitialStateParticlesP"})
             .Define("InitialStateMomentum", "InitialStateP4.Vect().Mag()")
             .Define("InitialStateEnergy",          RDFUtils::GetColumnSum, {"InitialStateParticlesE"})
             /* incoming neutrinos */
             .Define("InitialStateNuMu",            RDFUtils::GENIE::GetParticlesWithPDG<14>, {"InitialStateParticles"})
             .Define("InitialStateNuMu_P",          RDFUtils::GENIE::GetMomentum, {"InitialStateNuMu"})
             .Define("InitialStateNuMu_P4",         RDFUtils::GENIE::SumLorentzVectors, {"InitialStateNuMu_P"}) // just to obtain a TLVector instead of vector<TLVector>
             .Define("InitialStateAntiNuMu",        RDFUtils::GENIE::GetParticlesWithPDG<-14>, {"InitialStateParticles"})
             .Define("InitialStateAntiNuMu_P",      RDFUtils::GENIE::GetMomentum, {"InitialStateAntiNuMu"})
             .Define("InitialStateAntiNuMu_P4",     RDFUtils::GENIE::SumLorentzVectors, {"InitialStateAntiNuMu_P"})
             /*
                final state state particles _____________________________________________________________________________________
                generator-level final state: particles to be tracked by detector-level MC
             */
             .Define("StableFinalStateParticles",   RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStStableFinalState>, {"Particles"})
             .Define("StableFinalStateParticlesPDG",RDFUtils::GENIE::GetPDG, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesP",  RDFUtils::GENIE::GetMomentum, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesPx", RDFUtils::GetComponent<0>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesPy", RDFUtils::GetComponent<1>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesPz", RDFUtils::GetComponent<2>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesE",  RDFUtils::GetComponent<3>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateP4",          RDFUtils::GENIE::SumLorentzVectors, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateMomentum",    "StableFinalStateP4.Vect().Mag()")
             .Define("StableFinalStateEnergy",      RDFUtils::GetColumnSum, {"StableFinalStateParticlesE"})
             /*
                final state nuclear remnant _____________________________________________________________________________________
                low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle
             */
             .Define("FinalStateNuclearRemnant",    RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStFinalStateNuclearRemnant>, {"Particles"})
             .Define("FinalStateNuclearRemnantPDG", RDFUtils::GENIE::GetPDG, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearRemnantP",   RDFUtils::GENIE::GetMomentum, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearRemnantPx",  RDFUtils::GetComponent<0>, {"FinalStateNuclearRemnantP"})
             .Define("FinalStateNuclearRemnantPy",  RDFUtils::GetComponent<1>, {"FinalStateNuclearRemnantP"})
             .Define("FinalStateNuclearRemnantPz",  RDFUtils::GetComponent<2>, {"FinalStateNuclearRemnantP"})
             .Define("FinalStateNuclearRemnantE",   RDFUtils::GetComponent<3>, {"FinalStateNuclearRemnantP"})
             .Define("FinalStateNuclearP4",         RDFUtils::GENIE::SumLorentzVectors, {"FinalStateNuclearRemnantP"})
             .Define("FinalStateNuclearMomentum", "FinalStateNuclearP4.Vect().Mag()")
             .Define("FinalStateNuclearEnergy",     RDFUtils::GetColumnSum, {"FinalStateNuclearRemnantE"})
             /*
                EVENT TOPOLOGY __________________________________________________________________________________________________
                organize final states
             */
             .Define("FinalStateTopology",          RDFUtils::GENIE::GetFinalStateTopology, {"StableFinalStateParticlesPDG"})
             .Define("FinalStateTopologyName",      [](GenieUtils::event_topology t){return t.GetTopologyName();}, {"FinalStateTopology"})
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

//RDFUtils::SANDRECO_________________________________________________________________

ROOT::VecOps::RVec<particle> RDFUtils::SANDRECO::ParticleGoodTrackFit(const ROOT::VecOps::RVec<particle>& particles){
    
    ROOT::VecOps::RVec<particle> good_particles;

    for(auto& p : particles){
        if(p.tr.chi2_ln==0 && p.tr.chi2_cr==0) good_particles.push_back(p);	
    }

    return good_particles;
}

ROOT::VecOps::RVec<particle> RDFUtils::SANDRECO::FilterPrimaries(const ROOT::VecOps::RVec<particle>& particles){
    
    ROOT::VecOps::RVec<particle> primaries;

    for(auto& p : particles){
        if(p.primary==1) primaries.push_back(p);	
    }

    return primaries;
}

template<int PDG>
ROOT::VecOps::RVec<particle> RDFUtils::SANDRECO::GetParticlesWithPDG(const ROOT::VecOps::RVec<particle>& particles){

    ROOT::VecOps::RVec<particle> filtered_particles;

    for(auto& p : particles){
        if(p.pdg==PDG) filtered_particles.push_back(p);
    }

    return filtered_particles;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::SANDRECO::GetMomentum(const ROOT::VecOps::RVec<particle>& particles){

    ROOT::VecOps::RVec<TLorentzVector> momenta;

    for(auto& p : particles){
        TLorentzVector momentum = {p.pxreco, p.pyreco, p.pzreco, p.Ereco};
        momenta.push_back(momentum);
    }

    return momenta;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::SANDRECO::GetTrackVertex(const ROOT::VecOps::RVec<particle>& particles){

    ROOT::VecOps::RVec<TLorentzVector> vertices;

    for(auto& p : particles){
        TLorentzVector vertex = {p.xreco, p.yreco, p.zreco, p.treco};
        vertices.push_back(vertex);
    }

    return vertices;
}

ROOT::RDF::RNode RDFUtils::SANDRECO::AddColumnsFromSANDRECO(ROOT::RDF::RNode& df){
    return df.Define("ParticlesGoodTrackFit", RDFUtils::SANDRECO::ParticleGoodTrackFit, {"particles"})
             .Define("PrimariesGoodTrackFit", RDFUtils::SANDRECO::FilterPrimaries, {"ParticlesGoodTrackFit"})
             // particles
             /*muons with good trackfit*/
             .Define("Muons",                 RDFUtils::SANDRECO::GetParticlesWithPDG<13>, {"PrimariesGoodTrackFit"})
             .Define("MuonsVtx",              RDFUtils::SANDRECO::GetTrackVertex, {"PrimariesGoodTrackFit"})
             .Define("MuonsVtxX",             RDFUtils::GetComponent<0>, {"MuonsVtx"})
             .Define("MuonsVtxY",             RDFUtils::GetComponent<1>, {"MuonsVtx"})
             .Define("MuonsVtxZ",             RDFUtils::GetComponent<2>, {"MuonsVtx"})
             .Define("MuonsP4",               RDFUtils::SANDRECO::GetMomentum, {"Muons"})
             .Define("MuonsP",                "MuonsP4[0].Vect().Mag()")
             ;
}

//RDFUtils::FASTRECO_________________________________________________________________

template<int PDG>
ROOT::VecOps::RVec<fast::particle> RDFUtils::FASTRECO::GetParticlesWithPDG(const std::map<int, fast::particle>& particle_map){

    ROOT::VecOps::RVec<fast::particle> particles_filter;

    for(auto& entry : particle_map){
        if(entry.second.parentid==-1 && entry.second.pdg==PDG) particles_filter.push_back(entry.second);
    }

    // particles_filter can contain either 1 or 0 entries
    return particles_filter;
}

ROOT::VecOps::RVec<TVector3> RDFUtils::FASTRECO::GetTrackVertex(const ROOT::VecOps::RVec<fast::particle>& particles){

    ROOT::VecOps::RVec<TVector3> vertices;

    for(auto& p : particles){
        TVector3 vertex = {p.position[0],p.position[1],p.position[2]};
        vertices.push_back(vertex);
    }

    return vertices;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::FASTRECO::GetReco4Momentum(const ROOT::VecOps::RVec<fast::particle>& particles){
    
    ROOT::VecOps::RVec<TLorentzVector> momenta;

    auto const MeV_to_GeV = 1E-3;

    if(particles.size()==0){
        momenta.push_back(RECO::DEFAULT_TLV);
        return momenta;
    }

    for (auto& p : particles)
    {
        TLorentzVector momentum = {p.recoP4[0] * MeV_to_GeV,
                                   p.recoP4[1] * MeV_to_GeV,
                                   p.recoP4[2] * MeV_to_GeV,
                                   p.recoP4[3] * MeV_to_GeV};
        
        momenta.push_back(momentum);                           
    }

    return momenta;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::FASTRECO::GetTrue4Momentum(const ROOT::VecOps::RVec<fast::particle>& particles){
    
    ROOT::VecOps::RVec<TLorentzVector> momenta;

    auto const MeV_to_GeV = 1E-3;

    if(particles.size()==0){
        momenta.push_back(RECO::DEFAULT_TLV);
        return momenta;
    }

    for (auto& p : particles)
    {
        TLorentzVector momentum = {p.trueP4[0] * MeV_to_GeV,
                                   p.trueP4[1] * MeV_to_GeV,
                                   p.trueP4[2] * MeV_to_GeV,
                                   p.trueP4[3] * MeV_to_GeV};
        
        momenta.push_back(momentum);                           
    }

    return momenta;
}

template<int coord>
double RDFUtils::FASTRECO::GetEventVertex(const ROOT::VecOps::RVec<TVector3>& vertices){
    if(vertices.size()==RECO::DEFAULT_DOUBLE){
        return 0.;
    }else{
        return vertices[0][coord];
    }
}

ROOT::RDF::RNode RDFUtils::FASTRECO::AddColumnsFromFASTRECO(ROOT::RDF::RNode& df){
    return df
            /*muons*/
             .Define("Muons",                 RDFUtils::FASTRECO::GetParticlesWithPDG<13>, {"particles"})
             .Define("MuonsTrueP4",           RDFUtils::FASTRECO::GetTrue4Momentum, {"Muons"})
             .Define("MuonsRecoP4",           RDFUtils::FASTRECO::GetReco4Momentum, {"Muons"})
            /*protons*/
             .Define("Protons",               RDFUtils::FASTRECO::GetParticlesWithPDG<2212>, {"particles"})
            /*neutrons*/
             .Define("Neutrons",              RDFUtils::FASTRECO::GetParticlesWithPDG<2112>, {"particles"})
             //
             .Define("MuonsVtx",              RDFUtils::FASTRECO::GetTrackVertex, {"Muons"})
             .Define("MuonsVtxX",             RDFUtils::FASTRECO::GetEventVertex<0>, {"MuonsVtx"})
             .Define("MuonsVtxY",             RDFUtils::FASTRECO::GetEventVertex<1>, {"MuonsVtx"})
             .Define("MuonsVtxZ",             RDFUtils::FASTRECO::GetEventVertex<2>, {"MuonsVtx"})
             
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