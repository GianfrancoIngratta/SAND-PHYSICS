#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "GeoUtils.h"
#include "TChain.h"

TGeoManager* geo = nullptr;

//RDFUtils______________________________________________________________________

ROOT::RDataFrame RDFUtils::InitDF(const char* production, const char* tree_name){
    
    TChain* chain = new TChain(tree_name, tree_name);

    chain->Add(production);

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

//RDFUtils::GENIE_________________________________________________________________

std::string RDFUtils::GENIE::InteractionTarget(const ROOT::VecOps::RVec<int>& pdg){
    return GenieUtils::PDG2Name(pdg[1]);
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

ROOT::VecOps::RVec<int> RDFUtils::GENIE::GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    ROOT::VecOps::RVec<int> pdgs;
    for(auto& p : particles) pdgs.push_back(p.Pdg());
    return pdgs;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::GENIE::GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    // input : particles
    // output : particles momentum TLorentzVector
    ROOT::VecOps::RVec<TLorentzVector> momenta;
    for(auto& p : particles){
        TLorentzVector momentum = *p.GetP4();
        momenta.push_back(momentum);
        }
    return momenta;
}

TLorentzVector RDFUtils::GENIE::SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL){
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
             .Define("InteractionTarget",           RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             .Define("InteractionTargetFromGEO",    RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
             // all particles (hadrons laptons) and nuclei
             .Define("Particles",                   RDFUtils::GENIE::AllGenieParticles, {"StdHepPdg","StdHepStatus","StdHepP4"})
             //initial state
             /* generator-level initial state */
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
             //final state state particles
             /* generator-level final state: particles to be tracked by detector-level MC */
             .Define("StableFinalStateParticles",   RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStStableFinalState>, {"Particles"})
             .Define("StableFinalStateParticlesPDG",RDFUtils::GENIE::GetPDG, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesP",  RDFUtils::GENIE::GetMomentum, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesPx", RDFUtils::GetComponent<0>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesPy", RDFUtils::GetComponent<1>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesPz", RDFUtils::GetComponent<2>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateParticlesE",  RDFUtils::GetComponent<3>, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateP4",          RDFUtils::GENIE::SumLorentzVectors, {"StableFinalStateParticlesP"})
             .Define("StableFinalStateMomentum", "StableFinalStateP4.Vect().Mag()")
             .Define("StableFinalStateEnergy",      RDFUtils::GetColumnSum, {"StableFinalStateParticlesE"})
             // final state nuclear remnant
             /* low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle */
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
             // EVENT TOPOLOGY
             /* organize final states */
             .Define("FinalStateTopology",          RDFUtils::GENIE::GetFinalStateTopology, {"StableFinalStateParticlesPDG"})
             .Define("FinalStateTopologyName",      [](GenieUtils::event_topology t){return t.GetTopologyName();}, {"FinalStateTopology"})
             /*mu*/
             .Define("FinalStateMuons",             RDFUtils::GENIE::GetParticlesWithPDG<13>, {"StableFinalStateParticles"})
             .Define("FinalStateMuonsP",            RDFUtils::GENIE::GetMomentum, {"FinalStateMuons"})
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
             .Define("FinalStatePhotons",          RDFUtils::GENIE::GetParticlesWithPDG<22>, {"StableFinalStateParticles"})
             .Define("FinalStatePhotonsP",         RDFUtils::GENIE::GetMomentum, {"FinalStatePhotons"})
             /*exhotic mesons : K0,KL,KS,K+-...*/ 
             .Define("FinalStateExhoticMesons",    RDFUtils::GENIE::GetExhoticMesons, {"StableFinalStateParticles"})
             .Define("FinalStateExhoticMesonsP",   RDFUtils::GENIE::GetMomentum, {"FinalStateExhoticMesons"})
             /*exhotic hadrons : sigma, lambda, ...*/ 
             .Define("FinalStateExhoticHadrons",   RDFUtils::GENIE::GetExhoticHadrons, {"StableFinalStateParticles"})
             .Define("FinalStateExhoticHadronsP",  RDFUtils::GENIE::GetMomentum, {"FinalStateExhoticHadrons"})
             /*recoiled nuclei */ 
             .Define("FinalStateNuclei",           RDFUtils::GENIE::GetRecoiledNuclei, {"StableFinalStateParticles"})
             .Define("FinalStateNucleiP",          RDFUtils::GENIE::GetMomentum, {"FinalStateNuclei"})
             ;
}

//RDFUtils::GEO___________________________________________________________________

std::string RDFUtils::GEO::GetMaterialFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }
   
    return geo->FindNode(x*100., y*100., z*100.)->GetVolume()->GetMaterial()->GetName();
}