#include "RDataFrameUtils.h"

#include "TChain.h"
#include "TMath.h"

TGeoManager* geo = nullptr;

//RDFUtils______________________________________________________________________

TChain* RDFUtils::InitTChain(TString production,
                            const char* tree_name,
                            unsigned int file_index_start,
                            unsigned int file_index_stop){
    
    TChain* chain = new TChain(tree_name, tree_name);
    if(!production.Contains("*")){
        LOG("W", "Pass multiple files of a production using * ");
        throw "";
    }
    for(unsigned int i = file_index_start; i < file_index_stop; i++){
        auto file = TString::Format(production.ReplaceAll("*","%d"), i, i);
        if (!gSystem->AccessPathName(file.Data())){
            chain->Add(file);
        }else{
            LOG("W", TString::Format("Skipping file %s", file.Data()).Data());
        }
    }
    std::cout << "tree name : " << tree_name << " number of events in chain : " << chain ->GetEntries() << "\n";
    return chain;                        
}

ROOT::RDataFrame RDFUtils::InitDF(TChain* input_chain){
    ROOT::RDataFrame* df = new ROOT::RDataFrame(*input_chain);
    if(!df){
        std::cout<<"DF NOT initialized in file "<<__FILE__<<" line "<<__LINE__<<"\n";
        throw "";
    }
    return *df;
}

void RDFUtils::PrintColumns(ROOT::RDataFrame& df){

    auto colNames = df.GetColumnNames();
    
    for (auto &&colName : colNames){
        std::cout << "colName : " << colName << ", type :" << df.GetColumnType(colName) << "\n";
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

double RDFUtils::SumDuble(const ROOT::VecOps::RVec<double>& v){
    double sum = 0;
    for (const auto& elem : v) {
        sum += elem;
    }
    return sum;
}


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

ROOT::VecOps::RVec<double> RDFUtils::GENIE::GetKinE(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    ROOT::VecOps::RVec<double> KinE = {};
    for(const auto& particle : particles){ 
        KinE.push_back(particle.KinE());
    }
    return KinE;
}

TLorentzVector RDFUtils::GENIE::SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL){
    // return the sum of N TLorentzVectors
    TLorentzVector sum = {0.,0.,0.,0.};
    for(auto& TL : VTL) sum += TL;
    return sum;
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

int RDFUtils::GENIE::CountChargedParticles(const genie::GHepParticle& fs_lepton, 
                                           const ROOT::VecOps::RVec<genie::GHepParticle>& fs_HadronicSystem){
    // count nof final state charged particles
    int count = 0;
    int lepton_charge = GenieUtils::database->Find(fs_lepton.Pdg())->Charge();
    if(lepton_charge != 0) count++;
    for(const auto& hadron : fs_HadronicSystem){
            auto hadron_charge = GenieUtils::database->Find(hadron.Pdg())->Charge();
            if(hadron_charge != 0) count++;
    }
    return count;
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

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetInitialState(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    /*
        Intial State Particles relevant for the kineamtic analysis : neutrino + target nucleon
    */
   ROOT::VecOps::RVec<genie::GHepParticle> initial_state;
   for(const auto& particle : particles){
        if(particle.Status() == genie::kIStInitialState){
            if((particle.Pdg() == 14)| // numu
              (particle.Pdg() == -14)| // antinumu
              (particle.Pdg() == 12)|  // numu
              (particle.Pdg() == -12)| // antinumu
              (particle.Pdg()==2212)   // interacion on hydrogen (free proton)
              ){initial_state.push_back(particle);}
        }else if(particle.Status() == genie::kIStNucleonTarget){ // p or n in the target nucleus
            initial_state.push_back(particle);
        }
        // else if((particle.Pdg()==20000000200)|
        //          (particle.Pdg()==20000000201)|
        //          (particle.Pdg()==20000000202)){ // genie pseudo particles (n+n) (n+p) (p+p)
        //     initial_state.push_back(particle);
        // }
   }
   return initial_state;
}

genie::GHepParticle RDFUtils::GENIE::GetFinalStateLepton(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    genie::GHepParticle fs_lepton;
    for(auto& particle : particles){
        if((particle.Status() == genie::kIStStableFinalState) & (particle.FirstMother() == 0)){
            fs_lepton = particle;
            break;
        }
    }
    return fs_lepton;
}

genie::GHepParticle RDFUtils::GENIE::GetIncomingNeutrino(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    genie::GHepParticle nu;
    for(auto particle : particles){
        if(particle.Status() == genie::kIStInitialState){
            if((particle.Pdg() == 14)|(particle.Pdg() == -14)|(particle.Pdg() == 12)|(particle.Pdg() == -12)){
                nu =  particle;
                break;
            }
        }
    }
    return nu;
}

genie::GHepParticle RDFUtils::GENIE::GetNucleonTarget(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    genie::GHepParticle t;
    for(auto particle : particles){
        if(particle.Status() == genie::kIStNucleonTarget){
            return particle;
        }else if((particle.Status() == genie::kIStInitialState)&(particle.Pdg()==2212)){
            return particle;
        }else{
        }
    }
    return t;
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
    // int ist_store = -10; modified
    std::vector<int> ist_store = {};
    for(const auto& particle : particles){
        particle_index ++;
        int ist_comp  = particle.Status();
        if(ist_comp == STATUS){
            // ist_store = particle_index;    //store this mother modified
            ist_store.push_back(particle_index);    //store this mother
            continue;
        }
        // if(particle.FirstMother() == ist_store) prim_had_syst.push_back(particle);
        auto it = std::find(ist_store.begin(), ist_store.end(), particle.FirstMother());
        if(it != ist_store.end()) prim_had_syst.push_back(particle);
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

TVector3 RDFUtils::GENIE::GetTransverseComponent(const TLorentzVector& v1, const TLorentzVector v2){
    /*
        Get v2 component on the plane perpendiculare to v1
    */
   TVector3 v1_direction = v1.Vect().Unit();
   double v2_component_on_v1 = v2.Vect().Dot(v1_direction);
   TVector3 v2_parallel_to_v1 = v2_component_on_v1 * v1_direction;
   return v2.Vect() - v2_parallel_to_v1;
}

ROOT::VecOps::RVec<double> RDFUtils::GENIE::DotProductWithAxis(const TVector3& axis, const ROOT::VecOps::RVec<TLorentzVector>& V){
    ROOT::VecOps::RVec<double> products = {};
    for(const auto& v : V ){
        products.push_back(axis.Dot(v.Vect()));
    }
    return products;
}

TLorentzVector RDFUtils::GENIE::GetNup4FromMu(double proton_mass, double neutron_mass, double muon_mass, 
                                              const TLorentzVector& lepton, const TVector3 nu_versor){
    TVector3 leptonP3 = lepton.Vect();
    double mu_emission_angle = nu_versor.Angle(leptonP3);
    double neutron_mass2 = neutron_mass * neutron_mass;
    double proton_mass2 = proton_mass * proton_mass;
    double muon_mass2 = muon_mass * muon_mass;
    
    double numerator = neutron_mass2 - proton_mass2 - muon_mass2 + 2 * lepton.T() * proton_mass;
    double denominator = 2 * (proton_mass - lepton.T() + leptonP3.Mag() * cos(mu_emission_angle));
    double calculated_nu_E = numerator / denominator;
    
    TVector3 nu_P3 = calculated_nu_E * nu_versor;
    TLorentzVector nu_P4 = {nu_P3.X(), nu_P3.Y(), nu_P3.Z(), calculated_nu_E};
    return nu_P4;
}

double IntersectWithInfCylinder(const TVector3& vtx, const TVector3& neutron_direction) {
    /*
        Calculate the intersection btw the neutron that starts from vtx and travels
        in the direction neutron_direction (versor) and SAND.
        Take the intersection along the direction of travel, discard the other one.
    */
    TVector3 start_point(vtx.X() - GeoUtils::DRIFT::SAND_CENTER_X, 
                         vtx.Y() - GeoUtils::DRIFT::SAND_CENTER_Y, 
                         vtx.Z() - GeoUtils::DRIFT::SAND_CENTER_Z);
    double A = neutron_direction.Z() * neutron_direction.Z() + neutron_direction.Y() * neutron_direction.Y();
    double B = 2 * (start_point.Z() * neutron_direction.Z() + start_point.Y() * neutron_direction.Y());
    double C = start_point.Z() * start_point.Z() + start_point.Y() * start_point.Y() - 
               (GeoUtils::SAND_INNER_VOL_DIAMETER / 2.0) * (GeoUtils::SAND_INNER_VOL_DIAMETER / 2.0);
    double discriminant = B * B - 4 * A * C;

    if (discriminant < 0) {
        std::cerr << "No real intersection.\n";
        return 0.; 
    }

    double t1 = (-B + sqrt(discriminant)) / (2 * A);
    double t2 = (-B - sqrt(discriminant)) / (2 * A);

    if (t1 >= 0 && t2 >= 0) {
        return std::min(t1, t2); // Prendi la più piccola tra le soluzioni positive
    } else if (t1 >= 0) {
        return t1;
    } else if (t2 >= 0) {
        return t2;
    } else {
        std::cerr << "Both intersections are negative, no valid intersection.\n";
        return std::numeric_limits<double>::infinity();
    }

}

TVector3 RDFUtils::GENIE::NeutronArrivalPosECAL(std::string units, 
                                                double vtx_x, 
                                                double vtx_y, 
                                                double vtx_z, 
                                                const TVector3& neutron_momentum) {
    /*
        Given the neutrino vertex position and the neutron momentum
        predict where neutron is expected to release energy.
        !!! THE RETURNED POINT'S UNITS IS mm 
    */
    TVector3 start_point(vtx_x, vtx_y, vtx_z);
    // convert to mm
    if (units == "m") {
        start_point *= 1e3;
    } else if (units == "cm") {
        start_point *= 10.;
    } else if (units == "mm") {
        start_point *= 1.;
    } else {
        std::cerr << "Error: Unsupported unit \"" << units << "\". Supported units are m, cm, mm.\n";
        throw std::invalid_argument("Unsupported unit");
    }

    TVector3 neutron_direction = neutron_momentum.Unit();
    double t = IntersectWithInfCylinder(start_point, neutron_direction);

    if (fabs(start_point.X() + t * neutron_direction.X() - GeoUtils::DRIFT::SAND_CENTER_X) 
            >= GeoUtils::SAND_INNER_VOL_X_LENGTH / 2.0) {// neutron intersect ECAL in one of the 2 endcaps
        
        double t1 = (GeoUtils::SAND_CENTER_X - start_point.X()
                   + GeoUtils::SAND_INNER_VOL_X_LENGTH / 2.0) / neutron_direction.X();
        
        double t2 = (GeoUtils::SAND_CENTER_X - start_point.X()
                   - GeoUtils::SAND_INNER_VOL_X_LENGTH / 2.0 ) / neutron_direction.X();
        if(t1<0){
            t = t2;
        }else{
            t = t1;
        }
    }
    /*
        assume neutron release energy somewhere around half thickness ECAL module
    */
    t += GeoUtils::ECAL::ECAL_thickness/2.;
    return start_point + t * neutron_direction;
}

double RDFUtils::GENIE::GetNeutronTOF(std::string units_l,
                                      std::string units_E,
                                      const double vtx,
                                      const double vty,
                                      const double vtz,
                                      TVector3 X3,
                                      TVector3 P3
                                      ){
    /*
        Get Expected neutron tof in ECAL given the vertex poi and it energy
    */
    TVector3 vertex = {vtx, vty, vtz};
    if(units_l == "m"){
        vertex *= 1e3;
    }else if(units_l == "cm"){
        vertex *= 10.;
    }else{}
    
    if(units_E == "GeV") P3 *= 1e3;
    
    double c = 299.792; // [mm/ns]
    double dist = (X3 - vertex).Mag(); // [mm]
    double Ptot = P3.Mag(); // [MeV]
    double neutron_mass = GenieUtils::database->Find(2112)->Mass() * 1e3; // [MeV] 
    double E = sqrt(Ptot*Ptot + neutron_mass*neutron_mass); // [Mev]
    double gamma = E / neutron_mass;
    double beta = sqrt(1 - 1/(gamma*gamma));

    return dist / (beta * c); // [ns]
}

ROOT::RDF::RNode RDFUtils::AddConstantsToDF(ROOT::RDataFrame& df){
    return df
            // units
            .Define("GENIE_UNIT_LEGTH",   [](){std::string u = "m";return u;})
            .Define("GENIE_UNIT_ENERGY",  [](){std::string u = "GeV";return u;})
            .Define("EDEP_UNIT_LEGTH",    [](){std::string m = "mm";return m;})
            .Define("EDEP_UNIT_ENERGY",   [](){std::string m = "MeV";return m;})
            .Define("TGEO_UNIT_LEGTH",    [](){std::string m = "cm";return m;})
            //
            .Define("PROTON_MASS_GeV",    [](){return GenieUtils::GetMass(2212);})
            .Define("NEUTRON_MASS_GeV",   [](){return GenieUtils::GetMass(2112);})
            .Define("MUON_MASS_GeV",      [](){return GenieUtils::GetMass(13);})
            ;
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDF::RNode& df){
    return df
            // exclude some unphysical events under name unknown
             .Filter([](TObjString& ev){return !(ev.GetString().Contains("Unknown"));}, {"EvtCode"})
             .Define("NeutrinoFlavor",              "StdHepPdg[0]")
             .Define("isCCEvent",                   RDFUtils::GENIE::IsCC, {"EvtCode"})
             .Define("EventType",                   RDFUtils::GENIE::EventType, {"EvtCode"})
             .Define("InteractionTargetPDG",        "StdHepPdg[1]")
             // !!! CHANNEL OF INTEREST : CCQE events on Hydrogen
             .Define("CCQEonHydrogen", [](std::string event_type, int target_pdg){return (event_type=="QES")*(target_pdg==2212);}, {"EventType", "InteractionTargetPDG"})
             .Define("Interaction_vtxX","EvtVtx[0]")
             .Define("Interaction_vtxY","EvtVtx[1]")
             .Define("Interaction_vtxZ","EvtVtx[2]")
             .Define("Interaction_vtxT","EvtVtx[3]")
             .Define("InteractionVolume",           RDFUtils::GEO::GetVolumeFromCoordinates, {"GENIE_UNIT_LEGTH",
                                                                                              "Interaction_vtxX",
                                                                                              "Interaction_vtxY",
                                                                                              "Interaction_vtxZ"})
             // !!! FIDUCIAL VOLUME CUT !!!
             .Define("isInFiducialVolume",          GeoUtils::DRIFT::IsInFiducialVolume, {"InteractionVolume",
                                                                                          "GENIE_UNIT_LEGTH",
                                                                                          "Interaction_vtxX",
                                                                                          "Interaction_vtxY",
                                                                                          "Interaction_vtxZ"})
             .Filter("isInFiducialVolume")
             // !!!
             .Define("InteractionTarget",           RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             // there is some problem with this -> TGeoMananger doesn't localize proprly the volume if file is .gdml and not .root
             .Define("InteractionTargetFromGEO",    RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
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
                incoming neutrino
             */
             .Define("IncomingNeutrino",            RDFUtils::GENIE::GetIncomingNeutrino, {"Particles"})
             .Define("IncomingNeutrinoPDG",         [](const genie::GHepParticle& nu){return nu.Pdg();}, {"IncomingNeutrino"})
             .Define("IncomingNeutrinoName",        [](const int nu_pdg){return GenieUtils::PDG2Name(nu_pdg);}, {"IncomingNeutrinoPDG"})
             .Define("IncomingNeutrinoP4",          [](const genie::GHepParticle& nu){return *nu.GetP4();}, {"IncomingNeutrino"})
             .Define("NuDirection", "IncomingNeutrinoP4.Vect().Unit()")
             /* 
                nucleon target (free proton - i.e. interaction on H - is NOT included using genie::kIStNucleonTarget in the category by GENIE) 
             */
             .Define("NucleonTarget",               RDFUtils::GENIE::GetNucleonTarget, {"Particles"})
             .Define("NucleonTargetNamePDG",        [](const genie::GHepParticle& t){return t.Pdg();}, {"NucleonTarget"})
             .Define("NucleonTargetName",           [](const int t_pdg){return GenieUtils::PDG2Name(t_pdg);}, {"NucleonTargetNamePDG"})
             .Define("NucleonTargetP4",             [](const genie::GHepParticle& t){return *t.GetP4();}, {"NucleonTarget"})
             /*
                initial state : incoming neutrino + nucleon target (include events on Hydrogen by hands not included otherwise)
             */
             .Define("InitialState",                            RDFUtils::GENIE::GetInitialState, {"Particles"})
             .Define("InitialStatePDG",                         RDFUtils::GENIE::GetPDG, {"InitialState"})
             .Define("InitialStateNames",                       RDFUtils::GENIE::PDG2Name, {"InitialState"})
             .Define("InitialStateMomenta",                     RDFUtils::GENIE::GetMomentum, {"InitialState"})
             .Define("InitialStateTotal4Momentum",              RDFUtils::GENIE::SumLorentzVectors, {"InitialStateMomenta"})
             .Define("InitialStateEmissionAngle",               [](TVector3 nu, TLorentzVector v){return nu.Angle(v.Vect());},{"NuDirection","InitialStateTotal4Momentum"})
             /*
                final state lepton
             */
             .Define("FinalStateLepton",                        RDFUtils::GENIE::GetFinalStateLepton, {"Particles"})
             .Define("FinalStateLeptonPDG",                     [](const genie::GHepParticle& l){return l.Pdg();}, {"FinalStateLepton"})
             .Define("FinalStateLeptonNames",                   [](const int l_pdg){return GenieUtils::PDG2Name(l_pdg);}, {"FinalStateLeptonPDG"})
             .Define("FinalStateLepton4Momentum",               [](const genie::GHepParticle& l){return *l.GetP4();}, {"FinalStateLepton"})
             .Define("FinalStateLeptonEmissionAngle",           [](TVector3 nu, TLorentzVector v){return nu.Angle(v.Vect());},{"NuDirection","FinalStateLepton4Momentum"})
             /*
                primary state hadronic system
             */
             .Define("PrimaryStateHadronicSystem",              RDFUtils::GENIE::PrimaryStateHadronicSystem, {"Particles", "EventType"})
             .Define("PrimaryStateHadronicSystemPDG",           RDFUtils::GENIE::GetPDG, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemNames",         RDFUtils::GENIE::PDG2Name, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemMomenta",       RDFUtils::GENIE::GetMomentum, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemTotal4Momentum",RDFUtils::GENIE::SumLorentzVectors, {"PrimaryStateHadronicSystemMomenta"})
             .Define("PrimaryStateHadronicSystemEmissionAngle", [](TVector3 nu, TLorentzVector v){return nu.Angle(v.Vect());},{"NuDirection","PrimaryStateHadronicSystemTotal4Momentum"})
             .Define("PrimaryStateHadronicSystemKinE",          RDFUtils::GENIE::GetKinE, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemTotalKinE",     RDFUtils::SumDuble, {"PrimaryStateHadronicSystemKinE"})
             .Define("PrimaryStateHadronicSystemTopology",      RDFUtils::GENIE::GetInteractionTopology, {"PrimaryStateHadronicSystem"})
             .Define("PrimaryStateHadronicSystemTopology_code", [](GenieUtils::event_topology t){return t.unique_code();}, {"PrimaryStateHadronicSystemTopology"})
             .Define("PrimaryStateHadronicSystemTopology_name", [](GenieUtils::event_topology t){return t.Name();}, {"PrimaryStateHadronicSystemTopology"})
             /*
                final state hadronic system
             */
             .Define("FinalStateHadronicSystem",                RDFUtils::GENIE::FinalStateHadronicSystem, {"Particles"})
             .Define("FinalStateHadronicSystemPDG",             RDFUtils::GENIE::GetPDG, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemNames",           RDFUtils::GENIE::PDG2Name, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemMomenta",         RDFUtils::GENIE::GetMomentum, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemTotal4Momentum",  RDFUtils::GENIE::SumLorentzVectors, {"FinalStateHadronicSystemMomenta"})
             .Define("FinalStateHadronicSystemEmissionAngle",   [](TVector3 nu, TLorentzVector v){return nu.Angle(v.Vect());},{"NuDirection","FinalStateHadronicSystemTotal4Momentum"})
             .Define("FinalStateHadronicSystemKinE",            RDFUtils::GENIE::GetKinE, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemTotalKinE",       RDFUtils::SumDuble, {"FinalStateHadronicSystemKinE"})
             .Define("FinalStateHadronicSystemTopology",        RDFUtils::GENIE::GetInteractionTopology, {"FinalStateHadronicSystem"})
             .Define("FinalStateHadronicSystemTopology_code",   [](GenieUtils::event_topology t){return t.unique_code();}, {"FinalStateHadronicSystemTopology"})
             .Define("FinalStateHadronicSystemTopology_name",   [](GenieUtils::event_topology t){return t.Name();}, {"FinalStateHadronicSystemTopology"})
              /*
                final state nuclear remnant
                low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle
              */
             .Define("NuclearRemnant",                RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStFinalStateNuclearRemnant>, {"Particles"})
             .Define("NuclearRemnantPDG",             RDFUtils::GENIE::GetPDG, {"NuclearRemnant"})
             .Define("NuclearRemnantMomenta",         RDFUtils::GENIE::GetMomentum, {"NuclearRemnant"})
             .Define("NuclearTotal4Momentum",         RDFUtils::GENIE::SumLorentzVectors, {"NuclearRemnantMomenta"})
             .Define("NuclearRemnantKinE",            RDFUtils::GENIE::GetKinE, {"NuclearRemnant"})
             .Define("NuclearRemnantTotalKinE",       RDFUtils::SumDuble, {"NuclearRemnantKinE"})
              /*
                Overall variables to check overall energy conservation
              */
             .Define("PrimaryStateTotal4Momantum" , "FinalStateLepton4Momentum + PrimaryStateHadronicSystemTotal4Momentum")
             .Define("FinalStateTotal4Momantum" , "FinalStateLepton4Momentum + FinalStateHadronicSystemTotal4Momentum + NuclearTotal4Momentum")
              /*
                For channel selection
              */
             // charge multiplicity
             .Define("NofFinalStateChargedParticles",               RDFUtils::GENIE::CountChargedParticles, {"FinalStateLepton","FinalStateHadronicSystem"}) 
             // Prediction on neutrons
             .Define("ExpectedNeutrinoP4FromMuon",                  RDFUtils::GENIE::GetNup4FromMu, {"PROTON_MASS_GeV",
                                                                                                     "NEUTRON_MASS_GeV",
                                                                                                     "MUON_MASS_GeV",
                                                                                                     "FinalStateLepton4Momentum",
                                                                                                     "NuDirection",
                                                                                                     })
             .Define("ExpectedHadronSystP3",                        "ExpectedNeutrinoP4FromMuon.Vect() - FinalStateLepton4Momentum.Vect()")
             .Define("ExpectedHadronSystEnergy",                    "ExpectedNeutrinoP4FromMuon.T() + PROTON_MASS_GeV - FinalStateLepton4Momentum.T()")
             .Define("ExpectedNeutronArrivalPositionECAL",          RDFUtils::GENIE::NeutronArrivalPosECAL, {"GENIE_UNIT_LEGTH",
                                                                                                             "Interaction_vtxX",
                                                                                                             "Interaction_vtxY",
                                                                                                             "Interaction_vtxZ",
                                                                                                             "ExpectedHadronSystP3"})
             .Define("ExpectedNeutronTOF",                          RDFUtils::GENIE::GetNeutronTOF, {"GENIE_UNIT_LEGTH",
                                                                                                     "GENIE_UNIT_ENERGY",
                                                                                                     "Interaction_vtxX",
                                                                                                     "Interaction_vtxY",
                                                                                                     "Interaction_vtxZ",
                                                                                                     "ExpectedNeutronArrivalPositionECAL",
                                                                                                     "ExpectedHadronSystP3"})                                                                                                             
             .Define("ExpectedFiredModuleByNeutron",                GeoUtils::ECAL::GetModuleIdFromPoint, {"EDEP_UNIT_LEGTH",
                                                                                                           "ExpectedNeutronArrivalPositionECAL"})                                                                                                                                                                                                
             // A Novel Approach to Neutrino-Hydrogen Measurements
             .Define("FinalStateLeptonTransverseP",                 "FinalStateLepton4Momentum.Vect() - FinalStateLepton4Momentum.Vect().Dot(NuDirection) * NuDirection")
             .Define("FinalStateHadronicSystemTransverseP",         "FinalStateHadronicSystemTotal4Momentum.Vect() - FinalStateHadronicSystemTotal4Momentum.Vect().Dot(NuDirection) * NuDirection")
             .Define("MissingTransverseMomentum",                   "(FinalStateLeptonTransverseP + FinalStateHadronicSystemTransverseP).Mag()")
             .Define("RmH",                                         "(MissingTransverseMomentum - FinalStateHadronicSystemTransverseP.Mag())/(MissingTransverseMomentum + FinalStateHadronicSystemTransverseP.Mag())")
             // https://arxiv.org/pdf/1507.00967
             .Define("DoubleTransverseAxis",                        "(NuDirection.Cross(FinalStateLepton4Momentum.Vect())).Unit()")
             .Define("FinalStateHadronicSystem_DoubleTransverseP",  RDFUtils::GENIE::DotProductWithAxis, {"DoubleTransverseAxis","FinalStateHadronicSystemMomenta"})
             .Define("DoubleTransverseMomentumImbalance",           RDFUtils::SumDuble, {"FinalStateHadronicSystem_DoubleTransverseP"})
             ;
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(ROOT::RDF::RNode& df){
    /* These variables allows the separation of interacions occurring on free proton (H)
    from those occuring on Carbon. Variables are explained in arXiv:2102.03346v1 and arXiv:1809.08752v2[Petti] */
    //____________________________________________________________________________
    return df
            // FinalHadronicSystemP4_TT : projection of FinalHadronicSystemP4 onto DoubleTransverseAxis perpendicular to neutrino direction
             .Define("DoubleTransverseAxis",       RDFUtils::TLVectorCrossProduct, {"InitialStateNuMu_P4","FinalStateMuonsP4"})
             ;
}

//RDFUtils::GEO______________________________________________________________________

std::string RDFUtils::GEO::GetVolumeFromCoordinates(std::string units, double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }
    // // convert to cm
    // if(units == "mm"){
    //     x = x * 0.1;
    //     y = y * 0.1;
    //     z = z * 0.1;
    // }else if(units == "m"){
    //     x = x * 100.;
    //     y = y * 100.;
    //     z = z * 100.;
    // }
    // convert to mm
    if(units == "cm"){
        x = x * 10.;
        y = y * 10.;
        z = z * 10.;
    }else if(units == "m"){
        x = x * 1000.;
        y = y * 1000.;
        z = z * 1000.;
    }
    auto navigator = geo->GetCurrentNavigator();

    if(!navigator) navigator = geo->AddNavigator();

    // TGeoManager wants cm units
    auto volume = navigator->FindNode(x, y, z)->GetVolume()->GetName();
   
    // if(TString::Format(volume).Contains("C3H6")){
    //     std::cout << "volume : " << volume << "x, y, z " << x << " " << y << " " << z << "\n";
    // }
    return volume;
}

std::string RDFUtils::GEO::GetMaterialFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }

    auto navigator = geo->GetCurrentNavigator();

    if(!navigator) navigator = geo->AddNavigator();

    auto material = navigator->FindNode(x*100., y*100., z*100.)->GetVolume()->GetMaterial()->GetName();

    return material;
}

//RDFUtils::EDEPSIM_______________________________________________________________

double RDFUtils::EDEPSIM::NeutrinoEnergyFromCCQEonH(const TLorentzVector& muon, double muon_angle){
    double muon_mass = 0.105658;
    double proton_mass = 0.938272;
    double neutron_mass = 0.939565;
    double neutrino_energy = (neutron_mass*neutron_mass - muon_mass*muon_mass - proton_mass*proton_mass + 2*proton_mass*muon.T()) / (2*(proton_mass - muon.T() + muon.Vect().Mag()*cos(muon_angle)));
    return neutrino_energy;
}

std::string RDFUtils::EDEPSIM::NOSPILL::EventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V){
    TString s = V[0].GetReaction();
    if(s.Contains("DIS")){
        return "DIS";
    }else if(s.Contains("RES")){
        return "RES";
    }else if(s.Contains("QES")){
        return "QES";
    }else if(s.Contains("COH")){
        return "COH";
    }else if(s.Contains("MEC")){
        return "MEC";
    }else if(s.Contains("Unknown")){
        return "Unknown";
    }
    else{
        return "Other";
    }
}

bool RDFUtils::EDEPSIM::NOSPILL::IsCCQEonH(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V){
    TString s = V[0].GetReaction();
    if(s.Contains("nu:-14;tgt:1000010010;N:2212;proc:Weak[CC],QES;")){
        return true;
    }else{
        return false;
    }
}

std::string RDFUtils::EDEPSIM::NOSPILL::FileName(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V){
    return V[0].GetFilename();
}

int RDFUtils::EDEPSIM::NOSPILL::NofPrimaries(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    auto Primaries = v[0].Particles;
    return Primaries.size();
}

ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesPDG(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    ROOT::VecOps::RVec<int> pdgs;
    auto Primaries = v[0].Particles;
    for(const auto& primary : Primaries){
        pdgs.push_back(primary.GetPDGCode());
    }
    return pdgs;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesP4(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    ROOT::VecOps::RVec<TLorentzVector> P4;
    auto Primaries = v[0].Particles;
    for(const auto& primary : Primaries){
        P4.push_back(primary.GetMomentum());
    }
    return P4;
}

TLorentzVector RDFUtils::EDEPSIM::NOSPILL::GetPrimaryHadronicP4(const ROOT::VecOps::RVec<TLorentzVector>& momenta, 
                                                                const ROOT::VecOps::RVec<int> pdgs){
    TLorentzVector hadronic_P4 = {0.,0.,0.,0.};
    for (auto i = 0u; i < pdgs.size(); i++)
    {
        if((fabs(pdgs[i])!=13) & (fabs(pdgs[i])!=12)) hadronic_P4 += momenta[i];
    }
    return hadronic_P4;                                                                                            
}

TLorentzVector RDFUtils::EDEPSIM::NOSPILL::GetPrimaryLeptonP4(const ROOT::VecOps::RVec<TLorentzVector>& momenta, 
                                                              const ROOT::VecOps::RVec<int> pdgs){
    TLorentzVector lepton = {0.,0.,0.,0.};                                                            
    for (auto i = 0u; i < pdgs.size(); i++)
    {
        if((fabs(pdgs[i])==13) | (fabs(pdgs[i])==12)) lepton = momenta[i];
    }
    return lepton;                                                                  
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::NOSPILL::GetEmissionAngle(const TVector3 nuDirection, const ROOT::VecOps::RVec<TLorentzVector>& primariesP4){
    ROOT::VecOps::RVec<double> angles;
    for (const auto& primaryP4 : primariesP4)
    {
        angles.push_back(primaryP4.Vect().Angle(nuDirection));
    }
    return angles;
}

ROOT::VecOps::RVec<std::string> RDFUtils::EDEPSIM::NOSPILL::PDG2Name(const ROOT::VecOps::RVec<int>& pdgs){
    // convert particles pdg into names
    ROOT::VecOps::RVec<std::string> names;
    for(const auto& pdg : pdgs){
        names.push_back(GenieUtils::PDG2Name(pdg));
    }
    return names;
}

ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesTrackId(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    ROOT::VecOps::RVec<int> ids;
    auto Primaries = v[0].Particles;
    for(const auto& primary : Primaries){
        ids.push_back(primary.GetTrackId());
    }
    return ids;
}


ROOT::VecOps::RVec<EDepUtils::track_hits> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesHits(TG4Event& ev, 
                                                                                       const ROOT::VecOps::RVec<int>& ids,
                                                                                       const ROOT::VecOps::RVec<int>& pdgs){
    /*
        return a vector of trac_hits, that is group hits
        that belongs to a given primary
    */
    std::map<int, EDepUtils::track_hits> map_grouped_hits;

    // initialize the map, each key a primary id
    for (size_t i = 0; i < ids.size(); i++)
    {
        EDepUtils::track_hits t;
        t.track_id = ids[i];
        t.pdg = pdgs[i];
        map_grouped_hits[ids[i]] = t;
    }

    for(auto hit: ev.SegmentDetectors["EMCalSci"]){
        // retreive the track id of the most important particle associated with this hit
        auto hit_id = hit.GetPrimaryId();

        // check if the track id associated to the hit belongs to one of the primaries id
        auto it = std::find(ids.begin(), ids.end(), hit_id);

        if (it != ids.end()){
            // update info of EDepUtils::track_hits that belongs to the primary whose id is hit_id
            auto hit_middle = (hit.GetStop() + hit.GetStart())*0.5;
            map_grouped_hits[hit_id].hit_edep.push_back(hit.GetEnergyDeposit());
            map_grouped_hits[hit_id].hit_LorentzVector.push_back(hit_middle);
        }
    }

    ROOT::VecOps::RVec<EDepUtils::track_hits> T;

    for(const auto& pair : map_grouped_hits) T.push_back(pair.second);

    return T;
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesEDepECAL(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_track_hits){
    ROOT::VecOps::RVec<double> primaries_edep;
    for(const auto& track_hits : vector_track_hits){
        double de = 0.;
        for(auto& e : track_hits.hit_edep) de += e;
        primaries_edep.push_back(de);
    }
    return primaries_edep;
}

TLorentzVector FindEarliestHit(const std::vector<TLorentzVector>& hits){
    auto minElementIt = std::min_element(hits.begin(), hits.end(), [](const TLorentzVector& a, const TLorentzVector& b) {
        return a.T() < b.T();
    });
    return (minElementIt != hits.end()) ? *minElementIt : TLorentzVector();
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::EDEPSIM::NOSPILL::GetPrimariesFirstHitECAL(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_primary_hits){
    ROOT::VecOps::RVec<TLorentzVector> primary_hit_ECAL;
    for(const auto& primary_hits : vector_primary_hits){ // run over each primary
        if(primary_hits.hit_LorentzVector.size()==0){ // primary has no hits in the ecal
            primary_hit_ECAL.push_back({0.,0.,0.,0.});
        }else{
            // find the earlies hit of a primary
            primary_hit_ECAL.push_back(FindEarliestHit(primary_hits.hit_LorentzVector)); 
        }
    }
    return primary_hit_ECAL;
}

template<int coordinate>
ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_primary_hits){
    ROOT::VecOps::RVec<double> vertices;
    for(const auto& primary_hits : vector_primary_hits){ // run over each primary
        if(primary_hits.hit_LorentzVector.size()==0){ // primary has no hits in the ecal
            vertices.push_back(-999.);
        }else{
            vertices.push_back(primary_hits.hit_LorentzVector[0][coordinate]);
        }
    }
    return vertices;
}

ROOT::RDF::RNode RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df){
    return df
             .Define("FileName",                RDFUtils::EDEPSIM::NOSPILL::FileName, {"Primaries"})
             .Define("NofEvents",               RDFUtils::EDEPSIM::SPILL::NofEventsPerSpill, {"Primaries"})
             .Filter("isInFiducialVolume")                                                                                       
             .Define("NofPrimaries",            RDFUtils::EDEPSIM::NOSPILL::NofPrimaries, {"Primaries"})
             .Define("PrimariesPDG",            RDFUtils::EDEPSIM::NOSPILL::GetPrimariesPDG, {"Primaries"})
             .Define("PrimariesName",           RDFUtils::EDEPSIM::NOSPILL::PDG2Name, {"PrimariesPDG"})
             .Define("PrimariesTrackId",        RDFUtils::EDEPSIM::NOSPILL::GetPrimariesTrackId, {"Primaries"})
             .Define("PrimariesP4",             RDFUtils::EDEPSIM::NOSPILL::GetPrimariesP4, {"Primaries"})
              // ECAL track_hits
             .Define("PrimariesHitsECAL",       RDFUtils::EDEPSIM::NOSPILL::GetPrimariesHits, {"Event", "PrimariesTrackId", "PrimariesPDG"})
             .Define("PrimariesHitsX",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<0>, {"PrimariesHitsECAL"})
             .Define("PrimariesHitsY",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<1>, {"PrimariesHitsECAL"})
             .Define("PrimariesHitsZ",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<2>, {"PrimariesHitsECAL"})
             // ECAL hits
            //  .Define("PrimariesHitsECALVol",    RDFUtils::GEO::GetVolumeFromCoordinates, {"EDEP_UNIT_LEGTH",
            //                                                                               "PrimariesHitsX",
            //                                                                               "PrimariesHitsY",
            //                                                                               "PrimariesHitsZ"})
             .Define("PrimariesEmissionAngle",  RDFUtils::EDEPSIM::NOSPILL::GetEmissionAngle, {"NuDirection", "PrimariesP4"})
             .Define("PrimariesFirstHitECAL",   RDFUtils::EDEPSIM::NOSPILL::GetPrimariesFirstHitECAL, {"PrimariesHitsECAL"})
             .Define("PrimariesEDepECAL",       RDFUtils::EDEPSIM::NOSPILL::GetPrimariesEDepECAL, {"PrimariesHitsECAL"})
            //
            //  .Define("EventType",               RDFUtils::EDEPSIM::NOSPILL::EventType, {"Primaries"})
            //  .Define("PrimariesVertexX",        [](const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices){return vertices[0].GetPosition()[0];}, {"Primaries"})
            //  .Define("Interaction_vtxX","EvtVtx[0]")
            //  .Define("Interaction_vtxY","EvtVtx[1]")
            //  .Define("Interaction_vtxZ","EvtVtx[2]")
            //  .Define("Interaction_vtxT","EvtVtx[3]")
             // !!! FIDUCIAL VOLUME CUT !!!
            //  .Define("InteractionVolume",           RDFUtils::GEO::GetVolumeFromCoordinates, {"GENIE_UNIT_LEGTH",
            //                                                                                   "Interaction_vtxX",
            //                                                                                   "Interaction_vtxY",
            //                                                                                   "Interaction_vtxZ"})
            //  // !!! FIDUCIAL VOLUME CUT !!!
            //  .Define("isInFiducialVolume",          GeoUtils::DRIFT::IsInFiducialVolume, {"InteractionVolume",
            //                                                                               "GENIE_UNIT_LEGTH",
            //                                                                               "Interaction_vtxX",
            //                                                                               "Interaction_vtxY",
            //                                                                               "Interaction_vtxZ"})
             /*
                All primary particles
            //  */
            //  .Define("CCQEonHydrogen",          RDFUtils::EDEPSIM::NOSPILL::IsCCQEonH,  {"Primaries"})
            //  .Define("NeutrinoP4",              RDFUtils::GENIE::SumLorentzVectors, {"PrimariesP4"})
            //  .Define("NeutrinoE", "NeutrinoP4[3]")
            //  .Define("NeutrinoDirection",       "NeutrinoP4.Vect().Unit()")
            //  .Define("PrimaryLeptonP4",         RDFUtils::EDEPSIM::NOSPILL::GetPrimaryLeptonP4,   {"PrimariesP4", "PrimariesPDG"})
            //  .Define("PrimaryHadronicSystP4",   RDFUtils::EDEPSIM::NOSPILL::GetPrimaryHadronicP4, {"PrimariesP4", "PrimariesPDG"})
            //  .Define("PrimaryHadronicSystEmissionAngle", [](const TLorentzVector& v,const TLorentzVector& w){return v.Vect().Angle(w.Vect());}, {"NeutrinoP4", "PrimaryHadronicSystP4"})
            //  .Define("PrimaryLeptonEmissionAngle", "PrimariesEmissionAngle[0]")
             // EXPECTED PRIMARY NEUTRON FROM INTERACTION ON ARGO
             ;
}

// SPILL FUNCTIONS _________________________________________________________________________

int RDFUtils::EDEPSIM::SPILL::NofEventsPerSpill(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices){
    return vertices.size();
}

ROOT::VecOps::RVec<std::string> RDFUtils::EDEPSIM::SPILL::EventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V){
    ROOT::VecOps::RVec<std::string> types = {};
    for(const auto& v : V){
        TString s = v.GetReaction();
        if(s.Contains("DIS")){
            types.push_back("DIS");
        }else if(s.Contains("RES")){
            types.push_back("RES");
        }else if(s.Contains("QES")){
            types.push_back("QES");
        }else if(s.Contains("COH")){
            types.push_back("COH");
        }else if(s.Contains("MEC")){
            types.push_back("MEC");
        }else if(s.Contains("Unknown")){
            types.push_back("Unknown");
        }
        else{
            types.push_back("Other");
        }
    }
    return types;
}

ROOT::VecOps::RVec<std::string> RDFUtils::EDEPSIM::SPILL::FileName(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V){
    ROOT::VecOps::RVec<std::string> files = {};
    for(const auto& v : V) files.push_back(v.GetFilename());
    return files;
}

template<int coordinate>
ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::SPILL::PrimariesVertex(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices)
{
    ROOT::VecOps::RVec<double> vertices_positions;
    for(auto& vertex: vertices)
    {
        auto vertex_position = vertex.GetPosition();
        vertices_positions.emplace_back(vertex_position[coordinate]);
    }
    return vertices_positions;
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::SPILL::PrimariesVertexTimeDiff(ROOT::VecOps::RVec<double>& primaries_time)
{
    ROOT::VecOps::RVec<double> time_diff;

    for(auto i=1u; i<primaries_time.size(); i++)
    {
        time_diff.emplace_back(primaries_time[i]-primaries_time[i-1]);
    }
    return time_diff;
}

ROOT::RDF::RNode RDFUtils::EDEPSIM::SPILL::AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df){
    return df
             .Define("NofEventsPerSpill", RDFUtils::EDEPSIM::SPILL::NofEventsPerSpill, {"Primaries"})
             .Define("PrimariesVertexX", RDFUtils::EDEPSIM::SPILL::PrimariesVertex<0>,{"Primaries"})
             .Define("PrimariesVertexY", RDFUtils::EDEPSIM::SPILL::PrimariesVertex<1>,{"Primaries"})
             .Define("PrimariesVertexZ", RDFUtils::EDEPSIM::SPILL::PrimariesVertex<2>,{"Primaries"})
             .Define("PrimariesVertexT", RDFUtils::EDEPSIM::SPILL::PrimariesVertex<3>,{"Primaries"})
             .Define("PrimariesVertexTimeDiff", RDFUtils::EDEPSIM::SPILL::PrimariesVertexTimeDiff,{"PrimariesVertexT"})
             ;
}

//RDFUtils::DIGIT_______________________________________________________________

int RDFUtils::DIGIT::NofFiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules){
    std::map<int, int> fired_modules_map;
    for(const auto mod : fired_cells_modules){
        fired_modules_map[mod] += 1;
    }
    int nof_fired_modules = fired_modules_map.size();
    return nof_fired_modules;
}

std::vector<int> RDFUtils::DIGIT::FiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules){
    std::map<int, int> fired_modules_map;
    std::vector<int> fired_modules;
    for(const auto mod : fired_cells_modules){
        fired_modules_map[mod] += 1;
    }
    for (const auto& pair : fired_modules_map) {
        fired_modules.push_back(pair.first);
    }
    return fired_modules;
}

int RDFUtils::DIGIT::NofClusters(const ROOT::VecOps::RVec<cluster>& clusters){
    return clusters.size();
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::DIGIT::Cluster2Vertex4Distance(double x,
                                                                            double y,
                                                                            double z,
                                                                            double t,
                                                                            const ROOT::VecOps::RVec<cluster>& clusters){
    TLorentzVector vertex = {x*1e3, y*1e3, z*1e3, t*1e3};
    ROOT::VecOps::RVec<TLorentzVector> differences;
    for(const auto& cluster : clusters){
        TLorentzVector cluster_X4 = {cluster.x, cluster.y, cluster.z, cluster.t};
        differences.push_back(cluster_X4 - vertex);
    }
    return differences;                                                                       
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::DIGIT::GetClusterX4(const ROOT::VecOps::RVec<cluster>& clusters){
    ROOT::VecOps::RVec<TLorentzVector> X4;
    for(const auto& cluster : clusters) X4.push_back({cluster.x, cluster.y, cluster.z, cluster.t});
    return X4;
}

template<int side>
ROOT::VecOps::RVec<double> RDFUtils::DIGIT::FiredECALGetTDC(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<double> TDCs;
    for(const auto& cell : cells){
        if(side==1){ // cell side 1
            if(cell.ps1.size()!=0){
                TDCs.push_back(cell.ps1[0].tdc);
            }else{
                TDCs.push_back(0.);
            }
        }else{ // cell side 2
            if(cell.ps2.size()!=0){
                TDCs.push_back(cell.ps2[0].tdc);
            }else{
                TDCs.push_back(0.);
            }
        }
    }
    return TDCs;
}

ROOT::VecOps::RVec<TLorentzVector> RDFUtils::DIGIT::ReconstructHitFromCell(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<TLorentzVector> reco_hits;
    for (const auto& cell : cells)
    {
        int broken_side = 0;
        double tdc1 = 0., tdc2 = 0.;
        TLorentzVector hit_reco;

        if(!GeoUtils::ECAL::isCellComplete(cell, broken_side)){ // broken cell
            // what we do in case of broke cells ?
        } else{ // complete cell
            tdc1 = cell.ps1.at(0).tdc;
            tdc2 = cell.ps2.at(0).tdc;
        }
        // if(cell.id < 25000){ // Barrel
        if(cell.mod < 30){ // Barrel
            hit_reco = {cell.x - GeoUtils::ECAL::XfromTDC(tdc1, tdc2),
                        cell.y,
                        cell.z,
                        GeoUtils::ECAL::TfromTDC(tdc1, tdc2, cell.l)};
        }else{ // Endcap
            hit_reco = {cell.x,
                        cell.y - GeoUtils::ECAL::XfromTDC(tdc1, tdc2),
                        cell.z,
                        GeoUtils::ECAL::TfromTDC(tdc1, tdc2, cell.l)};
        }
        reco_hits.push_back(hit_reco);
    }
    return reco_hits;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::STDistance(const TVector3& expected_neutron_hit3,
                                                       const double exp_neutron_tof,
                                                       const ROOT::VecOps::RVec<TLorentzVector>& cells_reco_hits){
    /*
        Return the space time distance btw 
        reconstruted hits from cells' tdc 
        and the expected neutron hit
    */
   ROOT::VecOps::RVec<double> ds2; // space-time distances
   double c = 299.792; // ligth speed [mm/ns]
   TLorentzVector expected_neutron_hit = {expected_neutron_hit3, exp_neutron_tof};
   for (const auto& cell_reco_hit : cells_reco_hits)
   {
        TLorentzVector delta = expected_neutron_hit - cell_reco_hit;
        double d = (delta.T()*c) * (delta.T()*c) - delta.Vect().Mag2(); 
        ds2.push_back(d);
   }
   return ds2;
}

ROOT::VecOps::RVec<int> RDFUtils::DIGIT::IsCellComplete(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> isComplete;
    for(const auto& cell : cells){
        if((cell.ps1.size()!=0) & (cell.ps2.size()!=0)){
            isComplete.push_back(1);
        }else{
            isComplete.push_back(0);
        }
    }
    return isComplete;
}

ROOT::RDF::RNode RDFUtils::DIGIT::AddColumnsFromDigit(ROOT::RDF::RNode& df){
    return df
            .Define("NofEventFiredModules",         RDFUtils::DIGIT::NofFiredECALMods , {"dg_cell.mod"})
            .Define("EventFiredModules",            RDFUtils::DIGIT::FiredECALMods , {"dg_cell.mod"})
            .Define("Fired_Cells_mod",              "dg_cell.mod")
            .Define("Fired_Cells_id",               "dg_cell.id")
            .Define("Fired_Cells_x",                "dg_cell.x")
            .Define("Fired_Cells_y",                "dg_cell.y")
            .Define("Fired_Cells_z",                "dg_cell.z")
            .Define("isCellComplete",               RDFUtils::DIGIT::IsCellComplete, {"dg_cell"})
            .Define("Fired_Cells_tdc1",             RDFUtils::DIGIT::FiredECALGetTDC<1>, {"dg_cell"})
            .Define("Fired_Cells_tdc2",             RDFUtils::DIGIT::FiredECALGetTDC<2>, {"dg_cell"})
            .Define("Cell_Reconstructed_hits",      RDFUtils::DIGIT::ReconstructHitFromCell, {"dg_cell"})
            .Define("STDistToNeutronExpectedHit",   RDFUtils::DIGIT::STDistance, {"ExpectedNeutronArrivalPositionECAL",
                                                                                  "ExpectedNeutronTOF",
                                                                                  "Cell_Reconstructed_hits"})
            // ECAL CLUSTER INFO
            // .Define("NofEventClusters",             RDFUtils::DIGIT::NofClusters, {"cluster"})
            // .Define("ClusterX4",                    RDFUtils::DIGIT::GetClusterX4, {"cluster"})
            // .Define("Cluster2Vertex4Distance",      RDFUtils::DIGIT::Cluster2Vertex4Distance, {"Interaction_vtxX",
            //                                                                               "Interaction_vtxY",
            //                                                                               "Interaction_vtxZ",
            //                                                                               "Interaction_vtxT",
            //                                                                               "cluster"})
            ;
}
