#include "RDataFrameUtils.h"

#include "TChain.h"
#include "TMath.h"

TGeoManager* geo = nullptr;
std::vector<dg_wire>* RecoUtils::event_digits = nullptr;

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


ROOT::RDF::RNode RDFUtils::Filter(ROOT::RDF::RNode& df, const char* condition, bool report){
    auto df_filtered = df.Filter(condition);
    if(report)
        df_filtered.Report()->Print();
    return df_filtered;
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

template<int coord>
ROOT::VecOps::RVec<double> RDFUtils::GetComponent3(const ROOT::VecOps::RVec<TVector3>& vTV){

    ROOT::VecOps::RVec<double> v(vTV.size());

    for (size_t i = 0; i < vTV.size(); i++)
    {
        v[i] = vTV[i][coord];
    }
    
    return v;
}

ROOT::VecOps::RVec<int> RDFUtils::FindMinimum(const ROOT::VecOps::RVec<double> v){
    /**
     * @brief Return a vector of 0s and 1s, where 1 indicates the minimum value in the input vector.
     */
    ROOT::VecOps::RVec<int> isMinimum(v.size(), 0);

    ROOT::VecOps::RVec<double> valid_values = v[v != RDFUtils::DEFAULT_NO_DATA];

    if (!valid_values.empty()) {
        double min_v = ROOT::VecOps::Min(valid_values);

        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] == min_v && v[i] != RDFUtils::DEFAULT_NO_DATA){ // redundant check
                isMinimum[i] = 1;
            }
        }
    }
    return isMinimum;
}

ROOT::VecOps::RVec<int> RDFUtils::FindMinimum2(const ROOT::VecOps::RVec<double> v,
                                               const ROOT::VecOps::RVec<int> condition){
    /**
     * @brief same as RDFUtils::FindMinimum but with additional condition
     */
    ROOT::VecOps::RVec<int> isMinimum(v.size(), 0);

    if(v.size() != condition.size()){
        throw std::runtime_error("v and condition must have the same length");
    }

    ROOT::VecOps::RVec<double> valid_values = v[v != RDFUtils::DEFAULT_NO_DATA];

    if (!valid_values.empty()) {
        double min_v = ROOT::VecOps::Min(valid_values);

        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] == min_v && condition[i] == 1){ // redoundant check
                isMinimum[i] = 1;
            }
        }
    }
    return isMinimum;
}

double RDFUtils::FindMinimum3(const ROOT::VecOps::RVec<double> values,
                              const ROOT::VecOps::RVec<int> condition){
    /**
     * @brief same as RDFUtils::FindMinimum but with additional condition
     * return 1 value only.

     * INPUT: can be a vector of reconstructed times and the relative condition[i]
     * that the value[i] satisfies

     * OUTPUT: the lowest value that whose condition[i] == 1
     * 
    */                            
    if(values.size() != condition.size()){
        throw std::runtime_error("v and condition must have the same length");
    }

    ROOT::VecOps::RVec<double> values_condition = values[condition == 1];

    if(values_condition.size() == 0){
        return RDFUtils::DEFAULT_NO_DATA;
    }else{
        return ROOT::VecOps::Min(values_condition);
    }
}

ROOT::VecOps::RVec<TVector3> RDFUtils::GetVect(const ROOT::VecOps::RVec<TLorentzVector>& vTL){

    ROOT::VecOps::RVec<TVector3> v(vTL.size());
    
    for (size_t i = 0; i < vTL.size(); i++) v[i] = vTL[i].Vect();

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

ROOT::VecOps::RVec<int> RDFUtils::IsNULLTLV(const ROOT::VecOps::RVec<TLorentzVector>& V){
    ROOT::VecOps::RVec<bool> isNULL;
    for(auto const& v : V){
        if(v.Mag() == 0.){
            isNULL.push_back(1);
        }else{
            isNULL.push_back(0);
        }
    }
    return isNULL;
}

ROOT::VecOps::RVec<int> RDFUtils::CompareVectors(const ROOT::VecOps::RVec<double>& V1, 
                                                  const ROOT::VecOps::RVec<double>& V2){
    ROOT::VecOps::RVec<int> isEqual(V1.size());
    const double tolerance = 1e-6;
    for (size_t i = 0; i < V1.size(); i++) {
        isEqual[i] = (fabs(V1[i] - V2[i]) < tolerance) ? 1 : 0;
    }
    return isEqual;
}

ROOT::VecOps::RVec<int> RDFUtils::IsValidVector(const ROOT::VecOps::RVec<double>& v){
    ROOT::VecOps::RVec<int> isValid(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        if (v[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            isValid[i] = 0;
        }else{
            isValid[i] = 1;
        }
        
    }
    return isValid;
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

ROOT::VecOps::RVec<double> RDFUtils::GENIE::GetBeta(const ROOT::VecOps::RVec<TLorentzVector>& P4){
    ROOT::VecOps::RVec<double> Beta = {};
    for(const auto& p4 : P4){
        Beta.push_back(p4.Beta());
    }
    return Beta;
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

double GeoUtils::IntersectWithInfCylinder(const TVector3& vtx, const TVector3& direction, double diameter) {
    /*
        Calculate the intersection btw a line that starts from vtx and travels
        in the direction neutron_direction (versor) and an infinite cylinder.
        Take the intersection along the direction of travel, discard the other one.
    */
    TVector3 start_point(vtx.X() - GeoUtils::DRIFT::SAND_CENTER_X, 
                         vtx.Y() - GeoUtils::DRIFT::SAND_CENTER_Y, 
                         vtx.Z() - GeoUtils::DRIFT::SAND_CENTER_Z);
    
    double A = direction.Z() * direction.Z() + direction.Y() * direction.Y();
    double B = 2 * (start_point.Z() * direction.Z() + start_point.Y() * direction.Y());
    double C = start_point.Z() * start_point.Z() + start_point.Y() * start_point.Y() - (diameter / 2.0) * (diameter / 2.0);
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

double GeoUtils::IntersectWithTube(std::string units, 
                         const TVector3& vertex, 
                         const TVector3& direction,
                         double tube_diameter,
                         double tube_length) {
    /*
        Given the neutrino vertex position and the neutron momentum
        predict where neutron is expected to release energy.
        !!! THE RETURNED POINT'S UNITS IS mm 
    */
    // convert to mm
    auto vertex_ = vertex;
    if (units == "m") {
        vertex_ = vertex * 1e3;
    } else if (units == "cm") {
        vertex_ = vertex * 10.;
    } else if (units == "mm") {
        vertex_ = vertex;
    } else {
        std::cerr << "Error: Unsupported unit \"" << units << "\". Supported units are m, cm, mm.\n";
        throw std::invalid_argument("Unsupported unit");
    }
    TVector3 direction_versor = direction.Unit();
    double t = GeoUtils::IntersectWithInfCylinder(vertex_, direction_versor, tube_diameter);

    if (fabs(vertex_.X() + t * direction_versor.X() - GeoUtils::DRIFT::SAND_CENTER_X) >= tube_length / 2.0) {
        
        double t1 = (GeoUtils::SAND_CENTER_X + tube_length / 2.0 - vertex_.X()) / direction_versor.X();
        double t2 = (GeoUtils::SAND_CENTER_X - tube_length / 2.0 - vertex_.X()) / direction_versor.X();
      
        if(t1<0){
            t = t2;
        }else{
            t = t1;
        }
    }

    return t;
}

// double RDFUtils::GENIE::GetBeta(double Ptot, double E){
//     return Ptot / E;
// }

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
    double E = sqrt(Ptot * Ptot + neutron_mass * neutron_mass); // [Mev]
    // double gamma = E / neutron_mass;
    // double beta = sqrt(1 - 1/(gamma*gamma));
    double beta = Ptot / E;

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
            .Define("NEUTRON_MASS_MeV",   [](){return GenieUtils::GetMass(2112)*1e3;})
            .Define("MUON_MASS_GeV",      [](){return GenieUtils::GetMass(13);})
            .Define("MUON_MASS_MeV",      [](){return GenieUtils::GetMass(13)*1e3;})
            // offset
            .Define("TIME_CALIBRATION_OFFSET", [](){return -0.90;}) // ns calubrated on earlist cells 
            .Define("X_CALIBRATION_OFFSET",    [](){return +0.;}) // mm
            .Define("Y_CALIBRATION_OFFSET",    [](){return +0.;}) // mm
            // .Define("X_CALIBRATION_OFFSET",    [](){return +0.64;}) // mm
            // .Define("Y_CALIBRATION_OFFSET",    [](){return +0.41;}) // mm
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
             .Define("InteractionVolume_short",     GeoUtils::InteractionVolume_short, {"InteractionVolume"})                                                                                              
             .Define("isInFiducialVolume",          GeoUtils::DRIFT::IsInFiducialVolume, {"InteractionVolume",
                                                                                          "GENIE_UNIT_LEGTH",
                                                                                          "Interaction_vtxX",
                                                                                          "Interaction_vtxY",
                                                                                          "Interaction_vtxZ"})
            //  .Filter("isInFiducialVolume")
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
             .Define("FinalStateHadronicSystemBeta",            RDFUtils::GENIE::GetBeta, {"FinalStateHadronicSystemMomenta"})
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
            //  .Define("ExpectedHadronSystP4", [](const TVector3 P3, const double E){TLorentzVector v = {P3, E}; return v;}, 
            //                                                         {"ExpectedHadronSystP3", "ExpectedHadronSystEnergy"})
             .Define("ExpectedHadronSystBeta",                      "ExpectedHadronSystP3.Mag() / ExpectedHadronSystEnergy")                                                                                                             
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
    // convert to cm
    if(units == "mm"){
        x = x * 0.1;
        y = y * 0.1;
        z = z * 0.1;
    }else if(units == "m"){
        x = x * 100.;
        y = y * 100.;
        z = z * 100.;
    }
    // // convert to mm
    // if(units == "cm"){
    //     x = x * 10.;
    //     y = y * 10.;
    //     z = z * 10.;
    // }else if(units == "m"){
    //     x = x * 1000.;
    //     y = y * 1000.;
    //     z = z * 1000.;
    // }
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
    auto Primaries = v[0].Particles;
    ROOT::VecOps::RVec<TLorentzVector> P4(Primaries.size());
    for(size_t i = 0; i < Primaries.size(); ++i){
        P4[i] = Primaries[i].GetMomentum();
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
    ROOT::VecOps::RVec<double> angles(primariesP4.size());
    for (size_t i = 0; i < primariesP4.size(); i++)
    {
        auto initial_momentum =  primariesP4[i].Vect();
        double argument = initial_momentum.Dot(nuDirection) / ((initial_momentum.Mag())*(nuDirection.Mag()));
        angles[i] = acos(argument);
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

    int j = 0;
    for(auto hit: ev.SegmentDetectors["EMCalSci"]){
        // retreive the track id of the most important particle associated with this hit
        auto hit_primary_id = hit.GetPrimaryId();

        // check if the track id associated to the hit belongs to one of the primaries id
        auto it = std::find(ids.begin(), ids.end(), hit_primary_id);

        if (it != ids.end()){
            // update info of EDepUtils::track_hits that belongs to the primary whose id is hit_primary_id
            auto hit_middle = (hit.GetStop() + hit.GetStart())*0.5;
            map_grouped_hits[hit_primary_id].h_indices.push_back(j);
            map_grouped_hits[hit_primary_id].hit_edep.push_back(hit.GetEnergyDeposit());
            map_grouped_hits[hit_primary_id].hit_LorentzVector.push_back(hit_middle);
        }
        j++;
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

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::NOSPILL::GetDirectionInECAL(const ROOT::VecOps::RVec<TLorentzVector>& P4,
                                                                           const double vtxX,
                                                                           const double vtxY,
                                                                           const double vtxZ,
                                                                           const ROOT::VecOps::RVec<TLorentzVector>& primaries_first_hit){
    /**
     * @brief Calculates the angle between the direction of emission of the primary particle and
     * its direction when it is in ECAL.
     * For each primary particle, if a valid hit in the ECAL exists, the angle is computed and returned.
     * If no hit is found, a default value of RDFUtils::DEFAULT_NO_DATA is returned.
     * @param P4  primaries initial momenta (RVec<TLorentzVector>).
     * @param vtxX, vtxY, vtxZ  Interaction vertex coordinates (in meters).
     * @param primaries_first_hit  First hit positions of primary particles in ECAL (RVec<TLorentzVector>).
     * @return RVec<double>  Angles (in radians) or RDFUtils::DEFAULT_NO_DATA if no hit.
     * @note Ensure all input coordinates are in meters.
    */
     
    ROOT::VecOps::RVec<double> directions(primaries_first_hit.size());
    TVector3 vtx_mm = {vtxX * 1e3, vtxY * 1e3, vtxZ * 1e3}; // [mm]                                                             
    for (size_t i = 0; i < primaries_first_hit.size(); i++)
    {
        // if no hit in ECAL defaul is null TLV, return default directopn RDFUtils::DEFAULT_NO_DATA
        if(primaries_first_hit[i].Mag() == 0.){
            directions[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            TVector3 intial_momentum = P4[i].Vect();
            TVector3 first_hit_ = primaries_first_hit[i].Vect(); // [mm]
            TVector3 direction_in_ECAL = first_hit_ - vtx_mm;
            double argument = direction_in_ECAL.Dot(intial_momentum) / ((direction_in_ECAL.Mag())*(intial_momentum.Mag()));
            directions[i] = acos(argument);
        }
    }
    return directions;
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
            vertices.push_back(RDFUtils::DEFAULT_NO_DATA);
        }else{
            vertices.push_back(primary_hits.hit_LorentzVector[0][coordinate]);
        }
    }
    return vertices;
}

ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::GetHitTrajectoryId(const ROOT::VecOps::RVec<int>& pmts_hindex, 
                                                              TG4Event& ev){
    /**
    * @brief Retrieves the primary trajectory IDs for ECAL PMT hits.
    * 
    * Given a vector of PMT hit indices (`pmts_hindex`) and a TG4Event (`ev`), 
    * this function returns the primary trajectory IDs associated with each hit. 
    * If the hit index is invalid or there are no TDC data, it returns a default value (`RDFUtils::DEFAULT_NO_DATA`). 
    * Throws an exception if the hit index exceeds the available ECAL hits or if the vector sizes don't match.
    * 
    * @param pmts_hindex Vector of PMT hit indices.
    * @param ev TG4Event structure with event data.
    * @return ROOT::VecOps::RVec<int> Vector of primary trajectory IDs.
    */                                                    
    ROOT::VecOps::RVec<int> track_ids;
    auto n_hits = ev.SegmentDetectors["EMCalSci"].size();  

    for(const auto& pmt_hindex : pmts_hindex){  // index is the id of the hit that generated the tdc
        if(pmt_hindex == RDFUtils::DEFAULT_NO_DATA){  // pmt hasn't any tdc
            track_ids.push_back(RDFUtils::DEFAULT_NO_DATA);
        }else{
            // convert pmt_hindex to unsigned for safe comparison
            if(static_cast<unsigned>(pmt_hindex) >= n_hits){
                std::cout << "Hit index is out of bounds for ECAL hits.\n";
                throw std::out_of_range("Hit index exceeds available ECAL hits.");
            }

            auto hit_j = ev.SegmentDetectors["EMCalSci"][pmt_hindex];
            track_ids.push_back(hit_j.GetPrimaryId());
        }
    }

    if(pmts_hindex.size() != track_ids.size()){
        std::cout << "Mismatch between the size of PMT TDC h_index and track IDs.\n";
        throw std::runtime_error("1-1 correspondence between PMT TDC h_index and track IDs not found.");
    }

    return track_ids;
}

template<int PDG>
ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::CheckTrajId(const ROOT::VecOps::RVec<int>& track_ids,
                                                       TG4Event& ev){
    /**
    * @brief Retrieves 1, 0 if the track id belongs to a primary particle.
    * 
    * Given a vector of PMT track indices and a TG4Event (`ev`), 
    * this function check if the id is the one of the primary particle
    * 
    * @param track_ids Vector of track ids.
    * @param ev TG4Event structure with event data.
    * @return ROOT::VecOps::RVec<int> Vector of 1 or 0 (yes/no).
    */
   ROOT::VecOps::RVec<int> IsPrimaryTraj(track_ids.size());
   auto Trajectories = ev.Trajectories;
   for (size_t i = 0; i < track_ids.size(); i++)
   {
        if(track_ids[i] == RDFUtils::DEFAULT_NO_DATA){
            IsPrimaryTraj[i] = 0;
        }else{
            auto this_trj = Trajectories[track_ids[i]];
            if( (this_trj.GetParentId() == -1) & (this_trj.GetPDGCode() == PDG)){
                IsPrimaryTraj[i] = 1;
            }else{
                IsPrimaryTraj[i] = 0;
            }
        }
   }
   return IsPrimaryTraj;
}

// ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::GetNeutronEKin(const ROOT::VecOps::RVec<int>& track_ids,
//                                                            const ROOT::VecOps::RVec<int>& is_primary_neutron, 
//                                                            TG4Event& ev,
//                                                            const double neutron_mass){
//     auto Trajectories = ev.Trajectories;
//     ROOT::VecOps::RVec<double> E_kin(track_ids.size());
//     for (size_t i = 0; i < track_ids.size(); i++)
//     {
//         if(is_primary_neutron[i]==1){
//             auto this_trj = Trajectories[track_ids[i]];
//             auto this_trj_p = this_trj.GetInitialMomentum();
//             auto this_trj_K = sqrt(this_trj_p * this_trj_p - neutron_mass * neutron_mass) - neutron_mass;
//             E_kin[i] = this_trj_K;
//         }else{
//             E_kin[i] == RDFUtils::DEFAULT_NO_DATA;
//         }
//     }
//     return E_kin;                                                        
// }

ROOT::VecOps::RVec<TG4Trajectory> RDFUtils::EDEPSIM::GetTrajectories(const ROOT::VecOps::RVec<int>& track_ids,
                                                                  TG4Event& ev){
    /**
    * @brief Retrieves all the TG4Trajectory that have fired the cells
    *        of the event.
    * 
    * Given a vector of track indices and a TG4Event (`ev`), 
    * this function returns the unique trajectories associated
    * with the ids.
    * 
    * @param track_ids Vector of track ids.
    * @param ev TG4Event structure with event data.
    * @return ROOT::VecOps::RVec<TG4Trajectory> Vector of TG4Trajectory.
    */
    // Remove duplicate track ids
    std::set<int> unique_values(track_ids.begin(), track_ids.end());
    ROOT::VecOps::RVec<int> track_ids_unique(unique_values.begin(), unique_values.end());
    
    // Create a default trajectory for cases where no valid track is found
    TG4Trajectory Dummy;
    ROOT::VecOps::RVec<TG4Trajectory> trajectories;

    for (const auto& track_id : track_ids_unique) {
        if (track_id == RDFUtils::DEFAULT_NO_DATA || track_id < 0) {
            // Add the dummy trajectory for invalid track ids (-999 or negative)
            trajectories.push_back(Dummy);
        } else if (track_id < static_cast<int>(ev.Trajectories.size())) {
            // Add valid trajectories based on track_id
            trajectories.push_back(ev.Trajectories[track_id]);
        }
    }

    return trajectories;
}

ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::ExpandPDG(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    /**
     * @brief For each TG4Trajectory point append its pdg to the output vector
     * Useful later to be analyzed with Pandas DF
    */
   ROOT::VecOps::RVec<int> expanded_id;
   for (const auto& trajectory : trajectories)
   {
        if(trajectory.GetTrackId() == -1){ // is dummy trj
            expanded_id.push_back(RDFUtils::DEFAULT_NO_DATA);
        }else{
            for (size_t i = 0; i < trajectory.Points.size(); i++)
            {
                expanded_id.push_back(trajectory.GetPDGCode());
            }
        }
   }
   return expanded_id;
}

ROOT::VecOps::RVec<int> RDFUtils::EDEPSIM::ExpandIndex(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories){
    /**
     * @brief For each TG4Trajectory point append its id to the output vector
     * Useful later to be analyzed with Pandas DF
    */
   ROOT::VecOps::RVec<int> expanded_id;
   for (const auto& trajectory : trajectories)
   {
        if(trajectory.GetTrackId() == -1){ // is dummy trj
            expanded_id.push_back(RDFUtils::DEFAULT_NO_DATA);
        }else{
            for (size_t i = 0; i < trajectory.Points.size(); i++)
            {
                expanded_id.push_back(trajectory.GetTrackId());
            }
        }
   }
   return expanded_id;
}

ROOT::VecOps::RVec<TG4TrajectoryPoint> RDFUtils::EDEPSIM::GetTrjPointsAll(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories) {
    /**
    * @brief Collects trajectory points from a vector of TG4Trajectory objects.
    *
    * Iterates through each trajectory, adding a dummy point if the TrackId is -1,
    * or appending all valid points to the output vector.
    *
    * @param trajectories A vector of TG4Trajectory objects.
    * @return A vector of TG4TrajectoryPoint objects.
    */
    ROOT::VecOps::RVec<TG4TrajectoryPoint> points;
    TG4TrajectoryPoint dummy;

    for (const auto& trajectory : trajectories) {
        if (trajectory.GetTrackId() == -1) { // Dummy trajectory check
            points.push_back(dummy); // Add dummy point
        } else {
            for (const auto& point : trajectory.Points) {
                points.push_back(point); // Append valid points
            }
        }
    }
    return points;
}

template<int component>
ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::GetPointComponent(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points){
    
    ROOT::VecOps::RVec<double> v(points.size());
    for (size_t i = 0; i < points.size(); i++)
    {
        v[i] = points[i].GetPosition()[component];
    }
    return v;
}

template<int component>
ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::GetPointMomentum(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points){
    
    ROOT::VecOps::RVec<double> v(points.size());
    for (size_t i = 0; i < points.size(); i++)
    {
        v[i] = points[i].GetMomentum()[component];
    }
    return v;
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::GetPointProcess(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points){
    
    ROOT::VecOps::RVec<double> v(points.size());
    for (size_t i = 0; i < points.size(); i++)
    {
        v[i] = points[i].GetProcess();
    }
    return v;
}

// ROOT::VecOps::RVec<TLorentzVector> RDFUtils::EDEPSIM::ExpectedNeutronTjr(TVector3 vertex, 
//                                                                          const TVector3& momentum_vector, 
//                                                                          const ROOT::VecOps::RVec<double>& v){
//     // !! vertex [mm] e momentum_vector [MeV], length of v gives nof of points!!
//     ROOT::VecOps::RVec<TG4TrajectoryPoint> points(v.size());
//     auto direction = momentum_vector.Unit();
//     for (size_t i = 0; i < v.size(); i++)
//     {
//         x = vertex.X() +  ;
//     }
    
    
// }


ROOT::VecOps::RVec<TLorentzVector> RDFUtils::EDEPSIM::GetHitFromIndex(const ROOT::VecOps::RVec<int>& h_index, TG4Event& ev){
    /*
        Retreive the hit in ECAL cooresponding to index j
    */
    ROOT::VecOps::RVec<TLorentzVector> hits;
    for(const auto& j : h_index){
        if(j == RDFUtils::DEFAULT_NO_DATA){ // incomplete cell, no tdc recorded
            hits.push_back({0.,0.,0.,0.}); // dummy TLVector
        }else{
            auto hit = ev.SegmentDetectors["EMCalSci"][j];
            auto middle_hit = (hit.GetStart() + hit.GetStop())*0.5;
            hits.push_back(middle_hit);
        }
    }
    return hits;
}

ROOT::VecOps::RVec<double> RDFUtils::EDEPSIM::GetHitEdepFromIndex(const ROOT::VecOps::RVec<int>& h_index, TG4Event& ev){
    /*
        Retreive the hit in ECAL cooresponding to index j
    */
    ROOT::VecOps::RVec<double> hit_e;
    for(const auto& j : h_index){
        if(j == RDFUtils::DEFAULT_NO_DATA){ // incomplete cell, no tdc recorded
            hit_e.push_back(RDFUtils::DEFAULT_NO_DATA); // dummy TLVector
        }else{
            auto hit = ev.SegmentDetectors["EMCalSci"][j];
            hit_e.push_back(hit.GetEnergyDeposit());
        }
    }
    return hit_e;
}

ROOT::RDF::RNode RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df){
    return df
             .Define("FileName",                RDFUtils::EDEPSIM::NOSPILL::FileName, {"Primaries"})
             .Define("NofEvents",               RDFUtils::EDEPSIM::SPILL::NofEventsPerSpill, {"Primaries"})
             .Define("NofPrimaries",            RDFUtils::EDEPSIM::NOSPILL::NofPrimaries, {"Primaries"})
             .Define("PrimariesPDG",            RDFUtils::EDEPSIM::NOSPILL::GetPrimariesPDG, {"Primaries"})
             .Define("PrimariesName",           RDFUtils::EDEPSIM::NOSPILL::PDG2Name, {"PrimariesPDG"})
             .Define("PrimariesTrackId",        RDFUtils::EDEPSIM::NOSPILL::GetPrimariesTrackId, {"Primaries"})
             .Define("PrimariesP4",             RDFUtils::EDEPSIM::NOSPILL::GetPrimariesP4, {"Primaries"})
             .Define("PrimariesBeta",           RDFUtils::GENIE::GetBeta, {"PrimariesP4"})
            //  .Define("Primaries_beta",          RDFUtils::GENIE::GetBeta, {"PrimariesP4"})
             .Define("Primaries_vtxX",          [](const ROOT::VecOps::RVec<TLorentzVector>& vertices){return vertices[0].X();},{"Primaries.Position"})
             .Define("Primaries_vtxY",          [](const ROOT::VecOps::RVec<TLorentzVector>& vertices){return vertices[0].Y();},{"Primaries.Position"})
             .Define("Primaries_vtxZ",          [](const ROOT::VecOps::RVec<TLorentzVector>& vertices){return vertices[0].Z();},{"Primaries.Position"})
             .Define("Primaries_vtxT",          [](const ROOT::VecOps::RVec<TLorentzVector>& vertices){return vertices[0].T();},{"Primaries.Position"})
              // ECAL track_hits
             .Define("PrimariesHitsECAL",       RDFUtils::EDEPSIM::NOSPILL::GetPrimariesHits, {"Event", "PrimariesTrackId", "PrimariesPDG"})
             .Define("PrimariesHitsX",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<0>, {"PrimariesHitsECAL"})
             .Define("PrimariesHitsY",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<1>, {"PrimariesHitsECAL"})
             .Define("PrimariesHitsZ",          RDFUtils::EDEPSIM::NOSPILL::GetECALHitPos<2>, {"PrimariesHitsECAL"})
             // ECAL hits
             .Define("PrimariesEmissionAngle",  RDFUtils::EDEPSIM::NOSPILL::GetEmissionAngle, {"NuDirection", "PrimariesP4"})
             .Define("PrimariesFirstHitECAL",   RDFUtils::EDEPSIM::NOSPILL::GetPrimariesFirstHitECAL, {"PrimariesHitsECAL"})
            //  .Define("IsECALHitMissing",        RDFUtils::IsNULLTLV, {"PrimariesFirstHitECAL"})
             .Define("PrimariesEDepECAL",       RDFUtils::EDEPSIM::NOSPILL::GetPrimariesEDepECAL, {"PrimariesHitsECAL"})
             .Define("DeviationAngle",RDFUtils::EDEPSIM::NOSPILL::GetDirectionInECAL,  {"PrimariesP4",
                                                                                                  "Interaction_vtxX",
                                                                                                  "Interaction_vtxY",
                                                                                                  "Interaction_vtxZ",
                                                                                                  "PrimariesFirstHitECAL"})
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
    ROOT::VecOps::RVec<double> TDCs(cells.size());
    for (size_t i = 0; i < cells.size(); i++){
        if(side==1){ // cell side 1
            if(cells[i].ps1.size()!=0){
                TDCs[i] = cells[i].ps1[0].tdc;
            }else{
                TDCs[i] = 0.;
            }
        }else{ // cell side 2
            if(cells[i].ps2.size()!=0){
                TDCs[i] = cells[i].ps2[0].tdc;
            }else{
                TDCs[i] = 0.;
            }
        }
    }
    return TDCs;
}

template<int side>
ROOT::VecOps::RVec<double> RDFUtils::DIGIT::FiredECALGetADC(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<double> ADCs(cells.size());
    for (size_t i = 0; i < cells.size(); i++){
        if(side==1){ // cell side 1
            if(cells[i].ps1.size()!=0){
                ADCs[i] = cells[i].ps1[0].adc;
            }else{
                ADCs[i] = 0.;
            }
        }else{ // cell side 2
            if(cells[i].ps2.size()!=0){
                ADCs[i] = cells[i].ps2[0].adc;
            }else{
                ADCs[i] = 0.;
            }
        }
    }
    return ADCs;
}

int RDFUtils::DIGIT::Get_TDC_hindex(const dg_ps& photo_sensor){
    for (const auto& pe : photo_sensor.photo_el) {
        if(pe.time == photo_sensor.tdc){
            return pe.h_index;
        }
    }
    return RDFUtils::DEFAULT_NO_DATA;
}   

template<int side>
ROOT::VecOps::RVec<int> RDFUtils::DIGIT::GetHindexOfTDC(const ROOT::VecOps::RVec<dg_cell>& cells){
    /*
        Get hit index that corresponding to the observed TDC
    */
    ROOT::VecOps::RVec<int> h_indices;
    for(const auto& cell : cells){
        if(side == 1){ // cell side 1
            if(cell.ps1.size()!=0){
                h_indices.push_back(RDFUtils::DIGIT::Get_TDC_hindex(cell.ps1[0]));
            }else{// if cell has no TDC, set defaul index
                h_indices.push_back(RDFUtils::DEFAULT_NO_DATA);
            }
        }else if(side == 2){ // cell side 1
            if(cell.ps2.size()!=0){
                h_indices.push_back(RDFUtils::DIGIT::Get_TDC_hindex(cell.ps2[0]));
            }else{
                h_indices.push_back(RDFUtils::DEFAULT_NO_DATA);
            }
        }else{
            std::cout << "calorimeter has not side " << side << "\n";
            throw "";
        }
    }
    return h_indices;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::TfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                     const double calibration_const){
    ROOT::VecOps::RVec<double> t_hit_reco;
    for (const auto& cell : cells)
    {
        int broken_side = 0;
        double tdc1 = 0., tdc2 = 0.;
        double t_reco = RDFUtils::DEFAULT_NO_DATA;

        if(!GeoUtils::ECAL::isCellComplete(cell, broken_side)){ // broken cell
            // t_reco = RDFUtils::DEFAULT_NO_DATA;
        } else{ // complete cell
            tdc1 = cell.ps1.at(0).tdc;
            tdc2 = cell.ps2.at(0).tdc;
            t_reco = GeoUtils::ECAL::TfromTDC(tdc1, tdc2, cell.l) + calibration_const;
        }
        t_hit_reco.push_back(t_reco);
    }
    return t_hit_reco;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::EfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                     const double calibration_const,
                                                     const ROOT::VecOps::RVec<double>& e_true){
    ROOT::VecOps::RVec<double> e_reco;
    int i=0;
    for (const auto& cell : cells)
    {
        int broken_side = 0;
        double e = RDFUtils::DEFAULT_NO_DATA;

        if(!GeoUtils::ECAL::isCellComplete(cell, broken_side)){ // broken cell
            // t_reco = RDFUtils::DEFAULT_NO_DATA;
        } else{ // complete cell
            double d1 = 0.5 * cell.l + (GeoUtils::ECAL::XfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc) + calibration_const);
            double d2 = 0.5 * cell.l - (GeoUtils::ECAL::XfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc) + calibration_const);
            e = GeoUtils::ECAL::EfromADC(cell.ps1.at(0).adc, cell.ps2.at(0).adc, d1, d2, cell.lay);
            // std::cout << "e_reco :" << e << ", e_true :" << e_true[i] << "\n";
        }
        e_reco.push_back(e);
        i++;
    }
    return e_reco;
}

TVector3 RDFUtils::DIGIT::WeightedCell(const ROOT::VecOps::RVec<TVector3> x_reco,
                                             const ROOT::VecOps::RVec<double> e_reco,
                                             const ROOT::VecOps::RVec<int> fired_by_primary){
    double x_w = 0.;                                        
    double y_w = 0.;                                        
    double z_w = 0.;                                        
    double etot = 0.;

    for (size_t i = 0; i < fired_by_primary.size(); i++)
    {
        if(fired_by_primary[i] == 1){
            if((x_reco[i].Z() != RDFUtils::DEFAULT_NO_DATA) && 
               (e_reco[i] != RDFUtils::DEFAULT_NO_DATA) // redoudant but safe
               ){
                x_w += x_reco[i].X() * e_reco[i];
                y_w += x_reco[i].Y() * e_reco[i];
                z_w += x_reco[i].Z() * e_reco[i];
                etot += e_reco[i];
            }
        }
    }
    return {x_w/etot, y_w/etot, z_w/etot};                                        
}

ROOT::VecOps::RVec<TVector3> RDFUtils::DIGIT::XfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                       const double calibration_x,
                                                       const double calibration_y){
    ROOT::VecOps::RVec<TVector3> reco_hits;
    for (const auto& cell : cells)
    {
        int broken_side = 0;
        double tdc1 = 0., tdc2 = 0.;
        TVector3 hit_reco = {RDFUtils::DEFAULT_NO_DATA,
                             RDFUtils::DEFAULT_NO_DATA,
                             RDFUtils::DEFAULT_NO_DATA};

        if(!GeoUtils::ECAL::isCellComplete(cell, broken_side)){ // broken cell
            // what we do in case of broke cells ?
        } else{ // complete cell
            tdc1 = cell.ps1.at(0).tdc;
            tdc2 = cell.ps2.at(0).tdc;
        }
        // if(cell.id < 25000){ // Barrel
        if(cell.mod < 30){ // Barrel
            hit_reco = {cell.x - GeoUtils::ECAL::XfromTDC(tdc1, tdc2) + calibration_x, cell.y, cell.z};
        }else{ // Endcap
            hit_reco = {cell.x, cell.y - GeoUtils::ECAL::XfromTDC(tdc1, tdc2) + calibration_y, cell.z};
        }
        reco_hits.push_back(hit_reco);
    }
    return reco_hits;
}

double RDFUtils::DIGIT::WeightedCellTime(const ROOT::VecOps::RVec<double>& t_hit_reco,
                                         const ROOT::VecOps::RVec<double>& e_reco,
                                         const ROOT::VecOps::RVec<int>& isCompatible){
    double t_w = 0.;
    double e_tot = 0.;
    for (size_t i = 0; i < t_hit_reco.size(); i++)
    {
        if((t_hit_reco[i] != RDFUtils::DEFAULT_NO_DATA) | (isCompatible[i] == 1)){ // skip broken cells
            t_w += t_hit_reco[i] * e_reco[i];
            e_tot += e_reco[i];
        }
    }
    return t_w / e_tot;                                                
}

ROOT::VecOps::RVec<TVector3> RDFUtils::DIGIT::GetExpectedHitPosition(TVector3 vertex,
                                                                    const TVector3& momentum_vector,
                                                                    const ROOT::VecOps::RVec<dg_cell>& cells) {
    
    /**
     * @brief Computes expected hit positions in calorimeter cells.
     *
     * Given the particle's vertex and momentum, calculates expected hit positions in each calorimeter cell.
     * Chack that the intersection happpen inside the ECAL
     *
     * @param vertex The origin point of the particle (TVector3).
     * @param momentum_vector The momentum direction of the particle (TVector3).
     * @param cells A vector of calorimeter cells, each with center coordinates.
     * 
     * @return A vector of expected hit positions, or {-999, -999, -999} for cells where the hit is not valid.
     */                                                                        
    ROOT::VecOps::RVec<TVector3> Expected_hit_position(cells.size());

    // Calculate the minimum t required to intersect the calorimeter (ECAL)
    TVector3 neutron_direction = momentum_vector.Unit();
    
    double min_t = GeoUtils::IntersectWithTube("mm",
                                     vertex, 
                                     neutron_direction, 
                                     GeoUtils::SAND_INNER_VOL_DIAMETER, 
                                     GeoUtils::SAND_INNER_VOL_X_LENGTH
                                     ); // ecal entry point
    
    double max_t = GeoUtils::IntersectWithTube("mm",
                                     vertex, 
                                     neutron_direction, 
                                     GeoUtils::SAND_INNER_VOL_DIAMETER + 2 * GeoUtils::ECAL::thickness, 
                                     GeoUtils::SAND_INNER_VOL_X_LENGTH + 2 * GeoUtils::ECAL::thickness
                                     ); // ecal exit point

    for (size_t i = 0; i < cells.size(); i++) {
        // Calculate t for each cell based on the z-coordinate
        double t = (cells[i].z - vertex.Z()) / neutron_direction.Z();
        
        if ( t > min_t && t < max_t) { // point is in ECAL range
            Expected_hit_position[i] = vertex + t * neutron_direction;
        } else {
            // If t is too small, return a placeholder invalid position
            Expected_hit_position[i] = TVector3(RDFUtils::DEFAULT_NO_DATA, RDFUtils::DEFAULT_NO_DATA, RDFUtils::DEFAULT_NO_DATA);
        }
    }
    
    return Expected_hit_position;
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

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetFlightLength(const TVector3& vtx, 
                                                            ROOT::VecOps::RVec<TVector3>& hits){
    ROOT::VecOps::RVec<double> Flight_lenghts(hits.size());
    for (size_t i = 0; i < hits.size(); i++)
    {
        Flight_lenghts[i] = (hits[i] - vtx).Mag();
    }
    return Flight_lenghts;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetBeta(const ROOT::VecOps::RVec<double>& flight_length,
                                                   const ROOT::VecOps::RVec<double>& t_hit_reco,
                                                   double vertex_time){
    
    ROOT::VecOps::RVec<double> beta_reco(t_hit_reco.size());
    double c = 299.792; // mm / ns
    for (size_t i = 0; i < t_hit_reco.size(); i++)
    {
        if (t_hit_reco[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            beta_reco[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            beta_reco[i] = (flight_length[i]) / ((t_hit_reco[i] - vertex_time) * c);
        }
    }
    return beta_reco;                                                   
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetPFromBeta(const ROOT::VecOps::RVec<double>& beta,
                                                         const double neutron_mass){
    ROOT::VecOps::RVec<double> ptot_reco(beta.size());
    for (size_t i = 0; i < beta.size(); i++)
    {
        if (beta[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            ptot_reco[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            ptot_reco[i] = neutron_mass * beta[i] / sqrt(1 - beta[i]*beta[i]);
        }
    }
    return ptot_reco;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetKinEFromBeta(const ROOT::VecOps::RVec<double>& beta,
                                                         const double neutron_mass){
    ROOT::VecOps::RVec<double> K(beta.size());
    for (size_t i = 0; i < beta.size(); i++)
    {
        if (beta[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            K[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            K[i] = neutron_mass * ( 1 / sqrt(1 - beta[i]*beta[i]) - 1);
        }
    }
    return K;
}

ROOT::VecOps::RVec<TVector3> RDFUtils::DIGIT::GetRecoP3(const ROOT::VecOps::RVec<double>& ptot_reco,
                                                        const TVector3& expected_direction){
    ROOT::VecOps::RVec<TVector3> p_reco(ptot_reco.size());
    
    TVector3 dummy = {RDFUtils::DEFAULT_NO_DATA, 
                      RDFUtils::DEFAULT_NO_DATA, 
                      RDFUtils::DEFAULT_NO_DATA};
    
    for (size_t i = 0; i < ptot_reco.size(); i++)
    {
        if (ptot_reco[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            p_reco[i] = dummy;
        }else{
            p_reco[i] = ptot_reco[i] * expected_direction;
        }
    }
    return p_reco;                                     
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetMissingPT(const TLorentzVector& lepton_P4,
                                                         const ROOT::VecOps::RVec<TVector3>& neutron_reco_p3){
    
    ROOT::VecOps::RVec<double> missing_pt(neutron_reco_p3.size());

    TVector3 lepton_p3 = lepton_P4.Vect();

    for (size_t i = 0; i < neutron_reco_p3.size(); i++)
    {
        if (neutron_reco_p3[i].Z() == RDFUtils::DEFAULT_NO_DATA)
        {
            missing_pt[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            missing_pt[i] = (lepton_p3 - neutron_reco_p3[i]).Mag(); 
        }
    }
    return missing_pt;                                                
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::GetTOF(const ROOT::VecOps::RVec<double>& flight_length, 
                                                   const double beta){
    /*
        calculate time of flight of particle with beta and a given flight_length
    */
    double c = 299.792; // [mm/ns]
    ROOT::VecOps::RVec<double> TOF(flight_length.size());
    for (size_t i = 0; i < flight_length.size(); i++)
    {
        TOF[i] = flight_length[i] / (c * beta);
    }
    return TOF;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::TimeResiduals(const ROOT::VecOps::RVec<double>& expected_t_hit,
                                                          const ROOT::VecOps::RVec<double>& reconstruced_t_hit){
    /**
     * @brief Get time residuals between expected time and reconstructed time.
     * if cell is broken, return default data
     *  */
    if(expected_t_hit.size() != reconstruced_t_hit.size()){
        throw std::runtime_error("expected_t_hit and reconstruced_t_hit must have the same length");
    }
    
    ROOT::VecOps::RVec<double> time_residuals(expected_t_hit.size());
    
    for (size_t i = 0; i < expected_t_hit.size(); i++)
    {
        if (expected_t_hit[i] == RDFUtils::DEFAULT_NO_DATA)
        {
            time_residuals[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            time_residuals[i] = expected_t_hit[i] - reconstruced_t_hit[i];
        }
        
    }
    return time_residuals;
}

ROOT::VecOps::RVec<double> RDFUtils::DIGIT::SpaceResiduals(const ROOT::VecOps::RVec<TVector3>& expected_x_hit,
                                                           const ROOT::VecOps::RVec<TVector3>& reconstructed_x_hit){
    /**
     * @brief Get space residuals between expected and reconstructed hit position.
     * if cell is broken, return default data
     *  */
    if(expected_x_hit.size() != reconstructed_x_hit.size()){
        throw std::runtime_error("expected_x_hit and reconstruced_t_hit must have the same length");
    }
    
    ROOT::VecOps::RVec<double> space_residuals(expected_x_hit.size());
    
    for (size_t i = 0; i < expected_x_hit.size(); i++)
    {
        if ((expected_x_hit[i].Z() == RDFUtils::DEFAULT_NO_DATA) | (reconstructed_x_hit[i].Z() == RDFUtils::DEFAULT_NO_DATA))
        {
            space_residuals[i] = RDFUtils::DEFAULT_NO_DATA;
        }else{
            space_residuals[i] = (expected_x_hit[i] - reconstructed_x_hit[i]).Mag();
        }
    }
    return space_residuals;
}

ROOT::VecOps::RVec<int> RDFUtils::DIGIT::IsCompatible(const ROOT::VecOps::RVec<double>& space_residuals,
                                                      const ROOT::VecOps::RVec<double>& time_residuals){
    /*
      Input: time/space residuals between the expected and the reconstructed hit in the cell
      and look for those cell compatible in the space and in the time
    */
    ROOT::VecOps::RVec<int> isCompatible(space_residuals.size());
    
    for (size_t i = 0; i < space_residuals.size(); i++){
        if((space_residuals[i] == RDFUtils::DEFAULT_NO_DATA)){
            isCompatible[i] = 0;
        }else{
            // 10 cm precision in the determination of neutron position and 0.8 ns in time
            if((fabs(space_residuals[i]) <= 100) & (fabs(time_residuals[i]) <= 0.8)){
                isCompatible[i] = 1;
            }else{
                isCompatible[i] = 0;
            }
        }
    }
    return isCompatible;
}

ROOT::VecOps::RVec<int> RDFUtils::DIGIT::IsCandidateCell(const ROOT::VecOps::RVec<int>& isCompatible,
                                                         const ROOT::VecOps::RVec<double>& variable) {
    /**
    * @brief .
    * 
    * @param isCompatible ROOT::VecOps::RVec<int> - whether each cell is spacetime-compatible (1 or 0).
    * @param variable ROOT::VecOps::RVec<double>
    * 
    * @return ROOT::VecOps::RVec<int> - A vector where each element is 1 if the corresponding cell is a candidate, or 0 otherwise.
    */
    
    // Initialize the output vector with zeros
    ROOT::VecOps::RVec<int> IsCandidate(variable.size(), 0);

    // Filter out invalid values (i.e., DEFAULT_NO_DATA) from the input vector
    ROOT::VecOps::RVec<double> variable_valid = variable[variable != RDFUtils::DEFAULT_NO_DATA];

    // Check if there are any valid data points; if not, return the zero-initialized output
    if (variable_valid.empty()) {
        return IsCandidate;
    }

    // Find the minimum value among the valid hits
    double min_variable_valid = ROOT::VecOps::Min(variable_valid);

    // Loop through each element to determine if it's a candidate
    for (size_t i = 0; i < variable.size(); i++) {
        // Check if the cell is space-compatible and has a valid hit time
        if (isCompatible[i] == 1 && variable[i] != RDFUtils::DEFAULT_NO_DATA) {
            // Mark as candidate if the hit time matches the minimum valid hit time
            if (variable[i] == min_variable_valid) {
                IsCandidate[i] = 1;
            }
        }
    }
    // Return the final candidate vector
    return IsCandidate;
}

int RDFUtils::DIGIT::isCandidateSignal(ROOT::VecOps::RVec<int>& cell_has_coincidence){
    /*
        if event has candidate cell (a cell with a coincidence) is candidate signal 
    */
    for (size_t i = 0; i < cell_has_coincidence.size(); i++)
    {
        if(cell_has_coincidence[i] == 1){
            return 1;
        }
    }
    return 0;
}

ROOT::RDF::RNode RDFUtils::DIGIT::AddColumnsFromDigit(ROOT::RDF::RNode& df){
    return df
            .Define("NofEventFiredModules",         RDFUtils::DIGIT::NofFiredECALMods ,             {"dg_cell.mod"})
            .Define("EventFiredModules",            RDFUtils::DIGIT::FiredECALMods ,                {"dg_cell.mod"})
            .Define("Fired_Cells_mod",              "dg_cell.mod")
            .Define("Fired_Cells_id",               "dg_cell.id")
            .Define("Fired_Cells_x",                "dg_cell.x")
            .Define("Fired_Cells_y",                "dg_cell.y")
            .Define("Fired_Cells_z",                "dg_cell.z")
            .Define("isCellComplete",               RDFUtils::DIGIT::IsCellComplete,                {"dg_cell"})
            .Define("Fired_Cells_tdc1",             RDFUtils::DIGIT::FiredECALGetTDC<1>,            {"dg_cell"})
            .Define("Fired_Cells_adc1",             RDFUtils::DIGIT::FiredECALGetADC<1>,            {"dg_cell"})
            .Define("Fired_Cell_tdc1_hindex",       RDFUtils::DIGIT::GetHindexOfTDC<1>,             {"dg_cell"})
            .Define("trackid_producing_tdc1",       RDFUtils::EDEPSIM::GetHitTrajectoryId,          {"Fired_Cell_tdc1_hindex", "Event"})
            .Define("Fired_by_primary_neutron",     RDFUtils::EDEPSIM::CheckTrajId<2112>,           {"trackid_producing_tdc1", "Event"})
            .Define("Fired_by_primary_antimu",      RDFUtils::EDEPSIM::CheckTrajId<-13>,            {"trackid_producing_tdc1", "Event"})
            .Define("ECALactive_trajectories",      RDFUtils::EDEPSIM::GetTrajectories,             {"trackid_producing_tdc1", "Event"})
            .Define("Fired_Cells_tdc2",             RDFUtils::DIGIT::FiredECALGetTDC<2>,            {"dg_cell"})
            .Define("Fired_Cells_adc2",             RDFUtils::DIGIT::FiredECALGetADC<2>,            {"dg_cell"})
            .Define("Fired_Cell_tdc2_hindex",       RDFUtils::DIGIT::GetHindexOfTDC<2>,             {"dg_cell"})
            .Define("Fired_Cell_true_hit2",         RDFUtils::EDEPSIM::GetHitFromIndex,             {"Fired_Cell_tdc2_hindex", "Event"})
            .Define("trackid_producing_tdc2",       RDFUtils::EDEPSIM::GetHitTrajectoryId,          {"Fired_Cell_tdc2_hindex", "Event"})
            /*
                True
            */
            .Define("Primaries_vtx", [](const double vtx, const double vty, const double vtz){
                TVector3 v = {vtx, vty, vtz}; return v;},                                           {"Primaries_vtxX","Primaries_vtxY","Primaries_vtxZ"})
            .Define("Fired_Cell_true_hit1",         RDFUtils::EDEPSIM::GetHitFromIndex,             {"Fired_Cell_tdc1_hindex", "Event"})
            .Define("Fired_Cell_true_Hit_e",        RDFUtils::EDEPSIM::GetHitEdepFromIndex,         {"Fired_Cell_tdc1_hindex","Event"})
            .Define("Fired_Cell_true_Hit_x",        RDFUtils::GetComponent<0>,                      {"Fired_Cell_true_hit1"})
            .Define("Fired_Cell_true_Hit_y",        RDFUtils::GetComponent<1>,                      {"Fired_Cell_true_hit1"})
            .Define("Fired_Cell_true_Hit_z",        RDFUtils::GetComponent<2>,                      {"Fired_Cell_true_hit1"})
            .Define("Fired_Cell_true_Hit_t",        RDFUtils::GetComponent<3>,                      {"Fired_Cell_true_hit1"})
            .Define("Fired_Cell_true_HitPos",       RDFUtils::GetVect,                              {"Fired_Cell_true_hit1"})
            .Define("True_FlightLength",            RDFUtils::DIGIT::GetFlightLength,               {"Primaries_vtx", "Fired_Cell_true_HitPos"})
            /*
                Expected neutron hit in cells from MC truth
            */
            .Alias("ExpectedNeutron_Beta",          "ExpectedHadronSystBeta")
            .Define("ExpectedNeutron_HitPosition",  RDFUtils::DIGIT::GetExpectedHitPosition,        {"Primaries_vtx", "ExpectedHadronSystP3", "dg_cell"})
            .Define("ExpectedNeutron_HitPosition_x",RDFUtils::GetComponent3<0>,                     {"ExpectedNeutron_HitPosition"})
            .Define("ExpectedNeutron_HitPosition_y",RDFUtils::GetComponent3<1>,                     {"ExpectedNeutron_HitPosition"})
            .Define("ExpectedNeutron_HitPosition_z",RDFUtils::GetComponent3<2>,                     {"ExpectedNeutron_HitPosition"})
            .Define("ExpectedNeutron_FlightLength", RDFUtils::DIGIT::GetFlightLength,               {"Primaries_vtx", "ExpectedNeutron_HitPosition"})
            .Define("ExpectedNeutron_TOF",          RDFUtils::DIGIT::GetTOF,                        {"ExpectedNeutron_FlightLength", "ExpectedNeutron_Beta"})
            .Define("ExpectedNeutron_HitTime",      "ExpectedNeutron_TOF + Primaries_vtxT")
            /*
                Recontructed hit from cells TDC
            */
            .Define("Reconstructed_HitPosition",    RDFUtils::DIGIT::XfromTDC,                      {"dg_cell", "X_CALIBRATION_OFFSET", "Y_CALIBRATION_OFFSET"})
            .Define("Reconstructed_HitPosition_x",  RDFUtils::GetComponent3<0>,                     {"Reconstructed_HitPosition"})
            .Define("Reconstructed_HitPosition_y",  RDFUtils::GetComponent3<1>,                     {"Reconstructed_HitPosition"})
            .Define("Reconstructed_HitPosition_z",  RDFUtils::GetComponent3<2>,                     {"Reconstructed_HitPosition"})
            .Define("Reconstructed_HitTime",        RDFUtils::DIGIT::TfromTDC,                      {"dg_cell", "TIME_CALIBRATION_OFFSET"})
            .Define("Reconstructed_Energy",         RDFUtils::DIGIT::EfromTDC,                      {"dg_cell","X_CALIBRATION_OFFSET","Fired_Cell_true_Hit_e"})                                                                                    
            .Define("IsEarliestCell_neutron",       RDFUtils::FindMinimum2,                         {"Reconstructed_HitTime", "Fired_by_primary_neutron"})
            .Define("IsEarliestCell_antimu",        RDFUtils::FindMinimum2,                         {"Reconstructed_HitTime", "Fired_by_primary_antimu"})
            .Define("Reconstructed_FlightLength",   RDFUtils::DIGIT::GetFlightLength,               {"Primaries_vtx", "Reconstructed_HitPosition"})
        
            // reconstructed neutron
            .Define("Reconstructed_Beta",           RDFUtils::DIGIT::GetBeta,                       {"ExpectedNeutron_FlightLength", "Reconstructed_HitTime", "Primaries_vtxT"})
            .Define("Reconstructed_Ptot_GeV",       RDFUtils::DIGIT::GetPFromBeta,                  {"Reconstructed_Beta", "NEUTRON_MASS_GeV"})
            .Define("Reconstructed_neutronKin_MeV", RDFUtils::DIGIT::GetKinEFromBeta,               {"Reconstructed_Beta", "NEUTRON_MASS_MeV"})
            ;
}

ROOT::VecOps::RVec<dg_cell> RDFUtils::DIGIT::FitlerCells(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                        const ROOT::VecOps::RVec<int>& is_complete,
                                                        const ROOT::VecOps::RVec<int>& is_fired_by_neutron){
    /*
        Given all the fired cells of a given event, retreive only the complete cells fired by signal neutron
    */                                                                
    ROOT::VecOps::RVec<dg_cell> filtered_cells;

    for (size_t i = 0; i < is_fired_by_neutron.size(); i++)
    {
        if(is_fired_by_neutron[i] == 1 && is_complete[i] == 1){
            filtered_cells.push_back(cells[i]);
        }
    }
    return filtered_cells;
}

ROOT::VecOps::RVec<int>  RDFUtils::DIGIT::GetCellMod(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> mods(cells.size());
    for (size_t i = 0; i < cells.size(); i++)
    {
        mods[i] = cells[i].mod;
    }
    return mods;
}

ROOT::VecOps::RVec<int>  RDFUtils::DIGIT::GetCellId(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> ids(cells.size());
    for (size_t i = 0; i < cells.size(); i++)
    {
        ids[i] = cells[i].id;
    }
    return ids;
}

ROOT::VecOps::RVec<double>  RDFUtils::DIGIT::GetCellX(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> X(cells.size());
    for (size_t i = 0; i < cells.size(); i++)
    {
        X[i] = cells[i].x;
    }
    return X;
}

ROOT::VecOps::RVec<double>  RDFUtils::DIGIT::GetCellY(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> Y(cells.size());
    for (size_t i = 0; i < cells.size(); i++)
    {
        Y[i] = cells[i].y;
    }
    return Y;
}

ROOT::VecOps::RVec<double>  RDFUtils::DIGIT::GetCellZ(const ROOT::VecOps::RVec<dg_cell>& cells){
    ROOT::VecOps::RVec<int> Z(cells.size());
    for (size_t i = 0; i < cells.size(); i++)
    {
        Z[i] = cells[i].z;
    }
    return Z;
}

// // add
// ROOT::VecOps::RVec<TVector3> RDFUtils::DIGIT::GetDirection(const TVector3& vertex,
//                                                            ROOT::VecOps::RVec<TVector3>& hit_reco){
//     ROOT::VecOps::RVec<TVector3> directions(hit_reco.size());
//     for (size_t i = 0; i < hit_reco.size(); i++)
//     {
//         if(hit_reco[i].Z() == RDFUtils::DEFAULT_NO_DATA){
//             directions[i] == RDFUtils::DEFAULT_NO_DATA;
//         }else{
//             directions[i] ==  (hit_reco[i] - vertex).Unit();
//         }
//     }
//     return directions                                                    
// }

ROOT::RDF::RNode RDFUtils::DIGIT::GetInfoCellsFromSignal(ROOT::RDF::RNode& df){
    return df
            .Filter("CCQEonHydrogen")
            .Define("filtered_cells",               RDFUtils::DIGIT::FitlerCells,              {"dg_cell", "isCellComplete", "Fired_by_primary_neutron"})
            .Define("neutrons_cells_mod",           RDFUtils::DIGIT::GetCellMod,               {"filtered_cells"})
            .Define("neutrons_cells_id",            RDFUtils::DIGIT::GetCellId,                {"filtered_cells"})
            .Define("neutrons_cells_x",             RDFUtils::DIGIT::GetCellX,                 {"filtered_cells"})
            .Define("neutrons_cells_y",             RDFUtils::DIGIT::GetCellY,                 {"filtered_cells"})
            .Define("neutrons_cells_z",             RDFUtils::DIGIT::GetCellZ,                 {"filtered_cells"})
            /*
                TRUE HIT
            */
            .Define("neutrons_cells_hindex1",       RDFUtils::DIGIT::GetHindexOfTDC<1>,        {"filtered_cells"})                                                         
            .Define("true_hit",                     RDFUtils::EDEPSIM::GetHitFromIndex,        {"neutrons_cells_hindex1", "Event"})
            .Define("true_hit_e",                   RDFUtils::EDEPSIM::GetHitEdepFromIndex,    {"neutrons_cells_hindex1","Event"})
            .Define("true_hit_x",                   RDFUtils::GetComponent<0>,                 {"true_hit"})
            .Define("true_hit_y",                   RDFUtils::GetComponent<1>,                 {"true_hit"})
            .Define("true_hit_z",                   RDFUtils::GetComponent<2>,                 {"true_hit"})
            .Define("true_hit_t",                   RDFUtils::GetComponent<3>,                 {"true_hit"})
            .Define("true_hit_x3",                  RDFUtils::GetVect,                         {"true_hit"})
            .Define("true_hit_FlightLength",        RDFUtils::DIGIT::GetFlightLength,          {"Primaries_vtx", "true_hit_x3"})
            .Define("true_hit_is_earliest",         RDFUtils::FindMinimum,                     {"true_hit_t"})
            /*
                PREDICTED HIT
            */
           .Define("predicted_hit_x3",              RDFUtils::DIGIT::GetExpectedHitPosition,   {"Primaries_vtx", "PredictedNeutron_P3_GeV", "filtered_cells"})
           .Define("predicted_hit_x",               RDFUtils::GetComponent3<0>,                {"predicted_hit_x3"})
           .Define("predicted_hit_y",               RDFUtils::GetComponent3<1>,                {"predicted_hit_x3"})
           .Define("predicted_hit_z",               RDFUtils::GetComponent3<2>,                {"predicted_hit_x3"})
           .Define("predicted_hit_FlightLength",    RDFUtils::DIGIT::GetFlightLength,          {"Primaries_vtx", "predicted_hit_x3"})
           .Define("predicted_hit_TOF",             RDFUtils::DIGIT::GetTOF,                   {"predicted_hit_FlightLength", "PredictedNeutron_Beta"})
           .Define("predicted_hit_t",               "predicted_hit_TOF + Primaries_vtxT")
           /*
                RECONSTRUCTED HIT
           */
           .Define("reconstructed_hit",              RDFUtils::DIGIT::XfromTDC,              {"filtered_cells", "X_CALIBRATION_OFFSET", "Y_CALIBRATION_OFFSET"})
           .Define("reconstructed_hit_x",            RDFUtils::GetComponent3<0>,             {"reconstructed_hit"})
           .Define("reconstructed_hit_y",            RDFUtils::GetComponent3<1>,             {"reconstructed_hit"})
           .Define("reconstructed_hit_z",            RDFUtils::GetComponent3<2>,             {"reconstructed_hit"})
           .Define("reconstructed_hit_t",            RDFUtils::DIGIT::TfromTDC,              {"filtered_cells", "TIME_CALIBRATION_OFFSET"})
           .Define("reconstructed_hit_e",            RDFUtils::DIGIT::EfromTDC,              {"filtered_cells","X_CALIBRATION_OFFSET","true_hit_e"})                                                                                    
           .Define("reconstructed_hit_is_earliest",  RDFUtils::FindMinimum,                  {"reconstructed_hit_t"})
           .Define("reconstructed_hit_FlightLength", RDFUtils::DIGIT::GetFlightLength,       {"Primaries_vtx", "reconstructed_hit"})
           .Define("reconstructed_neutron_beta",     RDFUtils::DIGIT::GetBeta,               {"reconstructed_hit_FlightLength", "reconstructed_hit_t", "Primaries_vtxT"})
           /*
                NEUTRON
           */
            .Alias("neutron_P4_true",                "FinalStateHadronicSystemTotal4Momentum")
            .Define("neutron_beta_true",             "neutron_P4_true.Beta()")
            .Alias("predicted_neutron_beta",         "PredictedNeutron_Beta")
            ;
}

ROOT::RDF::RNode RDFUtils::DIGIT::GetFilteredTrajectories(ROOT::RDF::RNode& df){
    /**
     * @brief Take as input columns of AddColumnsFromDigit and retreive info
     * about trajecotries that have fired ECAL cells
     */
    return df
          .Define("trackid",                      RDFUtils::EDEPSIM::ExpandIndex,               {"ECALactive_trajectories"})
          .Define("pdg",                          RDFUtils::EDEPSIM::ExpandPDG,                 {"ECALactive_trajectories"})
          .Define("Points_all_trajectories",      RDFUtils::EDEPSIM::GetTrjPointsAll,           {"ECALactive_trajectories"})
          .Define("point_x",                      RDFUtils::EDEPSIM::GetPointComponent<0>,      {"Points_all_trajectories"})
          .Define("point_y",                      RDFUtils::EDEPSIM::GetPointComponent<1>,      {"Points_all_trajectories"})
          .Define("point_z",                      RDFUtils::EDEPSIM::GetPointComponent<2>,      {"Points_all_trajectories"})
          .Define("point_px",                     RDFUtils::EDEPSIM::GetPointMomentum<0>,       {"Points_all_trajectories"})
          .Define("point_py",                     RDFUtils::EDEPSIM::GetPointMomentum<1>,       {"Points_all_trajectories"})
          .Define("point_pz",                     RDFUtils::EDEPSIM::GetPointMomentum<2>,       {"Points_all_trajectories"})
          .Define("process",                      RDFUtils::EDEPSIM::GetPointProcess,           {"Points_all_trajectories"})
        ;
}

// ROOT::RDF::RNode RDFUtils::DIGIT::GetFilteredTrajectories_(ROOT::RDF::RJittedFilter& df){
//     /**
//      * @brief Take as input columns of AddColumnsFromDigit and retreive info
//      * about trajecotries that have fired ECAL cells
//      */
//     return df
//           .Define("trackid",                      RDFUtils::EDEPSIM::ExpandIndex,               {"ECALactive_trajectories"})
//           .Define("pdg",                          RDFUtils::EDEPSIM::ExpandPDG,                 {"ECALactive_trajectories"})
//           .Define("Points_all_trajectories",      RDFUtils::EDEPSIM::GetTrjPointsAll,           {"ECALactive_trajectories"})
//           .Define("point_x",                      RDFUtils::EDEPSIM::GetPointComponent<0>,      {"Points_all_trajectories"})
//           .Define("point_y",                      RDFUtils::EDEPSIM::GetPointComponent<1>,      {"Points_all_trajectories"})
//           .Define("point_z",                      RDFUtils::EDEPSIM::GetPointComponent<2>,      {"Points_all_trajectories"})
//           .Define("point_t",                      RDFUtils::EDEPSIM::GetPointComponent<3>,      {"Points_all_trajectories"})
//           .Define("point_px",                     RDFUtils::EDEPSIM::GetPointMomentum<0>,       {"Points_all_trajectories"})
//           .Define("point_py",                     RDFUtils::EDEPSIM::GetPointMomentum<1>,       {"Points_all_trajectories"})
//           .Define("point_pz",                     RDFUtils::EDEPSIM::GetPointMomentum<2>,       {"Points_all_trajectories"})
//           .Define("point_E",                      RDFUtils::EDEPSIM::GetPointMomentum<3>,       {"Points_all_trajectories"})
//           .Define("process",                      RDFUtils::EDEPSIM::GetPointProcess,           {"Points_all_trajectories"})
//         ;
// }

// RDFUtils::RECO_______________________________________________________________________________

// int RDFUtils::RECO::GetNofFiredWires(const ROOT::VecOps::RVec<dg_wire>& wires){
//     // STA FUNZIONE CRASHAAAA
//     return wires.size();
// }
int RDFUtils::RECO::GetNofFiredWires(const ROOT::VecOps::RVec<int>& wires_id){
    int nof_wires = static_cast<int>(wires_id.size());
    return nof_wires;
}

int RDFUtils::RECO::Test(bool b, const MinuitFitInfos& i){
    std::cout << "is good event : " << b << "\n";
    if(b){
        std::cout << "MinValue : " << i.MinValue << "\n";
        return 1;
    }else{
        return 0;
    }
}

int RDFUtils::RECO::isCandidateSignal(std::string interaction_volume, int charge_multiplicity, int nof_cell_candidates){
    /***
     * @brief events is candidate signa if:
     * - the interaction vertex is in one of C3H6 target
     * - the charge multiplicity is 1
     * - the nof cell with possible space-time coincidence is at lest 1
     */
    int is_candidate = 0;
    if((interaction_volume == "C3H6_Target") && (charge_multiplicity == 1) && (nof_cell_candidates > 0)){
        return 1;
    }
    return is_candidate;
}

ROOT::RDF::RNode RDFUtils::RECO::AddColumnsFromDriftReco(ROOT::RDF::RNode& df){
    return df
            // .Filter("KeepThisEvent==1")
            .Define("nof_fired_wires", RDFUtils::RECO::GetNofFiredWires, {"fired_wires.did"})
            .Define("pass_nof_wires_cut", "nof_fired_wires > 69")
            /*
                TRUE ANTIMUONS
            */
            .Alias("Antimuon_dip_true",                             "true_helix.dip_")
            .Alias("Antimuon_curvature_radius_true",                "true_helix.R_")
            .Alias("Antimuon_Phi0_true",                            "reco_helix.Phi0_")
            .Alias("Antimuon_x0_true",                              "true_helix.x0_")
            .Alias("Antimuon_pt_true",                              "pt_true")
            .Alias("Antimuon_p_true",                               "p_true")
            .Define("Antimuon_ptot_true",                           "sqrt(p_true.X()*p_true.X() + p_true.Y()*p_true.Y() + p_true.Z()*p_true.Z())")
            /*
                RECONSTRUCTED ANTIMUON
            */
            .Alias("chi2_fit_zy",                                   "fit_infos_zy.MinValue")
            .Alias("chi2_fit_xz",                                   "fit_infos_xz.MinValue")
            .Alias("Antimuon_curvature_radius_reco",                "reco_helix.R_")
            .Alias("Antimuon_dip_reco",                             "reco_helix.dip_")
            .Alias("Antimuon_Phi0_reco",                            "reco_helix.Phi0_")
            .Alias("Antimuon_x0_reco",                              "reco_helix.x0_")
            .Alias("Antimuon_pt_reco",                              "pt_reco")
            .Alias("Antimuon_p_reco",                               "p_reco")
            .Define("Antimuon_ptot_reco",                           "sqrt(p_reco.X()*p_reco.X() + p_reco.Y()*p_reco.Y() + p_reco.Z()*p_reco.Z())")
            .Define("Antimuon_reconstructed_P4",                    [](const TVector3& P3, const double m, const double ptot){TLorentzVector v = {P3, sqrt(m*m + ptot*ptot)}; return v;}, 
                                                                    {"Antimuon_p_reco", "MUON_MASS_MeV","Antimuon_ptot_reco"})
            /*
                PREDICTED NEUTRON
            */
            .Define("Antimuon_reconstructed_P4_GeV",                 "Antimuon_reconstructed_P4 * 1e-3")
            .Define("Neutrino_reconstructed_P4_GeV",                 RDFUtils::GENIE::GetNup4FromMu,            {"PROTON_MASS_GeV","NEUTRON_MASS_GeV", "MUON_MASS_GeV", "Antimuon_reconstructed_P4_GeV","NuDirection"})
            .Define("PredictedNeutron_P3_GeV",                       "Neutrino_reconstructed_P4_GeV.Vect() - Antimuon_reconstructed_P4_GeV.Vect()")
            .Define("PredictedNeutron_E_GeV",                        "sqrt(PredictedNeutron_P3_GeV.Mag() * PredictedNeutron_P3_GeV.Mag() + NEUTRON_MASS_GeV * NEUTRON_MASS_GeV)")
            .Define("PredictedNeutron_Beta",                         "PredictedNeutron_P3_GeV.Mag() / PredictedNeutron_E_GeV")
            .Define("PredictedNeutron_Angle",                        "Antimuon_reconstructed_P4_GeV.Vect().Angle(PredictedNeutron_P3_GeV)")
            /*
                PREDICTED HITS
            */
           .Define("ExpectedNeutron_HitPosition_",                   RDFUtils::DIGIT::GetExpectedHitPosition,   {"Primaries_vtx", "PredictedNeutron_P3_GeV", "dg_cell"})
           .Define("ExpectedNeutron_HitPosition_x_",                 RDFUtils::GetComponent3<0>,                {"ExpectedNeutron_HitPosition_"})
           .Define("ExpectedNeutron_HitPosition_y_",                 RDFUtils::GetComponent3<1>,                {"ExpectedNeutron_HitPosition_"})
           .Define("ExpectedNeutron_HitPosition_z_",                 RDFUtils::GetComponent3<2>,                {"ExpectedNeutron_HitPosition_"})
           .Define("ExpectedNeutron_FlightLength_",                  RDFUtils::DIGIT::GetFlightLength,          {"Primaries_vtx", "ExpectedNeutron_HitPosition_"})
           .Define("ExpectedNeutron_TOF_",                           RDFUtils::DIGIT::GetTOF,                   {"ExpectedNeutron_FlightLength_", "PredictedNeutron_Beta"})
           .Define("Expected_HitTime_",                              "ExpectedNeutron_TOF_ + Primaries_vtxT")
           /*
                CELLS WITH COINCIDENCES
           */
            .Define("Residuals_HitTime_",                           RDFUtils::DIGIT::TimeResiduals,             {"Expected_HitTime_", "Reconstructed_HitTime"})
            .Define("Residuals_HitSpace_",                          RDFUtils::DIGIT::SpaceResiduals,            {"ExpectedNeutron_HitPosition_", "Reconstructed_HitPosition"})
            .Define("IsCompatible",                                 RDFUtils::DIGIT::IsCompatible,              {"Residuals_HitSpace_", "Residuals_HitTime_"})
            .Define("earliest_compatible_cell",                     RDFUtils::FindMinimum3,                     {"Reconstructed_HitTime","IsCompatible"})
            .Define("nof_compatible_cells",                         [](const ROOT::VecOps::RVec<int>& compatibles){return static_cast<int>(compatibles[compatibles==1].size());}, {"IsCompatible"})
            .Define("candidate_signal_event",                       RDFUtils::RECO::isCandidateSignal,          {"InteractionVolume_short","NofFinalStateChargedParticles","nof_compatible_cells"})
    ;
}