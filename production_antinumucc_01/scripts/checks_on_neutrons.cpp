TGeoManager* geo = NULL;

ROOT::VecOps::RVec<int> GetPrimariesPDG(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    ROOT::VecOps::RVec<int> pdgs;
    auto Primaries = v[0].Particles;
    for(const auto& primary : Primaries){
        pdgs.push_back(primary.GetPDGCode());
    }
    return pdgs;
}

bool is_signal_event(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    auto reaction = v[0].GetReaction();
    // std::cout << "reaction " << reaction << "\n";
    return true;
}


bool IsCC(TObjString& s){
    if(s.GetString().Contains("[NC]")){
        return false;
    }else{
        return true;
    }
}


std::string EventType(TObjString& s){
    
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



double Getnu(ROOT::VecOps::RVec<double>& P4){
    TLorentzVector neutrino = {P4[0], P4[1], P4[2], P4[3]};
    TLorentzVector proton = {P4[4], P4[5], P4[6], P4[7]};
    TLorentzVector lepton = {P4[8], P4[9], P4[10], P4[11]};
    TLorentzVector neutron = {P4[12], P4[13], P4[14], P4[15]};

    double neutron_mass = 0.939;
    double neutron_momentum = neutron.Vect().Mag();
    double neutron_energy = sqrt(neutron_momentum*neutron_momentum + neutron_mass*neutron_mass);
    double nu = neutron_energy - neutron_mass;
    return nu;
}

double GetQ2(ROOT::VecOps::RVec<double>& P4){
    TLorentzVector neutrino = {P4[0], P4[1], P4[2], P4[3]};
    TLorentzVector proton = {P4[4], P4[5], P4[6], P4[7]};
    TLorentzVector lepton = {P4[8], P4[9], P4[10], P4[11]};
    TLorentzVector neutron = {P4[12], P4[13], P4[14], P4[15]};

    auto Q2 = (neutrino.Vect() - lepton.Vect()).Mag2();
    return Q2;
}

TLorentzVector get_nu_energy(ROOT::VecOps::RVec<double>& P4){
    TLorentzVector neutrino = {P4[0], P4[1], P4[2], P4[3]};
    return neutrino;
}

TLorentzVector GetProtonInitialMomentum(ROOT::VecOps::RVec<double>& P4){
    return {P4[4], P4[5], P4[6], P4[7]};
}

TLorentzVector GetMuonInitialMomentum(ROOT::VecOps::RVec<double>& P4){
    return {P4[8], P4[9], P4[10], P4[11]};
}

TLorentzVector GetNeutronInitialMomentum(ROOT::VecOps::RVec<double>& P4){
    return {P4[12], P4[13], P4[14], P4[15]};
}

std::vector<TG4TrajectoryPoint> GetNeutronTrajectoryPoints(ROOT::VecOps::RVec<TG4Trajectory>& trajectories) {
    std::vector<TG4TrajectoryPoint> neutron_points;

    for (size_t i = 1; i < trajectories.size(); i++) {
        int pdg = trajectories.at(i).GetPDGCode();
        int parent_id = trajectories.at(i).GetParentId();
        if ((pdg == 2112) && (parent_id == -1)) {
            return trajectories.at(i).Points;
        }
    }
    return neutron_points;
}

TVector3 GetVertex(ROOT::VecOps::RVec<double>& EvtVtx){
    return {EvtVtx[0]*1e3, EvtVtx[1]*1e3, EvtVtx[2]*1e3};
}

const double SAND_CENTER_X = 0.;
const double SAND_CENTER_Y = -2384.73;
const double SAND_CENTER_Z = 23910.;
const double SAND_INNER_VOL_X_LENGTH = 3380.0;
const double SAND_INNER_VOL_DIAMETER = 4000.0;

double IntersectWithInfCylinder(const TVector3& vtx, const TVector3& direction, double diameter) {
    /*
        Calculate the intersection btw a line that starts from vtx and travels
        in the direction neutron_direction (versor) and an infinite cylinder.
        Take the intersection along the direction of travel, discard the other one.
    */
    TVector3 start_point(vtx.X() - SAND_CENTER_X, 
                         vtx.Y() - SAND_CENTER_Y, 
                         vtx.Z() - SAND_CENTER_Z);
    
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

double IntersectWithTube(const TVector3& vertex, const TVector3& direction) {
    /*
        Given the neutrino vertex position and the neutron momentum
        predict where neutron is expected to release energy.
        !!! THE RETURNED POINT'S UNITS IS mm 
    */
    // convert to mm

    double t = IntersectWithInfCylinder(vertex, direction, SAND_INNER_VOL_DIAMETER);

    if (fabs(vertex.X() + t * direction.X() - SAND_CENTER_X) >= SAND_INNER_VOL_X_LENGTH / 2.0) {
        
        double t1 = (SAND_CENTER_X + SAND_INNER_VOL_X_LENGTH / 2.0 - vertex.X()) / direction.X();
        double t2 = (SAND_CENTER_X - SAND_INNER_VOL_X_LENGTH / 2.0 - vertex.X()) / direction.X();
      
        if(t1<0){
            t = t2;
        }else{
            t = t1;
        }
    }

    return t;
}


TLorentzVector GetExpected1stHit(TVector3 vertex, TLorentzVector neutron_initial_momentum){

    TVector3 neutron_direction = neutron_initial_momentum.Vect().Unit();

    float min_t = IntersectWithTube(vertex, neutron_direction);

    TVector3 entry_point = vertex + min_t * neutron_direction;

    float neutron_momentum = neutron_initial_momentum.Vect().Mag();

    float neutron_energy = sqrt(neutron_momentum*neutron_momentum + 0.939*0.939);

    double beta_neutron = (neutron_momentum) / neutron_energy;

    float flight_length = (vertex - entry_point).Mag();

    float c = 299; // mm/ns
    
    double expected_time = flight_length / (beta_neutron * c);

    return {entry_point, expected_time + 1.};
}

TLorentzVector Get1stHitECAL(const std::vector<TG4TrajectoryPoint>& points){

    TLorentzVector entry_point = {0.,0.,0.,0.};

    auto navigator = geo->GetCurrentNavigator();

    if(!navigator) navigator = geo->AddNavigator();

    for(const auto& point : points) 
    {
        auto pos = point.GetPosition();
        auto node = navigator->FindNode(pos.X()*0.1, pos.Y()*0.1, pos.Z()*0.1);
        TString volume = node->GetName();
        // std::cout << volume << "\n";
        if(volume.Contains("ECAL")) return pos;
    }
    return entry_point;
}

// ROOT::VecOps::RVec<double> GetTimeRes(const TLorentzVector& expected, const ROOT::VecOps::RVec<TLorentzVector>& truth){
//     ROOT::VecOps::RVec<double> res(expected.size());
//     for(size_t i = 0; i < expected.size(); i++)
//     {
//         res[i] = expected[i].T() - truth[i].T();
//     }
//     return res;
// }

// ROOT::VecOps::RVec<double> GetSpaceRes(const ROOT::VecOps::RVec<TLorentzVector>& expected, const ROOT::VecOps::RVec<TLorentzVector>& truth){
//     ROOT::VecOps::RVec<double> res(expected.size());
//     for(size_t i = 0; i < expected.size(); i++)
//     {
//         res[i] = (expected[i].Vect() - truth[i].Vect()).Mag();
//     }
//     return res;
// }

int NofHits(TG4Event& ev, int id)
{
    auto calo_hits = ev.SegmentDetectors["EMCalSci"];
    int hits = 0;
    for(size_t j = 0; j < calo_hits.size(); j++){
        
        int hit_primary_id = calo_hits[j].GetPrimaryId();
        
        if(hit_primary_id == id) hits ++;
    }
    return hits;
}

TLorentzVector GetNeutron1stECALhit(TG4Event& ev, int primary_neutron_id){

    auto calo_hits = ev.SegmentDetectors["EMCalSci"];

    for(size_t j = 0; j < calo_hits.size(); j++){
       int hit_primary_id = calo_hits[j].GetPrimaryId();
    //    std::cout << "hit_primary_id " << hit_primary_id << "\n";
       if(hit_primary_id == primary_neutron_id)
       {   
           auto middle = (calo_hits[j].GetStart() + calo_hits[j].GetStop())*0.5;
           return middle;
       }
    }
    return {0.,0.,0.,0.};
}

int GetPrimaryNeutronTrackId(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v){
    auto Primaries = v[0].Particles;
    // std::cout << "primaries size : " << Primaries.size() << "\n";
    // std::cout << "code 2 " << Primaries[1].GetTrackId() << "\n";
    for(const auto& primary : Primaries){
        if(primary.GetPDGCode() == 2112){
            return primary.GetTrackId();
        }
    }
    return -2;
}

int checks_on_neutrons(){

    auto PRODUCTION_FOLDER = TString::Format("/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc_01/production_0016/");
    auto DATA_PREUNFOLD = TString::Format("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc_01/data_preunfold/");
    auto FLUX = TString::Format("%sfluxes_merged.root", DATA_PREUNFOLD.Data());

    TFile *fFlux = TFile::Open(FLUX.Data());
    TH1F* hE = new TH1F("hE", "Distribuzione di E;E (GeV);Conteggi", 24, 0, 6);
    TTree* flux = (TTree*)fFlux->Get("flux");
    flux->Draw("E >> hE");


    TChain* chain_genie = new TChain("gRooTracker", "gRooTracker");
    TChain* chain_edep = new TChain("EDepSimEvents", "EDepSimEvents");

    TString fgenie("");
    TString fedep(""); 
    
    for(size_t i = 1; i < 10; i++)
    {
        int ii = i;
        fgenie = TString::Format("%srun_1600%d/events-in-SANDtracker.1600%d.gtrac.root", PRODUCTION_FOLDER.Data(), ii, ii);
        fedep = TString::Format("%srun_1600%d/events-in-SANDtracker.1600%d.edep-sim.root", PRODUCTION_FOLDER.Data(), ii, ii);
        chain_genie -> Add(fgenie.Data());
        chain_edep -> Add(fedep.Data());
    }

    // for(size_t i = 10; i < 99; i++)
    // {
    //     int ii = i;
    //     fgenie = TString::Format("%srun_160%d/events-in-SANDtracker.160%d.gtrac.root", PRODUCTION_FOLDER.Data(), ii, ii);
    //     fedep = TString::Format("%srun_160%d/events-in-SANDtracker.160%d.edep-sim.root", PRODUCTION_FOLDER.Data(), ii, ii);
    //     chain_genie -> Add(fgenie.Data());
    //     chain_edep -> Add(fedep.Data());
    // }

    // for(size_t i = 100; i < 248; i++)
    // {
    //     int ii = i;
    //     fgenie = TString::Format("%srun_16%d/events-in-SANDtracker.16%d.gtrac.root", PRODUCTION_FOLDER.Data(), ii, ii);
    //     fedep = TString::Format("%srun_16%d/events-in-SANDtracker.16%d.edep-sim.root", PRODUCTION_FOLDER.Data(), ii, ii);
    //     chain_genie -> Add(fgenie.Data());
    //     chain_edep -> Add(fedep.Data());
    // }

    // chain_genie -> Add(fgenie2.Data());
    // chain_genie -> Add(fgenie3.Data());
    
    // chain_edep -> Add(fedep2.Data());
    // chain_edep -> Add(fedep3.Data());

    chain_edep -> AddFriend(chain_genie, "genie");

    std::cout << "file edep input : " << fedep.Data() << " \n";
    std::cout << "nod events in chain : " << chain_edep -> GetEntries() << " \n";

    ROOT::RDataFrame df(*chain_edep);

    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    // ROOT::EnableImplicitMT();
    geo->SetMaxThreads(100);

    auto df_ = df
                .Define("PrimariesPDG", GetPrimariesPDG, {"Primaries"})
                .Define("EventType", EventType, {"EvtCode"})
                .Define("InteractionTargetPDG", "StdHepPdg[1]")
                .Define("vertex", GetVertex ,{"EvtVtx"})
                .Define("CCQEonHydrogen", [](std::string event_type, int target_pdg){return (event_type=="QES")*(target_pdg==2212);}, {"EventType", "InteractionTargetPDG"})
                .Filter("CCQEonHydrogen")
                .Define("nu", Getnu, {"StdHepP4"})
                .Define("Q2", GetQ2, {"StdHepP4"})
                .Define("neutrino_P4", get_nu_energy, {"StdHepP4"})
                .Define("neutrino_energy", "neutrino_P4.T()")
                .Define("proton_intial_momentum", GetProtonInitialMomentum, {"StdHepP4"})
                .Define("muon_intial_momentum", GetMuonInitialMomentum, {"StdHepP4"})
                .Define("muon_track_id", "0")
                .Define("neutron_intial_momentum", GetNeutronInitialMomentum, {"StdHepP4"})
                .Define("neutron_track_id", GetPrimaryNeutronTrackId, {"Primaries"})
                .Define("nof_hits_neutron", NofHits, {"Event", "neutron_track_id"})
                .Define("nof_hits_muon", NofHits, {"Event", "muon_track_id"})
                .Define("neutron_trajectories", GetNeutronTrajectoryPoints, {"Trajectories"})
                .Define("neutron_1stECAL_hit", GetNeutron1stECALhit, {"Event", "neutron_track_id"})
                .Define("true_ECAL_entry_point", Get1stHitECAL, {"neutron_trajectories"})
                .Define("expected_ECAL_entry_point", GetExpected1stHit, {"vertex", "neutron_intial_momentum"})
                .Define("time_residuals", "expected_ECAL_entry_point.T() - true_ECAL_entry_point.T()")
                .Define("space_residuals", "(expected_ECAL_entry_point.Vect() - true_ECAL_entry_point.Vect()).Mag()")
                ;


    df_.Snapshot("tCheckNeutrons", "check_neutrons_edep_.root", 
                                        {
                                        "PrimariesPDG",
                                        "EvtDXSec",
                                        "EvtXSec",
                                        "Q2",
                                        "nu",
                                        "neutrino_energy",
                                        "muon_intial_momentum",
                                        "proton_intial_momentum",
                                        "neutron_intial_momentum",
                                        "neutron_trajectories",
                                        "neutron_1stECAL_hit",
                                        "true_ECAL_entry_point",
                                        "expected_ECAL_entry_point",
                                        "time_residuals",
                                        "space_residuals",
                                        "nof_hits_neutron",
                                        "nof_hits_muon",
                                        });
    
    // auto antineutrino_flux = df_.Histo1D({"neutrino_energy", "interacting neutrino true Energy;[GeV]", 24u, 0., 6.}, "neutrino_energy");

    // TCanvas* canvas = new TCanvas("","",900,700);

    // antineutrino_flux->Draw("HIST");
    // flux->Draw("HIST SAME");
    // canvas->Draw();

    return 0;
}