import time
start = time.time()

import numpy as np
import sys
import pandas as pd
sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-physics/notebook/python_tools')
from ROOT_tools import ROOT_tools
from MultiPlotter import MultiPlotter
from ROOT2Pandas import Converter
from SampleManager import Sample, Manager
tool = ROOT_tools()
import uproot4 as upr
print("Import time:", time.time() - start)

columns_df = [
    # /*
    #     EVENT INFO
    # */
    "FileName",
    "CCQEonHydrogen",
    "EventType",
    "NuDirection",
    "Interaction_vtxX",
    "Interaction_vtxY",
    "Interaction_vtxZ",
    "Interaction_vtxT",
    "InteractionVolume_short",
    "NofFinalStateChargedParticles",
    "PrimaryStateHadronicSystemTopology_name",
    "InteractionTarget",
    "candidate_signal_event",
    "nof_fired_wires",
    "IncomingNeutrinoP4",
    "FinalStateHadronicSystemTotal4Momentum",
    "Antimuon_p_true",
    # /*
    #     RECONSTRUCTED ANTIMUON
    # */
    "Antimuon_reconstructed_P4",
    # /*
    #     PREDICTED NEUTRON
    # */
    "Neutrino_reconstructed_P4_GeV",
    "PredictedNeutron_P3_GeV",
    "PredictedNeutron_E_GeV",
    "PredictedNeutron_Beta",
    # /*
    #     EVENT FIRED CELL GENERAL INFO    
    # */
    "Fired_Cells_mod",
    "Fired_Cells_id",
    "Fired_Cells_x",
    "Fired_Cells_y",
    "Fired_Cells_z",
    "isCellComplete",
    "AreTDCsConsistent",
    "Fired_Cells_adc1",
    "Fired_Cells_tdc1",
    "Fired_Cells_adc2",
    "Fired_Cells_tdc2",
    "Fired_Cell_true_hit1",
    "Fired_Cell_true_hit2",
    "Fired_by_primary_neutron",
    "Fired_by_primary_antimu",  
    # /*
    #     TRUE NEUTRON HITS
    # */
    "Fired_Cell_true_Hit_x",
    "Fired_Cell_true_Hit_y",
    "Fired_Cell_true_Hit_z",
    "Fired_Cell_true_Hit_t",
    "Fired_Cell_true_Hit_e",
    "True_FlightLength",
    # /*
    #     PREDICTED NEUTRON HITS
    # */
    "ExpectedNeutron_HitPosition_x_",
    "ExpectedNeutron_HitPosition_y_",
    "ExpectedNeutron_HitPosition_z_",
    "ExpectedNeutron_FlightLength_",
    "ExpectedNeutron_TOF_",
    "Expected_HitTime_",
    # /*
    #     RECONSTRUCTED NEUTRON HITS
    # */
    "Reconstructed_HitPosition_x",
    "Reconstructed_HitPosition_y",
    "Reconstructed_HitPosition_z",
    "Reconstructed_HitTime",
    "Reconstructed_Energy",
    "IsEarliestCell_neutron",
    "Reconstructed_FlightLength",
    # /*
    #     CELLS WITH COINCIDENCES
    # */
   "Residuals_HitTime_",
   "Residuals_HitSpace_",
   "IsCompatible",
   "earliest_compatible_cell",
   "nof_compatible_cells",
#    /*
#         RECONSTRUCTED NEUTRON AND NEUTRINO ENERGY
#    */
    "reconstructed_neutron_KinE_MeV",
]

production_folder="/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/production_0000"

# for i in range(0,100):
#     start=i*10
#     end = start + 9
#     fInput  = f"{production_folder}/ecal_prediction/events-in-SANDtracker.{start}.to.{end}.ecal_prediction.root"
#     fOutput = f"{production_folder}/ecal_prediction/csv/events-in-SANDtracker.{start}.to.{end}.ecal_prediction.csv"
#     print(f"file start : {start}, file end : {end}")
#     converter = Converter([fInput], "ecal_prediction")

#     start = time.time()
#     print("converting dataframe")
#     df = converter.CreatePandas(columns=columns_df, rename=True, indices=['FileName'])
#     # df_signal = df[df.CCQEonHydrogen==1]

#     print("production converted, time to convert: ", time.time() - start)

#     start = time.time()
#     df.to_csv(fOutput)
#     print("production to csv, time to write: ", time.time() - start)
#     print("*"*20)

import glob 

files=glob.glob(production_folder+"/ecal_prediction/csv/*")
df0 = pd.read_csv(files[0])
df_signal = df0[df0.CCQEonHydrogen==1]

for i, f in enumerate(files[1:]):
    percentage = (i / len(files[1:])) * 100
    print(f"\r{percentage:.2f}% completed.", end="")
    print()
    df = pd.read_csv(f)
    df_signal = df_signal.append(df[df.CCQEonHydrogen==1])

df_signal.to_csv(f"{production_folder}/ecal_prediction/csv/events-in-SANDtracker.0.to.999.ecal_prediction.CCQEonH.csv")