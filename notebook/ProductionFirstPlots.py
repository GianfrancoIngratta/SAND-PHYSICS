import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# import particle
import ROOT as r
import math
import sys
import glob

sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-physics/notebook/python_tools')
from ROOT_tools import ROOT_tools
from MultiPlotter import MultiPlotter
from ROOT2Pandas import Converter
from SampleManager import Sample, Manager
tool = ROOT_tools()

production = glob.glob("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_reverse_current_volSAND/events-in-volSAND.*.to.*.edep.analysed.root")

print("list of files in the production : ", production)
print("")
print("Intializing converter")
converter = Converter([production[0]], "edep_extended")

columns = [
            "FileName",
            "EventId",
            "EventType",
            "CCQEonHydrogen",
            "NofEvents",
            "Interaction_vtxX",
            "Interaction_vtxY",
            "Interaction_vtxZ",
            "Interaction_vtxT",
            "IncomingNeutrinoP4",
            "InteractionVolume",
            "InteractionTarget",
            "FinalStateLeptonPDG",
            "FinalStateLeptonNames",
            "FinalStateLepton4Momentum",
            "FinalStateLeptonEmissionAngle",
            "NofPrimaries",
            "NofFinalStateChargedParticles",
            "FinalStateHadronicSystemPDG", 
            "FinalStateHadronicSystemNames", 
            "FinalStateHadronicSystemMomenta", 
            "FinalStateHadronicSystemTotal4Momentum", 
            "FinalStateHadronicSystemEmissionAngle", 
            "FinalStateHadronicSystemTotalKinE", 
]

df = converter.CreatePandas(
    columns = columns,
    rename = True,
    indices = ['FileName']
)

print("Number of Events in dataframe : ", len(df))