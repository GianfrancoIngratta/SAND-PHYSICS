#!/bin/bash
# PYTHIA6
FILE_INDEX_START=$1
export LD_LIBRARY_PATH=/opt/exp_software/neutrino/PYTHIA6/Pythia6Support/v6_424/lib:${LD_LIBRARY_PATH}

# ROOT
source /opt/exp_software/neutrino/ROOT/v6.20.00_py3/bin/thisroot.sh

# GEANT4
source /opt/exp_software/neutrino/GEANT4/4.10.05/bin/geant4.sh

# EDEP-SIM
source /opt/exp_software/neutrino/EDEPSIM/setup.sh

# In this PATH there are the libraries transferred during job submission
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/storage/gpfs_data/neutrino/users/gi/prod-scripts

# GENIE V2
source /opt/exp_software/neutrino/GENIEv2/Generator/setup.sh
source /storage/gpfs_data/neutrino/users/gi/sand-physics/setup.sh # to run executable in sand-physics

# SANDRECO
source /storage/gpfs_data/neutrino/users/gi/sand-reco/setup.sh

cd /storage/gpfs_data/neutrino/users/gi/sand-physics/

./build/bin/AnalyseEDepSim ${FILE_INDEX_START}