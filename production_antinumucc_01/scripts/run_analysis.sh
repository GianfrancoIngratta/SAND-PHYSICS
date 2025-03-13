#!/bin/bash

# Sorgente dello script di setup
source /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/setup.sh

# Lancia ROOT ed esegue i comandi specificati
root -l -i <<EOF
.L /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/src/RooUnfold.h
.L efficiency_plots_test.cpp
efficiency_plots_test()
EOF
