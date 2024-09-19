import uproot4 as upr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import particle
import ROOT
import math
import numpy as np
import sys

sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-physics/notebook/python_tools')
from ROOT_tools import ROOT_tools
from MultiPlotter import MultiPlotter
from ROOT2Pandas import Converter

tool = ROOT_tools()

file_name = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/events-in-SANDtracker.0.to.100.ecal-digit.analysed.root"

tree_name = "digit_extended"

converter = Converter(file_name, tree_name)

columns_df = ['FileName',
              'EventId',
              'EventType',
              'CCQEonHydrogen',
              'Interaction_vtxX',
              'Interaction_vtxY',
              'Interaction_vtxZ',
              'InteractionTarget',
              'InteractionVolume',
              'PrimaryStateHadronicSystemTopology_name',
              'NofFinalStateChargedParticles',
              'PrimaryStateHadronicSystemTotalKinE',
              'MissingTransverseMomentum',
              'ExpectedHadronSystP3',
              'ExpectedNeutronArrivalPositionECAL',
              'ExpectedFiredModuleByNeutron',
 ]

columns_primaries = ['FileName',
                      'EventId',
                      'EventType',
                      'CCQEonHydrogen',
                      'Interaction_vtxX',
                      'Interaction_vtxY',
                      'Interaction_vtxZ',
                      'Interaction_vtxT',
                      'PrimariesPDG',
                      'PrimariesTrackId',
                      'PrimariesP4',
                      'PrimariesFirstHitECAL',
                      'PrimariesEDepECAL',
                      'PrimariesEmissionAngle',
                      'IsECALHitMissing',
                      'DeviationAngle',
                      ]

columns_fired_cells = ['FileName',
                       'CCQEonHydrogen',
                       'EventId',
                       'Fired_Cells_mod',
                       'Fired_Cells_id',
                       'Fired_Cells_x',
                       'Fired_Cells_y',
                       'Fired_Cells_z',
                       'isCellComplete',
                       'Fired_Cells_tdc1',
                       'who_produced_tdc1',
                       'Fired_Cells_tdc2',
                       'who_produced_tdc2',
                       'Fired_Cell_true_hit1',
                       'Fired_Cell_true_hit2',
                       'Cell_Reconstructed_hit',
                       'ExpectedNeutronHit',
                        ]

primaries = converter.CreatePandas(
    columns = columns_primaries,
    rename = True,
    indices = ['FileName','EventId']
)

df = converter.CreatePandas(
    columns = columns_df,
    rename = False,
    indices = ['FileName','EventId']
)

fired_cells = converter.CreatePandas(
    columns = columns_fired_cells,
    rename = True,
    indices = ['FileName','EventId']
)

# ADD COLUMNS

primaries['dist_Vtx2ECAL'] = np.sqrt((primaries['PrimariesFirstHitECAL_x'] - primaries['Interaction_vtxX']*1e3)**2 + 
                                     (primaries['PrimariesFirstHitECAL_y'] - primaries['Interaction_vtxY']*1e3)**2 + 
                                     (primaries['PrimariesFirstHitECAL_z'] - primaries['Interaction_vtxZ']*1e3)**2)
def get_mass_from_pdgid(pdgid):
    try:
        particle_info = particle.Particle.from_pdgid(pdgid)
        if particle_info is not None:
            return particle_info.mass
        else:
            return np.NaN
    except Exception as e:
        print(f"Error occurred for PDG ID {pdgid}: {e}")
        return np.NaN


def InteractionVolume_short(target):
    if("C3H6Target" in target):
        return "C3H6_Target"
    elif("CTarget" in target):
        return "C_Target"
    else:
        return "Other"
print(primaries.columns)
primaries["mass"] = primaries['PrimariesPDG'].apply(get_mass_from_pdgid)
primaries["E_kin"] = primaries['PrimariesP4_t'] - primaries['mass']
primaries["gamma"] = primaries['PrimariesP4_t'] / (primaries["mass"])
primaries["beta"] =np.sqrt(1. - 1. / primaries["gamma"].values**2)

df["InteractionVolume_short"] = df.apply(lambda row: InteractionVolume_short(row['InteractionVolume']), axis=1)


# SIGNAL (MC TRUTH)

events_signal = df[df.CCQEonHydrogen==1]
neutron_signal = primaries[(primaries.CCQEonHydrogen==1)&(primaries.PrimariesPDG==2112)]
muon_signal = primaries[(primaries.CCQEonHydrogen==1)&(primaries.PrimariesPDG==-13)]

neutrons_signal_w_ECALhits = neutron_signal[neutron_signal.IsECALHitMissing == 0]
muons_signal_w_ECALhits = neutron_signal[muon_signal.IsECALHitMissing == 0]

neutrons_signal_w_ECALhits_index = neutrons_signal_w_ECALhits.index.get_level_values(0)
nof_neutrons_signal_at_least_1ECAL_hit = len(neutrons_signal_w_ECALhits)
tot_nof_neutrons_signal = len(neutron_signal)
print(f" nof of neutrons with at least 1 hit in ECAL {nof_neutrons_signal_at_least_1ECAL_hit}, {nof_neutrons_signal_at_least_1ECAL_hit / tot_nof_neutrons_signal *100} [%]")

fired_cells = fired_cells.join(df.ExpectedFiredModuleByNeutron, how = 'left')
neutron_signal_complete_fired_cells = fired_cells[(fired_cells.CCQEonHydrogen==1) 
                                                  & (fired_cells.isCellComplete==1)
                                                  & (fired_cells.who_produced_tdc1 == 1) # neutron has track id 1
                                                  & (fired_cells.who_produced_tdc2 == 1)
                                                  & (fired_cells.Fired_Cells_mod == fired_cells.ExpectedFiredModuleByNeutron)
                                                  ]
neutron_signal_complete_fired_cells_barrels = neutron_signal_complete_fired_cells[neutron_signal_complete_fired_cells.Fired_Cells_mod<30]
neutron_signal_complete_fired_cells_encaps = neutron_signal_complete_fired_cells[neutron_signal_complete_fired_cells.Fired_Cells_mod>=30]
neutron_signal_complete_fired_cells

muon_signal_complete_fired_cells = fired_cells[(fired_cells.CCQEonHydrogen==1) 
                                                  & (fired_cells.isCellComplete==1)
                                                  & (fired_cells.who_produced_tdc1 == 0) # muon has id track id 0
                                                  & (fired_cells.who_produced_tdc2 == 0)
                                                  ]
# for complete cell fired by muons add information time first hit (edepsim) in ECAL by muon
muon_signal_complete_fired_cells = muon_signal_complete_fired_cells.join(muons_signal_w_ECALhits.PrimariesFirstHitECAL_t, how='left')

# PLOTS
all_hists = []
print(f"{'-'*20}")
#___________ NEW PAGE
print("Neutron True Kinetic energy and tof (signal event)")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(20, 6))

# First histogram: Neutron True Kinetic Energy
plotter.plot_hist(
    data=neutron_signal['E_kin'],
    bins=np.arange(0, 800, 10),
    color='blue',
    xlabel="[MeV]",
    ylabel="Counts"
)
plotter.axes[0].set_title("Neutron True Kinetic energy (signal event)", fontsize=15)
plotter.next_plot()

# Second histogram: Neutron True Beta
plotter.plot_hist(
    data=neutron_signal['beta'],
    bins=np.arange(0, 1, 0.02),
    color='blue',
    xlabel="Beta",
    ylabel="Counts"
)
plotter.axes[1].set_title("Neutron True Beta (signal event)", fontsize=15)

all_hists.append(plotter)

#___________ NEW PAGE
print("Neutron angle of deviation wrt initial direction")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 6), suptitle = "Neutrons from signal")

plotter.plot_hist(
    data=neutron_signal["DeviationAngle"],
    bins=np.arange(0,np.pi, 0.1),
    color='blue',
    xlabel='Angle of Deviation [deg]',
    ylabel="counts"
)

plotter.axes[0].set_title("Neutron Angle of Deviation from Initial direction", fontsize=15)

plotter.next_plot()

plotter.plot_hist2d(
    x=neutron_signal['E_kin'],
    y=neutron_signal['DeviationAngle'],
    bins_x=np.arange(0, 800, 10),
    bins_y=np.arange(0,np.pi, 0.1),
    # log_scale=True,
    xlabel='Neutron Kinetic Energy [MeV]',
    ylabel='Angle of Deviation [deg]',
)

# plotter.next_plot()

# plotter.plot_hist2d(
#     x=neutron_signal.loc[neutron_signal_complete_fired_cells_barrels.index.get_level_values(0).unique()]['DeviationAngle'],
#     y=neutron_signal_complete_fired_cells_barrels['Fired_Cell_true_hit1_t'] - neutron_signal_complete_fired_cells_barrels['Cell_Reconstructed_hit_t'],
#     bins_x=np.arange(-np.pi,np.pi, 0.1),
#     bins_y= np.arange(-50, 50, 1),
#     # log_scale=True,
#     xlabel='Angle of Deviation [deg]',
#     ylabel=r'neutron $tof_{true}$ - $tof_{reco}$',
# )

all_hists.append(plotter)

#___________ NEW PAGE
print("Neutron ADC count, Kinetic Energy and Energy deposit comparison")

#___________ NEW PAGE
print("Comparison Time of flight to ECAL for neutrons and muons")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(20, 6))

# First histogram: Time of flight to ECAL for neutrons and muons
plotter.plot_hist(
    data=neutrons_signal_w_ECALhits.PrimariesFirstHitECAL_t,
    bins=np.arange(0, 200, 1),
    label='neutron',
    color='blue',
    xlabel="[ns]",
    ylabel="Counts",
    xlim=(0, 200),
    ticks=np.arange(0, 200, 10)
)
plotter.plot_hist(
    data=muon_signal.loc[neutrons_signal_w_ECALhits_index].PrimariesFirstHitECAL_t,
    bins=np.arange(0, 200, 1),
    label='muon',
    color='red',
)
plotter.set_labels(xlabel="[ns]", ylabel="Counts", fontsize=15)
plotter.axes[0].set_yscale('log')
plotter.axes[0].set_title("Time of flight to ECAL (signal event)", fontsize=15)
plotter.next_plot(plot_legend=True)

# Second histogram: Difference between muon and neutron time of flight
plotter.plot_hist(
    data=muon_signal.loc[neutrons_signal_w_ECALhits_index].PrimariesFirstHitECAL_t - neutrons_signal_w_ECALhits.PrimariesFirstHitECAL_t,
    bins=np.arange(-100, 10, 1),
    xlabel="[ns]",
    ylabel="Counts",
    xlim=(-100, 10),
    ticks=np.arange(-100, 10, 10)
)
plotter.axes[1].set_title("muon tof - neutron tof", fontsize=15)

all_hists.append(plotter)

#___________ NEW PAGE
print("Neutron ECAL hits comparison: true (EDEPSIM) vs reconstructed using cell TDCs")
plotter = MultiPlotter(nrows=1, ncols=3, figsize=(30, 8), suptitle="Neutron ECAL hits components: true (EDEPSIM) vs reconstructed using cell TDCs")

# First 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x=neutron_signal_complete_fired_cells_barrels['Fired_Cell_true_hit1_x'],
    y=neutron_signal_complete_fired_cells_barrels['Cell_Reconstructed_hit_x'],
    bins_x=np.arange(-2000, 2000, 25),
    bins_y=np.arange(-2000, 2000, 25),
    log_scale=True,
    xlabel=r'hit $x_{true}$ (edepsim)',
    ylabel=r'hit $x_{reco}$ from TDCs of cell',
)
plotter.next_plot()

# Second 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x=neutron_signal_complete_fired_cells['Fired_Cell_true_hit1_y'],
    y=neutron_signal_complete_fired_cells['Cell_Reconstructed_hit_y'],
    bins_x=np.arange(-4500, 100, 25),
    bins_y=np.arange(-4500, 100, 25),
    log_scale=True,
    xlabel=r'hit $y_{true}$ (edepsim)',
    ylabel=r'hit $y_{reco}$ from TDCs of cell',
)
plotter.next_plot()

# Third 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x=neutron_signal_complete_fired_cells['Fired_Cell_true_hit1_z'],
    y=neutron_signal_complete_fired_cells['Cell_Reconstructed_hit_z'],
    bins_x=np.arange(22000, 22800, 25),
    bins_y=np.arange(22000, 22800, 25),
    log_scale=True,
    xlabel=r'hit $z_{true}$ (edepsim)',
    ylabel=r'hit $z_{reco}$ from TDCs of cell',
)

all_hists.append(plotter)

#___________ NEW PAGE
print("Missing Transverse Momentum comparison C and C3H6")
ratio_C12_H = 311685 / 55008
ratio_mass_target_C3H6_mass_target_C = 3.2 / 0.7

plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 8), suptitle="Missing Transverse Momentum")

plotter.plot_hist(
    data=df[df.CCQEonHydrogen==0]['MissingTransverseMomentum'],
    bins=np.arange(-0.1,1.5,0.01),
    weights=np.ones(len(df[df.CCQEonHydrogen==0]['MissingTransverseMomentum'])) * (1/ratio_C12_H),
    color='blue',
    label=r'$\overline{\nu}_\mu$ CC on Carbon',
)

plotter.plot_hist(
    data=df[df.CCQEonHydrogen==1]['MissingTransverseMomentum'], 
    bins=np.arange(-0.1, 1.5, 0.01),
    color='red', 
    label=r'$\overline{\nu}_\mu$ CC on H',
    xlabel="[GeV]",
)

plotter.add_legend(labels=[r'$\overline{\nu}_\mu$ CC on Carbon', r'$\overline{\nu}_\mu$ CC on H'])
plotter.axes[plotter.current_ax].set_title(r"Missing Transverse Momentum : $p_T^m = |\vec{p_T}^\mu - \vec{p_T}^H|$", fontsize=16)
plotter.axes[plotter.current_ax].set_yscale("log")

plotter.next_plot()

plotter.plot_hist(
    data=df[df.InteractionVolume_short=='C_Target']['MissingTransverseMomentum'], 
    bins=np.arange(-0.1, 1.5, 0.01),
    color='blue', 
    label=r'$\overline{\nu}_\mu$ CC on Carbon Target'
)

plotter.plot_hist(
    data=df[df.InteractionVolume_short=='C3H6_Target']['MissingTransverseMomentum'], 
    bins=np.arange(-0.1, 1.5, 0.01),
    weights=np.ones(len(df[df.InteractionVolume_short=='C3H6_Target']['MissingTransverseMomentum'])) * (1 / ratio_mass_target_C3H6_mass_target_C),
    color='red', 
    label=r'$\overline{\nu}_\mu$ CC on C3H6 Target'
)

plotter.add_legend(labels=[r'$\overline{\nu}_\mu$ CC on Carbon Target', r'$\overline{\nu}_\mu$ CC on C3H6 Target'])
plotter.set_labels(xlabel="[GeV]")
plotter.set_limits(ylim=(1, None))  # For log scale y-axis
plotter.axes[plotter.current_ax].set_yscale("log")
plotter.axes[plotter.current_ax].set_title(r"Missing Transverse Momentum : $p_T^m = |\vec{p_T}^\mu - \vec{p_T}^H|$", fontsize=16)

all_hists.append(plotter)

#___________ NEW PAGE
print("Neutron tof to ECAL (signal, complete cells) comparison: true, prediction, reconstructed")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 8), suptitle="Neutron from signal, complete cells")

plotter.plot_hist2d(
    x=neutron_signal_complete_fired_cells['Fired_Cell_true_hit2_t'],
    y=neutron_signal_complete_fired_cells['ExpectedNeutronHit_t'],
    bins_x=np.arange(0,40,0.2),
    bins_y=np.arange(0,40,0.2),
    log_scale=True,
    # label='',
    xlabel='True neutron tof (edepsim) [ns]',
    ylabel=r'Predicted neutron tof from $\mu^+$ [ns]',
    xlim=[0,40],
    ylim=[0,40]
)

plotter.axes[plotter.current_ax].plot([0, 40], [0, 40], color='red', linestyle='--', linewidth=2)

plotter.next_plot()

plotter.plot_hist2d(
    x=neutron_signal_complete_fired_cells['Fired_Cell_true_hit1_t'],
    y=neutron_signal_complete_fired_cells['Cell_Reconstructed_hit_t'],
    bins_x=np.arange(0,40,0.2),
    bins_y=np.arange(0,40,0.2),
    log_scale=True,
    xlabel='True neutron neutron tof (edepsim) [ns]',
    ylabel='Reconstructed neutron tof from cell [ns]',
    xlim=[0,40],
    ylim=[0,40]
)

plotter.axes[plotter.current_ax].plot([0, 40], [0, 40], color='red', linestyle='--', linewidth=2)

all_hists.append(plotter)

#___________ NEW PAGE
print("Muon tof to ECAL (signal, complete cells) comparison: true, reconstructed")
plotter = MultiPlotter(nrows=1, ncols=1, figsize=(10, 8), suptitle=r"$\mu^+$ from signal , complete cells")

plotter.plot_hist2d(
    x=muon_signal_complete_fired_cells['Fired_Cell_true_hit1_t'],
    y=muon_signal_complete_fired_cells['Cell_Reconstructed_hit_t'],
    bins_x=np.arange(0, 50, 0.1),
    bins_y=np.arange(0, 50, 0.1),
    log_scale=True,
    xlabel=r'True tof $\mu^+$ (edepsim) [ns]',
    ylabel=r'Reconstructed $\mu^+$ tof from fired cells [ns]',
    xlim=[0,20],
    ylim=[0,20],
)

plotter.axes[plotter.current_ax].plot([0, 20], [0, 20], color='red', linestyle='--', linewidth=2)

all_hists.append(plotter)

#_________ WRITE PDF

MultiPlotter().save_multiple_figures_to_pdf('/storage/gpfs_data/neutrino/users/gi/sand-physics/scratch/multiple_plots.pdf', all_hists)

