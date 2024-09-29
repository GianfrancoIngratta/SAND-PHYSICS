import matplotlib.pyplot as plt
import numpy as np
import particle
import ROOT
import numpy as np
import sys
import pandas as pd
import glob

sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-physics/notebook/python_tools')
from ROOT_tools import ROOT_tools
from MultiPlotter import MultiPlotter
from ROOT2Pandas import Converter
from SampleManager import Sample, Manager
tool = ROOT_tools()

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

def calculate_tof_ns(distance_mm, beta):
    c = 3e8  # Velocità della luce in m/s
    distance_m = distance_mm / 1000  # Convertire la distanza da mm a metri
    tof_s = distance_m / (beta * c)  # Calcolare il tempo di volo in secondi
    tof_ns = tof_s * 1e9  # Convertire il tempo di volo in nanosecondi
    return tof_ns

def GetLinePoints(vx, vy, vz, vtx, vty, vtz):
    x, y, z = ([] for i in range(3))
    
    # Normalize the direction vector
    norm = np.sqrt(vx**2 + vy**2 + vz**2)
    direction = np.array([vx / norm, vy / norm, vz / norm])
    
    start = np.array([vtx, vty, vtz])
    
    for i in np.arange(0, 4, 0.1):
        # Compute the point on the line at step i
        point = start + i * direction
        # Append the coordinates to the respective lists
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    
    return {"x":x, "y":y,"z" : z}

def InteractionVolume_short(target):
    if("C3H6Target" in target):
        return "C3H6_Target"
    elif("CTarget" in target):
        return "C_Target"
    else:
        return "Other"
    
production = glob.glob("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/events-in-SANDtracker.*.to.*.ecal-digit.analysed.root")

converter = Converter(production, "digit_extended")

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
                       'Fired_Cells_adc1',
                       'Fired_Cells_tdc1',
                       'who_produced_tdc1',
                       'Fired_Cells_adc2',
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
    indices = ['FileName', 'PrimariesTrackId']
)
print(primaries)

df = converter.CreatePandas(
    columns = columns_df,
    rename = False,
    indices = ['FileName']
)
print(df)

fired_cells = converter.CreatePandas(
    columns = columns_fired_cells,
    rename = True,
    indices = ['FileName','Fired_Cells_id']
)
print(fired_cells)

primaries['dist_Vtx2ECAL'] = np.sqrt((primaries['PrimariesFirstHitECAL_x'] - primaries['Interaction_vtxX']*1e3)**2 + 
                                     (primaries['PrimariesFirstHitECAL_y'] - primaries['Interaction_vtxY']*1e3)**2 + 
                                     (primaries['PrimariesFirstHitECAL_z'] - primaries['Interaction_vtxZ']*1e3)**2)
primaries["mass"] = primaries['PrimariesPDG'].apply(get_mass_from_pdgid)
primaries["E_kin"] = primaries['PrimariesP4_t'] - primaries['mass']
primaries["gamma"] = primaries['PrimariesP4_t'] / (primaries["mass"])
primaries["beta"] =np.sqrt(1. - 1. / primaries["gamma"].values**2)

df["InteractionVolume_short"] = df.apply(lambda row: InteractionVolume_short(row['InteractionVolume']), axis=1)

primaries['HasChangedDirection'] = primaries.apply(
    lambda row: 1 if abs(row['DeviationAngle']) < 0.01 else (np.nan if row['DeviationAngle'] == -999 else 0),
    axis=1
)

df_manager = Manager(df, "df", reference_index="FileName")
primaries_manager = Manager(primaries, "primaries", reference_index="FileName")
fired_cells_manager = Manager(fired_cells, "fired_cells", reference_index="FileName")

df_manager.DefineSample("signal", "CCQEonHydrogen==1")

primaries_manager.DefineSample("signal", "CCQEonHydrogen==1")
primaries_manager.DefineSample("signal_mu+", "CCQEonHydrogen==1 & PrimariesPDG==-13")
primaries_manager.DefineSample("signal_neutrons", "CCQEonHydrogen==1 & PrimariesPDG==2112")
primaries_manager.DefineSample("signal_neutrons_w_hits", "CCQEonHydrogen==1 & PrimariesPDG==2112 & IsECALHitMissing == 0")
primaries_manager.DefineSample("signal_neutrons_no_direction_change", "CCQEonHydrogen==1 & PrimariesPDG==2112 & HasChangedDirection == 0")
primaries_manager.DefineSample("signal_neutrons_w_direction_change", "CCQEonHydrogen==1 & PrimariesPDG==2112 & HasChangedDirection == 1")
primaries_manager.DefineSample("signal_neutrons_reconstructable", "CCQEonHydrogen==1 & PrimariesPDG==2112 & HasChangedDirection == 0 & IsECALHitMissing == 0")

fired_cells_manager.DefineSample("complete_fired_by_track_1", "isCellComplete==1 & who_produced_tdc1==1")
fired_cells_manager.DefineSample("complete_fired_by_track_0", "isCellComplete==1 & who_produced_tdc1==0")

primaries_manager.CombineSamples("signal_neutrons_reconstructable", 
                                 fired_cells_manager,
                                 "complete_fired_by_track_1",
                                 "signal_neutrons_reconstructable_complete")

primaries_manager.CombineSamples("signal_mu+",
                                 fired_cells_manager,
                                 "complete_fired_by_track_0",
                                 "signal_antimu_complete")

fired_cells_manager.CombineSamples("complete_fired_by_track_1",
                                   primaries_manager,
                                   "signal_neutrons_reconstructable",
                                   "fired_by_reconstructable_neutron")

fired_cells_manager.CombineSamples("complete_fired_by_track_0",
                                   primaries_manager,
                                   "signal_antimu_complete",
                                   "fired_by_antimuon")

signal_antimuons = primaries_manager.GetSample("signal_mu+")
signal_antimuons_complete = primaries_manager.GetSample("signal_antimu_complete")

signal_neutrons = primaries_manager.GetSample("signal_neutrons")
signal_neutrons_no_direction_change = primaries_manager.GetSample("signal_neutrons_no_direction_change")
signal_neutrons_w_hits = primaries_manager.GetSample("signal_neutrons_w_hits")
signal_neutrons_reconstructable = primaries_manager.GetSample("signal_neutrons_reconstructable")
signal_neutrons_reconstructable_complete = primaries_manager.GetSample("signal_neutrons_reconstructable_complete")

fired_by_reconstructable_neutron = fired_cells_manager.GetSample("fired_by_reconstructable_neutron")
fired_by_antimuons = fired_cells_manager.GetSample("fired_by_antimuon")

print(f'total number of events in fiducial volume {len(df)}, CCQE {df_manager.GetSample("signal").size()/len(df)*1e2} [%]')
print(f'total simulated neutrons from signal {signal_neutrons.size()}')
print(f'with no directions before ECAL {signal_neutrons_no_direction_change.size()}, {signal_neutrons_no_direction_change.size()/signal_neutrons.size() * 1e2} [%]')
print(f'... and at least 1 hit in ECAL {signal_neutrons_reconstructable.size()}, {signal_neutrons_reconstructable.size()/signal_neutrons.size() * 1e2} [%]')
print(f'... that produce complete cells {signal_neutrons_reconstructable_complete.size()}, {signal_neutrons_reconstructable_complete.size()/signal_neutrons.size() * 1e2} [%]')

all_hists = []

print("Neutron True Kinetic energy and tof (signal event)")
plotter = MultiPlotter(nrows=1, ncols=3, figsize=(30, 10))
# First histogram: Neutron True Kinetic Energy
plotter.plot_hist(
    data = signal_neutrons.dataframe['E_kin'],
    bins = np.arange(0, 800, 10),
    color = 'blue',
    xlabel = "[MeV]",
    ylabel = "Counts",
    log_scale = True
)
plotter.axes[0].set_title("Neutron True Kinetic energy (signal event)", fontsize=15)
plotter.next_plot()

# Second histogram: Neutron True Beta
plotter.plot_hist(
    data = signal_neutrons.dataframe['beta'],
    bins = np.arange(0, 1, 0.02),
    color = 'blue',
    xlabel = "Beta",
    ylabel = "Counts"
)
plotter.axes[1].set_title("Neutron True Beta (signal event)", fontsize=15)

plotter.next_plot()

plotter.plot_hist2d(
    x = signal_neutrons.dataframe['E_kin'],
    y = signal_neutrons.dataframe['beta'],
    bins_x = np.arange(0, 800, 10),
    bins_y = np.arange(0, 1, 0.01),
    xlabel = "Neutron Kinetic Energy [MeV]",
    ylabel = "Neutron beta",
    # log_scale = True
)

all_hists.append(plotter)

print("Neutron angle of deviation wrt initial direction")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 8), suptitle = "Neutrons from signal")

plotter.plot_hist(
    data = signal_neutrons.dataframe["HasChangedDirection"],
    bins = [0,1,2],
    color = 'blue',
    xlabel = 'Has Changed Direction',
    ylabel = "counts",
)


plotter.next_plot()

plotter.plot_hist(
    data = primaries_manager.GetSample("signal_neutrons_no_direction_change").dataframe["E_kin"],
    bins = np.arange(0, 800, 20),
    color = 'red',
    label = 'no interaction before ECAL',
)

plotter.plot_hist(
    data = primaries_manager.GetSample("signal_neutrons_w_direction_change").dataframe["E_kin"],
    bins = np.arange(0, 800, 20),
    color = 'blue',
    label = 'interaction before ECAL',
)

plotter.axes[1].set_yscale("log")
plotter.axes[1].set_title("Neutron Kinetic energy (signal)", fontsize=15)
plotter.set_labels(xlabel="[MeV]", ylabel="Counts", fontsize=15)

plotter.next_plot(plot_legend=True)

all_hists.append(plotter)

plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 10), )

plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.AddInfo(signal_neutrons_reconstructable_complete, 'E_kin'),
    y = fired_by_reconstructable_neutron.dataframe['Fired_Cells_adc1'] + fired_by_reconstructable_neutron.dataframe['Fired_Cells_adc2'],
    bins_x = np.arange(0, 500, 10),
    bins_y = np.arange(0, 200, 5),
    xlabel = "Neutron Kinetic Energy [MeV]",
    ylabel = "ADC1 + ADC2 count on cell PMTs",
    log_scale = True
)

plotter.next_plot()

plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.AddInfo(signal_neutrons_reconstructable_complete, 'E_kin'),
    y = fired_by_reconstructable_neutron.AddInfo(signal_neutrons_reconstructable_complete, 'PrimariesEDepECAL'),
    bins_x = np.arange(0, 500, 10),
    bins_y = np.arange(0, 100, 2.5),
    xlabel = "Neutron Kinetic Energy [MeV]",
    ylabel = "Energy Deposit in Cell Fibers (edepsim) [MeV]",
    log_scale = True
)

all_hists.append(plotter)

print("Comparison Time of flight to ECAL for neutrons and muons")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(20, 6))

# First histogram: Time of flight to ECAL for neutrons and muons
plotter.plot_hist(
    data = signal_neutrons_w_hits.dataframe['PrimariesFirstHitECAL_t'],
    bins = np.arange(0, 200, 1),
    label = 'neutron',
    color = 'blue',
    xlabel = "[ns]",
    ylabel = "Counts",
    xlim = (0, 200),
    ticks = np.arange(0, 200, 10)
)
plotter.plot_hist(
    data=signal_antimuons.dataframe.loc[signal_neutrons_w_hits.SampleIndex()]['PrimariesFirstHitECAL_t'],
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
    data = signal_antimuons.dataframe.loc[signal_neutrons_w_hits.SampleIndex()]['PrimariesFirstHitECAL_t'].values - signal_neutrons_w_hits.dataframe['PrimariesFirstHitECAL_t'].values,
    bins = np.arange(-100, 10, 1),
    xlabel = "[ns]",
    ylabel = "Counts",
    xlim = (-100, 10),
    ticks = np.arange(-100, 10, 10)
)
plotter.axes[1].set_title("muon tof - neutron tof", fontsize=15)

all_hists.append(plotter)

print("Neutron ECAL hits comparison: true (EDEPSIM) vs reconstructed using cell TDCs")
plotter = MultiPlotter(nrows=1, ncols=3, figsize=(30, 10), suptitle="Neutron ECAL hits components: true (EDEPSIM) vs reconstructed using cell TDCs")

# First 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_x'],
    y = fired_by_reconstructable_neutron.dataframe['Cell_Reconstructed_hit_x'],
    bins_x = np.arange(-2000, 2000, 25),
    bins_y = np.arange(-2000, 2000, 25),
    log_scale = True,
    xlabel = r'hit $x_{true}$ (edepsim)',
    ylabel = r'hit $x_{reco}$ from TDCs of cell',
)
plotter.next_plot()

# Second 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_y'],
    y = fired_by_reconstructable_neutron.dataframe['Cell_Reconstructed_hit_y'],
    bins_x = np.arange(-4500, 100, 25),
    bins_y = np.arange(-4500, 100, 25),
    log_scale = True,
    xlabel = r'hit $y_{true}$ (edepsim)',
    ylabel = r'hit $y_{reco}$ from TDCs of cell',
)
plotter.next_plot()

# Third 2D histogram plot with logarithmic scale
plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_z'],
    y = fired_by_reconstructable_neutron.dataframe['Cell_Reconstructed_hit_z'],
    bins_x = np.arange(22000, 22800, 25),
    bins_y = np.arange(22000, 22800, 25),
    log_scale = True,
    xlabel = r'hit $z_{true}$ (edepsim)',
    ylabel = r'hit $z_{reco}$ from TDCs of cell',
)

all_hists.append(plotter)

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

print("Neutron tof to ECAL (signal, complete cells) comparison: true, prediction, reconstructed")
plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 8), suptitle="Neutron from signal, complete cells")

plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_t'],
    y = fired_by_reconstructable_neutron.dataframe['ExpectedNeutronHit_t'],
    bins_x = np.arange(0,40,0.2),
    bins_y = np.arange(0,40,0.2),
    log_scale = True,
    # label='',
    xlabel = 'True neutron tof (edepsim) [ns]',
    ylabel = r'Predicted neutron tof from $\mu^+$ [ns]',
    xlim = [0,40],
    ylim = [0,40]
)

plotter.axes[plotter.current_ax].plot([0, 40], [0, 40], color='red', linestyle='--', linewidth=2)

plotter.next_plot()

plotter.plot_hist2d(
    x = fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_t'],
    y = fired_by_reconstructable_neutron.dataframe['Cell_Reconstructed_hit_t'],
    bins_x = np.arange(0,40,0.2),
    bins_y = np.arange(0,40,0.2),
    log_scale = True,
    xlabel = 'True neutron neutron tof (edepsim) [ns]',
    ylabel = 'Reconstructed neutron tof from cell [ns]',
    xlim = [0,40],
    ylim = [0,40]
)

plotter.axes[plotter.current_ax].plot([0, 40], [0, 40], color='red', linestyle='--', linewidth=2)

all_hists.append(plotter)

plotter = MultiPlotter(nrows=1, ncols=2, figsize=(25, 8), suptitle="Neutron from signal, complete cells")

plotter.plot_hist(
    fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_t'] - fired_by_reconstructable_neutron.dataframe['ExpectedNeutronHit_t'],
    bins = np.arange(-30,30,0.3),
    xlabel = "[ns]",
)

plotter.axes[plotter.current_ax].set_title(r"$\mu^+$ tof true - predicted", fontsize = 20)

plotter.next_plot()

plotter.plot_hist(
    fired_by_reconstructable_neutron.dataframe['Fired_Cell_true_hit1_t'] - fired_by_reconstructable_neutron.dataframe['Cell_Reconstructed_hit_t'],
    bins = np.arange(-30,30,0.3),
    xlabel = "[ns]",
)
plotter.axes[plotter.current_ax].set_title(r"$\mu^+$ tof true - reco", fontsize = 20)

all_hists.append(plotter)

print("Muon tof to ECAL (signal, complete cells) comparison: true, reconstructed")
plotter = MultiPlotter(nrows=1, ncols=1, figsize=(10, 8), suptitle=r"$\mu^+$ from signal , complete cells")

plotter.plot_hist2d(
    x=fired_by_antimuons.dataframe['Fired_Cell_true_hit1_t'],
    y=fired_by_antimuons.dataframe['Cell_Reconstructed_hit_t'],
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
                 
MultiPlotter().save_multiple_figures_to_pdf('/storage/gpfs_data/neutrino/users/gi/sand-physics/scratch/multiple_plots.pdf', all_hists)