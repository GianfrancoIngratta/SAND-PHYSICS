# import seaborn as sns
import time
start = time.time()
import matplotlib.pyplot as plt
import numpy as np
# import particle
# import ROOT
# import sys
# import pandas as pd
# import glob
# sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-physics/notebook/python_tools')
# from ROOT_tools import ROOT_tools
# from MultiPlotter import MultiPlotter
# from ROOT2Pandas import Converter
# from SampleManager import Sample, Manager
# tool = ROOT_tools()
import plotly.graph_objects as go
print("Import time:", time.time() - start)


def plot_SAND3D(fig, opacity=0.4):
    # Parametri del cilindro
    sand_radius = 2000  # Raggio del cilindro
    sand_length = 3300  # Lunghezza del cilindro lungo X
    sand_center = [0, -2384.73, 23910]  # Il centro del cilindro nel piano YZ
    
    # Genera i punti per la superficie laterale del cilindro
    theta = np.linspace(0, 2 * np.pi, 100)  # Angolo attorno al cilindro
    x_sand = np.linspace(-sand_length / 2, sand_length / 2, 50) + sand_center[0]  # Lunghezza cilindro (asse X)
    theta, x = np.meshgrid(theta, x_sand)  # Griglia 2D per la superficie laterale

    # Coordinate nel piano ZY per la superficie laterale
    z_sand = sand_radius * np.cos(theta) + sand_center[2]  # Coordinata Z della circonferenza del cilindro
    y_sand = sand_radius * np.sin(theta) + sand_center[1]  # Coordinata Y della circonferenza del cilindro

    # Genera i punti per i tappi (dischi) alle estremità
    theta_cap = np.linspace(0, 2 * np.pi, 100)  # Angolo per i tappi
    r_cap = np.linspace(0, sand_radius, 50)  # Raggio dei tappi
    theta_cap, r = np.meshgrid(theta_cap, r_cap)  # Griglia 2D per i tappi

    # Coordinate dei tappi nel piano ZY, con X fisso a +sand_length/2 e -sand_length/2
    x_cap_top = np.full_like(theta_cap, sand_length / 2) + sand_center[0]  # Tappo superiore, X fisso
    x_cap_bottom = np.full_like(theta_cap, -sand_length / 2) + sand_center[0]  # Tappo inferiore, X fisso

    # Coordinate ZY per i tappi (dischi)
    z_cap = r * np.cos(theta_cap) + sand_center[2]  # Coordinata Z per i tappi
    y_cap = r * np.sin(theta_cap) + sand_center[1]  # Coordinata Y per i tappi

    # Plotta la superficie laterale del cilindro
    fig.add_surface(x=x, y=z_sand, z=y_sand, colorscale='Blues', opacity=opacity, showscale=False)

    # Plotta il tappo superiore (a X = +sand_length/2)
    fig.add_surface(x=x_cap_top, y=z_cap, z=y_cap, colorscale='Blues', opacity=opacity, showscale=False)

    # Plotta il tappo inferiore (a X = -sand_length/2)
    fig.add_surface(x=x_cap_bottom, y=z_cap, z=y_cap, colorscale='Blues', opacity=opacity, showscale=False)

    # Imposta le etichette degli assi e il range, mantenendo il rapporto proporzionale tra gli assi
    fig.update_layout(
        scene=dict(
            xaxis_title='X',  # L'asse X rappresenta la lunghezza del cilindro
            yaxis_title='Z',  # L'asse Z rappresenta la circonferenza del cilindro nel piano ZY
            zaxis_title='Y',  # L'asse Y rappresenta la circonferenza del cilindro nel piano ZY
            aspectmode='manual',  # Mantenere il rapporto proporzionale tra gli assi
            aspectratio=dict(x=0.8, y=1, z=1),  # Regola il rapporto per rendere l'asse X più corto
            xaxis=dict(range=[-sand_length / 2 - 500, sand_length / 2 + 500]),  # Range asse X (lunghezza del cilindro)
            yaxis=dict(range=[sand_center[2] - sand_radius - 500, sand_center[2] + sand_radius + 500]),  # Range asse Z
            zaxis=dict(range=[sand_center[1] - sand_radius - 500, sand_center[1] + sand_radius + 500])  # Range asse Y
        ),
        width=800,
        height=600
    )

def GetExpectedTrajectory(vtx, direction, nof_points = 5000):
    direction = direction / np.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)
    x,y,z = [],[],[]
    for i in range(nof_points):
        x.append(vtx[0] + i*direction[0])
        y.append(vtx[1] + i*direction[1])
        z.append(vtx[2] + i*direction[2])
    return {"x":x,"y":y,"z":z}

def plot_track3D(fig, vector_x, vector_y, vector_z, label = '', color = 'blue', style = 'dash', line_width=3):
    fig.add_trace(go.Scatter3d(
        x = vector_x, 
        y = vector_z, 
        z = vector_y,
        mode='lines', 
        name=label, 
        line=dict(color=color, 
                  dash=style, 
                  width=line_width)  # Aggiunto il parametro width
    ))

def plot_points3D(fig, vector_x, vector_y, vector_z, label = '', color = 'blue', symbol = 'cross', size = 5):
    fig.add_trace(go.Scatter3d(
    x=vector_x, 
    y=vector_y, 
    z=vector_z,
    mode='markers',  # Cambia da 'lines' a 'markers'
    name=label, 
    marker=dict(
        color=color,
        symbol=symbol,  # Imposta il simbolo come croce
        size=size,  # Imposta la dimensione dei marker
        line=dict(width=2, color=color),
        opacity=0.8
    )
))

def plot_event3D(df, trajectories, fired_cells, event_name:str, plot_hit_reco = True, plot_hit_true = False):
    event = df.loc[event_name]
    print("Calculate expected neutron")
    predicted_neutron = GetExpectedTrajectory([event.Interaction_vtxX * 1e3, event.Interaction_vtxY * 1e3, event.Interaction_vtxZ * 1e3],
                      [event.PredictedNeutron_P3_GeVfX, event.PredictedNeutron_P3_GeVfY, event.PredictedNeutron_P3_GeVfZ])
    print("read event trajectories")
    event_trajectory = trajectories.loc[event_name]
    print("read fired wires")
    event_fired_cells = fired_cells.loc[event_name]
    compatible_cell = event_fired_cells.query("IsCompatible==1")

    track_ids_event = event_trajectory.trackid.unique()
    labels = {2112:"neutron", 22:"gamma", -13:"antimuons", -11:"positron", 11:"electron"}
    colors = {2112:"blue", 22:"violet", -13:"green", -11:"violet", 11:"violet"}

    print(event)

    # Crea il grafico interattivo con Plotly
    fig = go.Figure()

    for track_id in track_ids_event[track_ids_event!=-999]:

        plot_track3D(fig, event_trajectory.query(f"trackid == {track_id}").point_x, 
                          event_trajectory.query(f"trackid == {track_id}").point_y, 
                          event_trajectory.query(f"trackid == {track_id}").point_z, 
                          label=labels[event_trajectory.query(f"trackid == {track_id}").pdg.unique()[0]],
                          color=colors[event_trajectory.query(f"trackid == {track_id}").pdg.unique()[0]],
                          style='solid',
                          line_width=5)

    plot_track3D(fig, predicted_neutron['x'], 
                      predicted_neutron['y'], 
                      predicted_neutron['z'], 
                      label="predicted neutron",
                      color='blue',
                      style='dash',
                      line_width=5)

    # plot_track3D(fig, antimu_trj.point_x, 
    #                   antimu_trj.point_y, 
    #                   antimu_trj.point_z, 
    #                   label="true antimu",
    #                   color='green',
    #                   style='solid',
    #                   line_width=5)

    plot_points3D(fig, event_fired_cells.ExpectedNeutron_HitPosition_x_,
                       event_fired_cells.ExpectedNeutron_HitPosition_z_,
                       event_fired_cells.ExpectedNeutron_HitPosition_y_,
                       label='predicted hits',
                       color='blue',
                       size = 8)
    if(plot_hit_reco==1):
            plot_points3D(fig, event_fired_cells.Reconstructed_HitPosition_x,
                           event_fired_cells.Reconstructed_HitPosition_z,
                           event_fired_cells.Reconstructed_HitPosition_y,
                           label='ECAL hit reco',
                           color='orange',
                           size = 8)
    
    if(plot_hit_true == True):
            plot_points3D(fig, event_fired_cells.Fired_Cell_true_Hit_x,
                           event_fired_cells.Fired_Cell_true_Hit_z,
                           event_fired_cells.Fired_Cell_true_Hit_y,
                           label='ECAL hit true',
                           color='violet',
                           size = 8)

    plot_points3D(fig, compatible_cell.Reconstructed_HitPosition_x,
                       compatible_cell.Reconstructed_HitPosition_z,
                       compatible_cell.Reconstructed_HitPosition_y,
                       label='compatible cell',
                       color='red',
                       size = 8)

    plot_SAND3D(fig)

    # fig.write_html('grafico_interattivo.html')

    fig.update_layout(
        width=1e3,  # Imposta la larghezza
        height=1e3  # Imposta l'altezza
    )
    fig.show()

fig = go.Figure()

plot_SAND3D(fig)

fig.write_html('SAND3D.html')

