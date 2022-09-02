import os
import warnings
import random

import numpy as np
import plotly.graph_objs as go

from lobsterpy.structuregraph.graph import LobsterGraph
from itertools import product, permutations
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from sklearn.neighbors import NearestNeighbors



def create_plot(structuregraph: LobsterGraph):
    atom_number = []

    node_x = []
    node_y = []
    node_z = []
    coords = []

    axis_x = []
    axis_y = []
    axis_z = []
    axes = []

    edge_x = []
    edge_y = []
    edge_z = []
    edges = []

    sga = SpacegroupAnalyzer(structuregraph.sg.structure)
    structure = sga.get_conventional_standard_structure()

    a = structure.lattice.a
    b = structure.lattice.b
    c = structure.lattice.c
    alpha = structure.lattice.alpha
    beta = structure.lattice.beta
    gamma = structure.lattice.gamma
    alpha_rad = np.deg2rad(alpha)
    beta_rad = np.deg2rad(beta)
    gamma_rad = np.deg2rad(gamma)
    origin = np.array([0, 0, 0])

    # x-axis
    x = np.array([a, 0, 0])
    axes.append((origin, x))
    # y-axis
    y = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    axes.append((origin, y))
    # z-axis
    z = np.array([c * np.cos(beta_rad), c * np.cos(alpha_rad), c * np.sin(gamma_rad)])
    axes.append((origin, z))
    # x-axis parallel in y-direction
    axes.append((y, y + x))
    # x-axis parallel in z-direction
    axes.append((z, z + x))
    # x-axis parallel in yz-direction
    axes.append((y + z, y + z + x))
    # y-axis parallel in x-direction
    axes.append((x, x + y))
    # y-axis parallel in z-direction
    axes.append((z, z + y))
    # y-axis parallel in xz-direction
    axes.append((x + z, x + z + y))
    # z-axis parallel in x-direction
    axes.append((x, x + z))
    # z-axis parallel in y-direction
    axes.append((y, y + z))
    # z-axis parallel in xy-direction
    axes.append((x + y, x + y + z))

    for axis in axes:
        start, end = axis
        axis_x += [start[0], end[0], None]
        axis_y += [start[1], end[1], None]
        axis_z += [start[2], end[2], None]

    frac_coords = list(structure.frac_coords.copy())
    number_individual_atoms = len(frac_coords)
    cart_crystal_axis_mat = np.stack((x, y, z), axis=-1)
    atoms = dict()
    new_coords = []
    frac_coords = list(structure.frac_coords.copy())
    new_frac_coords = frac_coords.copy()
    for i, coord in enumerate(frac_coords):
        coords.append(np.dot(cart_crystal_axis_mat, coord))
        atoms[i] = dict()
        atoms[i]["element"] = structure[i].specie
        atoms[i]["number"] = structure[i].specie.number
        indices0 = []
        indices1 = []
        for j, c in enumerate(coord):
            if c < 0.01:
                indices0.append(j)
            elif c > 0.99:
                indices1.append(j)
        n0 = len(indices0)
        n1 = len(indices1)
        add0 = [np.zeros(shape=3)]
        add1 = [np.zeros(shape=3)]
        for j in range(1, n0 + 1):
            v = [1] * j + [0] * (n0 - n1 - j)
            ps0 = set(permutations(v))
            for p0 in ps0:
                a0 = np.zeros(shape=3)
                for index, permutation in zip(indices0, p0):
                    a0[index] = permutation
                add0.append(a0)
        for j in range(1, n1 + 1):
            v = [-1] * j + [0] * (n1 - n0 - j)
            ps1 = set(permutations(v))
            for p1 in ps1:
                a1 = np.zeros(shape=3)
                for index, permutation in zip(indices1, p1):
                    a1[index] = permutation
                add1.append(a1)
        adds = product(add0, add1)
        for a0, a1 in adds:
            if np.linalg.norm(a0 + a1) > 0:
                new_frac_coord = coord.copy() + a0 + a1
                new_frac_coords.append(new_frac_coord)
                new_coord = np.array(coord.copy()) + a0 + a1
                new_coord = np.dot(cart_crystal_axis_mat, new_coord)
                new_coords.append(new_coord)
                m = number_individual_atoms - 1 + len(new_coords)
                atoms[m] = dict()
                atoms[m]["element"] = structure[i].specie
                atoms[m]["number"] = structure[i].specie.number
    coords += new_coords
    frac_coords += new_frac_coords

    #***#
    nodes = []
    elements = dict()
    edges = []
    for node in structuregraph.sg.graph.nodes.data():
        element = node[1]["specie"]
        coord_env = node[1]["properties"]["env"]
        nodes.append({"element": element, "coordination_environment": coord_env,})
        elements[element] = []
    for i, coord in enumerate(coords):
        element = atoms[i]["element"]
        element = str(element)
        elements[element].append(coord)
    for node1, node2, data in structuregraph.sg.graph.edges.data():
        element1 = atoms[node1]["element"]
        element1 = str(element1)
        env1 = nodes[node1]["coordination_environment"]
        env1 = int(env1[-1]) # das ist nicht im allgemeinen nötig
        for coord in elements[element1]:
            element2 = nodes[node2]["element"]
            #env2 = atoms[node2]["coordination_environment"]
            #env2 = int(env2[-1]) # das ist nicht im allgemeinen nötig
            vec = elements[element2].copy()
            vec.append(coord)
            vec = np.stack(vec)
            knn = NearestNeighbors(n_neighbors=env1+1)
            knn.fit(vec)
            distance_mat, neighbours_mat = knn.kneighbors(vec)
            for neighbor in neighbours_mat[-1][1:]:
                edges.append((coord, elements[element2][neighbor]))
    #***#

    for i, coord in enumerate(coords):
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_z.append(coord[2])
        atom_number.append(atoms[i]["number"])

    for start, end in edges:
        edge_x += [start[0], end[0], None]
        edge_y += [start[1], end[1], None]
        edge_z += [start[2], end[2], None]

    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
        showspikes=False
    )

    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(axis),
            yaxis=dict(axis),
            zaxis=dict(axis),
        ),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        hovermode="closest",
        height=850
    )

    fig = go.Figure(layout=layout)


    for start, end in edges:
        fig.add_trace(
            go.Scatter3d(
                x=[start[0], end[0], None],
                y=[start[1], end[1], None],
                z=[start[2], end[2], None],
                mode="lines",
                line={
                    "width": 2,
                    "color": "black"
                },
                hoverinfo="none",
                customdata=[1, 1, None]
            )
        )


    for axis in axes:
        start, end = axis
        axis_x += [start[0], end[0], None]
        axis_y += [start[1], end[1], None]
        axis_z += [start[2], end[2], None]

    axes_trace = go.Scatter3d(
        x=axis_x,
        y=axis_y,
        z=axis_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="grey", width=1)
    )

    fig.add_trace(axes_trace)

    for i, coord in enumerate(coords):
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_z.append(coord[2])
        atom_number.append(atoms[i]["number"])

    node_trace = go.Scatter3d(
        x=node_x,
        y=node_y,
        z=node_z,
        mode="markers",
        hoverinfo="none",
        marker=dict(
            symbol="circle",
            size=6,
            color=atom_number,
            colorscale="Viridis",
            line=dict(color="rgb(50,50,50)", width=0.5),
        ),
    )

    fig.add_trace(node_trace)

    #return fig
    fig.show()



warnings.filterwarnings(action='ignore')

#dir = "mp-10143/"
#crystalsystem = "cubic"
filepath = os.path.expanduser("~/automationresults")
mps = os.listdir(filepath)
rand_index = random.randrange(len(mps))
dir = mps[rand_index]
path = os.path.join(filepath, dir)
#path = os.path.join(os.path.expanduser("~/RAW_files_phonon/"), crystalsystem, dir)
path_to_poscar = os.path.join(path, "POSCAR")
path_to_charge = os.path.join(path, "CHARGE.lobster")
path_to_icobilist = os.path.join(path, "ICOBILIST.lobster")
path_to_icooplist = os.path.join(path, "ICOOPLIST.lobster")
path_to_icohplist = os.path.join(path, "ICOHPLIST.lobster")
path_to_cohpcar = os.path.join(path, "COHPCAR.lobster")
path_to_madelung = os.path.join(path, "MadelungEnergies.lobster")

print("hier0")
testgraph = LobsterGraph(
    path_to_poscar=path_to_poscar,
    path_to_charge=path_to_charge,
    path_to_icobilist=path_to_icobilist,
    path_to_icooplist=path_to_icooplist,
    path_to_icohplist=path_to_icohplist,
    path_to_cohpcar=path_to_cohpcar,
    path_to_madelung=path_to_madelung,
    which_bonds="all",
    #start=-2,
    add_additional_data_sg=True
)

print("hier1")
plotfig = create_plot(testgraph)
print("hier2")

#for item in testgraph.sg.graph.edges.data():
#    print(item)
#    print()

#print(testgraph.sg.structure.frac_coords)

#print()

#sga = SpacegroupAnalyzer(testgraph.sg.structure)
#s = sga.get_conventional_standard_structure()
#print(s.frac_coords)

for item in testgraph.sg.graph.nodes.data():
    element = item[1]["specie"]
    coordination_environment = item[1]["properties"]["env"]
    #print(f"Element {element}, Coordination Environment {coordination_environment}")

structure = testgraph.sg.structure
sga = SpacegroupAnalyzer(structure)
structure = sga.get_conventional_standard_structure()

