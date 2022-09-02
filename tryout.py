import os
import warnings
import numpy as np
import plotly.graph_objs as go
from lobsterpy.structuregraph.graph import LobsterGraph
from itertools import product, permutations
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from main import UnitCell



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
    eq_atoms = dict()
    new_coords = []
    frac_coords = list(structure.frac_coords.copy())
    new_frac_coords = frac_coords.copy()
    for i, coord in enumerate(frac_coords):
        coords.append(np.dot(cart_crystal_axis_mat, coord))
        atoms[i] = dict()
        atoms[i]["element"] = structure[i].specie
        atoms[i]["number"] = structure[i].specie.number
        eq_atoms[i] = []
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
                eq_atoms[i].append(m)
    coords += new_coords
    frac_coords += new_frac_coords

    for i, coord in enumerate(coords):
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_z.append(coord[2])
        atom_number.append(atoms[i]["number"])

    limit = 0.01
    for node1, node2, data in structuregraph.sg.graph.edges.data():
        #start = frac_coords[node1]
        #end = frac_coords[node2] + data["to_jimage"]
        start = structuregraph.sg.structure.frac_coords[node1]
        end = structuregraph.sg.structure.frac_coords[node2] + data["to_jimage"]
        start_coord = np.dot(cart_crystal_axis_mat, start)
        end_coord = np.dot(cart_crystal_axis_mat, end)
        d0 = np.linalg.norm(start_coord - end_coord)
        if (-limit <= end[0] <= limit + 1) and \
                (-limit <= end[1] <= limit + 1) and \
                (-limit <= end[2] <= limit + 1):
            edges.append((start_coord, end_coord))
        for eq_atom1 in eq_atoms[node1]:
            shift = frac_coords[eq_atom1] - frac_coords[node1]
            new_end = end + shift
            if (-limit <= new_end[0] <= limit + 1) and \
                    (-limit <= new_end[1] <= limit + 1) and \
                    (-limit <= new_end[2] <= limit + 1):
                start_coord = coords[eq_atom1]
                end_coord = np.dot(cart_crystal_axis_mat, new_end)
                d1 = np.linalg.norm(start_coord - end_coord)
                if abs(d0 - d1) <= limit:
                    edges.append((start_coord, end_coord))

    """
    for start, end in self.edges:
        start_coord = self.coords[start]
        end_coord = self.coords[end]
        edge_x += [start_coord[0], end_coord[0], None]
        edge_y += [start_coord[1], end_coord[1], None]
        edge_z += [start_coord[2], end_coord[2], None]
    """

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

    return fig



warnings.filterwarnings(action='ignore')

dir = "mp-10143/"
path = os.path.join(os.path.expanduser("~/automationresults"), dir)
path_to_poscar = os.path.join(path, "POSCAR")
path_to_charge = os.path.join(path, "CHARGE.lobster")
path_to_icobilist = os.path.join(path, "ICOBILIST.lobster")
path_to_icooplist = os.path.join(path, "ICOOPLIST.lobster")
path_to_icohplist = os.path.join(path, "ICOHPLIST.lobster")
path_to_cohpcar = os.path.join(path, "COHPCAR.lobster")
path_to_madelung = os.path.join(path, "MadelungEnergies.lobster")

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

plotfig = create_plot(testgraph)

#for item in testgraph.sg.graph.edges.data():
#    print(item)
#    print()

#print(testgraph.sg.structure.frac_coords)

#print()

#sga = SpacegroupAnalyzer(testgraph.sg.structure)
#s = sga.get_conventional_standard_structure()
#print(s.frac_coords)

structure = testgraph.sg.structure
sga = SpacegroupAnalyzer(structure)
cell = UnitCell(sga.get_conventional_standard_structure())
cell.plot_extended_cell([
    #[1,1,0,0],
    #[1,-1,0,0],
    #[1,0,1,0],
    #[1,0,-1,0],
    #[1,0,0,1],
    #[1,0,0,-1]
])