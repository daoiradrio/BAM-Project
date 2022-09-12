import os
import random
import warnings

import numpy as np
import plotly.graph_objs as go

from lobsterpy.structuregraph.graph import LobsterGraph
from itertools import product, permutations

from transformator import Transformator
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer



def create_plot(lobstergraph: LobsterGraph):
    cells = []

    atom_number = []

    node_x = []
    node_y = []
    node_z = []
    coords = []

    axis_x = []
    axis_y = []
    axis_z = []
    axes = []

    edges = []

    structure = lobstergraph.sg.structure

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

    cart_crystal_axis_mat = np.stack((x, y, z), axis=-1)

    """
    tol = 0.01
    for node1, node2, data in lobstergraph.sg.graph.edges.data():
        frac_coord1 = lobstergraph.sg.structure.frac_coords[node1]
        #coord1 = np.dot(cart_crystal_axis_mat, frac_coord1)
        frac_coord2 = lobstergraph.sg.structure.frac_coords[node2] + data["to_jimage"]
        #coord2 = np.dot(cart_crystal_axis_mat, frac_coord2)
        if (-tol <= frac_coord2[0] <= 1+tol) and \
           (-tol <= frac_coord2[1] <= 1+tol) and \
           (-tol <= frac_coord2[2] <= 1+tol):
            coord1 = np.dot(cart_crystal_axis_mat, frac_coord1)
            coord2 = np.dot(cart_crystal_axis_mat, frac_coord2)
            edges.append((coord1, coord2))
    """

    new_coords = []
    frac_coords = []
    new_frac_coords = []
    atoms = dict()
    eq_atoms = dict()
    k = len(structure.frac_coords)
    for i, frac_coord in enumerate(structure.frac_coords.copy()):
        frac_coords.append(frac_coord)
        coord = np.dot(cart_crystal_axis_mat, frac_coord)
        coords.append(coord)
        atoms[i] = {
            "element": str(structure[i].specie),
            "number": structure[i].specie.number,
        }
        eq_atoms[i] = []
        indices0 = []
        indices1 = []
        for j, c in enumerate(frac_coord):
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
                new_frac_coord = frac_coord + a0 + a1
                new_frac_coords.append(new_frac_coord)
                new_coord = np.dot(cart_crystal_axis_mat, new_frac_coord)
                new_coords.append(new_coord)
                atoms[k] = {
                    "element": str(structure[i].specie),
                    "number": structure[i].specie.number
                }
                eq_atoms[i].append(k)
                k += 1
    coords += new_coords
    frac_coords += new_frac_coords
    cells.append(coords)

    tol = 0.01
    for node1, node2, data in lobstergraph.sg.graph.edges.data():
        frac_coord1 = lobstergraph.sg.structure.frac_coords[node1]
        frac_coord2 = lobstergraph.sg.structure.frac_coords[node2] + data["to_jimage"]
        if (-tol <= frac_coord2[0] <= 1+tol) and \
           (-tol <= frac_coord2[1] <= 1+tol) and \
           (-tol <= frac_coord2[2] <= 1+tol):
            coord1 = np.dot(cart_crystal_axis_mat, frac_coord1)
            coord2 = np.dot(cart_crystal_axis_mat, frac_coord2)
            edges.append((coord1, coord2))
        for eq_atom1 in eq_atoms[node1]:
            start = frac_coords[eq_atom1]
            shift = start - frac_coord1
            end = frac_coord2 + shift
            if (-tol <= end[0] <= 1+tol) and \
               (-tol <= end[1] <= 1+tol) and \
               (-tol <= end[2] <= 1+tol):
                start = np.dot(cart_crystal_axis_mat, start)
                end = np.dot(cart_crystal_axis_mat, end)
                edges.append((start, end))

    #"""
    data = list(lobstergraph.sg.graph.edges.data())

    xs = [vec["to_jimage"][0] for _, _, vec in data]
    x_min = min(xs)
    x_max = max(xs)

    ys = [vec["to_jimage"][1] for _, _, vec in data]
    y_min = min(ys)
    y_max = max(ys)

    zs = [vec["to_jimage"][2] for _, _, vec in data]
    z_min = min(zs)
    z_max = max(zs)

    counter = len(atoms)
    for i, (dim_min, dim_max) in enumerate([(x_min, x_max), (y_min, y_max), (z_min, z_max)]):
        shift = np.array([0, 0, 0])
        shift[i] = 1.0
        new_axes = []
        new_edges = []
        new_cells = []
        for j in range(dim_min, 0):
            for start, end in axes:
                new_start = start + np.dot(cart_crystal_axis_mat, j * shift)
                new_end = end + np.dot(cart_crystal_axis_mat, j * shift)
                new_axes.append((new_start, new_end))
            for start, end in edges:
                new_start = start + np.dot(cart_crystal_axis_mat, j * shift)
                new_end = end + np.dot(cart_crystal_axis_mat, j * shift)
                new_edges.append((new_start, new_end))
            for cell in cells:
                new_coords = []
                for k, coord in enumerate(cell):
                    new_coord = coord + np.dot(cart_crystal_axis_mat, j * shift)
                    new_coords.append(new_coord)
                    atoms[counter] = {
                        "element": atoms[k]["element"],
                        "number": atoms[k]["number"]
                    }
                    counter += 1
                new_cells.append(new_coords)
        for j in range(1, dim_max+1):
            for start, end in axes:
                new_start = start + np.dot(cart_crystal_axis_mat, j * shift)
                new_end = end + np.dot(cart_crystal_axis_mat, j * shift)
                new_axes.append((new_start, new_end))
            for start, end in edges:
                new_start = start + np.dot(cart_crystal_axis_mat, j * shift)
                new_end = end + np.dot(cart_crystal_axis_mat, j * shift)
                new_edges.append((new_start, new_end))
            for cell in cells:
                new_coords = []
                for k, coord in enumerate(cell):
                    new_coord = coord + np.dot(cart_crystal_axis_mat, j * shift)
                    new_coords.append(new_coord)
                    atoms[counter] = {
                        "element": atoms[k]["element"],
                        "number": atoms[k]["number"]
                    }
                    counter += 1
                new_cells.append(new_coords)
        axes += new_axes
        edges += new_edges
        cells += new_cells
    #"""

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
        #start = np.dot(cart_crystal_axis_mat, new_edges[i])
        #end = np.dot(cart_crystal_axis_mat, new_edges[i+1])
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

    for cell in cells:
        for i, coord in enumerate(cell):
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



def get_random_structure() -> str:
    filepath = os.path.expanduser("~/automationresults")
    mps = os.listdir(filepath)
    rand_index = random.randrange(len(mps))
    dir = mps[rand_index]
    path = os.path.join(filepath, dir)
    print(f"RANDOMLY CHOSEN STRUCTURE: {dir}")
    return path



def get_chosen_structure(file: str = None) -> str:
    if file is None:
        path = get_random_structure()
    else:
        path = os.path.join(os.path.expanduser("~/automationresults"), file)
    return path



warnings.filterwarnings(action='ignore')

#dir = "mp-10143/"
dir = "mp-510401"
#dir = "mp-2384"
path = get_chosen_structure(dir)
#path = get_random_structure()
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



"""
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.plotter import CohpPlotter

dir = get_random_structure()
COHPCAR_path = os.path.join(dir, "COHPCAR.lobster")
POSCAR_path = os.path.join(dir, "POSCAR")

completecohp = CompleteCohp.from_file(
    fmt="LOBSTER", filename=COHPCAR_path, structure_file=POSCAR_path
)
#"""