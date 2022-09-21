import os
import random
import warnings

import numpy as np
import plotly.graph_objs as go

from lobsterpy.structuregraph.graph import LobsterGraph
from pymatgen.electronic_structure.cohp import CompleteCohp
from itertools import product, permutations

warnings.filterwarnings(action='ignore')



def get_primitive_cell(lobstergraph: LobsterGraph, completecohp: CompleteCohp) -> (list, dict, dict):
    """
    function for building first full primitive cell with edges

    :param lobstergraph: graph object containing information about nodes (atoms) and edges (bonds)
    :return cells: list of primitive cells, contains one cell after this function is finished
    :return atoms: dictionary of atoms in primitive cell, keys are enumerated atom numbers, values element symbol and
                   element number
    :return eq_atoms: dictionary of atoms equivalent the already existing ones in the graph object, keys are
                      enumerated atoms numbers, values are enumerated atom numbers of equivalent new atoms
    """

    # initialization block
    cell = {
        "atoms": {},
        "axes": [],
        "edges": [],
        "edge_properties": {
            "cohp_plot": [],
            "bond_length": [],
            "icobi": [],
            "icoop": [],
            "icohp": [],
            "icohp_bonding_perc": []
        },
    }
    structure = lobstergraph.sg.structure
    num_atoms = len(structure.frac_coords)
    eq_atoms = dict()
    tol = 0.01

    # iterate over all already existing nodes, save cartesian coordinates, fractional coordinates, element and
    # element number
    for i, frac_coord in enumerate(structure.frac_coords.copy()):
        cell["atoms"][i] = {
            "element": str(structure[i].specie),
            "number": structure[i].specie.number,
            "frac_coord": frac_coord,
        }
        eq_atoms[i] = []
        # for every node find fractional coordinates with entries 0 or 1 in some dimension, for these nodes equivalent
        # nodes/coordinates (due to translational symmetry) will be created
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
        # if there is more than one entry with 0 or 1 create all possible permutations for the translational shift
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
        # create new nodes/coordinates by adding translational shift vectors
        # save new fractional coordinates, cartesian coordinates, element and element number
        for a0, a1 in adds:
            if np.linalg.norm(a0 + a1) > 0:
                cell["atoms"][num_atoms] = {
                    "element": str(structure[i].specie),
                    "number": structure[i].specie.number,
                    "frac_coord": frac_coord + a0 + a1,
                }
                eq_atoms[i].append(num_atoms)
                num_atoms += 1

    # iterate over all already existing nodes
    for i, (node1, node2, data) in enumerate(lobstergraph.sg.graph.edges.data()):
        # extract edge properties
        cohp_data = completecohp.get_cohp_by_label(label=data["bond_label"]).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x = [spinup_cohps[i] + spindown_cohps[i] for i, _ in enumerate(spinup_cohps)]
        y = [energies[i] - fermi_energy for i, _ in enumerate(energies)]
        frac_coord1 = cell["atoms"][node1]["frac_coord"]
        frac_coord2 = cell["atoms"][node2]["frac_coord"] + data["to_jimage"]
        # check if edges lies within cell
        if (-tol <= frac_coord2[0] <= 1+tol) and \
           (-tol <= frac_coord2[1] <= 1+tol) and \
           (-tol <= frac_coord2[2] <= 1+tol):
            # add edge by fractional coordinates of connected nodes to cell
            cell["edges"].append((frac_coord1, frac_coord2))
            # add edge properties to cell
            cell["edge_properties"]["cohp_plot"].append((x, y))
            cell["edge_properties"]["bond_length"].append(data["bond_length"])
            cell["edge_properties"]["icobi"].append(data["ICOBI"])
            cell["edge_properties"]["icoop"].append(data["ICOOP"])
            cell["edge_properties"]["icohp"].append(data["ICOHP"])
            cell["edge_properties"]["icohp_bonding_perc"].append(data["ICOHP_bonding_perc"])
        # iterate over all edges that are equivalent to the current one due to translational symmetry
        for eq_atom in eq_atoms[node1]:
            start = cell["atoms"][eq_atom]["frac_coord"]
            shift = start - frac_coord1
            end = frac_coord2 + shift
            # check if edges lies within cell
            if (-tol <= end[0] <= 1+tol) and \
               (-tol <= end[1] <= 1+tol) and \
               (-tol <= end[2] <= 1+tol):
                # add edge by fractional coordinates of connected nodes to cell
                cell["edges"].append((start, end))
                # add edge properties to cell
                cell["edge_properties"]["cohp_plot"].append((x, y))
                cell["edge_properties"]["bond_length"].append(data["bond_length"])
                cell["edge_properties"]["icobi"].append(data["ICOBI"])
                cell["edge_properties"]["icoop"].append(data["ICOOP"])
                cell["edge_properties"]["icohp"].append(data["ICOHP"])
                cell["edge_properties"]["icohp_bonding_perc"].append(data["ICOHP_bonding_perc"])

    return cell



def get_primitive_supercell(
        lobstergraph: LobsterGraph, cell: dict, cart_crystal_axis_matrix: np.ndarray
) -> list:
    """
    function to build a primitive supercell from a primitive cell based on connectivity information in graph
    object ("to_jimage" vectors)

    :param lobstergraph: graph object containing information about nodes (atoms) and edges (bonds)
    :param cells: list containing primitive cell to duplicate in order to build supercell
    :param atoms: dictionary of atoms in primitive cell, keys are enumerated atom numbers, values element symbol and
                  element number
    :param cart_crystal_axis_matrix:
    :return: cells: list containing all primitive cells that make up primitive supercell
    """

    data = list(lobstergraph.sg.graph.edges.data())

    xs = [vec["to_jimage"][0] for _, _, vec in data]
    # maximum x-component of all "to_jimage" vectors determines size of supercell in +x direction
    x_min = min(xs)
    # minimum x-component of all "to_jimage" vectors determines size of supercell in -x direction
    x_max = max(xs)

    ys = [vec["to_jimage"][1] for _, _, vec in data]
    # maximum y-component of all "to_jimage" vectors determines size of supercell in +y direction
    y_min = min(ys)
    # minimum y-component of all "to_jimage" vectors determines size of supercell in -y direction
    y_max = max(ys)

    zs = [vec["to_jimage"][2] for _, _, vec in data]
    # maximum z-component of all "to_jimage" vectors determines size of supercell in +z direction
    z_min = min(zs)
    # minimum z-component of all "to_jimage" vectors determines size of supercell in -z direction
    z_max = max(zs)

    # iterate over x, y, z dimension
    num_atoms = len(cell["atoms"])
    cells = [cell]
    for i, (dim_min, dim_max) in enumerate([(x_min, x_max), (y_min, y_max), (z_min, z_max)]):
        # create new cell by shifting existing one in x, y or z direction
        shift = np.array([0, 0, 0])
        shift[i] = 1.0
        new_cells = []
        for cell in cells:
            # repeat every cell x_min/y_min/z_min times in negative direction and x_max/y_max/z_max times in
            # positive direction
            for j in [k for k in range(dim_min, dim_max+1) if k != 0]:
                new_cell = {
                    "atoms": {},
                    "axes": [],
                    "edges": [],
                    "edge_properties": {
                        "cohp_plot": [],
                        "bond_length": [],
                        "icobi": [],
                        "icoop": [],
                        "icohp": [],
                        "icohp_bonding_perc": []
                    },
                }
                # add axes to new cell
                for start, end in cell["axes"]:
                    new_start = start + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_end = end + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_cell["axes"].append((new_start, new_end))
                # add edges (bonds) and edge properties (bond properties) to new cell
                # (the latter based on equivalence)
                for l, (start, end) in enumerate(cell["edges"]):
                    new_start = start + j * shift
                    new_end = end + j * shift
                    new_cell["edges"].append((new_start, new_end))
                    new_cell["edge_properties"]["cohp_plot"].append(cell["edge_properties"]["cohp_plot"][l])
                    new_cell["edge_properties"]["bond_length"].append(cell["edge_properties"]["bond_length"][l])
                    new_cell["edge_properties"]["icobi"].append(cell["edge_properties"]["icobi"][l])
                    new_cell["edge_properties"]["icoop"].append(cell["edge_properties"]["icoop"][l])
                    new_cell["edge_properties"]["icohp"].append(cell["edge_properties"]["icohp"][l])
                    new_cell["edge_properties"]["icohp_bonding_perc"].append(
                        cell["edge_properties"]["icohp_bonding_perc"][l]
                    )
                # add nodes (atoms) to new cell as coordinates
                for atom in cell["atoms"].values():
                    new_cell["atoms"][num_atoms] = {
                        "element": atom["element"],
                        "number": atom["number"],
                        "frac_coord": atom["frac_coord"] + j * shift,
                    }
                    num_atoms += 1
                new_cells.append(new_cell)
        cells += new_cells
    return cells



def create_plot(lobstergraph: LobsterGraph, completecohp: CompleteCohp) -> go.Figure:
    """
    Creation of an interactive 3D plot of a compound's primitive supercell, containing information about site and
    bond properties.

    :param lobstergraph: LobsterGraph object, contains information about connectivity, site and bond properties in a
                         graph-like manner
    :return: fig: visualization of primitive supercell by 3D scatter plot
    """

    # initialization block
    atom_number = []

    node_x = []
    node_y = []
    node_z = []

    axis_x = []
    axis_y = []
    axis_z = []
    axes = []

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

    # cartesian x-axis
    x = np.array([a, 0, 0])
    axes.append((origin, x))
    # cartesian y-axis
    y = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    axes.append((origin, y))
    # cartesian z-axis
    z = np.array([c * np.cos(beta_rad), c * np.cos(alpha_rad), c * np.sin(gamma_rad)])
    axes.append((origin, z))
    # cartesian x-axis parallel in y-direction
    axes.append((y, y + x))
    # cartesian x-axis parallel in z-direction
    axes.append((z, z + x))
    # cartesian x-axis parallel in yz-direction
    axes.append((y + z, y + z + x))
    # cartesian y-axis parallel in x-direction
    axes.append((x, x + y))
    # cartesian y-axis parallel in z-direction
    axes.append((z, z + y))
    # cartesian y-axis parallel in xz-direction
    axes.append((x + z, x + z + y))
    # cartesian z-axis parallel in x-direction
    axes.append((x, x + z))
    # cartesian z-axis parallel in y-direction
    axes.append((y, y + z))
    # cartesian z-axis parallel in xy-direction
    axes.append((x + y, x + y + z))

    # matrix to transform fractional to cartesian coordinates
    cart_crystal_axis_matrix = np.stack((x, y, z), axis=-1)

    cell = get_primitive_cell(lobstergraph, completecohp)
    cell["axes"] = axes
    cells = get_primitive_supercell(lobstergraph, cell, cart_crystal_axis_matrix)

    # layout of plot axes
    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
        showspikes=False
    )

    # overall layout of plot
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
        height=820,
        width=1195
    )

    fig = go.Figure(layout=layout)

    # build up structure plot by including the primitive cells in plot
    for cell in cells:
        # collect node data for plot
        for atom in cell["atoms"].values():
            coord = np.dot(cart_crystal_axis_matrix, atom["frac_coord"])
            node_x.append(coord[0])
            node_y.append(coord[1])
            node_z.append(coord[2])
            atom_number.append(atom["number"])
        # add edges (bonds) to plot, the bond properties are saved as custom hover data
        # every edge is included as separate trace to make them separately accessible for the hover events
        for j, (start, end) in enumerate(cell["edges"]):
            start = np.dot(cart_crystal_axis_matrix, start)
            end = np.dot(cart_crystal_axis_matrix, end)
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
                    customdata=[
                        [
                            cell["edge_properties"]["cohp_plot"][j],
                            cell["edge_properties"]["bond_length"][j],
                            cell["edge_properties"]["icobi"][j],
                            cell["edge_properties"]["icoop"][j],
                            cell["edge_properties"]["icohp"][j],
                            cell["edge_properties"]["icohp_bonding_perc"][j]
                        ],
                        [
                            cell["edge_properties"]["cohp_plot"][j],
                            cell["edge_properties"]["bond_length"][j],
                            cell["edge_properties"]["icobi"][j],
                            cell["edge_properties"]["icoop"][j],
                            cell["edge_properties"]["icohp"][j],
                            cell["edge_properties"]["icohp_bonding_perc"][j]
                        ],
                        None
                    ]
                )
            )
        # collect axes data for the plot
        for start, end in cell["axes"]:
            axis_x += [start[0], end[0], None]
            axis_y += [start[1], end[1], None]
            axis_z += [start[2], end[2], None]

    # add axes of primitive cells to plot
    axes_trace = go.Scatter3d(
        x=axis_x,
        y=axis_y,
        z=axis_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="grey", width=1)
    )
    fig.add_trace(axes_trace)

    # add nodes (atoms) with color depending on element to plot
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



def get_structure_plot(
        path_to_poscar: str,
        path_to_charge: str,
        path_to_icobilist: str,
        path_to_icooplist: str,
        path_to_icohplist: str,
        path_to_cohpcar: str,
        path_to_madelung: str
):
    lobstergraph = LobsterGraph(
        path_to_poscar=path_to_poscar,
        path_to_charge=path_to_charge,
        path_to_icobilist=path_to_icobilist,
        path_to_icooplist=path_to_icooplist,
        path_to_icohplist=path_to_icohplist,
        path_to_cohpcar=path_to_cohpcar,
        path_to_madelung=path_to_madelung,
        which_bonds="all",
        # start=-2,
        add_additional_data_sg=True
    )

    completecohp = CompleteCohp.from_file(
        fmt="LOBSTER", filename=path_to_cohpcar, structure_file=path_to_poscar
    )

    return create_plot(lobstergraph, completecohp)


def get_dummy_cohp_plot() -> go.Figure:
    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        visible=False,
        title="",
        showspikes=False
    )

    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=axis,
            yaxis=axis,
        ),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        height=440,
        width=600,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False)
    )

    return go.Figure(layout=layout)



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