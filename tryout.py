import os
import random
import warnings

import numpy as np
import plotly.graph_objs as go

from lobsterpy.structuregraph.graph import LobsterGraph
from pymatgen.electronic_structure.cohp import CompleteCohp
from itertools import product, permutations



def get_primitive_cell(lobstergraph: LobsterGraph, completecohp: CompleteCohp) -> (list, dict, dict):
    """
    function for building first full primitive cell with edges

    :param lobstergraph:
    :param cart_crystal_axis_matrix:
    :return cells:
    :return frac_coords:
    :return atoms:
    :return eq_atoms:
    """

    # initialization block
    cells = [
        {
            "frac_coords": [],
            "edges": [],
            "axes": [],
            "edge_properties": {
                "cohp_plot": [],
                "bond_length": [],
                "icobi": [],
                "icoop": [],
                "icohp": [],
                "icohp_bonding_perc": []
            },
        }
    ]
    new_frac_coords = []
    atoms = dict()
    eq_atoms = dict()
    structure = lobstergraph.sg.structure
    num_atoms = len(structure.frac_coords)
    tol = 0.01

    # iterate over all already existing nodes, save cartesian coordinates, fractional coordinates, element and
    # element number
    for i, frac_coord in enumerate(structure.frac_coords.copy()):
        cells[0]["frac_coords"].append(frac_coord)
        atoms[i] = {
            "element": str(structure[i].specie),
            "number": structure[i].specie.number,
        }
        # for every node find fractional coordinates with entries 0 or 1 in some dimension, for these nodes equivalent
        # nodes/coordinates (due to translational symmetry) will be created
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
                new_frac_coord = frac_coord + a0 + a1
                new_frac_coords.append(new_frac_coord)
                atoms[num_atoms] = {
                    "element": str(structure[i].specie),
                    "number": structure[i].specie.number
                }
                eq_atoms[i].append(num_atoms)
                num_atoms += 1

    # add new cartesian coordinates to create first full primitive cell
    cells[0]["frac_coords"] += new_frac_coords

    # iterate over all already existing nodes
    for i, (node1, node2, data) in enumerate(lobstergraph.sg.graph.edges.data()):
        cohp_data = completecohp.get_cohp_by_label(label=data["bond_label"]).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x = [spinup_cohps[i] + spindown_cohps[i] for i, _ in enumerate(spinup_cohps)]
        y = [energies[i] - fermi_energy for i, _ in enumerate(energies)]
        frac_coord1 = cells[0]["frac_coords"][node1]
        frac_coord2 = cells[0]["frac_coords"][node2] + data["to_jimage"]
        # check if edges lies within cell
        if (-tol <= frac_coord2[0] <= 1+tol) and \
           (-tol <= frac_coord2[1] <= 1+tol) and \
           (-tol <= frac_coord2[2] <= 1+tol):
            cells[0]["edges"].append((frac_coord1, frac_coord2))
            cells[0]["edge_properties"]["cohp_plot"].append((x, y))
            cells[0]["edge_properties"]["bond_length"].append(data["bond_length"])
            cells[0]["edge_properties"]["icobi"].append(data["ICOBI"])
            cells[0]["edge_properties"]["icoop"].append(data["ICOOP"])
            cells[0]["edge_properties"]["icohp"].append(data["ICOHP"])
            cells[0]["edge_properties"]["icohp_bonding_perc"].append(data["ICOHP_bonding_perc"])
        # iterate over all edges that are equivalent to the current one due to translational symmetry
        for eq_atom1 in eq_atoms[node1]:
            start = cells[0]["frac_coords"][eq_atom1]
            shift = start - frac_coord1
            end = frac_coord2 + shift
            # check if edges lies within cell
            if (-tol <= end[0] <= 1+tol) and \
               (-tol <= end[1] <= 1+tol) and \
               (-tol <= end[2] <= 1+tol):
                cells[0]["edges"].append((start, end))
                cells[0]["edge_properties"]["cohp_plot"].append((x, y))
                cells[0]["edge_properties"]["bond_length"].append(data["bond_length"])
                cells[0]["edge_properties"]["icobi"].append(data["ICOBI"])
                cells[0]["edge_properties"]["icoop"].append(data["ICOOP"])
                cells[0]["edge_properties"]["icohp"].append(data["ICOHP"])
                cells[0]["edge_properties"]["icohp_bonding_perc"].append(data["ICOHP_bonding_perc"])

    return cells, atoms, eq_atoms



def get_primitive_supercell(
        lobstergraph: LobsterGraph, cells: list, atoms: dict, cart_crystal_axis_matrix: np.ndarray
) -> list:
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

    # block for building primitive supercell
    # iterate over x, y, z dimension
    # x_max/y_max/z_max times in positive direction
    num_atoms = len(atoms)
    for i, (dim_min, dim_max) in enumerate([(x_min, x_max), (y_min, y_max), (z_min, z_max)]):
        shift = np.array([0, 0, 0])
        shift[i] = 1.0
        new_cells = []
        # repeat axes, nodes and edges x_min/y_min/z_min times in negative direction
        for cell in cells:
            for j in range(dim_min, 0):
                new_cell = {
                    "frac_coords": [],
                    "edges": [],
                    "axes": [],
                    "edge_properties": {
                        "cohp_plot": [],
                        "bond_length": [],
                        "icobi": [],
                        "icoop": [],
                        "icohp": [],
                        "icohp_bonding_perc": []
                    }
                }
                for start, end in cell["axes"]:
                    new_start = start + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_end = end + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_cell["axes"].append((new_start, new_end))
                for k, (start, end) in enumerate(cell["edges"]):
                    new_start = start + j * shift
                    new_end = end + j * shift
                    new_cell["edges"].append((new_start, new_end))
                    new_cell["edge_properties"]["cohp_plot"].append(cell["edge_properties"]["cohp_plot"][k])
                    new_cell["edge_properties"]["bond_length"].append(cell["edge_properties"]["bond_length"][k])
                    new_cell["edge_properties"]["icobi"].append(cell["edge_properties"]["icobi"][k])
                    new_cell["edge_properties"]["icoop"].append(cell["edge_properties"]["icoop"][k])
                    new_cell["edge_properties"]["icohp"].append(cell["edge_properties"]["icohp"][k])
                    new_cell["edge_properties"]["icohp_bonding_perc"].append(
                        cell["edge_properties"]["icohp_bonding_perc"][k]
                    )
                for k, frac_coord in enumerate(cell["frac_coord"]):
                    new_cell["frac_coords"].append(frac_coord + j * shift)
                    new_cell["edge_properties"].append(cell["edge_properties"][k])
                    atoms[num_atoms] = {
                        "element": atoms[k]["element"],
                        "number": atoms[k]["number"]
                    }
                    num_atoms += 1
                new_cells.append(new_cell)
            for j in range(1, dim_max + 1):
                new_cell = {
                    "frac_coords": [],
                    "edges": [],
                    "axes": [],
                    "edge_properties": {
                        "cohp_plot": [],
                        "bond_length": [],
                        "icobi": [],
                        "icoop": [],
                        "icohp": [],
                        "icohp_bonding_perc": []
                    }
                }
                for start, end in cell["axes"]:
                    new_start = start + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_end = end + np.dot(cart_crystal_axis_matrix, j * shift)
                    new_cell["axes"].append((new_start, new_end))
                for k, (start, end) in enumerate(cell["edges"]):
                    new_start = start + j * shift
                    new_end = end + j * shift
                    new_cell["edges"].append((new_start, new_end))
                    new_cell["edge_properties"]["cohp_plot"].append(cell["edge_properties"]["cohp_plot"][k])
                    new_cell["edge_properties"]["bond_length"].append(cell["edge_properties"]["bond_length"][k])
                    new_cell["edge_properties"]["icobi"].append(cell["edge_properties"]["icobi"][k])
                    new_cell["edge_properties"]["icoop"].append(cell["edge_properties"]["icoop"][k])
                    new_cell["edge_properties"]["icohp"].append(cell["edge_properties"]["icohp"][k])
                    new_cell["edge_properties"]["icohp_bonding_perc"].append(
                        cell["edge_properties"]["icohp_bonding_perc"][k]
                    )
                for k, frac_coord in enumerate(cell["frac_coords"]):
                    new_cell["frac_coords"].append(frac_coord + j * shift)
                    atoms[num_atoms] = {
                        "element": atoms[k]["element"],
                        "number": atoms[k]["number"]
                    }
                    num_atoms += 1
                new_cells.append(new_cell)
        cells += new_cells
    return cells



def create_plot(lobstergraph: LobsterGraph, completecohp: CompleteCohp):
    """
    Creation of an interactive 3D plot of a compound's primitive supercell, containing information about site and
    bond properties.

    :param lobstergraph: LobsterGraph object, contains information about connectivity, site and bond properties in a
                         graph-like manner
    :return: plotly Figure, visualization of primitive supercell by 3D scatter plot
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

    cells, atoms, eq_atoms = get_primitive_cell(lobstergraph, completecohp)

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

    cells[0]["axes"] = axes

    # matrix to transform fractional to cartesian coordinates
    cart_crystal_axis_matrix = np.stack((x, y, z), axis=-1)

    cells = get_primitive_supercell(lobstergraph, cells, atoms, cart_crystal_axis_matrix)

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
        height=850
    )

    fig = go.Figure(layout=layout)

    i = 0
    for cell in cells:
        for frac_coord in cell["frac_coords"]:
            coord = np.dot(cart_crystal_axis_matrix, frac_coord)
            node_x.append(coord[0])
            node_y.append(coord[1])
            node_z.append(coord[2])
            atom_number.append(atoms[i]["number"])
            i += 1
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
        for start, end in cell["axes"]:
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
    #fig.show()



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

dir = "mp-10143/"
#dir = "mp-510401"
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

COHPCAR_path = os.path.join(os.path.expanduser("~/automationresults"), dir, "COHPCAR.lobster")
POSCAR_path = os.path.join(os.path.expanduser("~/automationresults"), dir, "POSCAR")

completecohp = CompleteCohp.from_file(
    fmt="LOBSTER", filename=COHPCAR_path, structure_file=POSCAR_path
)

plotfig = create_plot(testgraph, completecohp)