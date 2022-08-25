import os
import json
import itertools

import numpy as np

from pymatgen.analysis.graphs import StructureGraph, Structure, Lattice
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from structuregraph_helpers.create import get_structure_graph
from structuregraph_helpers.plotting import plotly_plot_structure_graph
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer



def structgraph_pymatgen(directory: str) -> StructureGraph:
    ext = os.path.expanduser("~/automationresults/")
    filepath = ext + directory
    chemenvlobster = LobsterNeighbors(
        are_coops=False,
        filename_ICOHP=filepath + "ICOHPLIST.lobster",
        structure=Structure.from_file(filepath + "POSCAR"),
        additional_condition=0
    )
    return chemenvlobster.get_bonded_structure(structure=chemenvlobster.structure)


def structgraph_structuregraph_helpers(directory: str) -> StructureGraph:
    ext = os.path.expanduser("~/automationresults/")
    filepath = ext + directory + "POSCAR"
    return get_structure_graph(Structure.from_file(filepath), method="minimaldistance")


def plot_json(sg: StructureGraph, plotfile: str = "plot.eps") -> None:
    with open('data.json', 'w') as f:
        json.dump(sg.to_json(), f)
    sg.draw_graph_to_file(plotfile, hide_image_edges=False, algo="neato")


def plot_structuregraph_helpers(sg: StructureGraph) -> None:
    figure = plotly_plot_structure_graph(sg)
    figure.show()


def build_supercell(structure: Structure, a_scale: float = 1.0, b_scale: float = 1.0, c_scale: float = 1.0) -> Structure:
    scale_matrix = [
        [a_scale, 0.0, 0.0],
        [0.0, b_scale, 0.0],
        [0.0, 0.0, c_scale]
    ]

    sga = SpacegroupAnalyzer(structure)
    conventional_structure = sga.get_conventional_standard_structure()
    conventional_structure.make_supercell(scaling_matrix=scale_matrix)

    a = conventional_structure.lattice.a
    b = conventional_structure.lattice.b
    c = conventional_structure.lattice.c
    alpha = conventional_structure.lattice.alpha
    beta = conventional_structure.lattice.beta
    gamma = conventional_structure.lattice.gamma

    if a_scale != 1.0:
        a = a / a_scale
    if b_scale != 1.0:
        b = b / b_scale
    if c_scale != 1.0:
        c = c / c_scale

    new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    conventional_structure.lattice = new_lattice

    return conventional_structure



combinations = itertools.combinations_with_replacement([-1, 0, 1], 3)
for comb in combinations:
    permutations = set(itertools.permutations(comb))
    for perm in permutations:
        print(perm)
        print(np.linalg.norm(perm))
        print("__________________")