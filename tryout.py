import os
import warnings

import numpy as np

from lobsterpy.structuregraph.graph import LobsterGraph
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from itertools import product, permutations
from helper import create_graph_plot

warnings.filterwarnings(action='ignore')



dir = "mp-10143/"
#dir = "mp-510401"
#dir = "mp-2384"

path = os.path.join(os.path.expanduser("~/automationresults"), dir)

path_to_poscar = os.path.join(path, "POSCAR")
path_to_charge = os.path.join(path, "CHARGE.lobster")
path_to_icobilist = os.path.join(path, "ICOBILIST.lobster")
path_to_icooplist = os.path.join(path, "ICOOPLIST.lobster")
path_to_icohplist = os.path.join(path, "ICOHPLIST.lobster")
path_to_cohpcar = os.path.join(path, "COHPCAR.lobster")
path_to_madelung = os.path.join(path, "MadelungEnergies.lobster")

#"""
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

#"""
completecohp = CompleteCohp.from_file(
    fmt="LOBSTER", filename=path_to_cohpcar, structure_file=path_to_poscar
)
#"""

create_graph_plot(lobstergraph)