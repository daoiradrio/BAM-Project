import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_primitive_cell
from helpers_for_tests import lobstergraph, completecohp



def test_cohp_plot():
    bond_label = -1
    cohp_data = completecohp.get_cohp_by_label(label=bond_label).as_dict()
    spinup_cohps = cohp_data["COHP"]["1"]
    spindown_cohps = cohp_data["COHP"]["-1"]
    energies = cohp_data["energies"]
    fermi_energy = cohp_data["efermi"]
    x = [spinup_cohps[j] + spindown_cohps[j] for j, _ in enumerate(spinup_cohps)]
    y = [energies[j] - fermi_energy for j, _ in enumerate(energies)]



def test_edge_properties():
    data = get_primitive_cell(lobstergraph, completecohp)
    assert (data != None)