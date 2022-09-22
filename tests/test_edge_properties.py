import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_primitive_cell
from helpers_for_tests import lobstergraph, completecohp



def test_edge_properties():
    cell = get_primitive_cell(lobstergraph, completecohp)

    for test_edge, (_, _, true_edge) in zip(cell["edges"].values(), lobstergraph.sg.graph.edges.data()):
        if test_edge["bond_length"] != true_edge["bond_length"]:
            assert False
        elif test_edge["icobi"] != true_edge["ICOBI"]:
            assert False
        elif test_edge["icoop"] != true_edge["ICOOP"]:
            assert False
        elif test_edge["icohp"] != true_edge["ICOHP"]:
            assert False
        elif test_edge["icohp_bonding_perc"] != true_edge["ICOHP_bonding_perc"]:
            assert False

    assert True



def test_cohp_plot():
    cell = get_primitive_cell(lobstergraph, completecohp)
    
    for bond_label, edge in cell["edges"].items():
        cohp_data = completecohp.get_cohp_by_label(label=bond_label).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x_true = [spinup_cohps[j] + spindown_cohps[j] for j, _ in enumerate(spinup_cohps)]
        y_true = [energies[j] - fermi_energy for j, _ in enumerate(energies)]
        
        x_test = edge["cohp_plot"][0]
        y_test = edge["cohp_plot"][1]

        for test_point, true_point in zip(x_test, x_true):
            if test_point != true_point:
                assert False
        for test_point, true_point in zip(y_test, y_true):
            if test_point != true_point:
                assert False

    assert True