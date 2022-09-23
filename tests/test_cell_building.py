import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_primitive_cell
from helpers_for_tests import lobstergraph, completecohp
from itertools import product



def test_primitive_cell_building():
    cell = get_primitive_cell(lobstergraph, completecohp)

    true_vecs = [vec for vec in product([0, 1], repeat=3)]
    test_vecs = [cell["atoms"][bond_label]["frac_coord"] for bond_label in cell["atoms"].keys()]

    tol = 0.01

    for true_vec in true_vecs:
        flag = False
        for test_vec in test_vecs:
            diff_a = abs(true_vec[0] - test_vec[0])
            diff_b = abs(true_vec[1] - test_vec[1])
            diff_c = abs(true_vec[2] - test_vec[2])
            if (diff_a <= tol) and (diff_b <= tol) and (diff_c <= tol):
                flag = True
                break
        if not flag:
            assert False

    assert True