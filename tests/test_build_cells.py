import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_primitive_cell
from helpers_for_tests import lobstergraph, completecohp



def test_primitive_cell_building():
    data = get_primitive_cell(lobstergraph, completecohp)
    assert (data != None)