from collections import namedtuple
from copy import deepcopy
import numpy as np


def powers2(num):
    '''List (descending) of powers of two in the number'''
    powers = [int(power)
              for power, value in enumerate(reversed(format(num, 'b')))
              if value != '0']

    return powers[::-1]

# assert all(sum(2**p for p in powers2(n)) == n for n in (9, 13, 425, 123))

Tile = namedtuple('Tile', ('coords', 'cells', 'master_vertices', 'slave_vertices', 'mappings'))


def make_tile(x, cells, master_vertices, slave_vertices, vertex_mappings):
    '''Freeze the tile from data'''
    # The tile consists of coordinates x and cells as indices in to the 
    # cells array. master/slave vertices define a map for gluing in the
    # direction. The periodic maps for the remaining dirs are in vertex_
    # mappings
    return Tile(deepcopy(x), deepcopy(cells),
                np.copy(master_vertices), np.copy(slave_vertices),
                [vm.copy() for vm in vertex_mappings])
