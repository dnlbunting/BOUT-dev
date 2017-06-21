
from grid import rectangular_grid
import numpy as np

def test_getPoloidalGrid():
    # Create a simple non-periodic grid
    grid = rectangular_grid(10,10,10, yperiodic=False)

    # Check number of y points
    assert grid.numberOfPoloidalGrids() == 10
    
    # Try to get points outside the domain
    p, y = grid.getPoloidalGrid(-1)
    assert p == None

    p, y = grid.getPoloidalGrid(10)
    assert p == None

def test_getPoloidalGrid_periodic():
    # Create a periodic grid
    grid = rectangular_grid(10,10,10, yperiodic=True)

    assert grid.numberOfPoloidalGrids() == 10

    p_last, y_last = grid.getPoloidalGrid(9) # Last in domain

    assert p_last is not None
    
    p, y = grid.getPoloidalGrid(-1)

    assert p is p_last
    assert y < y_last
