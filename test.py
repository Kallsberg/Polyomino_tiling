# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple MIP model
# using the matrix API:
#  maximize
#        x +   y + 2 z
#  subject to
#        x + 2 y + 3 z <= 4
#        x +   y       >= 1
#        x, y, z binary

import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from Polyomino import Constraints,Variants

try:
    R = np.ones((16,16))

    p = [ np.array( [ [1,0,0] ,
                   [1,0,0] ,
                   [1,1,1] ,
                   [1,1,1] ] ) ] # L-octomino
    
    # Construct orientations and compute number of fixed polyominoes    
    polyominos, NUM = Variants( p )

    # Create a new model
    m = gp.Model("matrix1")

    

    # Build (sparse) constraint matrix
    # data = np.array([1.0, 2.0, 3.0, -1.0, -1.0])
    # row = np.array([0, 0, 0, 1, 1])
    # col = np.array([0, 1, 2, 0, 1])

    # A = sp.csr_matrix((data, (row, col)), shape=(2, 3))
    A, var, var_k = Constraints(R, polyominos)

    # Create variables
    x = m.addMVar(shape = var, vtype = GRB.BINARY, name = "x")

    # Set objective
    obj = np.ones(var)
    
    print(var)
    
    print(obj.T @ x)    
    m.setObjective(obj @ x, GRB.MAXIMIZE)

    # Build rhs vector
    rhs = np.ones(16*16)

    # Add constraints
    m.addConstr(A @ x <= rhs, name="c")

    # Optimize model
    m.optimize()

    print(x.X)
    print('Obj: %g' % m.objVal)

except gp.GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError:
    print('Encountered an attribute error')