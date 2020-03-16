from Polyomino import Variants,Constraints, Solution, PlotSol, Solve
from numpy     import ones, array

# Size of grid to tile
M, N = 16,16

# Generating the rectangle
R   = ones( (M,N) )

# Test case
case = 4

# Free Polyominoes, 4 test cases
if case == 0:
    p = [ array( [ [1,0,0] ,
                   [1,0,0] ,
                   [1,1,1] ,
                   [1,1,1] ] ) ] # L-octomino
    
    # Construct orientations and compute number of fixed polyominoes    
    polyominos, NUM = Variants( p )
elif case == 1:
    p = [ array( [ [1,1] ,
                    [1,0] ] ) ] # L-trimino
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p )                    
elif case == 2:
    p = [ array( [ [1,1] ] ) ] # Domino
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p )                    
elif case == 3:
    # Set of polyominoes
    p = [ array( [ [1] , [1] , [1] ] ) ,        # Triomino
          array( [ [1,0] , [1,0] , [1,1] ] ) ,  # L-tetromino
          array( [ [1] ] ) ]                    # Monomino
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p )
elif case == 4:
    # Set of polyominoes
    p = [ array( [ [1,1,0] , [1,1,0] , [1,1,1], [1,0,0] ] ) ]        # Octomino
          
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p )

# Solve the problem 
A_sol, obj = Solve( R, polyominos )

# Print solution
PlotSol( A_sol, polyominos, obj )