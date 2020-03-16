from numpy     import ones, array
from Polyomino import Solve, PlotSol
# Size of grid to tile
M, N = 11,11
print(M)
# Generating the rectangle
R   = ones((M,N))

# Fixed polyominos used in Yet another math programming consultant article
polyominos = [ array( [ [1] , [1] , [1] , [1] ] ) ,
               array( [ [0,1] , [1,1] , [0,1] ] ) ,  
               array( [ [1,1] , [1,1] ] ) ,
               array( [ [1,1,1] , [1,0,0] ] ) ,
               array( [ [0,1] , [1,1] , [1,0] ] ) ]

# Solving the layout problem
A_sol, obj = Solve( R , polyominos , Output = True )

# Plotting the solution
PlotSol( A_sol, polyominos, obj )