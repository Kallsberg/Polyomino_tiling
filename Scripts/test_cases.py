from Polyomino import Variants,Constraints, Solution, PlotSol, Solve
from numpy     import ones, array


# Size of grid to tile
M, N = 40,40

# Test case
case = 1

# Free Polyominoes, 4 test cases
if case == 0:
    # Size of grid to tile
    M, N = 40,40
    # running time was 446.77 seconds and 0 holes 
    p = [ array( [ [1,0,0] ,
                   [1,0,0] ,
                   [1,1,1] ,
                   [1,1,1] ] ) ] # L-octomino
    
    # Construct orientations and compute number of fixed polyominoes    
    polyominos, NUM = Variants( p )
elif case == 1:
    # Size of grid to tile
    M, N = 40,40

    p = [ array( [ [1,1,0,0] ,
                   [1,1,1,1] ,
                   [1,1,1,1] ] ) ] # L-trimino
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p ) 

    print(polyominos)                   
elif case == 2:
    p = [ array( [ [1,1],
                   [1,0] ] ) ] # L-triomino
    
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
    # Size of grid to tile
    M, N = 40,40
    # running time was 202.29 seconds and 0 holes 
    # Set of polyominoes
    p = [ array( [ [1,0,0] ,
                   [1,1,1] ] ), # L-tetromino
          array( [ [1,1,0,0] ,
                   [1,1,1,1] ,
                   [1,1,1,1] ] )      
                                ] 
          
    
    # Construct orientations and compute number of fixed polyominoes
    polyominos, NUM = Variants( p )



# Generating the rectangle
R   = ones( (M,N) )

# Solve the problem 
K_ = 0
A_sol, obj, status = Solve( R , polyominos , K = K_, Output = True , UpperBound = 1 )#, Zeros=[0,41,82,123,164],Ones=[1])
if status != 2:
    while status != 2:
        K_ += sum(sum(sum(p)))
        print('K: {:6}'.format(K_))

        A_sol, obj, status = Solve( R, polyominos, Output = True, K = K_ ,UpperBound = 1 )

# Print solution
PlotSol( A_sol, polyominos, obj, Title = False, Colorbar = False )