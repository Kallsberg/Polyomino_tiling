from numpy import array, flip, transpose, rot90, shape, append, zeros,max, where
from cvxopt import matrix, spmatrix, solvers, normal


def ExchangeMatrix(order):
    from cvxopt import matrix
    J = matrix(0.0,(order,order))
    for i in range(order):
        J[i,order-i-1] = 1
    return J

def Equal( p , q ):
    # The purpose of this function is to determine whether two polyominos are equal
    # INPUT:
    #     p : user input polyomino
    #     q : rotated/flipped polyomino
    # OUTPUT
    #     true/false

    # Storing sizes of the polyominos
    m_p, n_p = p.shape
    m_q, n_q = q.shape

    # First check if the dimensions are equal
    if m_p != m_q:
        return False
    elif n_p != n_q:
        return False
    else:
        # Next check if the number of squares are equal
        sum_p = sum( sum( p ) )
        sum_q = sum( sum( q ) )
        if sum_p != sum_q:
            return False
        else:
            # Check if elements are the same (elementwise product)
            sum_prod = sum( sum( p*q ) ) 
            if sum_prod != sum_p:
                return False
            else:
                return True

def Variants( p , reflect = True):
    # The purpose of this funtion is to compute number of orienations and construct the orientaitons of a 
    # given polyomino p
    #
    # INPUT:
    #     p : userdefined polyomino
    #
    # OUTPUT:
    #     polyominos: list of all the orientations of the polyomino
    #     number    : number of orientations of a given polyomino
    
    # Initialize counter and list
    number     = 0
    polyominos = []

    # If reflection is considered(one-sided polyominoes do not reflect)
    if reflect == True:
        ref = 2
    else:
        ref = 1
     

    # Main computation
    for k in range( len( p ) ):
        for reflect in range( ref ):
            for rotate in range( 4 ):
                
                if reflect == 1:
                    # Rotation with reflection
                    q = flip( rot90( p[k] , rotate ), 1 )
                else:
                    # Rotation without reflection
                    q = rot90( p[k], rotate )

                different = True
                for j in range( int( number ) ):
                    # Checking if the polyomino q is a new polyomino
                    if Equal( q,polyominos[j] ):
                        # If q is equal to a previous polyomino break for loop
                        different = False
                        break

                if different:
                    # If q is a new polyomino add to list
                    number += 1
                    polyominos.append( q )
    
    return polyominos, number

def Index( R ):
    # The purpose of this function is to index every cell in the reactangle in the tiling problem
    #
    # INPUT:
    #     R : matrix of zeros, the shape of tiling area
    # OUTOUT:
    #     R_index : matrix with enumerated cells
    
    # Get size of matrix R
    M,N = R.shape
    
    # Initialize counter
    l = 0

    # Allocate memory
    R_index = matrix( 0, (M,N) )

    for i in range( M ):
        for j in range( N ):
            l += 1 # update counter
        
            # Assign value to cell
            R_index[i,j] = int( l )

    num_eq =  int( max( R_index ) ) # number of equations
    
    return num_eq, R_index

def Constraints( R , polyominos ):
    '''
    The purpose for this funtion is to generate the constraint matrix
    for the problem
        max c*x
        s.t. 
        A*x <= y
        for x,y binary
    
    INPUT:
        R          : Rectangle to tile
        polyominoes: list of polyominoes
    
    OUTPUT:
        A    : Constraint matrix
        var  : Total number of variables
        var_k: Number of variables assigned to each polyomino
    '''

    M, N            = R.shape
    num_eq, R_index = Index( R ) 
    
    # Initialize counters
    var_tot = 0
    var_k   = zeros( len( polyominos ) )

    # Main computation
    for k in range( len( polyominos ) ):               # for each polyomino
        m_poly,n_poly = polyominos[k].shape            # size
        m_R,n_R       = R.shape
        
        # Number of unit squares in the polyomino
        squares = sum( sum( polyominos[k] ) )
        
        # Main computation
        for i in range( m_R - m_poly + 1 ):            # row
            for j in range( n_R - n_poly + 1 ):        # column 
                r,c = i + m_poly, j + n_poly           # placement
                # R[i:r,j:c] += polyominos[k]                 
                if sum( sum( R[i:r,j:c] * polyominos[k] ) ) == squares:
                    var_tot += 1                            
                    
                    # variables for each polymino
                    var_k[k] = var_tot                      
                    if k > 0:                               
                        var_k[k] = var_tot - sum( var_k[0:k] )
                        
    # Initialize matrix for constraints
    A   = zeros( ( num_eq , var_tot ) )
    var = 0
    
    # Main construction of A
    for k in range( len( polyominos ) ):
        m_poly,n_poly = polyominos[k].shape                     # size
        for i in range( R.shape[0] - m_poly + 1 ):              # row
            for j in range( R.shape[1] - n_poly + 1 ):          # column
                
                r , c = i + m_poly, j + n_poly                  # placement
                eqn = R_index[i:r , j:c] * polyominos[k]        # variable numbers
                
                # Only use positve values
                mg0 , ng0 = where(eqn > 0)
                
                # Update vector
                eqn = eqn[mg0,ng0] 
                
                # Update Constraint matrix
                A[eqn.astype('int') - 1,var] = 1

                # Advance to next variable
                var += 1
    return A, var, var_k

def Solution( x , R , polyominos ):
    '''
    The purpose of this function is to compute the solution to the problem
        max c*x
        s.t.
        A*x <= y
        for x,y binary
    INPUTS:
        x          : Solution
        R          : Rectangle
        polyominos : set of polyominoes
    OUTPUT:
        A_sol : Matrix with polyomino number, for plotting
    '''

    # Construciton of solution
    M,N = R.shape
    A_sol = zeros( (M,N) )
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    var = 0 # initialize counter

    # Main Construciton of solution
    # Looping through all the fixed polyominoes
    for k in range( len( polyominos ) ):
        m_poly,n_poly = polyominos[k].shape                     # size
        
        # Running through each cell in R
        for i in range( R.shape[0] - m_poly + 1 ):              # row
            for j in range( R.shape[1] - n_poly + 1) :          # column
                    
                r, c = i + m_poly, j + n_poly                   # placement
                if x[var] == 1:
                    A_sol[i:r, j:c] += polyominos[k]*( k + 1 )  # Enumerate cells by polymino number
                var += 1  # update conter
    return A_sol

def PlotSol( sol_matrix , polyominos , obj_val , Psize = (8,8) , Title = True, Colorbar = True):
    '''
    The purpose of this function is to plot the solution for the minimization
    problem
    
    INPUT:
        sol_matrix: Solution matrix
        polyominos: Fixed polyominos
        obj_val   : objective value for optimal solution
    
    OUTPUT:
        Plot of solution
    '''

    import matplotlib.pyplot as plt
    import numpy             as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Get matrix size
    M,N = sol_matrix.shape

    plt.figure(figsize = Psize)
    if Title == True:
        plt.title('Size of rectangle is %d$\\times$%d \n%d of %d squares are filled' % ( M, N, obj_val, M*N ) )
    ax = plt.gca()
    
    # Get discrete colormap
    cmap         = plt.get_cmap( 'jet', len( polyominos )  )
    value        = 0
    masked_array = np.ma.masked_where( sol_matrix == value, sol_matrix )

    cmap.set_bad( color = 'w' )
    # Plot and set limits .5 outside true range
    mat = ax.imshow( masked_array, cmap = cmap, vmin = 1 - .5, vmax = len( polyominos ) + .5)

    # Tell the colorbar to tick at integers
    if Colorbar == True:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cax = plt.colorbar( mat, ticks = np.arange( 0, len( polyominos ) + 1 ),cax = cax )

    

    ax.set_xticks( np.arange( -.5, N, 1 ) )
    ax.set_yticks( np.arange( -.5, M, 1 ) )
    ax.grid( which = 'major', color = 'k', linestyle = '-', linewidth = 2 )
    ax.tick_params( colors = 'w' )                                            # removes tick numbering
    
    # Title and show plot
    
    plt.show(ax)   

def KnapSack(poly_area,M,N):
    from gurobipy import GRB
    import gurobipy as gp
    
    # number of variables for problem
    n_var = len(poly_area)
    
    # Defining model
    m = gp.Model('KnapSack')
    m.setParam( 'OutputFlag', False )
    # Defining variables
    x = m.addMVar(shape = n_var,vtype = GRB.INTEGER,name = 'x')
    
    # Setting objective and constraint
    m.setObjective(poly_area @ x, GRB.MAXIMIZE)
    
    m.addConstr(poly_area @ x <= M*N , name = "c")

    # Solving the problem
    m.optimize() 
    
    if GRB.OPTIMAL == m.status:
        return m.objVal
    
def Solve( R, polyominos , Zeros = [] , Ones = [] , Output = True, MaxIter = None, K = None, UpperBound = None):
    '''
    This function is made to solve the tiling problem given a set of fixed polyominoes.

    INPUTS:
        R           : Rectangle to tile
        polyominoes : array of fixed polyominoes
        Zeros       : array of cell numbers for equality constraints with 0 on the RHS(optional)
        Ones        : array of cell numbers for equality constraints with 1 on the RHS(optional)
        Output      : boolean variable(default True), set to false if outputflag is not wanted
        MaxIter     : integer, if the user wants a number of max iterations. The number given is then 
                      multiplied with the  number of variables in the model
        K           : integer, K is used if the user knows there are going to be a specific number of holes
        UpperBound  : integer or boolean, default None. If True a knapsack problem is solved to find upper
                      bound, if integer the given number is the upper bound.
    '''

    from cvxopt   import matrix, solvers
    from gurobipy import GRB
    from numpy    import ones,zeros,unique
    import gurobipy as gp

    m = gp.Model("Polyomino")
    m.setParam( 'OutputFlag', Output )

    # Create variables
    M,N = R.shape
    

    
    # Constraints, RHS, and cost vector
    ineq_cons,var_tot, var_k = Constraints( R, polyominos ) 
    num_eq                   = ineq_cons.shape[0]
    # ineq_cons                = matrix( ineq_cons, ( num_eq,var_tot ), 'd' )
    obj                      = ones( ( var_tot ) ) 
    
    
    x = m.addMVar(shape = var_tot, vtype = GRB.BINARY, name = "x")
    y = m.addMVar(shape = M*N, vtype = GRB.BINARY, name = "y")
    
    
    # Adding weight corriponding to the number of squared in the polyomino
    tmp = 0
    i   = 0
    for k in range( len( polyominos ) ):
        squares = int( sum( sum( polyominos[k] ) ) ) 
    
        for i in range( int( var_k[k] ) ): 
            obj[i + tmp] = squares
        tmp += i + 1 
    
    
    # Solving the binary problem
    # Using Gurobi to solve the MIP
    
    # Set objective
    m.setObjective(ones(M*N) @ y, GRB.MAXIMIZE)

    if MaxIter != None:
        # The number of max interations is based on the number of varaibles
        m.setParam('IterationLimit',var_tot*MaxIter)
    
    if type(UpperBound) == bool:
        if UpperBound == True:
            P_sum = zeros(len(polyominos))
            for i in range(len(polyominos)):
                P_sum[i] = sum(sum(polyominos[i]))

            # Upperbound is found by solving a Knapsack problem
            UpperBound = KnapSack(unique(P_sum),M,N) - len(Zeros)
            
            if K == None:
                K = 0

            m.addConstr(sum(y) <= UpperBound - K)
            
    else:
        
        UpperBound = (M*N - len( Zeros ))//sum(sum(polyominos[0]))*sum(sum(polyominos[0]))
        m.addConstr(sum(y) == UpperBound - K)
        
    # Add constraints
    m.addConstr(ineq_cons @ x == y, name = "c")
    if len(Zeros) != 0:
        for i in range(len(Zeros)):
            m.addConstr(ineq_cons[Zeros[i],:] @ x == 0)
    if len(Ones) != 0:
        for i in range(len(Ones)):
            m.addConstr(ineq_cons[Ones[i],:] @ x == 1)
    # if K != None:
    #     m.addConstr(ones(M*N) @ y  == UpperBound - K, name = "feas")

    # Optimize model
    m.optimize() 

    if GRB.OPTIMAL == m.status:
        
        # Generate a solution matrix
        A_sol = Solution( x.X, R, polyominos )
        
        return A_sol , m.objVal, m.status
    else:
        return zeros((M,N)), 0, m.status
    
    