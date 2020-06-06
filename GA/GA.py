import random
from numpy import array, ones, zeros
import matplotlib.pyplot as plt


def printlist(rect):
    len(rect)
    for i in range(len(rect)):
        print(rect[i])
        
def add_poly(rect , poly , m = 0 , n = 0 , orr = 1):
    tmp_ = rect

    fit = True
    
    Mp,Np = poly.shape
    M,N   = rect.shape
    for i in range(Mp):
        for j in range(Np):
            if poly[i,j] == 1:
                if ((m + i ) >= M or (n + j ) >=  N):
                    fit = False

                elif rect[m+i,n+j] != 0 and fit == True:
                    fit = False
    if fit == True:
        for i in range(Mp):
            for j in range(Np):
                if ((m + i ) < M and (n + j ) < N):
                    tmp_[m+i,n+j] += poly[i,j]*orr
                
    return fit,tmp_

def counts(rect,poly,m = 0,n = 0):
    # Initialize
    CE_ = 0
    BE_ = 0
    Mp,Np = poly.shape
    M,N   = rect.shape
    # 
    for i in range(Mp):
        for j in range(Np):
            if poly[i,j] == 1 and ((m+i+1) < M and (n+j+1) < N):
                if rect[m+i,n+j] != 0:
                    CE_ = 0
                    BE_ = 0
                    return CE_,BE_
                else:
                    if (m+i) > M:
                        break
                    if (n+j) > N:
                        break
                    if (m + i + 1 ) == M:
                        BE_ += 1
                    if (n + j + 1) == N:
                        BE_ +=1
                    if (m + i) == 0:
                        BE_ += 1
                    if (n + j) == 0:
                        BE_ += 1    
                    if (n+j) < N:
                        if rect[m + i - 1,n + j] != 0:
                            CE_ += 1
                    if (n+j) < N and (m+i+1)< M:
                        if rect[m + i + 1,n + j] != 0:
                            CE_ +=1
                    if rect[m + i,n + j - 1] != 0:     
                        CE_ += 1
                    if n+j+1 < N:
                        if rect[m + i,n + j + 1] != 0:
                            CE_ += 1

    return CE_,BE_


def fitness(CommomEgdes,BoundaryEdges,weights = [1.0,0.0]):
    return weights[0]*CommomEgdes + weights[1]*BoundaryEdges
    
def FIND(number,rect):
    
    M,N = rect.shape
    for i in range(M):
        for j in range(N):
            if rect[i,j] == number:
                return i,j    

            
def ExchangeMatrix(order):
    from cvxopt import matrix
    J = matrix(0.0,(order,order))
    for i in range(order):
        J[i,order-i-1] = 1
    return J

def Orientation(p,label):
    M,N = p.shape
    if label == 0 or label == '000':
        return p
    if label == 1 or label == '001':
        p = ExchangeMatrix(M)@p
        return p
    if label == 2 or label == '010':
        p = p@ExchangeMatrix(N)
        return p
    if label == 3 or label == '011':
        p = ExchangeMatrix(M)@p@ExchangeMatrix(N)
        return p
    if label == 4 or label == '100':
        p = p.T
        return p
    if label == 5 or label == '101':
        p = ExchangeMatrix(N)@p.T
        return p
    if label == 6 or label == '110':
        p = p.T@ExchangeMatrix(M)
        return p
    if label == 7 or label == '111':
        p = ExchangeMatrix(N)@p.T@ExchangeMatrix(M)
        return p

def PriorityMatrix(rect):
    M,N = rect.shape

    Xc,Yc = M//2,N//2
    

    
    R_prior = zeros((M,N))
    if M%2 == 0:
        number = 1
        for k in range(Xc+1):
            for i in range(2*k):
                for j in range(2*k):
                    if R_prior[Xc-k+i,Yc-k+j] == 0 and ((Xc - k + i) != (Xc-1) or (Yc - k + j) != (Yc-1)):
                        R_prior[Xc-k+i,Yc-k+j] = number
                        number += 1
    else:
        number = 1
        for k in range(Xc+1):
            for i in range(2*k+1):
                for j in range(2*k+1):
                    if R_prior[Xc-k+i,Yc-k+j] == 0 and ((Xc - k + i) != Xc  or (Yc - k + j) != Yc):

                        R_prior[Xc-k+i,Yc-k+j] = number
                        number += 1

    return R_prior



def PlotSol( sol_matrix  , obj_val , Psize = (8,8) ):
    # The purpose of this function is to plot the solution for the minimization
    # problem
    #
    # INPUT:
    #     sol_matrix: Solution matrix
    #     polyominos: Fixed polyominos
    #     obj_val   : objective value for optimal solution

    import matplotlib.pyplot as plt
    import numpy             as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Get matrix size
    M,N = sol_matrix.shape

    plt.figure(figsize = Psize)
    plt.title('Size of rectangle is %d$\\times$%d \n%d of %d squares are filled' % ( M, N, obj_val, M*N ) )
    ax = plt.gca()
    # Get discrete colormap
    polyominos   = np.max(sol_matrix)
    cmap         = plt.get_cmap( 'jet',  polyominos )  
    value        = 0
    masked_array = np.ma.masked_where( sol_matrix == value, sol_matrix )

    cmap.set_bad( color = 'w' )
    # Plot and set limits .5 outside true range
    mat = ax.imshow( masked_array, cmap = cmap, vmin = 1 - .5, vmax = ( polyominos ) + .5)

    # Tell the colorbar to tick at integers
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax = plt.colorbar( mat, ticks = np.arange( 0, ( polyominos ) + 1 ),cax = cax )

    

    ax.set_xticks( np.arange( -.5, N, 1 ) )
    ax.set_yticks( np.arange( -.5, M, 1 ) )
    ax.grid( which = 'major', color = 'k', linestyle = '-', linewidth = 2 )
    ax.tick_params( colors = 'w' )                                            # removes tick numbering
    
    # Title and show plot
    
    plt.show(ax) 