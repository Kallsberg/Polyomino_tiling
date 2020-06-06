# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple MIP model
# using the matrix API:
#  maximize
#        x +   y + 2 z
#  subject to
#        x + 2 y + 3 z <= 4
#        x +   y       >= 1
#        x, y, z binary

from Polyomino import Variants,Constraints, Solution, PlotSol, Solve
from numpy     import ones, array, zeros, mean, savetxt,linspace

import time              as time
import matplotlib.pyplot as plt
# Size of grid to tile
Mmin  = 14
Mmax  = 14
M_    =  linspace(Mmin,Mmax,Mmax-Mmin+1)
k     = len(M_)
Tests = 1
T     = zeros((k,Tests))
FN    = zeros((k,Tests))

for i in range(k):
    for j in range(Tests):
        print('-------------------')
        print('M = N =',M_[i])
        K_    = 0
        t = time.time()
        N = int(M_[i])
        M = N
        # Generating the rectangle
        R   = ones( (M,N) )

        # Free Polyominoes, 4 test cases
        p = [ array( [ [1,0,0] ,
                       [1,1,1] ] )]#,
            #   array( [ [0,1,0] ,
            #            [1,1,1]]) ]
            
        # Construct orientations and compute number of fixed polyominoes    
        polyominos, NUM = Variants( p )
        
        # Solve the problem
        
        A_sol, obj, status = Solve( R , polyominos , Output = True , UpperBound = False , MaxIter = None)
        
        if status != 2:
            while status != 2:
                K_ += 1
                print('K: {:6}'.format(K_))
    
                A_sol, obj, status = Solve( R, polyominos, Output = True, K = K_ ,UpperBound = True , MaxIter=1e3)
        
        print('Status: {:}'.format(status))
        print('Obj: {:7}'.format(obj))
        T[i,j]  = time.time() - t
        FN[i,j] = obj

FN_avg = zeros((k,1))
T_avg  = zeros((k,1))        

for i in range(k):
    FN_avg[i] = M_[i]**2 - mean(FN[i,:])
    T_avg[i] = mean(T[i,:])

# Plotting results
plt.semilogy(M_,T_avg,'-o')
plt.grid(which = 'both')
plt.xlabel('M = N')
plt.ylabel('Average time [sec]')
plt.show()

plt.plot(M_,FN_avg,'-o')
plt.grid(which = 'both')
plt.xlabel('M=N')
plt.ylabel('Number of holes')
plt.show()


# # Save data
# savetxt('FN_feasi.txt',FN, delimiter = ',') 
# savetxt('FN_avg_feasi.txt',FN_avg, delimiter = ',') 
# savetxt('Time_feasi.txt',T, delimiter = ',') 
# savetxt('Time_avg_feasi.txt',T_avg, delimiter = ',')     