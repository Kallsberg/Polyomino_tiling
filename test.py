import matplotlib.pyplot as plt

from Polyomino import Variants,ExchangeMatrix
from numpy     import array, zeros
from cvxopt    import matrix

# polys = []

# p = [array([[1],[1],[1],[1]]),
#      array([[1,0],[1,0],[1,1]])]



# p = [array([[1,1],[1,0],[1,0]])]                     # L-trimino
# polyominos,number = Variants(p, reflect= False)
# for i in range(len(polyominos)):
#     print(polyominos[i])


# m_p,n_p = polyominos[0].shape

# polys = zeros((m_p+4,int((m_p+n_p)*number/2) + 20))
# k = 0
# # for k in range(polyominos):


# # polys[(0+2):(m_p+2),(0+2):(n_p+2)] += polyominos[k]*(k+1)
# # k = 1
# # m_p,n_p = polyominos[k].shape

# # polys[(0+2):(m_p+2),(0+2+n_p+2):(n_p+2+n_p+2)] += polyominos[k]*(k+1)
# n_old = 2
# for k in range(number):
#     m_p,n_p = polyominos[k].shape    
#     polys[2:(m_p+2),(n_old):(n_old+n_p)] += polyominos[k]*(k+1)
#     n_old = (4 + n_old)
# plt.imshow(polys,cmap = 'gist_rainbow')
# plt.show()    

J = ExchangeMatrix(5)
x = matrix(1.0,(5,1))

print(J*x)
# ## Bobyqa