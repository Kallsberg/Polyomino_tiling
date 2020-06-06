import random
from numpy import array, ones, zeros, sum, mean, savetxt, linspace
import matplotlib.pyplot as plt
import time as time
from GA import *

import multiprocessing as mp
from joblib import Parallel, delayed










class Individual:
    def __init__(self, chromosome): 
        self.chromosome = chromosome  
        self.fitness = self.cal_fitness() 

    @classmethod
    def mutated_genes(self): 
        ''' 
        create random genes for mutation 
        '''
        global GENES 
        gene = random.choice(GENES) 
        return gene 

    @classmethod
    def create_gnome(self): 
        ''' 
        create chromosome or string of genes 
        '''
        global gnome_len
        return [self.mutated_genes() for _ in range(int(gnome_len))] 

    def mate(self, par2): 
        ''' 
        Perform mating and produce new offspring 
        '''
  
        # chromosome for offspring 
        child_chromosome = [] 
        for gp1, gp2 in zip(self.chromosome, par2.chromosome):     
  
            # random probability   
            prob = random.random() 
  
            # if prob is less than 0.45, insert gene 
            # from parent 1  
            if prob < 0.45: 
                child_chromosome.append(gp1) 
  
            # if prob is between 0.45 and 0.90, insert 
            # gene from parent 2 
            elif prob < 0.95: 
                child_chromosome.append(gp2) 
  
            # otherwise insert random gene(mutate),  
            # for maintaining diversity 
            else: 
                child_chromosome.append(self.mutated_genes()) 
  
        # create new Individual(offspring) using  
        # generated chromosome for offspring 
        return Individual(child_chromosome) 
        
    def cal_fitness(self):
        from numpy import array
        ''' 
        Calculate fittness score, it is the number of 
        characters in string which differ from target 
        string. 
        '''
        global M,N
        global P,R0
        global inds
        
        R      = zeros((M,N))
        fit    = 0
        m,n    = 0,0
        C      = 0 
        pols   = 0
        for k in range(len(self.chromosome)):
            
            P_ = Orientation(P,int(self.chromosome[k]))
            Mp,Np = P_.shape
            fits = False
            fit = 0.0
            m,n = 0,0
            
            if pols != 0:    
                
                for ind in range(M*N):  
                    i,j = inds[ind][1:]
                    CE,BE = counts(R,P_,i,j)
                    fit_ = fitness(CE,BE)
                    if fit_ > fit:
                        fit = fit_
                        m,n = i,j
                        
                fits,R = add_poly(R,P_,m,n,int(self.chromosome[k])+1)
            else:
                fits,R = add_poly(R,P_,(M-Mp)//2,(N-Np)//2,int(self.chromosome[k])+1)
            
            if fits == True:
                pols += 1
                C += fit

        
        C = 0.4*pols + 0.6*C
            
        
        return -C 










GENES           = '01234567'
POPULATION_SIZE = 75
MaxGen          = 50

global M,N
global P,q
global R0
global gnome_len
global inds

M_min = 8
M_max = 8
M_    = linspace( M_min , M_max , M_max - M_min + 1)
k     = len(M_)
Tests = 2
FN    = zeros((MaxGen,k*Tests))
T     = zeros((k*Tests,1))
chrom = []
inds  = []
def FITNESS(i):
    return population[i].fitness
for i in range(k):
    for ii in range(Tests):
        # print(i)
        t = time.time()
        M = int(M_[i])
        N = M


        

        R = zeros((M,N))
        P = array([[1,0,0],[1,1,1]])
        q = sum(P)
        gnome_len = (M*N)//q

        R0 = PriorityMatrix(R)

        for j in range(M*N):
            jj,jjj = FIND(j,R0)
            inds.append([j,jj,jjj])
            
        #current generation 
        generation = 1

        found = False
        population = [] 
        
        # create initial population 
        for _ in range(POPULATION_SIZE): 
                    gnome = Individual.create_gnome() 
                    population.append(Individual(gnome)) 
        # for i in range(len(population)):
                # print("".join(population[i].chromosome))
                # print(population[i].fitness)
        while not found and generation < MaxGen: 
            
            # sort the population in increasing order of fitness score 
            population = sorted(population, key = lambda x:x.fitness) 
            
            
            print(population[0])
            # Otherwise generate new offsprings for new generation 
            new_generation = [] 

            # Perform Elitism, that mean 10% of fittest population 
            # goes to the next generation 
            s = int((5*POPULATION_SIZE)/100) 
            new_generation.extend(population[:s]) 

            # From 50% of fittest population, Individuals  
            # will mate to produce offspring 
            s = (POPULATION_SIZE)//2
            for _ in range(s): 
                parent1 = random.choice(population[:s]) 
                parent2 = random.choice(population[:s]) 
                child   = parent1.mate(parent2) 
                new_generation.append(child) 

            population = new_generation 

            # print("Generation: {}\tString: {}\tFitness: {}".format(generation,"".join(population[0].chromosome), population[0].fitness)) 
            FN[generation-1,i*Tests + ii] = -population[0].fitness    
            generation += 1
            
            
        # print("Generation: {}\tString: {}\tFitness: {}".format(generation,"".join(population[0].chromosome),population[0].fitness)) 
        FN[generation-1,i*Tests + ii] = -population[0].fitness
        chrom.append(population[0].chromosome)

        T[i*Tests + ii] = time.time() - t
        print('--------------------',i*Tests + ii)

        # R0 = zeros((M,N))
        # R  = zeros((M,N))

        # R0 = PriorityMatrix(R0)

        # chromo = population[0].chromosome

        # fit    = 0
        # m,n    = 0,0
        # C      = 0 
        # pols   = 0
        # for k in range(len(chromo)):
        #     P_ = Orientation(P,int(chromo[k]))
        #     fits = False
        #     fit = 0.0
        #     m,n = 0,0
            
        #     if pols != 0:    
        #         for ind in range(M*N):  
        #             i,j = FIND(ind,R0)
        #             CE,BE = counts(R,P_,i,j)
        #             fit_ = fitness(CE,BE)
        #             if fit_ > fit:
        #                 fit = fit_
        #                 m,n = i,j



        #         fits,R = add_poly(R,P_,m,n,int(chromo[k])+1)
        #     else:
        #         fits,R = add_poly(R,P_,(M-P_.shape[0])//2,(N-P_.shape[1])//2,int(chromo[k])+1)
            
        #     if fits == True:
        #         pols += 1
        #         C += fit

        # C = 0.4*pols + 0.6*C

        # print('Fitness:     {:}'.format(C))
        # print('Succesfully: {:}'.format(pols))

        # PlotSol( R  , pols*q , Psize = (8,8) )

        # PlotSol( R[2:(M-2),2:(N-2)]  , pols*q , Psize = (8,8) )# No border

# Computing averages
FN_avg = zeros((MaxGen,k))
T_avg  = zeros((k,1))
for i in range(k):
    FN_avg[:,i] = mean(FN[:,i*Tests:((i+1)*Tests-1)],axis=1)
    T_avg[i]    = mean(T[i*Tests:((i+1)*Tests-1)])


# Save data
# savetxt('FN.txt',FN, delimiter = ',') 
# savetxt('FN_avg.txt',FN_avg, delimiter = ',') 
# savetxt('Time.txt',T, delimiter = ',') 
# savetxt('Time_avg.txt',T_avg, delimiter = ',') 
# with open("file.txt", "w") as output:  # save the best chromosomes from each test
#     output.write(str(chrom))
# print('Data has been saved')

# Plotting results
plt.semilogy(M_,T_avg,'-o')
plt.grid(which = 'both')
plt.xlabel('M = N')
plt.ylabel('Average time [sec]')
plt.show()

plt.plot(FN_avg)
plt.grid(which = 'both')
plt.xlabel('Generation')
plt.ylabel('Average Fitness')
plt.legend(M_,bbox_to_anchor=(0.3, 1.0),ncol = 3)
plt.show()



# # Read data
# test = np.loadtxt('test.txt',delimiter = ',')



