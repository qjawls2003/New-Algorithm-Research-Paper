import random
import numpy
import math

def LNA(lb, ub, dim, iters):

    PopSize = (dim*2)+1
    particles = numpy.zeros((PopSize,dim))
    fit = numpy.zeros(PopSize)
    best = numpy.zeros(dim)
    bestfit = float("inf")
    
    pos = numpy.random.uniform(0,1,(dim))*(ub-lb)+lb
    convergence_curve=numpy.zeros(iters)
    restart = 0
    
    for l in range(iters):
        i = 0
        j = 1.618/((l*2)+1)
        ScopeFactor = numpy.random.uniform(j,i)*ub
        #print(ScopeFactor)
        #print(ScopeFactor)
        #print("iter: ", l)
        
        for i in range(dim):
            block = 0
            switch = 0
            count = 0
            for j in range(0, PopSize):

                if restart == 0:
                    particles[0,i] = pos[i].item()
                    restart = 1
                elif j == 0:
                    continue
                    
                else:
                  
                    if block == i:
          
                        if switch == 0 :
                            particles[j,i] = particles[0,i].item() + ScopeFactor
                            switch = 1
                        elif switch == 1:
                            particles[j,i] = particles[0,i].item() - ScopeFactor
                            switch = 0
                        count += 1
                        if (count % 2) == 0:
                            block += 1
                    else:
                        particles[j,i] = particles[0,i].item()
                        count += 1
                        if (count % 2) == 0:
                            block += 1
        
        best, fit = calcBestFitness(particles, PopSize, dim, bestfit)
        bestfit = min(fit)
        
        for i in range(dim):
            particles[0,i] = best[i]
    print("Best Solution: ", best, " Value: ", bestfit)
        

        
    
    
def calcBestFitness(particles, PopSize, dim, bestfit):
    fit = numpy.zeros(PopSize)
    best = numpy.zeros(dim)

    for i,pos in enumerate(particles):
        
        fitness = function1(pos)
       
        fit[i] = fitness
        
    best = particles[numpy.argmin(fit)]
 
    return best, fit
        

def function1(x):
    total=0
    for i in range(len(x)):
        total+= x[i]**2
    return total

LNA(-10,10,3,500)
                    
                    
                
            
