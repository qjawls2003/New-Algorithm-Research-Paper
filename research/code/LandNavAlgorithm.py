import random
import numpy
import math

                        # 3 dimensional position (used to solve a 4-dimensional problem space) example
                        #ScopeFactor = 2
                        #         Matrix: (numpy array)
                        #explorer: [x   y     z]
                        #scope 1 : [x+2 y     z]-block 1 |
                        #scope 2 : [x-2 y     z]-------- |
                        #scope 3 : [x   y+2   z]-block 2 |
                        #scope 4 : [x   y-2   z]-------- |
                        #scope 5 : [x   y   z+2]-block 3 |
                        #scope 6 : [x   y   z-2]-------- |

                        # s= scope particles
                        # E = explorer particle
                        #in a 3-D search space, the explorer searches with scope that looks like a cross.
                        #randomizing the ScopeFactor will make "crosses" of all different size and variations.
                        #3D Search Space (2 dimensional positions) simplified:
                        #y
                        #|      s
                        #|      |
                        #|      |
                        #|   s--E--s
                        #|      |
                        #|      |
                        #|      S
                        #x-----------------
                        
def LNA(lb, ub, dim, iters, func):

    PopSize = (dim*2)+1 #every dimension will have two scope particles, +1 is for the explorer particle
    particles = numpy.zeros((PopSize,dim)) #create a matrix of the particles, ex. 2 dimensions will have one explorer particles and 4 scope particles
                                            #(oriented East, West, North, South, respectively)
    fit = numpy.zeros(PopSize) #array of the fitness of each particle
    best = numpy.zeros(dim) #position of the best particle (with the miniumum fitness)
    bestfit = float("inf") #best value (minimum)
    
    pos = numpy.random.uniform(0,1,(dim))*(ub-lb)+lb #randomly places the explorer particle within the upper and lower bounds
    #convergence_curve=numpy.zeros(iters)
    restart = 0 #to avoid resetting the explorer again after one iteration
    repeated = 0 #used track repeats, to tighten the scope if the scope is too big
    prob = 0.1 #future use
    golden = 1.61803398875
    for l in range(iters): #for every iteration
        parity = l
        if repeated < iters/golden: #if the repeated bestfit values are more than iters/1.618, arbitrary
            i = 0
            j = golden/((l**2)+1) #arbitrary, will need to justify it through experimentation
        else:
            i = 0
            j = golden/((l*repeated)+1) #arbitrary
        
        
        
        for i in range(dim): #for each dimensions, ex. x coord -> y coord -> z coord
            block = 0 #used to designate which coordinate will be "scoping"
            switch = 0 #binary, since the scope can only to either subtracted or added from the explorer's coordinate
            count = 0 #used to assign blocks
            for j in range(0, PopSize): #for every particle
                ScopeFactor = numpy.random.uniform(i,j)*(ub/((repeated+1)**2)) #each particle will get a random scopefactor depending on iterations and repetition
                if restart == 0: #if it is the first iteration, assign the initial position to the explorer
                    particles[0,i] = pos[i].item()
                    restart = 1
                elif j == 0: #if not, keep going
                    continue
                    
                else:
                  
                    if block == i: #if the block matches the dimension
                        
                        if switch == 0 : #first scope particle gets the positive movement from the explorer's current position
                            positivestep = particles[0,i].item() + ScopeFactor #move in a positive direction (North or East in a 2-dimensional position)
                            if positivestep <= ub: #the new step is within the upper bound of the problem
                                particles[j,i] = positivestep
                            else:
                                particles[j,i] = ub
                            switch = 1
                        elif switch == 1: #second scope particle gets the negative movement from the explorer's current position
                            negativestep = particles[0,i].item() - ScopeFactor #move in a negative direction (South or West in a 2-dimensional position)
                            if negativestep >= lb: #the new step is within the lower bound of the problem
                                particles[j,i] = negativestep
                            else:
                                particles[j,i] = lb
                            switch = 0
                        count += 1
                        if (count % 2) == 0: #goto the next block when done with the scope pair
                            block += 1
                    else:
                        particles[j,i] = particles[0,i].item() #if the values are remaining the same as the explorer's position
                        count += 1
                        if (count % 2) == 0:
                            block += 1
        
        best, fit = calcBestFitness(particles, PopSize, dim, bestfit, func) #calcuate the best (position of the best particle), fit (array of all particles' fitness)
        #oldbestfit = bestfit
        bestfit = min(fit) #get the minimum fitness
        if numpy.array_equal(best,particles[0]): #if the bestfit values repeats (used to tigthen the search scope)
            repeated += 1

        newprob = 0 #work in progress to add in tunneling/jumping but may not need it
        if prob < newprob: #placeholder, does not jump at all yet
            pos = numpy.random.uniform(i,j,(dim)) #assign random position to jump to
            for i in range(dim):
                particles[0,i] = pos[i]
        else:
            for i in range(dim):
                particles[0,i] = best[i]
    print("Best Solution: ", best, " Value: ", bestfit)
        

        
    
    
def calcBestFitness(particles, PopSize, dim, bestfit, func):
    fit = numpy.zeros(PopSize)
    best = numpy.zeros(dim)

    for i,pos in enumerate(particles):
        
        fitness = func(pos)
       
        fit[i] = fitness
        
    best = particles[numpy.argmin(fit)] #get the particle with the lowest fitness
 
    return best, fit
        

def function1(x): #f(0,0,...0) = 0
    total=0
    for i in range(len(x)):
        total+= (x[i])**2
    return total

def function2(coord): #Beale Function: f(3,0.5) = 0
    x = coord[0]
    y = coord[1]

    f = (1.5-x+(x*y))**2+(2.25-x+(x*(y**2)))**2+(2.625-x+(x*(y**3)))**2
    return f

def function3(coord): #Levi Function: f(1,1) = 0
    x = coord[0]
    y = coord[1]
    pi = math.pi
    f = ((math.sin(3*pi*x))**2)+((x-1)**2)*(1+(math.sin(3*pi*y))**2)+((y-1)**2)*(1+(math.sin(2*pi*y))**2)
    return f

#LNA(lowerbound, upperbound, positional dim, iterations, function)
#For positional dimension, it is one dimension less than the actual function search space.
for i in range(10): #run the test 10 times
    LNA(-10,10,2,500, function3)

#For the Levi function, since there are so many local minimas, search is difficult. This can be
#mitigated by figuring out the optimal "Scope" retraction and contraction function

          
                      
            
