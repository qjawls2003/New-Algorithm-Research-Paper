import random
import numpy
import math


class Particle:
    def __init__(self, dim, init_pos, num, bounds):
        self.id = num
        self.position = []
        self.velocity = []
        self.mass_p = 0
        self.mass_a = 0
        self.mass_i = 0
        self.fitness = 0
        self.dim = dim
        bound = bounds[1]-bounds[0]
        for i in range(dim):
            self.velocity.append(random.uniform(-1,1))
            self.position.append(random.uniform(0,1)*init_pos[i]*bound + bounds[0])


    def calculateFitness(self,costFunc):
        self.fitness = costFunc(self.position)
        return self.fitness
    

    def updatePosition(self,accel,dim):
        for j in range(dim):
            rand = random.random()
            self.velocity[j] = rand*self.velocity[j] + accel[self.id, j]
            self.position[j] = self.position[j] + self.velocity[j]



def G_Constant(curTime,totTime):
    alpha = 20
    G_init = 100
    e_val = numpy.exp(-alpha*float(curTime)/totTime)
    G = G_init*e_val
    return G

def euclideanDis(particle_i, particle_j,dim):
    distSum = 0
    for i in range(dim):
        distSum += (particle_j.position[i]-particle_i.position[i])**2
    return math.sqrt(distSum)


def forceCalc(i,t,iters,PopSize, dim,particles):
    G = G_Constant(t,iters)
    e = 0.0001
    M_pi = particles[i].mass_p
    kbest = int(PopSize*0.8)
    Force = numpy.zeros((PopSize,dim))
    accel = numpy.zeros((PopSize,dim))
    for y in range(PopSize):
        j = particles[y]
        for x in range(kbest):
            if j != i:
                R = euclideanDis(particles[i],j,dim)
                M_aj = j.mass_a
                for d in range(dim):
                    rand = random.random()
                    Force[y,d] = Force[y,d] + rand*(((M_pi*M_aj)*(j.position[d]-particles[i].position[d]))/(R+e))
            

    return Force*G
    
def massCalc(i,t,PopSize,dim, M_i):
    fitmax = max(fitness)
    fitmin = min(fitness)
    #fitsum = sum(fitness)
    #fitmean = fitsum/PopSize

    if fitmax == fitmin:
        M_i = numpy.zeros(PopSize)
    else:
        best = [fitmin,i]
        worst = [fitmax,i]
        bestfit = best[0]
        worstfit = worst[0]
        for p in range(PopSize):
            M_i[p] = (fitness[p] - worstfit)/(bestfit-worstfit)
    Msum = sum(M_i)
    for q in range(PopSize):
        M_i[q] = M_i[q]/Msum
    return M_i[i]

def accelCalc(i,t,PopSize, dim,force):
    accel = numpy.zeros((PopSize,dim))
    Force = force
    for x in range (PopSize):
        for y in range(dim):
            accel[x,y] = Force[x,y]
    return accel


class GSA():
    def __init__(self, PopSize,dim,costFunc,bounds,iterations):
        global best
        global worst
        global fitness
        global particles
        global M_i
        init_pos = [1,1]
        best = [999999999, 0]
        worst = [-999999999,0]
        fitness = numpy.zeros(PopSize)
        M_i = numpy.zeros(PopSize)
        particles = []
        for i in range(PopSize):
            particles.append(Particle(dim,init_pos,i,bounds))
            fitness[i] = particles[i].calculateFitness(costFunc)
            if fitness[i] < best[0]:
                    best = [fitness[i],i]
            elif fitness[i] > worst[0]:
                    worst = [fitness[i],i]

        counter = 0
        while counter < iterations:
            for idnum in range(PopSize):
                particles[idnum].mass_i = massCalc(idnum,counter, PopSize, dim, M_i)
                force = forceCalc(idnum,counter,iterations,PopSize, dim,particles)
                accel = accelCalc(idnum,counter,PopSize, dim,force)
                particles[idnum].updatePosition(accel,dim)
            counter += 1
        print("Solution: ")
        print("At: ", particles[best[1]].position)
        print("Result: ", best[0])

def func1(x):
    total=0
    for i in range(len(x)):
        total+=x[i]**2
    return total

if __name__ == "__GSA__":
    main()
bounds = [-10,10]
GSA(20,2,func1,bounds,100)
                    
                
