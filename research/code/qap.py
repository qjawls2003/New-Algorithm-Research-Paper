'''
Created on Aug 30, 2012
Modified for AI class on 9 Sep 2018

@author: alexander.mentis

Solve QAP using simulated annealing.
'''
#CDT Beomjin An

import copy
import math
import random

def init_flow(infile):
    """
    Initialize and return the flow matrix.
    
    Reads the file pointed to by infile. The file must have the following
    format:
    
    Line x should contain the list of flow values from department x+1 to 
    departments 1 through x-1, each flow separated by whitespace. A blank line
    terminates the file, so file comments can be inserted after one or more
    blank lines.
    """
    
    flows = []
    for line in infile:
        if line.strip():
            flows.append([int(flow) for flow in line.split()])
        else:
            break
    
    return flows

def init_locations(flows):
    """
    Set initial department locations randomly.
    """
    num_departments = len(flows) + 1 # flows doesn't have row for 1st department
    
    # assume rectangular layouts
    rows = math.floor(math.sqrt(num_departments))
    cols = math.ceil(num_departments / rows)
    
    dept_iter = iter(random.sample(range(num_departments), num_departments))

    return [[next(dept_iter) for col in range(cols)] for row in range(rows)]

def cost(locations, flows):
    """
    Calculates the cost based on the rectilinear distance between the source
    and destination times the flow.
    """

    total_cost = 0
    
    # flow is symmetrical, so to avoid double-counting flow, we only count flow
    # from locations below each current location and exit the loop as soon as
    # it reaches the current location
    for r1, r1_depts in enumerate(locations):
        for c1, dept1 in enumerate(r1_depts):
            try:
                for r2, r2_depts in enumerate(locations):
                    for c2, dept2 in enumerate(r2_depts):
                        if r2 == r1 and c2 == c1:
                            # break out of two inner loops
                            raise StopIteration
                        else:
                            # the flows lookup table is a half-matrix, so
                            # we have to make sure we use the largest department
                            # for the row and the smallest for the column
                            lo, hi = ((dept1, dept2) if dept1 < dept2 
                                                            else (dept2, dept1))
                            dist = abs(r2-r1) + abs(c2-c1)
                    
                            # the half-matrix has no row for the first 
                            # department, so we subtract 1 from the dept number
                            # to get the correct row; we never have to worry
                            # about 0 being the hi_dept, since another
                            # department will always be higher and force 0 to
                            # the the lo_dept
                            total_cost += flows[hi-1][lo] * dist
            except StopIteration:
                continue
    
    return total_cost

def swap(locations, r1, c1, r2, c2):
    """
    Swaps the departments at the specified x, y coordinates in the locations
    grid.
    """
    
    locations[r1][c1], locations[r2][c2] = locations[r2][c2], locations[r1][c1]
            
def move(locations):
    """
    Perturb the department arrangement by swapping two department locations. 
    Returns a tuple containing the locations swapped for use with undo swap, if
    necessary.
    """
    
    r1 = random.choice(range(len(locations)))
    c1 = random.choice(range(len(locations[r1])))

    r2 = random.choice(range(len(locations)))
    c2 = random.choice(range(len(locations[r2])))

    while r1 == r2 and c1 == c2:
        r2 = random.choice(range(len(locations)))
        c2 = random.choice(range(len(locations[r2])))
    
    swap(locations, r1, c1, r2, c2)
    
    return (r1, c1, r2, c2)

def init_temperature(locations, flows, init_accept_rate):
    """
    Calculate the initial annealing temperature.
    
    Following Dreo, et al. (2006), calculate the average energy change over 100
    random moves. Derive init_temp from exp(-avg_change/init_temp) = tau_0, 
    where tau_0 is provided by the user. A tau_0 value of 0.50 represents an 
    assumed poor initial configuration, whereas a tau_0 value of 0.20 represents 
    an assumed good one.
    """
    
    delta_E = []
    for trial in range(100):
        start_cost = cost(locations, flows)
        move(locations)
        end_cost = cost(locations, flows)
        delta_E.append(abs(end_cost - start_cost))
        
        avg_delta_E = sum(delta_E) / len(delta_E)
    
    return -(avg_delta_E) / math.log(init_accept_rate)

def simulated_annealing(locations, temp, flows):
    minTemp = 1
    a= 0.95
    current = (locations,cost(locations,flows))
    while temp > minTemp:
        i = 1
        while i <= 55:
            x = move(locations)
            new = (locations, cost(locations,flows))
            delta = new[1] - current[1]
            ap = math.exp(delta/temp)
            #print(temp , "%%%%%%%%%%" ,current[1], "&&&&&&&", ap)
            if delta < 0 or ap < random.uniform(1,50):
                current = new
            else:
                swap(locations,x[0],x[1],x[2],x[3])
                current = (locations, cost(locations,flows))
            i += 1
        temp = temp*a
    return current[1]

def find_solution(locations,flows):
    #temp = init_temperature(locations, flows, 0.2) not useful
    result = simulated_annealing(locations,100,flows)
    x = 1
    while result > 576: 
        if x > 10: #run the SA code 9 times
            return -1
        #temp = init_temperature(locations, flows, 0.2)
        result = simulated_annealing(locations,100,flows)
        print("finding solution: ", result)
        if result <= 576:
            print("Good Job code!")
        elif result >= 600:
            print("bad...")
        else:
            print("better...")
        x += 1
    return result

def main():
    """
    Program entry point. Parses command line arguments and contains the main
    simulated annealing loop.
    """
    
    # Read flow data and generate initial department locations
    with open("input.txt") as infile:
        flows = init_flow(infile)
    x = 1
    num_departments = len(flows) + 1
    locations = init_locations(flows)

    # Implement SA algorithm here
    x = find_solution(locations,flows)
    while x == -1: #run until solution is found
        x = find_solution(locations,flows)
    return x
        
    
if __name__ == '__main__':
    main()
