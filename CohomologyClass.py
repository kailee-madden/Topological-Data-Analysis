import numpy as np
import itertools
import matplotlib.pyplot as plt

class Cohomology(object):
    '''Calculates the persistence cohomology of a simplicial complex by using cochains, 
    cocycles, and coboundaries. Cohomology is isomorphic to homology.'''

    def __init__(self):
        self.simplices = self.simplices = {"0" : set(), "1" : set(), "2" : set(), "3" : set()}
        self.all_simplices = []
        self.Z = []
        self.I = {0 : {}, 1 : {}, 2 : {}} 
        self.barcodes = {0 : {}, 1 : {}, 2 : {}} 

    def calculate(self, highest_simplices):
        dim = len(highest_simplices[0])-1
        self.simplices[str(dim)] = highest_simplices

        while dim > 0:
            dim -= 1
            self.get_other_simplices(dim)
        
        self.cohomology()
        print(self.I)

    def get_other_simplices(self, dimension):
        '''Gets lower simplices'''
        for simp in self.simplices[str(dimension+1)]:
            for temp in itertools.combinations(simp, dimension+1): #get all the combinations of each tuple
                self.simplices[str(dimension)].add(temp)
    
    def cohomology(self):
        '''Tracks cocycles and the intervals in which they live.
        Creates barcodes for these intervals for later visualization.'''
        for value in self.simplices.values(): #loop through the dictionary of total simplices (ie 0-simplices comes first then 1-simplices etc.)
            v = sorted(list(value)) #sort the set from small to large
            for simplex in v: #loop through the ordered list of simplices to get each individual simplex
                self.all_simplices.append(simplex) #add all individual simplices to a total list of simplices 
        for simp_location, simplex in enumerate(self.all_simplices): #loop through all the simplices
            dim = len(simplex) - 1 #get the dimension of the simplex
            if dim == 0: #if it is a 0-simplex then it will not be the coboundary of anything
                self.Z.append(simplex) #add it to Z
                self.I[dim][simplex] = "infinity" #add it to I with no end on its interval
                self.barcodes[dim][simp_location] = "infinity" #add it to the barcodes (location is the relevant piece)
            else: #if not a 0-simplex then have to check if it kills any cocycles
                boundary = []
                for combination in itertools.combinations(simplex, dim): #get the lower simplices possible to make the boundary
                    boundary.append(combination)
                found_youngest = False
                for cocycle_location, cocycle in reversed(list(enumerate(self.Z))): #loop in reverse but have normal indices
                    if cocycle in boundary and found_youngest == False: #check that this is the coboundary of the cocycle and that it is the youngest cocycle
                        temp = cocycle
                        del self.Z[cocycle_location] #kill the youngest cocycle
                        self.Z.append(simplex) #add the dead simplex to Z
                        new_interval = {cocycle : simplex}
                        non_interval = {simplex : simplex}
                        self.I[dim-1].update(new_interval) #update the cocycle's interval with its death location
                        self.I[dim].update(non_interval) #add the immediately killed interval
                        for cocy_loc, cocy in enumerate(self.all_simplices): #need to have number for the cocycle
                            if cocy == cocycle:
                                self.barcodes[dim-1][cocy_loc] = simp_location #update the barcodes with the cocycle's associated number
                                break #break out of the for loop since accomplished goal
                        found_youngest = True
                    elif cocycle in boundary and found_youngest == True: #check that this is the coboundary of the cocycle and that we already found the youngest one
                        self.Z[cocycle_location] = ("{}+{}".format(cocycle, temp)) #add the killed cocycle to any other cocycles in the boundary
                if found_youngest == False: #we had no cocycles in the boundary
                    self.Z.append(simplex) #add the simplex to the list of cocycles
                    self.I[dim][simplex] = "infinity" #add an associated interval going to infinity
                    self.barcodes[dim][simp_location] = "infinity"

    def visualize_barcodes(self):
        '''Visualizes barcodes using matplotlib'''
        x_max = len(self.all_simplices) + 1
        print("If line goes to {}, indicates cycle is never killed (goes to infinity)". format(x_max))
        for dim, value in self.barcodes.items(): #creates a new plot for each dimension
            plt.figure()
            plt.title("Barcode for {} Homology".format(dim))
            y = 1
            for start, end in value.items(): #loops through the intervals for the dimension
                if end == "infinity":
                    end = x_max
                x = [start, end]
                y_vals = [y,y]
                plt.plot(x, y_vals) #add the interval to the plot
                y+= 1
            plt.ylim(0, y)
            plt.gca().axes.get_yaxis().set_visible(False) #hides the y axis since it is irrelevant for our purposes
            plt.xlim(0, x_max)
        plt.show()