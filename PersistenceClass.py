import numpy as np
import itertools
import matplotlib.pyplot as plt

class Persistence(object):
    '''Calculates the persistence homology of a given simplicial complex'''

    def __init__(self):
        self.simplices = {"0" : set(), "1" : set(), "2" : set(), "3" : set()}
        self.all_simplices = []
        self.boundary_matrix = 0
        self.V = 0
        self.R = 0
        self.barcodes = {0 : {}, 1 : {}, 2 : {}} 
    
    def calculate(self, highest_simplices):
        dim = len(highest_simplices[0])-1
        self.simplices[str(dim)] = highest_simplices

        while dim > 0:
            dim -= 1
            self.get_other_simplices(dim)

        self.make_boundary_matrix()
        self.combined_boundary_matrix()
        self.get_barcodes()
        print(self.barcodes)

    def get_other_simplices(self, dimension):
        '''Gets lower simplices'''
        for simp in self.simplices[str(dimension+1)]:
            for temp in itertools.combinations(simp, dimension+1): #get all the combinations of each tuple
                self.simplices[str(dimension)].add(temp)
        
    
    def make_boundary_matrix(self): 
        '''Create the basic boundary matrix using mod 2'''
        for value in self.simplices.values(): #loop through the dictionary of total simplices (ie 0-simplices comes first then 1-simplices etc.)
            v = sorted(list(value)) #sort the set from small to large
            for simplex in v: #loop through the ordered list of simplices to get each individual simplex
                self.all_simplices.append(simplex) #add all individual simplices to a total list of sets 

        if len(self.all_simplices) == 0: #check to make sure there are simplices
            return 0

        self.boundary_matrix = np.zeros((len(self.all_simplices), len(self.all_simplices))) #create a boundary matrix of all zeros

        for c, column in enumerate(self.all_simplices): #loop through columns and indices
            for r, row in enumerate(self.all_simplices):
                if set(row).issubset(column) and len(row)+1 == len(column): #check to see if the row tuple is contained in the column tuple 
                    #and check that the dimension of the row simplex is one less than the dimension of the column simplex
                    self.boundary_matrix[r,c] = 1

    def combined_boundary_matrix(self):
        '''Create the modified boundary matrix'''
        self.V = np.identity(len(self.all_simplices)) #create the identity matrix
        self.R = self.boundary_matrix #create the matrix we will be modifying

        column_wise_sum = np.sum(self.R, axis=0) #get the column sums
        for c in range(0, len(self.all_simplices)): #loop through the columns' indices
            if column_wise_sum[c] != 0: #if it is not a cycle then we will check for pivots in common with other columns
                for r in range(len(self.all_simplices)-1, -1, -1): #loop backwards so first 1 found is a pivot
                    if self.R[r, c] == 1: #if this row has the pivot
                        for previous_col in range(0, c): #loop through all columns less than the current one
                            if self.R[r, previous_col] == 1: #check if this is a 1
                                pivot = True
                                #now we check if it is a pivot ie are there lower rows that have a 1 in this new column
                                for i in range(r+1, len(self.all_simplices)): 
                                    if self.R[i, previous_col] == 1:
                                        pivot = False
                                        break #if so then we skip onward
                                if pivot == True: #otherwise this is a pivot meaning we must add the two matrix columns together
                                    for row in range(0, len(self.all_simplices)):
                                        temp = self.R[row, c] + self.R[row, previous_col]
                                        self.R[row, c] = temp % 2 #add the columns using mod 2
                                        temp2 = self.V[row, c] + self.V[row,previous_col]
                                        self.V[row,c] = temp2 % 2 #add the columns of the identity matrix using mod 2
                                break #break out of the smaller loop because we need to restart to search for a new pivot in this modified column



                            

    def get_barcodes(self):
        '''Gets the barcodes by searching for cycles and when they are killed.
        Using the simplex index in the matrix instead of the actual simplex itself
        in order to make visualization simpler.But can be easily modified to give 
        the actual simplices if desired.'''
        cws = np.sum(self.boundary_matrix, axis=0) #get the column sums
        for index, col in enumerate(self.all_simplices): #loop through columns and indices
            if cws[index] == 0: #check if the column (ie simplex) is a cycle
                dim = len(col) -1 #then calculate the associated dimension
                self.barcodes[dim][index] = "infinity" #add the cycle to the list of barcodes
            else: #if not a cycle then its pivot is killing a cycle
                for row_index, row in reversed(list(enumerate(self.all_simplices))): #loop through to find the pivot, so go in reverse but keep normal indices
                    if self.R[row_index, index] == 1:
                        dim = len(row) -1 #get the dimension of the row simplex because this is the cycle being killed
                        self.barcodes[dim][row_index] = index #modify the barcode to have the cycle killed in this column
                        break #break out of the loop since we have found our pivot and modified the barcode as needed
    
    def visualize_barcodes(self):
        x_max = len(self.all_simplices) + 1
        print("If line goes to {}, indicates cycle is never killed (goes to infinity)". format(x_max))
        for dim, value in self.barcodes.items():
            plt.figure() #creates a new plot for each dimension
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