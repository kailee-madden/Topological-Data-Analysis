import numpy as np
from numpy.linalg import matrix_rank
import itertools

class Homology(object):

    def __init__(self):
        self.simplices = {"0" : set(), "1" : set(), "2" : set(), "3" : set()}
        self.boundary_matrices = {"0": 0, "1" : 0, "2" : 0, "3" : 0}
        self.homology_dimensions = {"0": 0, "1" : 0, "2" : 0}
    
    def calculations(self, highest_simplices):
        '''Finds the dimension of the first three homologies of a simplicial complex given the highest existing simplices'''
        dim = len(highest_simplices[0])-1
        self.simplices[str(dim)] = highest_simplices

        while dim > 0:
            dim -= 1
            self.get_other_simplices(dim)

        for key in list(self.boundary_matrices):
            if key == "0":
                self.generalized_matrices(self.simplices[key], 1, key)
            else:
                self.generalized_matrices(self.simplices[key], self.simplices[str(int(key)-1)], key)
            
        for key in list(self.homology_dimensions):
            self.find_dim_homology(self.boundary_matrices[key], self.boundary_matrices[str(int(key)+1)], key)
        
        return
    

    def get_other_simplices(self, dimension):
        '''Gets lower simplices'''
        for simp in self.simplices[str(dimension+1)]:
            for temp in itertools.combinations(simp, dimension+1): #get all the combinations of each tuple
                self.simplices[str(dimension)].add(temp)
        
    def generalized_matrices(self, columns, rows, key):
        #create boundary matrix with all zeros
        if len(columns) > 0 and type(rows) == int:
            boundary_matrix = np.zeros((rows, len(columns)))
            boundary = False
        elif len(columns) > 0 and len(rows) > 0:
            boundary_matrix = np.zeros((len(rows), len(columns)))
            boundary = True
        else:
            boundary = False
            boundary_matrix = 0 #if we don't have the necessary simplices to construct a boundary matrix then we just set it equal to zero

        #fill in the correct spots with 1s and -1s
        if boundary == True:
            n = 0
            for row in rows:
                m = 0
                for column in columns:
                    if set(row).issubset(column): #check to see if the row tuple is contained in the column tuple
                        for element in column:
                            if element not in row: #find the element that is in the column tuple but not the row tuple
                                missing = column.index(element) #find the missing element's location
                                if missing % 2 == 0: #determine if the element's location is an even or odd spot then use the boundary formula to add a 1 or -1
                                    boundary_matrix[n,m] = 1
                                else:
                                    boundary_matrix[n,m] = -1
                    m += 1
                n += 1
        self.boundary_matrices[key] = boundary_matrix

    def find_dim_homology(self, boundary_matrix, higher_boundary_matrix, key):
        try: #even if the higher_boundary_matrix is 0 this will work since matrix_rank(0) = 0
            dim_ker = np.shape(boundary_matrix)[1] - matrix_rank(boundary_matrix) #use rank nullity theorem so don't have to find actual nullspace
            dim_im = matrix_rank(higher_boundary_matrix)
            self.homology_dimensions[key] = dim_ker - dim_im
        except ValueError: 
            print("There are no simplices of the necessary dimension or they were entered improperly.")
            return
        
    def give_dimension(self):
        print("The dimension of the zero homology is {}.\nThe dimension of the one homology is {}.\nThe dimension of the two homology is {}.".format(self.homology_dimensions["0"], self.homology_dimensions["1"], self.homology_dimensions["2"]))
        return