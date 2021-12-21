import itertools
from random import randint
from compiled_surface_code_error_subsets import XCtrlError, DephasingError

def generate_location_list():
    list_X = [[0,1],[0,1,2,3],[0,1,2],[0,1]] # X-entangling possible locations
    location_list_X = list(itertools.product(*list_X))
    list_Z = [[0,1],[4,5,6,7],[0,1,2],[0,1,2,3]] # Z-entangling possible locations
    location_list_Z = list(itertools.product(*list_Z))
    location_list = location_list_X + location_list_Z
    return location_list

def generate_random_location():
    location_list = generate_location_list()
    i = randint(0,len(location_list)-1)
    return location_list[i]

def generate_error_subset(errors_to_generate):
    error_subset = []
    for error in errors_to_generate.keys():
        for i in range(errors_to_generate[error]):
            location = list(generate_random_location())
            error_i = eval(error)(name=error, location=location)
            error_subset.append(error_i)
    return error_subset