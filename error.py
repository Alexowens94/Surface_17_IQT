from random import randint
from projectq.ops import X, Z, Y

class Error():
    """Error object"""
    def __init__(self, name="Unnamed error", location=None):
        self.name = name # Name of the error type
        self.location=location

    def insert_error(self):
        raise Exception("This error does not have a defined insert_error function.")

    def generate_location_list(self):
        raise Exception("This error does not have a defined generate_location_list function.")

    def generate_random_location(self, location_list):
        i = randint(0,len(location_list)-1)
        return location_list[i]

    @classmethod
    def generate_error_list(cls, error_subset, error_locations):
        # error_locations is a dict: e.g. { XCtrlError : [0,1,3,4] }
        # e.g. errors_to_generate = { "XCtrlError" : 5, "DephasingError" : 6 }
        error_list = []
        for error in error_subset.keys():
            # e.g. "XCtrlError"
            for i in range(error_subset[error]):
                # e.g. for i in range(5):
                error_i = eval(error)(name=error)
                location = list(error_i.generate_random_location(error_locations[type(error_i)]))
                error_i.location = location
                error_list.append(error_i)
        return error_list

class XCtrlError(Error):

    def __init__(self, name="XCtrlError", location=None):
        super().__init__(name=name, location=location)

    def insert_error(self, *args):
        for qubit in args:
            X | qubit
            #print("X error inserted.")

class DephasingError(Error):

    def __init__(self, name="DephasingError", location=None):
        super().__init__(name=name, location=location)

    def insert_error(self, *args):
        for qubit in args:
            Z | qubit
            #print("Z error inserted.")

class YCtrlError(Error):

    def __init__(self, name="YCtrlError", location=None):
        super().__init__(name=name, location=location)

    def insert_error(self, *args):
        for qubit in args:
            Y | qubit
            #print("Y error inserted.")

def create_error_subset_list():
    error_subset_list = []
    for i in range(6):
        error_subset = {"DephasingError" : i }
        error_subset_list.append(error_subset)
    return error_subset_list