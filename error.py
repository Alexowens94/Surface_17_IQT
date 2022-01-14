from random import randint
import itertools
from projectq.ops import X, Z

class Error():
    """Error object"""
    def __init__(self, name="Unnamed error", probability=None, location=None, gates = None):
        self.name = name # Name of the error type
        self.probability = probability # Probability of this error (not currently used)
        self.location = location # Location of the error - used in error subsets
        self.gates = gates # Gates which instigate this error, e.g. Rx gate can cause an XCtrlError (not currently used)
        self.location_list = self.generate_location_list()

    def insert_error(self):
        raise Exception("This error does not have a defined insert_error function.")

    def generate_location_list(self):
        raise Exception("This error does not have a defined generate_location_list function.")

    def generate_random_location(self):
        location_list = self.location_list
        i = randint(0,len(location_list)-1)
        return location_list[i]

    @classmethod
    def generate_error_list(cls, error_subset):
        # e.g. errors_to_generate = { "XCtrlError" : 5, "DephasingError" : 6 }
        error_list = []
        for error in error_subset.keys():
            # e.g. "XCtrlError"
            for i in range(error_subset[error]):
                # e.g. for i in range(5):
                error_i = eval(error)(name=error)
                location = list(error_i.generate_random_location())
                error_i.location = location
                error_list.append(error_i)
        return error_list

class XCtrlError(Error):

    def __init__(self, name="XCtrlError", probability=None, location=None, gates=None):
        super().__init__(name=name, probability=probability, location=location, gates=gates)

    def insert_error(self, *args):
        for qubit in args:
            X | qubit
            print("X error inserted.")

    @classmethod
    def generate_location_list(cls):
        #-----------------------------------------------------------------------#
        #                   LOCATION GENERATOR FLOWCHART                        #
        #-----------------------------------------------------------------------#
        #                                 Syndrome                              #
        #                           0                   1                       #
        #                       Timestep                                        #
        #               X-type              Z-type                              #
        #            0 1 2 3                    4 5 6 7                         #
        #           EntanglingOp                EntanglingOp                    #
        #         0       1       2             0       1       2               #
        #       cancel_data_rx                 cancel_data_rx                   #
        #   True          False             True          False                 #
        #   0             0   1               1           1   2                 #
        #                                   (only locations after Rx or Rxx)    #
        #-----------------------------------------------------------------------#

        position_list = [[0,1],[0,1,2,3,4,5,6,7],[0,1,2]] # X-entangling possible locations
        location_list = list((list(tup) for tup in itertools.product(*position_list)))
        # Example location_list = [[0,0,0], [0,0,1], [0,0,2], [0,1,0], [0,1,1]...]
        cancel_data_rx = [[True, False, False],
                          [False, True, True],
                          [True, False, True],
                          [False, True, False],
                          [True, False, True],
                          [False, True, False],
                          [True, False, False],
                          [False, True, True]]
        new_location_list = []
        for l in location_list:
            if 0 <= l[1] <= 3:
                complete_location_0 = l + [0]
                new_location_list.append(complete_location_0)
                # cancel_data_rx[stabilizer_timestep_position][entangling_op_position]
                if not cancel_data_rx[l[1]][l[2]]:
                    complete_location_1 = l + [1]
                    new_location_list.append(complete_location_1)
            if 4 <= l[1] <= 7:
                # cancel_data_rx[stabilizer_timestep_position][entangling_op_position]
                complete_location_1 = l + [1]
                new_location_list.append(complete_location_1)
                if not cancel_data_rx[l[1]][l[2]]:
                    complete_location_2 = l + [2]
                    new_location_list.append(complete_location_2)
        return new_location_list

class DephasingError(Error):

    def __init__(self, name="DephasingError", probability=None, location=None, gates=None):
        super().__init__(name=name, probability=probability, location=location, gates=gates)

    def insert_error(self, *args):
        for qubit in args:
            Z | qubit
            print("Z error inserted.")

    @classmethod
    def generate_location_list(cls):
        #-----------------------------------------------------------------------#
        #                   LOCATION GENERATOR FLOWCHART                        #
        #-----------------------------------------------------------------------#
        #                                 Syndrome                              #
        #                           0                   1                       #
        #                       Timestep                                        #
        #               X-type              Z-type                              #
        #            0 1 2 3                    4 5 6 7                         #
        #           EntanglingOp                EntanglingOp                    #
        #         0       1       2             0       1       2               #
        #       cancel_data_rx                 cancel_data_rx                   #
        #   True          False             True          False                 #
        #   0             0   1             0 1           0   1   2   3         #
        #                                       (locations after any gate)      #
        #-----------------------------------------------------------------------#

        position_list = [[0,1],[0,1,2,3,4,5,6,7],[0,1,2]] # X-entangling possible locations
        location_list = list((list(tup) for tup in itertools.product(*position_list)))
        # Example location_list = [[0,0,0], [0,0,1], [0,0,2], [0,1,0], [0,1,1]...]
        cancel_data_rx = [[True, False, False],
                          [False, True, True],
                          [True, False, True],
                          [False, True, False],
                          [True, False, True],
                          [False, True, False],
                          [True, False, False],
                          [False, True, True]]
        new_location_list = []
        for l in location_list:
            if 0 <= l[1] <= 3:
                complete_location_0 = l + [0]
                new_location_list.append(complete_location_0)
                # cancel_data_rx[stabilizer_timestep_position][entangling_op_position]
                if not cancel_data_rx[l[1]][l[2]]:
                    complete_location_1 = l + [1]
                    new_location_list.append(complete_location_1)
            if 4 <= l[1] <= 7:
                # cancel_data_rx[stabilizer_timestep_position][entangling_op_position]
                complete_location_0 = l + [0]
                complete_location_1 = l + [1]
                new_location_list.append(complete_location_0)
                new_location_list.append(complete_location_1)
                if not cancel_data_rx[l[1]][l[2]]:
                    complete_location_2 = l + [2]
                    complete_location_3 = l + [3]
                    new_location_list.append(complete_location_2)
                    new_location_list.append(complete_location_3)
        return new_location_list

def create_error_subset_list():
    error_subset_list = []
    for i in range(6):
        error_subset = {"DephasingError" : i }
        error_subset_list.append(error_subset)
    return error_subset_list