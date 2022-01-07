from random import randint
import numpy as np
import json
import itertools
from math import sqrt, pi
from projectq.ops import All, Measure, X, Y, Z, H, Rx, Ry, Rxx

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
    def generate_error_subset(cls, errors_to_generate):
        # e.g. errors_to_generate = { "XCtrlError" : 5, "DephasingError" : 6 }
        error_subset = []
        for error in errors_to_generate.keys():
            # e.g. "XCtrlError"
            for i in range(errors_to_generate[error]):
                # e.g. for i in range(5):
                error_i = eval(error)(name=error)
                location = list(error_i.generate_random_location())
                error_i.location = location
                error_subset.append(error_i)
        return error_subset

class XCtrlError(Error):

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

#
#                 Z ancilla0
#           data0-----------data1-----------data2
#             |               |               |
#             |    X ancilla1 |    Z ancilla2 |  X ancilla3
#             |               |               |
#           data3-----------data4-----------data5
#             |               |               |
# X ancilla4  |   Z ancilla5  |    X ancilla6 |
#             |               |               |
#           data6-----------data7------------data8
#
#                                 Z ancilla7

class EntanglingOperation():

    def __init__(self, location=None):
        self.location = location

    def x_type_entangling(self,dataq,ancillaq,cancel_data_rx=None,s=1,error_subset=None):
        """
        CNOT gate, as it appears in x-stabilizer compiled to native operations. All ancilla single qubits ops
        canceled 'by hand'.
        cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
        on the same data qubit (depedent on choices of s and v)
        """
        Rxx(s*pi/2) | (dataq,ancillaq)
        this_location = self.location + [0]
        for error in error_subset:
            if this_location == error.location:
                error.insert_error(dataq,ancillaq)
                print("error {} inserted at {}".format(error.name, this_location))
        if not cancel_data_rx:
            Rx(-s*pi/2) | dataq
            this_location = self.location + [1]
            for error in error_subset:
                if this_location == error.location:
                    error.insert_error(dataq)
                    print("error {} inserted at {}".format(error.name,this_location))

    def z_type_entangling(self,dataq,ancillaq,cancel_data_rx=None,s=1,v=1,error_subset=None):
        """
        CNOT gate, as it appears in z-stabilizer compiled to native operations. All ancilla single qubits ops
        canceled 'by hand'.
        cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
        on the same data qubit (depedent on choices of s and v)
        """
        Ry(v*pi/2) | dataq
        this_location = self.location + [0]
        for error in error_subset:
            if this_location == error.location:
                error.insert_error(dataq)
                print("error {} inserted at {}".format(error.name, this_location))
        Rxx(s * pi / 2) | (dataq, ancillaq)
        this_location = self.location + [1]
        for error in error_subset:
            if this_location == error.location:
                error.insert_error(dataq, ancillaq)
                print("error {} inserted at {}".format(error.name, this_location))
        if not cancel_data_rx:
            Rx(-s * pi / 2) | dataq
            this_location = self.location + [2]
            for error in error_subset:
                if this_location == error.location:
                    error.insert_error(dataq)
                    print("error {} inserted at {}".format(error.name,this_location))
        Ry(v * pi / 2) | dataq
        this_location = self.location + [3]
        for error in error_subset:
            if this_location == error.location:
                error.insert_error(dataq)
                print("error {} inserted at {}".format(error.name, this_location))

class StabiliserTimestep():

    def __init__(self, location=None, qu_ind=None, entangling_type = "X"):
        """
        attr:
            self.location (e.g. [0,0]) [Syndrome measurement location, stabilizer location]
            self.qu_ind (e.g. [[4,1], [8,6], [6,4]] ) The qubit indices needed for stabiliser_timesteps 1-8
                i.e. stabiliser_timestep_1 entangles data[4] & ancilla[1], then data[8] & ancilla[6] etc.
            self.ent_type ("X" or "Z") specifies whether this stabilizer timestep will do X-type or Z-type entangling
        """
        self.location = location
        self.qu_ind = qu_ind
        self.ent_type = entangling_type

    def run(self, data, ancilla, error_subset, cancel_data_rx=None):
        if self.ent_type == "X":
            loc_0 = self.location + [0]
            ent_op_0 = EntanglingOperation(location=loc_0)
            ent_op_0.x_type_entangling(data[self.qu_ind[0][0]], ancilla[self.qu_ind[0][1]],cancel_data_rx=cancel_data_rx[0],s=1,error_subset=error_subset)
            loc_1 = self.location + [1]
            ent_op_1 = EntanglingOperation(location=loc_1)
            ent_op_1.x_type_entangling(data[self.qu_ind[1][0]], ancilla[self.qu_ind[1][1]],cancel_data_rx=cancel_data_rx[1],s=1,error_subset=error_subset)
            loc_2 = self.location + [2]
            ent_op_2 = EntanglingOperation(location=loc_2)
            ent_op_2.x_type_entangling(data[self.qu_ind[2][0]], ancilla[self.qu_ind[2][1]],cancel_data_rx=cancel_data_rx[2],s=1,error_subset=error_subset)
        elif self.ent_type == "Z":
            loc_0 = self.location + [0]
            ent_op_0 = EntanglingOperation(location=loc_0)
            ent_op_0.z_type_entangling(data[self.qu_ind[0][0]], ancilla[self.qu_ind[0][1]],cancel_data_rx=cancel_data_rx[0],s=1,error_subset=error_subset)
            loc_1 = self.location + [1]
            ent_op_1 = EntanglingOperation(location=loc_1)
            ent_op_1.z_type_entangling(data[self.qu_ind[1][0]], ancilla[self.qu_ind[1][1]],cancel_data_rx=cancel_data_rx[1],s=1,error_subset=error_subset)
            loc_2 = self.location + [2]
            ent_op_2 = EntanglingOperation(location=loc_2)
            ent_op_2.z_type_entangling(data[self.qu_ind[2][0]], ancilla[self.qu_ind[2][1]],cancel_data_rx=cancel_data_rx[2],s=1,error_subset=error_subset)
        else:
            raise Exception("Entangling type not recognized, should be \"X\" or \"Z\". Please check stabilizer_timestep instantiation.")

class StabiliserCycle():
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    def __init__(self, location = None):
        self.location = location

    def run(self, data, ancilla, eng, reset=True, error_subset=None):
        if len(data)!=9:
            raise Exception('data qubit register does not correspond to the surface 17 QEC code')
        loc_0 = self.location + [0]
        loc_1 = self.location + [1]
        loc_2 = self.location + [2]
        loc_3 = self.location + [3]
        loc_4 = self.location + [4]
        loc_5 = self.location + [5]
        loc_6 = self.location + [6]
        loc_7 = self.location + [7]

        # The qubit indices needed for stabiliser_timesteps 1-8
        # i.e. stabiliser_timestep_1 entangles data[4] & ancilla[1], then data[8] & ancilla[6] etc.
        qu_ind_1 = [[4, 1], [8, 6], [6, 4]]
        qu_ind_2 = [[1, 1], [5, 6], [3, 4]]
        qu_ind_3 = [[3, 1], [7, 6], [5, 3]]
        qu_ind_4 = [[0, 1], [4, 6], [2, 3]]
        qu_ind_5 = [[1, 2], [3, 5], [7, 7]]
        qu_ind_6 = [[2, 2], [4, 5], [8, 7]]
        qu_ind_7 = [[4, 2], [6, 5], [0, 0]]
        qu_ind_8 = [[5, 2], [7, 5], [1, 0]]

        cancel_data_rx = [[True, False, False],
                          [False, True, True],
                          [True, False, True],
                          [False, True, False],
                          [True, False, True],
                          [False, True, False],
                          [True, False, False],
                          [False, True, True]]

        stabiliser_timestep_1 = StabiliserTimestep(location=loc_0, qu_ind=qu_ind_1, entangling_type="X")
        stabiliser_timestep_2 = StabiliserTimestep(location=loc_1, qu_ind=qu_ind_2, entangling_type="X")
        stabiliser_timestep_3 = StabiliserTimestep(location=loc_2, qu_ind=qu_ind_3, entangling_type="X")
        stabiliser_timestep_4 = StabiliserTimestep(location=loc_3, qu_ind=qu_ind_4, entangling_type="X")
        stabiliser_timestep_5 = StabiliserTimestep(location=loc_4, qu_ind=qu_ind_5, entangling_type="Z")
        stabiliser_timestep_6 = StabiliserTimestep(location=loc_5, qu_ind=qu_ind_6, entangling_type="Z")
        stabiliser_timestep_7 = StabiliserTimestep(location=loc_6, qu_ind=qu_ind_7, entangling_type="Z")
        stabiliser_timestep_8 = StabiliserTimestep(location=loc_7, qu_ind=qu_ind_8, entangling_type="Z")

        stabiliser_timestep_1.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[0])
        stabiliser_timestep_2.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[1])
        stabiliser_timestep_3.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[2])
        stabiliser_timestep_4.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[3])
        stabiliser_timestep_5.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[4])
        stabiliser_timestep_6.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[5])
        stabiliser_timestep_7.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[6])
        stabiliser_timestep_8.run(data,ancilla,error_subset, cancel_data_rx=cancel_data_rx[7])

        All(Measure) | ancilla
        eng.flush()
        syndrome_t = [int(q) for q in ancilla]
        if reset:
            for a in ancilla: #reset the ancillas to 0 at end of stab round (allow for repeat rounds)
                if int(a) == 1:
                    X | a

        return syndrome_t

