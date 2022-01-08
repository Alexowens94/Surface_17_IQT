import json
from math import pi

from projectq import MainEngine
from projectq.backends import Simulator
from projectq.ops import All, Measure, X, Z, H, Rx, Ry, Rxx


class EntanglingOperation():

    def __init__(self, location=None):
        self.location = location

    def x_type_entangling(self, dataq, ancillaq, cancel_data_rx=None, s=1, error_list=None):
        """
        CNOT gate, as it appears in x-stabilizer compiled to native operations. All ancilla single qubits ops
        canceled 'by hand'.
        cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
        on the same data qubit (depedent on choices of s and v)
        """
        Rxx(s * pi / 2) | (dataq, ancillaq)
        this_location = self.location + [0]
        for error in error_list:
            if this_location == error.location:
                error.insert_error(dataq, ancillaq)
                print("error {} inserted at {}".format(error.name, this_location))
        if not cancel_data_rx:
            Rx(-s * pi / 2) | dataq
            this_location = self.location + [1]
            for error in error_list:
                if this_location == error.location:
                    error.insert_error(dataq)
                    print("error {} inserted at {}".format(error.name, this_location))

    def z_type_entangling(self, dataq, ancillaq, cancel_data_rx=None, s=1, v=1, error_list=None):
        """
        CNOT gate, as it appears in z-stabilizer compiled to native operations. All ancilla single qubits ops
        canceled 'by hand'.
        cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
        on the same data qubit (depedent on choices of s and v)
        """
        Ry(v * pi / 2) | dataq
        this_location = self.location + [0]
        for error in error_list:
            if this_location == error.location:
                error.insert_error(dataq)
                print("error {} inserted at {}".format(error.name, this_location))
        Rxx(s * pi / 2) | (dataq, ancillaq)
        this_location = self.location + [1]
        for error in error_list:
            if this_location == error.location:
                error.insert_error(dataq, ancillaq)
                print("error {} inserted at {}".format(error.name, this_location))
        if not cancel_data_rx:
            Rx(-s * pi / 2) | dataq
            this_location = self.location + [2]
            for error in error_list:
                if this_location == error.location:
                    error.insert_error(dataq)
                    print("error {} inserted at {}".format(error.name, this_location))
        Ry(v * pi / 2) | dataq
        this_location = self.location + [3]
        for error in error_list:
            if this_location == error.location:
                error.insert_error(dataq)
                print("error {} inserted at {}".format(error.name, this_location))


class StabiliserTimestep():

    def __init__(self, location=None, qu_ind=None, entangling_type="X"):
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

    def run(self, data, ancilla, error_list, cancel_data_rx=None):
        if self.ent_type == "X":
            loc_0 = self.location + [0]
            ent_op_0 = EntanglingOperation(location=loc_0)
            ent_op_0.x_type_entangling(data[self.qu_ind[0][0]], ancilla[self.qu_ind[0][1]],
                                       cancel_data_rx=cancel_data_rx[0], s=1,
                                       error_list=error_list)
            loc_1 = self.location + [1]
            ent_op_1 = EntanglingOperation(location=loc_1)
            ent_op_1.x_type_entangling(data[self.qu_ind[1][0]], ancilla[self.qu_ind[1][1]],
                                       cancel_data_rx=cancel_data_rx[1], s=1,
                                       error_list=error_list)
            loc_2 = self.location + [2]
            ent_op_2 = EntanglingOperation(location=loc_2)
            ent_op_2.x_type_entangling(data[self.qu_ind[2][0]], ancilla[self.qu_ind[2][1]],
                                       cancel_data_rx=cancel_data_rx[2], s=1,
                                       error_list=error_list)
        elif self.ent_type == "Z":
            loc_0 = self.location + [0]
            ent_op_0 = EntanglingOperation(location=loc_0)
            ent_op_0.z_type_entangling(data[self.qu_ind[0][0]], ancilla[self.qu_ind[0][1]],
                                       cancel_data_rx=cancel_data_rx[0], s=1,
                                       error_list=error_list)
            loc_1 = self.location + [1]
            ent_op_1 = EntanglingOperation(location=loc_1)
            ent_op_1.z_type_entangling(data[self.qu_ind[1][0]], ancilla[self.qu_ind[1][1]],
                                       cancel_data_rx=cancel_data_rx[1], s=1,
                                       error_list=error_list)
            loc_2 = self.location + [2]
            ent_op_2 = EntanglingOperation(location=loc_2)
            ent_op_2.z_type_entangling(data[self.qu_ind[2][0]], ancilla[self.qu_ind[2][1]],
                                       cancel_data_rx=cancel_data_rx[2], s=1,
                                       error_list=error_list)
        else:
            raise Exception(
                "Entangling type not recognized, should be \"X\" or \"Z\". Please check stabilizer_timestep instantiation.")


class StabiliserCycle():
    def __init__(self, location=None):
        '''
        attr:
            location = [0]
        '''
        self.location = location

    def run(self, data, ancilla, eng, reset=True, error_list=[]):
        if len(data) != 9:
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

        stabiliser_timestep_1.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[0])
        stabiliser_timestep_2.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[1])
        stabiliser_timestep_3.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[2])
        stabiliser_timestep_4.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[3])
        stabiliser_timestep_5.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[4])
        stabiliser_timestep_6.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[5])
        stabiliser_timestep_7.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[6])
        stabiliser_timestep_8.run(data, ancilla, error_list, cancel_data_rx=cancel_data_rx[7])

        All(Measure) | ancilla
        eng.flush()
        syndrome_t = [int(q) for q in ancilla]
        if reset:
            for a in ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
                if int(a) == 1:
                    X | a

        return syndrome_t


class LogicalQubit():

    def __init__(self, correction_table=None):
        self.eng = MainEngine(Simulator())
        self.data = self.eng.allocate_qureg(9)
        self.ancilla = self.eng.allocate_qureg(8)
        self.state = 1
        self.basis = 'X'
        quiescent_cycle = StabiliserCycle(location=[0])
        self.quiescent_state = quiescent_cycle.run(self.data, self.ancilla, self.eng)
        self.leaked_q_reg = 17 * [0]
        self.correction_table = load_lookup_table("correction_table_depolarising.json")


    def measure_syndrome(self, location, error_list):
        stabiliser_cycle = StabiliserCycle(location=location)
        syndrome = stabiliser_cycle.run(self.data, self.ancilla, self.eng, error_list=error_list)
        return syndrome

    def lookup(self, syndrome, display=False):
        """

        Args:
            syndrome (np.array): The measured syndrome.
            table (string): The filename of containing
                            the syndrome to error look-up
                            table.
            display (bool): Choice whether to print the
                            fault syndrome and error vector.

        Returns:
            error_vec (list): The errors to correct.
        """
        key = str(syndrome).strip('[,]')
        error_vec = self.correction_table[key][0]
        if display:
            print('ft syndrome: {}'.format(syndrome))
            print('error vector: {}'.format(error_vec))
        return error_vec

    def apply_correction(self, error_vec):
        for i in range(9):
            if error_vec[i] == 1:
                X | self.data[i]
            if error_vec[i + 9] == 1:
                Z | self.data[i]
        return

    def measure_qubit(self):
        if self.basis == 'X':
            All(H) | self.data  # change Z -> X basis

        All(Measure) | self.data
        self.eng.flush()  # flush all gates (and execute measurements)
        data_meas = [int(q) for q in self.data]
        logic_measurement = sum(data_meas) % 2

        return logic_measurement


def load_lookup_table(filename):
    with open(filename, 'r') as infile:
        correction_table = json.load(infile)
    return correction_table
