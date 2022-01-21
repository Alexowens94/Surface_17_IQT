import json
from math import pi

from projectq import MainEngine
from projectq.backends import Simulator
from projectq.ops import All, Measure, X, Z, H, Rx, Ry, Rxx, Rzz
from error import *


pdd_error_model = {Rxx : [XCtrlError, DephasingError], Rx : [XCtrlError, DephasingError], Ry : [YCtrlError, DephasingError]}

#dressed_state_error_model = {Rxx : [XCtrlError, LeakageError], Rx : [XCtrlError, LeakageError], Ry : [YCtrlError, LeakageError]}

class NoisyGate():
    """
    NoisyGate objects hold the information of the location of
    a gate, and when applied they insert the appropriate errors
    for that gate. Depending on the circuit-scheme and the
    error model in use.
    """
    def __init__(self, gate, qubits, location=None, error_model={}):
        self.gate = gate
        self.qubits = qubits
        self.location = location
        self.error_model = error_model
        self.possible_errors = self.get_errors_from_error_model()

    def __str__(self):
        return "NoisyGate: \n Gate: {}, \n Qubits : {}, \n Location: {}, \n Error model: {}, \n " \
               "Errors, {}".format(self.gate, self.qubits, self.location, self.error_model, self.possible_errors)

    def get_errors_from_error_model(self):
        errors = []
        for gate in self.error_model.keys():
            if type(self.gate) == gate:
                for error in self.error_model[gate]:
                    errors.append(error())
        return errors

    def apply(self, error_list=[]):
        self.gate | (self.qubits)
        for error in error_list:
            if error.location == self.location:
                error.insert_error(*self.qubits)
                print("Error: {} inserted at {}".format(error.name, self.location))


class EntanglingOperation():

    """An entangling operation entangles one data qubit and one ancilla qubit.
    This class is not meant to be used directly but is inherited by classes such
    as XEntanglingOp."""

    def __init__(self, location=None, dataq=None, ancillaq=None, cancel_data_rx=None, s=None, v=None, error_model=None):
        self.location = location
        self.dataq = dataq
        self.ancillaq = ancillaq
        self.error_model = error_model
        self.cancel_data_rx = cancel_data_rx
        self.s=s
        self.v=v
        self.noisy_gates = self.build_noisy_gates()

    def build_noisy_gates(self):
        raise Exception("No build_noisy_gates function defined.")

    def run(self, error_list=[]):
        for noisy_gate in self.noisy_gates:
            noisy_gate.apply(error_list=error_list)


class XTypeEntanglingOp(EntanglingOperation):
    """A type of entangling operation which is used by the first four stabiliser steps in the MS-circuit-scheme."""

    def build_noisy_gates(self):
        noisy_gates = []
        rxx_gate = NoisyGate(Rxx(self.s * pi / 2), [self.dataq, self.ancillaq], self.location + [0], error_model=self.error_model)
        noisy_gates.append(rxx_gate)
        if not self.cancel_data_rx:
            rx_gate = NoisyGate(Rx(-self.s * pi / 2), [self.dataq], self.location + [1], error_model=self.error_model)
            noisy_gates.append(rx_gate)
        return noisy_gates


class ZTypeEntanglingOp(EntanglingOperation):
    """A type of entangling operation which is used by the last four stabiliser steps in the MS-circuit-scheme."""
    def build_noisy_gates(self):
        noisy_gates = []
        ry_gate = NoisyGate(Ry(self.v*pi/2), [self.dataq], self.location + [0], error_model=self.error_model)
        noisy_gates.append(ry_gate)
        rxx_gate = NoisyGate(Rxx(self.s*pi/2), [self.dataq, self.ancillaq], self.location + [1], error_model=self.error_model)
        noisy_gates.append(rxx_gate)
        if not self.cancel_data_rx:
            rx_gate = NoisyGate(Rx(-self.s*pi/2), [self.dataq], self.location + [2], error_model=self.error_model)
            noisy_gates.append(rx_gate)
            ry_gate = NoisyGate(Ry(self.v * pi / 2), [self.dataq], self.location + [3], error_model=self.error_model)
            noisy_gates.append(ry_gate)
        else:
            ry_gate = NoisyGate(Ry(self.v * pi / 2), [self.dataq], self.location + [2], error_model=self.error_model)
            noisy_gates.append(ry_gate)
        return noisy_gates


class CZTypeEntangling(EntanglingOperation):

    def run(self, dataq, ancillaq, error_list=None):
        rzz_gate = NoisyGate(Rzz(pi), self.location + [0], error_list)
        rzz_gate.apply(dataq, ancillaq)
        z_gate_a = NoisyGate(Z, self.location + [1], error_list)
        z_gate_a.apply(ancillaq)
        z_gate_d = NoisyGate(Z, self.location + [2], error_list)
        z_gate_d.apply(dataq)


class StabiliserTimestep():
    """During one StabiliserTimestep, 3 entangling operations are done in parallel."""

    def __init__(self, data=None, ancilla=None, location=None, qu_ind=None, entangling_type="X", cancel_data_rx=None, error_model=None):
        """
        attr:
            self.location (e.g. [0,0]) [Syndrome measurement location, stabilizer location]
            self.qu_ind (e.g. [[4,1], [8,6], [6,4]] ) The qubit indices needed for stabiliser_timesteps 1-8
                i.e. stabiliser_timestep_1 entangles data[4] & ancilla[1], then data[8] & ancilla[6] etc.
            self.ent_type ("X" or "Z") specifies whether this stabilizer timestep will do X-type or Z-type entangling
        """
        self.data = data
        self.ancilla = ancilla
        self.location = location
        self.error_model=error_model
        self.qu_ind = qu_ind
        self.ent_type = entangling_type
        self.cancel_data_rx = cancel_data_rx
        self.entangling_operations = self.build_entangling_operations()

    def build_entangling_operations(self):
        if self.ent_type == "X":
            ent_op_0 = XTypeEntanglingOp(location=self.location + [0], dataq=self.data[self.qu_ind[0][0]], ancillaq=self.ancilla[self.qu_ind[0][1]], cancel_data_rx=self.cancel_data_rx[0],s=1,error_model=self.error_model)
            ent_op_1 = XTypeEntanglingOp(location=self.location + [1], dataq=self.data[self.qu_ind[1][0]], ancillaq=self.ancilla[self.qu_ind[1][1]],
                                       cancel_data_rx=self.cancel_data_rx[1], s=1, error_model=self.error_model)
            ent_op_2 = XTypeEntanglingOp(location=self.location + [2], dataq=self.data[self.qu_ind[2][0]], ancillaq=self.ancilla[self.qu_ind[2][1]],
                                       cancel_data_rx=self.cancel_data_rx[2], s=1, error_model=self.error_model)
        elif self.ent_type == "Z":
            ent_op_0 = ZTypeEntanglingOp(location=self.location + [0], dataq=self.data[self.qu_ind[0][0]], ancillaq=self.ancilla[self.qu_ind[0][1]],
                                       cancel_data_rx=self.cancel_data_rx[0], s=1, v=1, error_model=self.error_model)
            ent_op_1 = ZTypeEntanglingOp(location=self.location + [1], dataq=self.data[self.qu_ind[1][0]], ancillaq=self.ancilla[self.qu_ind[1][1]],
                                       cancel_data_rx=self.cancel_data_rx[1], s=1, v=1, error_model=self.error_model)
            ent_op_2 = ZTypeEntanglingOp(location=self.location + [2], dataq=self.data[self.qu_ind[2][0]], ancillaq=self.ancilla[self.qu_ind[2][1]],
                                       cancel_data_rx=self.cancel_data_rx[2], s=1, v=1, error_model=self.error_model)
        else:
            raise Exception(
                "Entangling type not recognized, should be \"X\" or \"Z\". Please check stabilizer_timestep instantiation.")
        return [ent_op_0, ent_op_1, ent_op_2]

    def run(self, error_list=[]):
        for entangling_operation in self.entangling_operations:
            entangling_operation.run(error_list=error_list)


class StabiliserCycle():
    """During one StabiliserCycle, all the data qubits get entangled with all the appropriate ancillas."""
    def __init__(self, location=None, data=None, ancilla=None, circuit_type="MS", error_model=None, error_subset={}):
        '''
        attr:
            location = [0]
        '''
        self.data=data
        self.ancilla=ancilla
        self.location = location
        self.circuit_type = circuit_type
        self.error_model = error_model
        self.stabiliser_timesteps = self.build_stabiliser_timesteps(data, ancilla, self.circuit_type)
        self.all_gate_locations, self.error_locations, self.gate_locations = self.generate_gate_and_error_locations()
        self.errors_to_insert = Error.generate_error_list(error_subset, self.error_locations)

    def __str__(self):
        string = ""
        for timestep in self.stabiliser_timesteps:
            string+=str(timestep)
        return string

    def generate_gate_and_error_locations(self):
        gate_locations = { Rx : [], Ry : [], Rxx : [] }
        all_gate_locations = []
        error_locations = { XCtrlError : [], YCtrlError : [], DephasingError : []}
        for timestep in self.stabiliser_timesteps:
            for entangling_op in timestep.entangling_operations:
                for noisy_gate in entangling_op.noisy_gates:
                    all_gate_locations.append(noisy_gate.location)
                    gate_locations[type(noisy_gate.gate)].append(noisy_gate.location)
                    for error in noisy_gate.possible_errors:
                        error_locations[type(error)].append(noisy_gate.location)
        return all_gate_locations, error_locations, gate_locations

    def build_stabiliser_timesteps(self, data, ancilla, circuit_type):
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
        # i.e. stabiliser_timestep_0 entangles data[4] & ancilla[1], then data[8] & ancilla[6] etc.
        # The same for all circuit types
        qu_ind_0 = [[4, 1], [8, 6], [6, 4]]
        qu_ind_1 = [[1, 1], [5, 6], [3, 4]]
        qu_ind_2 = [[3, 1], [7, 6], [5, 3]]
        qu_ind_3 = [[0, 1], [4, 6], [2, 3]]
        qu_ind_4 = [[1, 2], [3, 5], [7, 7]]
        qu_ind_5 = [[2, 2], [4, 5], [8, 7]]
        qu_ind_6 = [[4, 2], [6, 5], [0, 0]]
        qu_ind_7 = [[5, 2], [7, 5], [1, 0]]

        if circuit_type == "MS":
            # Only relevant for MS circuit type

            cancel_data_rx = [[True, False, False],
                              [False, True, True],
                              [True, False, True],
                              [False, True, False],
                              [True, False, True],
                              [False, True, False],
                              [True, False, False],
                              [False, True, True]]

            stabiliser_timestep_0 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_0, qu_ind=qu_ind_0, entangling_type="X", cancel_data_rx=cancel_data_rx[0], error_model=self.error_model)
            stabiliser_timestep_1 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_1, qu_ind=qu_ind_1, entangling_type="X", cancel_data_rx=cancel_data_rx[1], error_model=self.error_model)
            stabiliser_timestep_2 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_2, qu_ind=qu_ind_2, entangling_type="X", cancel_data_rx=cancel_data_rx[2], error_model=self.error_model)
            stabiliser_timestep_3 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_3, qu_ind=qu_ind_3, entangling_type="X", cancel_data_rx=cancel_data_rx[3], error_model=self.error_model)
            stabiliser_timestep_4 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_4, qu_ind=qu_ind_4, entangling_type="Z", cancel_data_rx=cancel_data_rx[4], error_model=self.error_model)
            stabiliser_timestep_5 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_5, qu_ind=qu_ind_5, entangling_type="Z", cancel_data_rx=cancel_data_rx[5], error_model=self.error_model)
            stabiliser_timestep_6 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_6, qu_ind=qu_ind_6, entangling_type="Z", cancel_data_rx=cancel_data_rx[6], error_model=self.error_model)
            stabiliser_timestep_7 = StabiliserTimestep(data=self.data, ancilla=self.ancilla, location=loc_7, qu_ind=qu_ind_7, entangling_type="Z", cancel_data_rx=cancel_data_rx[7], error_model=self.error_model)
            stabiliser_timesteps = [stabiliser_timestep_0, stabiliser_timestep_1, stabiliser_timestep_2,
                                    stabiliser_timestep_3, stabiliser_timestep_4, stabiliser_timestep_5,
                                    stabiliser_timestep_6, stabiliser_timestep_7]
        if circuit_type == "CZ":
            # DO NOTHING
            print("Nothing done.")
        return stabiliser_timesteps

    def run(self, eng, reset=True):

        for timestep in self.stabiliser_timesteps:
            timestep.run(error_list=self.errors_to_insert)

        All(Measure) | self.ancilla
        eng.flush()
        syndrome_t = [int(q) for q in self.ancilla]
        if reset:
            for a in self.ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
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
        syndrome = stabiliser_cycle.run(self.data, self.ancilla, self.eng)
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
