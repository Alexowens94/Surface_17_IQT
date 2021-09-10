import random
import numpy as np
import pypartitions
import json
from math import sqrt, pi
from projectq.ops import All, CNOT, H, Measure, X, Z, Y, Rx, Ry, Rxx, Deallocate, Command, MatrixGate,\
    BasicGate, get_inverse

class RXXpi4Gate(BasicGate):
    """fixed angle (pi/4) 2q XX  gate class.
        Decomposing the CNOT to native gateset using standard Ry/Rx/Rxx gates resulted in errors"""

    @property
    def matrix(self):
        """Access to the matrix property of this gate."""
        return np.matrix([[1 / sqrt(2), 0, 0, -1j / sqrt(2)],
                    [0, 1 / sqrt(2), -1j / sqrt(2), 0],
                    [0, -1j / sqrt(2), 1 / sqrt(2), 0],
                    [-1j / sqrt(2), 0, 0, 1 / sqrt(2)]
                    ])

    def __str__(self):
        """Return a string representation of the object."""
        return "Rxx(pi/4)"
#: Shortcut (instance of) :class:`projectq.ops.RXXpi4Gate`
XXpi4 = RXXpi4Gate()
#: Inverse (and shortcut) of :class:`projectq.ops.RXXpi4Gate`
XXpi4dag = XXpi4dagger = get_inverse(XXpi4)
class RXpi2Gate(BasicGate):
    """fixed angle (pi/2) X rotation gate class.
        Decomposing the CNOT to native gateset using standard Ry/Rx/Rxx gates resulted in errors"""

    @property
    def matrix(self):
        """Access to the matrix property of this gate."""
        return np.matrix([[1 / sqrt(2), -1j / sqrt(2)], [-1j / sqrt(2), 1 / sqrt(2)]])

    def __str__(self):
        """Return a string representation of the object."""
        return "Rx(pi/2)"
#: Shortcut (instance of) :class:`projectq.ops.RXpi2Gate`
Xpi2 = RXpi2Gate()
#: Inverse (and shortcut) of :class:`projectq.ops.RXpi2Gate`
Xpi2dag = Xpi2dagger = get_inverse(Xpi2)
class RXminpi2Gate(BasicGate):
    """fixed angle (-pi/2) X rotation gate class.
        Decomposing the CNOT to native gateset using standard Ry/Rx/Rxx gates resulted in errors"""

    @property
    def matrix(self):
        """Access to the matrix property of this gate."""
        return np.matrix([[1 / sqrt(2), 1j / sqrt(2)], [1j / sqrt(2), 1 / sqrt(2)]])

    def __str__(self):
        """Return a string representation of the object."""
        return "Rx(-pi/2)"
#: Shortcut (instance of) :class:`projectq.ops.RYpi2Gate`
Xminpi2 = RXminpi2Gate()
class RYpi2Gate(BasicGate):
    """fixed angle (pi/2) Y rotation gate class.
        Decomposing the CNOT to native gateset using standard Ry/Rx/Rxx gates resulted in errors"""

    @property
    def matrix(self):
        """Access to the matrix property of this gate."""
        return np.matrix([[1 / sqrt(2), -1 / sqrt(2)], [1 / sqrt(2), 1 / sqrt(2)]])

    def __str__(self):
        """Return a string representation of the object."""
        return "Ry(pi/2)"
#: Shortcut (instance of) :class:`projectq.ops.RYpi2Gate`
Ypi2 = RYpi2Gate()
#: Inverse (and shortcut) of :class:`projectq.ops.RXXpi4Gate`
Ypi2dag = Ypi2dagger = get_inverse(Ypi2)
class RYminpi2Gate(BasicGate):
    """fixed angle (-pi/2) Y rotation gate class.
        Decomposing the CNOT to native gateset using standard Ry/Rx/Rxx gates resulted in errors"""

    @property
    def matrix(self):
        """Access to the matrix property of this gate."""
        return np.matrix([[1 / sqrt(2), 1 / sqrt(2)], [-1 / sqrt(2), 1 / sqrt(2)]])

    def __str__(self):
        """Return a string representation of the object."""
        return "Ry(-pi/2)"
#: Shortcut (instance of) :class:`projectq.ops.RYpi2Gate`
Yminpi2 = RYminpi2Gate()




def load_lookup_table(filename):
    with open(filename,'r') as infile:
        correction_table=json.load(infile)
    return correction_table

def insert_1q_error(qubit):
    rand = random.random()
    if rand<1/3:
        X|qubit
    if 1/3<rand<(2*1/3): #depolarizing noise model
        Y|qubit
    if (2*1/3)<rand<(3*1/3):
        Z|qubit

def insert_meas_error(qubit):
    X | qubit

def insert_2q_error(qubit1,qubit2):
    rand = random.random()
    if rand<1/15:
        #I|qubit1
        X|qubit2
    if 1/15<rand<(2*1/15):
        # I|qubit1
        Y|qubit2
    if (2*1/15)<rand<(3*1/15):
        # I|qubit1
        Z|qubit2
    if (3*1/15)<rand<(4*1/15):
        X|qubit1
        #I|qubit2
    if (4*1/15)<rand<(5*1/15):
        X|qubit1
        X|qubit2
    if (5*1/15)<rand<(6*1/15):
        X|qubit1
        Y | qubit2
    if (6*1/15)<rand<(7*1/15):
        X|qubit1
        Z | qubit2
    if (7*1/15)<rand<(8*1/15):
        Y|qubit1
        #I|qubit2
    if (8*1/15)<rand<(9*1/15):
        Y|qubit1
        X|qubit2
    if (9*1/15)<rand<(10*1/15):
        Y|qubit1
        Y|qubit2
    if (10*1/15)<rand<(11*1/15):
        Y|qubit1
        Z|qubit2
    if (11*1/15)<rand<(12*1/15):
        Z|qubit1
        # I|qubit2
    if (12*1/15)<rand<(13*1/15):
        Z|qubit1
        X|qubit2
    if (13*1/15)<rand<(14*1/15):
        Z|qubit1
        Y|qubit2
    if (14*1/15)<rand<(15*1/15):
        Z|qubit1
        Z|qubit2

def noisy_stabiliser_step_1(data,ancilla,num_1q_error,num_2q_error,num_gates_1q,num_gates_2q,display=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    # 4 single qubit gates, so 4 possible noise locations


    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    # change measurement basis to X for half the ancillae
    H | ancilla[1]
    if error_locations_1q[0]==1:
        if display:
            print('h1')
        insert_1q_error(ancilla[1])

    H | ancilla[3]
    if error_locations_1q[1]==1:
        if display:
            print('h2')
        insert_1q_error(ancilla[3])

    H | ancilla[4]
    if error_locations_1q[2]==1:
        if display:
            print('h3')
        insert_1q_error(ancilla[4])

    H | ancilla[6]
    if error_locations_1q[3]==1:
        if display:
            print('h4')
        insert_1q_error(ancilla[6])

    # entangle

    CNOT | (ancilla[1], data[1])
    if error_locations_2q[0]==1:
        if display:
            print('cn1')
        insert_2q_error(ancilla[1], data[1])

    CNOT | (ancilla[4], data[3])
    if error_locations_2q[1]==1:
        if display:
            print('cn2')
        insert_2q_error(ancilla[4], data[3])

    CNOT | (ancilla[6], data[5])
    if error_locations_2q[2]==1:
        if display:
            print('cn3')
        insert_2q_error(ancilla[6], data[5])


    CNOT | (data[2], ancilla[2])
    if error_locations_2q[3]==1:
        if display:
            print('cn4')
        insert_2q_error(data[2], ancilla[2])

    CNOT | (data[4], ancilla[5])
    if error_locations_2q[4]==1:
        if display:
            print('cn5')
        insert_2q_error(data[4], ancilla[5])

    CNOT | (data[8], ancilla[7])
    if error_locations_2q[5]==1:
        if display:
            print('cn6')
        insert_2q_error(data[8], ancilla[7])

    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
def noisy_stabiliser_step_2(data,ancilla,num_1q_error,num_2q_error,num_gates_1q,num_gates_2q,display=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    # entangle
    CNOT | (ancilla[1], data[4])
    if error_locations_2q[0]==1:
        if display:
            print('cn1')
        insert_2q_error(ancilla[1], data[4])

    CNOT | (ancilla[4], data[6])
    if error_locations_2q[1]==1:
        if display:
            print('cn2')
        insert_2q_error(ancilla[4], data[6])

    CNOT | (ancilla[6], data[8])
    if error_locations_2q[2]==1:
        if display:
            print('cn3')
        insert_2q_error(ancilla[6], data[8])


    CNOT | (data[1], ancilla[2])
    if error_locations_2q[3]==1:
        if display:
            print('cn4')
        insert_2q_error(data[1], ancilla[2])

    CNOT | (data[3], ancilla[5])
    if error_locations_2q[4]==1:
        if display:
            print('cn5')
        insert_2q_error(data[3], ancilla[5])

    CNOT | (data[7], ancilla[7])
    if error_locations_2q[5]==1:
        if display:
            print('cn6')
        insert_2q_error(data[7], ancilla[7])

    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
def noisy_stabiliser_step_3(data,ancilla,num_1q_error,num_2q_error,num_gates_1q,num_gates_2q,display=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    # entangle
    CNOT | (ancilla[1], data[0])
    if error_locations_2q[0]==1:
        if display:
            print('cn1')
        insert_2q_error(ancilla[1], data[0])

    CNOT | (ancilla[6], data[4])
    if error_locations_2q[1]==1:
        if display:
            print('cn2')
        insert_2q_error(ancilla[6], data[4])

    CNOT | (ancilla[3], data[2])
    if error_locations_2q[2]==1:
        if display:
            print('cn3')
        insert_2q_error(ancilla[3], data[2])


    CNOT | (data[1], ancilla[0])
    if error_locations_2q[3]==1:
        if display:
            print('cn4')
        insert_2q_error(data[1], ancilla[0])

    CNOT | (data[5], ancilla[2])
    if error_locations_2q[4]==1:
        if display:
            print('cn5')
        insert_2q_error(data[5], ancilla[2])

    CNOT | (data[7], ancilla[5])
    if error_locations_2q[5]==1:
        if display:
            print('cn6')
        insert_2q_error(data[7], ancilla[5])

    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
def noisy_stabiliser_step_4(data,ancilla,num_1q_error,num_2q_error,num_gates_1q,num_gates_2q,display=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    # 4 single qubit gates, so 4 possible noise locations
    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)
    # entangle

    CNOT | (ancilla[1], data[3])
    if error_locations_2q[0]==1:
        if display:
            print('cn1')
        insert_2q_error(ancilla[1], data[3])

    CNOT | (ancilla[6], data[7])
    if error_locations_2q[1]==1:
        if display:
            print('cn2')
        insert_2q_error(ancilla[6], data[7])

    CNOT | (ancilla[3], data[5])
    if error_locations_2q[2]==1:
        if display:
            print('cn3')
        insert_2q_error(ancilla[3], data[5])


    CNOT | (data[0], ancilla[0])
    if error_locations_2q[3]==1:
        if display:
            print('cn4')
        insert_2q_error(data[0], ancilla[0])

    CNOT | (data[4], ancilla[2])
    if error_locations_2q[4]==1:
        if display:
            print('cn5')
        insert_2q_error(data[4], ancilla[2])

    CNOT | (data[6], ancilla[5])
    if error_locations_2q[5]==1:
        if display:
            print('cn6')
        insert_2q_error(data[6], ancilla[5])
    # change measurement basis to X for half the ancillae
    H | ancilla[1]
    if error_locations_1q[0] == 1:
        if display:
            print('h1')
        insert_1q_error(ancilla[1])

    H | ancilla[3]
    if error_locations_1q[1] == 1:
        if display:
            print('h2')
        insert_1q_error(ancilla[3])

    H | ancilla[4]
    if error_locations_1q[2] == 1:
        if display:
            print('h3')
        insert_1q_error(ancilla[4])

    H | ancilla[6]
    if error_locations_1q[3] == 1:
        if display:
            print('h4')
        insert_1q_error(ancilla[6])

    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))

def native_noisy_hadamard(ancilla, error_locations_1q,display=False,draw=False):
    X | ancilla
    if error_locations_1q[0]==1:
        if display:
            print('X (hadamard) ancilla ')
        insert_1q_error(ancilla)
    if draw:
        Ry(pi/2)
    Ypi2dag | ancilla
    if error_locations_1q[1]==1:
        if display:
            print('0.5Y (hadamard) ancilla ')
        insert_1q_error(ancilla)

def native_noisy_ancilla_to_x_basis(ancilla,error_locations_1q,
                                    display=False,draw=False): # just pass the right error_locations
    # change measurement basis to X for half the ancillae
    native_noisy_hadamard(ancilla[1], error_locations_1q[0:2],draw=draw)
    native_noisy_hadamard(ancilla[3], error_locations_1q[2:4],draw=draw)
    native_noisy_hadamard(ancilla[4], error_locations_1q[4:6],draw=draw)
    native_noisy_hadamard(ancilla[6], error_locations_1q[6:8],draw=draw)

def native_noisy_CNOT(c,t,error_locations_1q=4*[0],error_locations_2q=0,v=1,s=1,
                      display=False,draw=False):
    """CNOT decomposition as per figure 1 of Dmitri Maslov 2017 New J. Phys. 19 023035"""
# Gate 1 (1q control)
    if v == 1:
        Ypi2 | c
    if v == -1:
        Ypi2dag | c
    if error_locations_1q[0] == 1:
        if display:
            print('0.5Yc1')
        insert_1q_error(c)
#Gate 2 (2q)
    if draw:
        Rxx(pi/4) | (c, t)
    else:
        if s == 1:
            XXpi4 | (c, t)
        if s == -1:
            XXpi4dag | (c, t)
    if error_locations_2q == 1:
        if display:
            print('XX1')
        insert_2q_error(c,t)
#Gate 3 (1q control)
    if s == 1:
        Xpi2dag | c
    if s == -1:
        Xpi2 | c
    if error_locations_1q[1] == 1:
        if display:
            print('0.5Xc1')
        insert_1q_error(c)
#Gate 4 (1q target)
    if v * s == 1:
        Xpi2dag | t
    if v * s == -1:
        Xpi2 | t
    if error_locations_1q[2] == 1:
        if display:
            print('0.5Xt1')
        insert_1q_error(t)
#Gate 5 (1q control)
    if v == 1:
        Ypi2dag | c
    if v == -1:
        Ypi2 | c
    if error_locations_1q[3] == 1:
        if display:
            print('0.5Yc1')
        insert_1q_error(c)

def nat_CNOT(c,t,v=1,s=1):
    """CNOT decomposition as per figure 1 of Dmitri Maslov 2017 New J. Phys. 19 023035"""

    if v==1:
        Ypi2 | c
    if v==-1:
        Ypi2dag | c

    if s==1:
        XXpi4 | (c,t)
    if s==-1:
        XXpi4dag | (c,t)

    if s==1:
        Xpi2dag | c
    if s==-1:
        Xpi2 | c

    if v*s == 1:
        Xpi2dag | t
    if v*s == -1:
        Xpi2 | t

    if v == 1:
        Ypi2dag | c
    if v == -1:
        Ypi2 | c
def native_noisy_entangle(data, d_ind, ancilla, a_ind, error_locations_1q,
                          error_locations_2q, display=False,draw=False):
    """

    :param data:
    :param d_ind: list of indices of data qubits involve in this entangling step - in order of appearance
    :param ancilla:
    :param a_ind: list of indices of ancilla qubits involved in this entangling step - in order of appearance
    :param error_locations_1q:
    :param error_locations_2q:
    :param display:
    :return:
    """

    native_noisy_CNOT(ancilla[a_ind[0]],data[d_ind[0]],error_locations_1q[0:4],
                      error_locations_2q[0], display= display)
    native_noisy_CNOT(ancilla[a_ind[1]],data[d_ind[1]],error_locations_1q[4:8],
                      error_locations_2q[1], display= display)
    native_noisy_CNOT(ancilla[a_ind[2]],data[d_ind[2]],error_locations_1q[8:12],
                      error_locations_2q[2], display= display)

    native_noisy_CNOT(data[d_ind[3]],ancilla[a_ind[3]], error_locations_1q[12:16],
                      error_locations_2q[3], display= display)
    native_noisy_CNOT(data[d_ind[4]],ancilla[a_ind[4]], error_locations_1q[16:20],
                      error_locations_2q[4], display= display)
    native_noisy_CNOT(data[d_ind[5]],ancilla[a_ind[5]], error_locations_1q[20:24],
                      error_locations_2q[5], display= display)


def native_noisy_stabiliser_step_1(data,ancilla,num_1q_error,num_2q_error,
                                   num_gates_1q,num_gates_2q,display=False,draw=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    error_locations_1q = generate_noise_locations(num_1q_error, num_gates_1q,[1] * num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    # switch half the ancillas to the x basis (requires 8 single qubit gates)
    native_noisy_ancilla_to_x_basis(ancilla,error_locations_1q[0:8],draw=draw)
    # entangle 6x 2 qubit, 24 (4*6) 1 qubit gates
    data_inds = [1, 3, 5, 2, 4, 8]
    ancilla_inds = [1, 4, 6, 2, 5, 7]
    native_noisy_entangle(data, data_inds, ancilla, ancilla_inds,
                         error_locations_1q, error_locations_2q, display,draw=draw)
    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
def native_noisy_stabiliser_step_2(data,ancilla,num_1q_error,num_2q_error,
                                   num_gates_1q,num_gates_2q,display=False,draw=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    error_locations_1q = generate_noise_locations(num_1q_error, num_gates_1q, [1] * num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
    # entangle 6x 2 qubit, 24 (4*6) 1 qubit gates
    ancilla_inds = [1, 4, 6, 2, 5, 7]
    data_inds = [4, 6, 8, 1, 3, 7]
    native_noisy_entangle(data, data_inds, ancilla, ancilla_inds,
                          error_locations_1q, error_locations_2q, display)


def native_noisy_stabiliser_step_3(data,ancilla,num_1q_error,num_2q_error,
                                   num_gates_1q,num_gates_2q,display=False,draw=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''

    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)

    # entangle 6x 2 qubit, 24 (4*6) 1 qubit gates
    data_inds = [0, 4, 2, 1, 5, 7]
    ancilla_inds = [1, 6, 3, 0, 2, 5]
    native_noisy_entangle(data, data_inds, ancilla, ancilla_inds,
                          error_locations_1q, error_locations_2q, display)
    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))
def native_noisy_stabiliser_step_4(data,ancilla,num_1q_error,num_2q_error,num_gates_1q,
                                   num_gates_2q,display=False,draw=False):
    '''

    :param data: list of data qubits
    :param ancilla: list of ancilla qubits
    :param num_1q_error: number of 1 qubit errors to introduce in this step, for efficient logical error estimation
    see Importance Sampling https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341
    :param num_2q_error: number of 2 qubit errors to introduce in this step
    :return: -
    '''
    error_locations_1q=generate_noise_locations(num_1q_error,num_gates_1q,[1]*num_gates_1q)
    error_locations_2q = generate_noise_locations(num_2q_error, num_gates_2q, [1] * num_gates_2q)
    # entangle 6x 2 qubit, 24 (4*6) 1 qubit gates
    data_inds = [3, 7, 5, 0, 4, 6]
    ancilla_inds = [1, 6, 3, 0, 2, 5]
    native_noisy_entangle(data, data_inds, ancilla, ancilla_inds,
                          error_locations_1q, error_locations_2q, display)
    # switch half the ancillas to the x basis (requires 8 single qubit gates)
    native_noisy_ancilla_to_x_basis(ancilla, error_locations_1q[24:])
    if display:
        print('step1 1q error locations {}'.format(error_locations_1q))
        print('step1 2q error locations {}'.format(error_locations_2q))

def stabiliser_step_1(data,ancilla):
    #change measurement basis to X for half the ancillae
    H | ancilla[1]
    H | ancilla[3]
    H | ancilla[4]
    H | ancilla[6]
    #entangle
    CNOT | (ancilla[1],data[1])
    CNOT | (ancilla[4],data[3])
    CNOT | (ancilla[6],data[5])


    CNOT | (data[2],ancilla[2])
    CNOT | (data[4],ancilla[5])
    CNOT | (data[8],ancilla[7])

def stabiliser_step_2(data,ancilla):
    CNOT | (ancilla[1],data[4])
    CNOT | (ancilla[4],data[6])
    CNOT | (ancilla[6],data[8])

    CNOT | (data[1],ancilla[2])
    CNOT | (data[3],ancilla[5])
    CNOT | (data[7],ancilla[7])
def stabiliser_step_3(data,ancilla):
    CNOT | (ancilla[1],data[0])
    CNOT | (ancilla[6],data[4])
    CNOT | (ancilla[3],data[2])

    CNOT | (data[1],ancilla[0])
    CNOT | (data[5],ancilla[2])
    CNOT | (data[7],ancilla[5])
def stabiliser_step_4(data, ancilla):
    CNOT | (ancilla[1],data[3])
    CNOT | (ancilla[6],data[7])
    CNOT | (ancilla[3],data[5])

    CNOT | (data[0],ancilla[0])
    CNOT | (data[4],ancilla[2])
    CNOT | (data[6],ancilla[5])
    # prepend Z basis measurement of half the ancillae with H to measure in the x basis
    H | ancilla[1]
    H | ancilla[3]
    H | ancilla[4]
    H | ancilla[6]
def stabilize(eng,ancilla,data,errors_per_func_1q,errors_per_func_2q,
              gates_per_func_1q,gates_per_func_2q,
              cycle,number_stab_steps_1cycle=4,gateset='native',draw=False):
    """
    Stabilizer cycle for the 17 surface code
    index repartition:

                Z ancilla0
          data0-----------data1-----------data2
            |               |               |
            |    X ancilla1 |    Z ancilla2 |  X ancilla3
            |               |               |
          data3-----------data4-----------data5
            |               |               |
X ancilla4  |   Z ancilla5  |    X ancilla6 |
            |               |               |
          data6-----------data7------------data8

                                Z ancilla7

    """
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    ind0=cycle*number_stab_steps_1cycle
    if gateset == 'circuit':
        noisy_stabiliser_step_1(data,ancilla,errors_per_func_1q[ind0],
                                errors_per_func_2q[ind0],gates_per_func_1q[ind0],gates_per_func_2q[ind0])
        noisy_stabiliser_step_2(data,ancilla,errors_per_func_1q[ind0+1],
                                errors_per_func_2q[ind0+1],gates_per_func_1q[ind0+1],gates_per_func_2q[ind0+1])
        noisy_stabiliser_step_3(data,ancilla,errors_per_func_1q[ind0+2],
                                errors_per_func_2q[ind0+2],gates_per_func_1q[ind0+2],gates_per_func_2q[ind0+2])
        noisy_stabiliser_step_4(data,ancilla,errors_per_func_1q[ind0+3],
                                errors_per_func_2q[ind0+3],gates_per_func_1q[ind0+3],gates_per_func_2q[ind0+3])
    if gateset == 'native':

        native_noisy_stabiliser_step_1(data, ancilla, errors_per_func_1q[ind0],
                                errors_per_func_2q[ind0], gates_per_func_1q[ind0],
                                gates_per_func_2q[ind0],draw=draw)
        native_noisy_stabiliser_step_2(data, ancilla, errors_per_func_1q[ind0 + 1],
                                errors_per_func_2q[ind0 + 1], gates_per_func_1q[ind0 + 1],
                                gates_per_func_2q[ind0 + 1],draw=draw)
        native_noisy_stabiliser_step_3(data, ancilla, errors_per_func_1q[ind0 + 2],
                                errors_per_func_2q[ind0 + 2], gates_per_func_1q[ind0 + 2],
                                gates_per_func_2q[ind0 + 2],draw=draw)
        native_noisy_stabiliser_step_4(data, ancilla, errors_per_func_1q[ind0 + 3],
                                errors_per_func_2q[ind0 + 3], gates_per_func_1q[ind0 + 3],
                                gates_per_func_2q[ind0 + 3],draw=draw)

    All(Measure) | ancilla
    eng.flush()
    syndrome_t = [int(q) for q in ancilla]
    if draw:
        return syndrome_t
    for a in ancilla: #reset the ancillas to 0 at end of stab round (allow for repeat rounds)
        if int(a) ==1:
            X|a


    return syndrome_t


def stabilize_noisy_ancilla_meas(eng, ancilla, data, error_locations_a_meas,
                                 errors_per_func_1q, errors_per_func_2q,
                                 gates_per_func_1q, gates_per_func_2q,
                                 cycle, number_stab_steps_1cycle=4, gateset='native'):
    """
    Stabilizer cycle for the 17 surface code
    index repartition:

                Z ancilla0
          data0-----------data1-----------data2
            |               |               |
            |    X ancilla1 |    Z ancilla2 |  X ancilla3
            |               |               |
          data3-----------data4-----------data5
            |               |               |
X ancilla4  |   Z ancilla5  |    X ancilla6 |
            |               |               |
          data6-----------data7------------data8

                                Z ancilla7

    """
    if len(data) != 9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    ind0 = cycle * number_stab_steps_1cycle
    if gateset == 'circuit':
        noisy_stabiliser_step_1(data, ancilla, errors_per_func_1q[ind0],
                                errors_per_func_2q[ind0], gates_per_func_1q[ind0], gates_per_func_2q[ind0])
        noisy_stabiliser_step_2(data, ancilla, errors_per_func_1q[ind0 + 1],
                                errors_per_func_2q[ind0 + 1], gates_per_func_1q[ind0 + 1], gates_per_func_2q[ind0 + 1])
        noisy_stabiliser_step_3(data, ancilla, errors_per_func_1q[ind0 + 2],
                                errors_per_func_2q[ind0 + 2], gates_per_func_1q[ind0 + 2], gates_per_func_2q[ind0 + 2])
        noisy_stabiliser_step_4(data, ancilla, errors_per_func_1q[ind0 + 3],
                                errors_per_func_2q[ind0 + 3], gates_per_func_1q[ind0 + 3], gates_per_func_2q[ind0 + 3])
    if gateset == 'native':
        native_noisy_stabiliser_step_1(data, ancilla, errors_per_func_1q[ind0],
                                       errors_per_func_2q[ind0], gates_per_func_1q[ind0], gates_per_func_2q[ind0])
        native_noisy_stabiliser_step_2(data, ancilla, errors_per_func_1q[ind0 + 1],
                                       errors_per_func_2q[ind0 + 1], gates_per_func_1q[ind0 + 1],
                                       gates_per_func_2q[ind0 + 1])
        native_noisy_stabiliser_step_3(data, ancilla, errors_per_func_1q[ind0 + 2],
                                       errors_per_func_2q[ind0 + 2], gates_per_func_1q[ind0 + 2],
                                       gates_per_func_2q[ind0 + 2])
        native_noisy_stabiliser_step_4(data, ancilla, errors_per_func_1q[ind0 + 3],
                                       errors_per_func_2q[ind0 + 3], gates_per_func_1q[ind0 + 3],
                                       gates_per_func_2q[ind0 + 3])
    ind_m = cycle*len(ancilla)
    for i,a in enumerate(ancilla):
        if error_locations_a_meas[ind_m+i] ==1:
            insert_meas_error(a)

    if gateset == 'compiled':
        compiled_noisy_stabiliser_step_1(data, ancilla, errors_per_func_1q[ind0],
                                       errors_per_func_2q[ind0], gates_per_func_1q[ind0], gates_per_func_2q[ind0])
        compiled_noisy_stabiliser_step_2(data, ancilla, errors_per_func_1q[ind0 + 1],
                                       errors_per_func_2q[ind0 + 1], gates_per_func_1q[ind0 + 1],
                                       gates_per_func_2q[ind0 + 1])
        compiled_noisy_stabiliser_step_3(data, ancilla, errors_per_func_1q[ind0 + 2],
                                       errors_per_func_2q[ind0 + 2], gates_per_func_1q[ind0 + 2],
                                       gates_per_func_2q[ind0 + 2])
        compiled_noisy_stabiliser_step_4(data, ancilla, errors_per_func_1q[ind0 + 3],
                                       errors_per_func_2q[ind0 + 3], gates_per_func_1q[ind0 + 3],
                                       gates_per_func_2q[ind0 + 3])
    ind_m = cycle*len(ancilla)
    for i,a in enumerate(ancilla):
        if error_locations_a_meas[ind_m+i] ==1:
            insert_meas_error(a)

    All(Measure) | ancilla

    eng.flush()
    syndrome_t = [int(q) for q in ancilla]
    for a in ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
        if int(a) == 1:
            X | a

    return syndrome_t

def stabilize_no_noise(eng,ancilla,data,reset=True):
    """
    Stabilizer cycle for the 17 surface code
    index repartition:

                Z ancilla0
          data0-----------data1-----------data2
            |               |               |
            |    X ancilla1 |    Z ancilla2 |  X ancilla3
            |               |               |
          data3-----------data4-----------data5
            |               |               |
X ancilla4  |   Z ancilla5  |    X ancilla6 |
            |               |               |
          data6-----------data7------------data8

                                Z ancilla7

    """
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    stabiliser_step_1(data,ancilla)
    stabiliser_step_2(data,ancilla)
    stabiliser_step_3(data,ancilla)
    stabiliser_step_4(data,ancilla)
    All(Measure) | ancilla
    eng.flush()
    syndrome_t = [int(q) for q in ancilla]
    # print('single cycle syndrome {}'.format(syndrome_t))
    if reset:
        for a in ancilla: #reset the ancillas to 0 at end of stab round (allow for repeat rounds)
            if int(a) ==1:
                X|a


    return syndrome_t

def generate_noise_locations(num_errors,num_locations,max_errors_per_location,display=False):
    '''for the importance sampling simulation method (explained in the Appendix of
    https://journals-aps-org.ezproxy.sussex.ac.uk/pra/pdf/10.1103/PhysRevA.96.032341)
    we need to insert N errors in random places, and assess the logical error rate given there are
     N errors.
      we achieve random placement of N errors by randomly partitioning N into
       the number of possible error locations'''
    if sum(max_errors_per_location) == num_errors:
        return max_errors_per_location
    if sum(max_errors_per_location) < num_errors:
        raise ValueError('{} errors cant fit into bins of max sizes {} '.format(num_errors,max_errors_per_location))
    bad_part = True
    while bad_part == True: #find partitions where no bin has more than allowed num errors,
        #e.g. for gates if you want x unique errors, max per location is 1
        flag=False
        # first generate a partition (allow zeros)
        part = pypartitions.rand_partitions(num_errors, num_locations, 1, zeros=True)[0]
        #print(part)
        #shuffle the partition in place, otherwise always have trailing 0s
        random.shuffle(part)
        for l1, l2 in zip(part, max_errors_per_location):
            if l1 > l2:
                flag=True
        bad_part = flag
    if display==True:
        print('founds suitable error partition {}'.format(part))


    return part

def decode(quiescent,syndrome,num_cycles,table,data):
    '''determine required corrections from measured syndromes, first process
     multiple syndrome measurments by the algorithm of Tomita,Svore in
    https://arxiv.org/pdf/1404.3747.pdf
    then use a lookup table as per https://iopscience.iop.org/article/10.1088/1367-2630/aab341

    for the sake of importance sampling we want to keep track of how many instances would stop
    after 1 error correction round
    '''
    # print('quiescent {}'.format(quiescent))
    # print('syndrome volume {}'.format(syndrome))
    if num_cycles == 1:
        ft_syndrome = (quiescent - syndrome[0]) % 2
        a_all_0 = np.all((ft_syndrome == 0))
        if a_all_0:
            a0_increment = 1
        else:
            a0_increment = 0
    else:
        ft_syndrome, a0_increment = process_syndromes(quiescent, syndrome, num_cycles)
    error_vec = lookup(ft_syndrome,table)
    apply_correction(error_vec,data)
    return a0_increment

def process_syndromes(quiescent,syndrome,num_cycles,display=False):
    flips = np.zeros((num_cycles, len(syndrome[0])))
    for i, s in enumerate(syndrome):
        if i == 0:
            flips[i, :] =  (quiescent - syndrome[i]) % 2
        else:
            flips[i, :] =  (syndrome[i - 1] - syndrome[i]) % 2
    # for i in range(num_cycles):
    #     print('flips {}'.format(flips[i]))
    ft_syndrome = np.zeros(len(syndrome[0]),dtype=np.int8)
    a_all_0 = np.all((flips[0,:] == 0))
    if display:
        print('flips synd1: {}'.format(flips[0, :]))
        print('flips synd2: {}'.format(flips[1, :]))
    if a_all_0:
        a0_increment = 1
    else:
        for i, a in enumerate(flips[0, :]):  # slight concern might be ignoring rule 3 of synd processing as per tomita
            if a == 1:
                ft_syndrome[i] = (a + flips[1, i]) % 2
        a0_increment = 0

    return ft_syndrome,a0_increment
def lookup(syndrome, table, display = False):

    key = str(syndrome).strip('[,]')
    error_vec = table[key][0]
    if display:
        print('ft syndrome: {}'.format(syndrome))
        print('error vector: {}'.format(error_vec))
    return error_vec

def apply_correction(error_vec,data):
    for i in range(9):
        if error_vec[i] == 1:
            X | data[i]
        if error_vec[i+9] == 1:
            Z | data[i]
    return

def get_results(cycles, g1q, e1q, g2q, e2q, runs, a0, a0_log_error, a_not0_log_error):
    res = {
        "runs": runs,
        "cycles": cycles,
        "1q_gates": g1q,
        "1q_errors": e1q,
        "2q_gates": g2q,
        "2q_errors": e2q,
        "syndrome_a_all_0": a0,
        "L_error_a_all_0": a0_log_error,
        "L_error_a_not_0": a_not0_log_error
          }

    if cycles == 1:
        if a0 != 0:
            res["L_e_rate_given_a0"] = a0_log_error/a0
        else:
            print('no runs would have stopped after cycle 1')
    if cycles > 1:
        res["a0_rate"] = a0/runs
        if (runs-a0) != 0:
            res["L_e_rate_given_a_not_0"] = a_not0_log_error/(runs-a0)
        else:
            print('all runs would have stopped after cycle 1')
    return res