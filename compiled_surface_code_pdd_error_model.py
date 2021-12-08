import random
import numpy as np
import json
from math import sqrt, pi
from projectq.ops import All, Measure, X, Y, Z, Rx, Ry, Rxx
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

def insert_2q_ctrl_error(qubit1, qubit2):
    X | qubit1
    X | qubit2

def insert_2q_motional_error(qubit1, qubit2):
    X | qubit1
    X | qubit2

def insert_1q_x_ctrl_error(qubit):
    X | qubit

def insert_1q_y_ctrl_error(qubit):
    Y | qubit

def insert_1q_dephasing_error(qubit):
    Z | qubit

def x_type_entangling(dataq,ancillaq,pxctrl=0,pxxctrl=0,pmot=0,pd=0,cancel_data_rx=False,s=1):
    """
    CNOT gate, as it appears in x-stabilizer compiled to native operations. All ancilla single qubits ops
    canceled 'by hand'.
    cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
    on the same data qubit (depedent on choices of s and v)
    """
    Rxx(s*pi/2) | (dataq,ancillaq)
    if random.random()<pxxctrl:
        insert_2q_ctrl_error(dataq, ancillaq)
    if random.random()<pmot:
        insert_2q_motional_error(dataq, ancillaq)
    if random.random()<pd:
        insert_1q_dephasing_error(dataq)
    if random.random()<pd:
        insert_1q_dephasing_error(ancillaq)
    if not cancel_data_rx:
        Rx(-s*pi/2) | dataq
        if random.random()<pxctrl:
            insert_1q_x_ctrl_error(dataq)
        if random.random()<pd:
            insert_1q_dephasing_error(dataq)

def z_type_entangling(dataq,ancillaq,pxctrl=0,pyctrl=0,pxxctrl=0,pmot=0,pd=0,cancel_data_rx=False,s=1,v=1):
    """
    CNOT gate, as it appears in z-stabilizer compiled to native operations. All ancilla single qubits ops
    canceled 'by hand'.
    cancel_data_rx: If True, 'by hand' remove Rx gate to cancel with the Rx in other entangling steps
    on the same data qubit (depedent on choices of s and v)
    """
    Ry(v*pi/2) | dataq
    if random.random()<pyctrl:
        insert_1q_y_ctrl_error(dataq)
    if random.random()<pd:
        insert_1q_dephasing_error(dataq)
    Rxx(s * pi / 2) | (dataq, ancillaq)
    if random.random() < pxxctrl:
        insert_2q_ctrl_error(dataq, ancillaq)
    if random.random() < pmot:
        insert_2q_motional_error(dataq, ancillaq)
    if random.random() < pd:
        insert_1q_dephasing_error(dataq)
    if random.random() < pd:
        insert_1q_dephasing_error(ancillaq)
    if not cancel_data_rx:
        Rx(-s * pi / 2) | dataq
        if random.random()<pxctrl:
            insert_1q_x_ctrl_error(dataq)
        if random.random()<pd:
            insert_1q_dephasing_error(dataq)
    Ry(v * pi / 2) | dataq
    if random.random() < pyctrl:
        insert_1q_y_ctrl_error(dataq)
    if random.random() < pd:
        insert_1q_dephasing_error(dataq)

def stabiliser_timestep_1(data,ancilla,pxxctrl,pxctrl,pmot,pd):
    """
    s, v and timestep ordering chosen to minimize single qubit gates
    """
    x_type_entangling(data[4], ancilla[1], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=1)
    x_type_entangling(data[8], ancilla[6], pxxctrl,pxctrl,pmot,pd,s=1)
    x_type_entangling(data[6], ancilla[4], pxxctrl,pxctrl,pmot,pd,s=1)

def stabiliser_timestep_2(data,ancilla,pxxctrl,pxctrl,pmot,pd):
    x_type_entangling(data[1], ancilla[1], pxxctrl,pxctrl,pmot,pd,s=-1)
    x_type_entangling(data[5], ancilla[6], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=-1)
    x_type_entangling(data[3], ancilla[4], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=-1)

def stabiliser_timestep_3(data,ancilla,pxxctrl,pxctrl,pmot,pd):
    x_type_entangling(data[3], ancilla[1], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=1)
    x_type_entangling(data[7], ancilla[6], pxxctrl,pxctrl,pmot,pd,s=1)
    x_type_entangling(data[5], ancilla[3], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=1)

def stabiliser_timestep_4(data,ancilla,pxxctrl,pxctrl,pmot,pd):
    x_type_entangling(data[0], ancilla[1], pxxctrl,pxctrl,pmot,pd,s=-1)
    x_type_entangling(data[4], ancilla[6], pxxctrl,pxctrl,pmot,pd,cancel_data_rx=True,s=-1)
    x_type_entangling(data[2], ancilla[3], pxxctrl,pxctrl,pmot,pd,s=-1)

def stabiliser_timestep_5(data,ancilla,pxxctrl,pxctrl,pyctrl,pmot,pd):
    z_type_entangling(data[1], ancilla[2], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=1)
    z_type_entangling(data[3], ancilla[5], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=1)
    z_type_entangling(data[7], ancilla[7], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=1)

def stabiliser_timestep_6(data,ancilla,pxxctrl,pxctrl,pyctrl,pmot,pd):
    z_type_entangling(data[2], ancilla[2], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=-1)
    z_type_entangling(data[4], ancilla[5], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=-1)
    z_type_entangling(data[8], ancilla[7], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=-1)

def stabiliser_timestep_7(data,ancilla,pxxctrl,pxctrl,pyctrl,pmot,pd):
    z_type_entangling(data[4], ancilla[2], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=1)
    z_type_entangling(data[6], ancilla[5], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=1)
    z_type_entangling(data[0], ancilla[0], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=1)

def stabiliser_timestep_8(data,ancilla,pxxctrl,pxctrl,pyctrl,pmot,pd):
    z_type_entangling(data[5], ancilla[2], pxxctrl,pxctrl,pyctrl,pmot,pd,v=1,s=-1)
    z_type_entangling(data[7], ancilla[5], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=-1)
    z_type_entangling(data[1], ancilla[0], pxxctrl,pxctrl,pyctrl,pmot,pd,cancel_data_rx=True,v=1,s=-1)

def stabilizer_cycle(data, ancilla, eng, reset=True, pxctrl=0,pxxctrl=0,pyctrl=0,pmot=0,pd=0):
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    stabiliser_timestep_1(data,ancilla,pxctrl,pxxctrl,pmot,pd)
    stabiliser_timestep_2(data,ancilla,pxctrl,pxxctrl,pmot,pd)
    stabiliser_timestep_3(data,ancilla,pxctrl,pxxctrl,pmot,pd)
    stabiliser_timestep_4(data,ancilla,pxctrl,pxxctrl,pmot,pd)
    stabiliser_timestep_5(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_6(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_7(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_8(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)


    All(Measure) | ancilla
    eng.flush()
    syndrome_t = [int(q) for q in ancilla]
    if reset:
        for a in ancilla: #reset the ancillas to 0 at end of stab round (allow for repeat rounds)
            if int(a) == 1:
                X | a

    return syndrome_t

def get_results_log_e(rounds, runs, incorrect_count, time):
    res = {
        "runs": runs,
        "rounds": rounds,
        "incorrect_count": incorrect_count,
        "time_taken": time,
          }
    return res

def get_results_log_e_run_til_fail(rounds, time, p1q, p2q):
    res = {
        "rounds_til_fail": rounds,
        "p1q_error": p1q,
        "p2q_error": p2q,
        "time_taken": time,
          }
    return res

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

def load_lookup_table(filename):
    with open(filename,'r') as infile:
        correction_table=json.load(infile)
    return correction_table