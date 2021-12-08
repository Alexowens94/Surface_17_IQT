import random
import numpy as np
import json
from math import sqrt, pi
from projectq.ops import All, Measure, X, Y, Z, Rx, Ry, Rxx

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

def x_type_entangling(dataq, ancillaq, leaked_reg, d_leak_ind, a_leak_ind, p_leak=0, pxctrl=0, pxxctrl=0, pmot=0, cancel_data_rx=False, s=1):
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s*pi/2) | (dataq, ancillaq)
        if random.random() < p_leak:
            leaked_reg[d_leak_ind] = 1 # Insert a leakage error
        if random.random() < p_leak:
            leaked_reg[a_leak_ind] = 1
        if random.random() < pxxctrl:
            insert_2q_ctrl_error(dataq, ancillaq) # Insert control error
        if random.random() < pmot:
            insert_2q_motional_error(dataq, ancillaq) # Insert motional error
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    if not cancel_data_rx:
        Rx(-s*pi/2) | dataq
        if random.random() < pxctrl:
            insert_1q_x_ctrl_error(dataq)

def z_type_entangling(dataq, ancillaq, leaked_reg, d_leak_ind, a_leak_ind, p_leak=0, pxctrl=0, pyctrl=0, pxxctrl=0, pmot=0, cancel_data_rx=False, s=1, v=1):
    Ry(v*pi/2) | dataq
    if random.random()<pyctrl:
        insert_1q_y_ctrl_error(dataq)
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s * pi / 2) | (dataq, ancillaq)
        if random.random() < p_leak:
            leaked_reg[d_leak_ind] = 1
        if random.random() < p_leak:  # currently roll for both qubits to leak independent of each other
            leaked_reg[a_leak_ind] = 1
        if random.random() < pxxctrl:
            insert_2q_ctrl_error(dataq, ancillaq)
        if random.random() < pmot:
            insert_2q_motional_error(dataq, ancillaq)
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    if not cancel_data_rx:
        Rx(-s * pi / 2) | dataq
        if random.random() < pxctrl:
            insert_1q_x_ctrl_error(dataq)
    Ry(-v * pi / 2) | dataq
    if random.random()<pyctrl:
        insert_1q_y_ctrl_error(dataq)

def stabiliser_timestep_1(data, ancilla, leaked_reg, p_leak, pxctrl, pxxctrl, pmot):
    x_type_entangling(data[4], ancilla[1], leaked_reg, d_leak_ind=4, a_leak_ind=1+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=1)
    x_type_entangling(data[8], ancilla[6], leaked_reg, d_leak_ind=8, a_leak_ind=6+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=1)
    x_type_entangling(data[6], ancilla[4], leaked_reg, d_leak_ind=6, a_leak_ind=4+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=1)

def stabiliser_timestep_2(data, ancilla, leaked_reg, p_leak, pxctrl, pxxctrl, pmot):
    x_type_entangling(data[1], ancilla[1], leaked_reg, d_leak_ind=1, a_leak_ind=1+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=-1)
    x_type_entangling(data[5], ancilla[6], leaked_reg, d_leak_ind=5, a_leak_ind=6+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=-1)
    x_type_entangling(data[3], ancilla[4], leaked_reg, d_leak_ind=3, a_leak_ind=4+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=-1)

def stabiliser_timestep_3(data, ancilla, leaked_reg, p_leak, pxctrl, pxxctrl, pmot):
    x_type_entangling(data[3], ancilla[1], leaked_reg, d_leak_ind=3, a_leak_ind=1+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=1)
    x_type_entangling(data[7], ancilla[6], leaked_reg, d_leak_ind=7, a_leak_ind=6+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=1)
    x_type_entangling(data[5], ancilla[3], leaked_reg, d_leak_ind=5, a_leak_ind=3+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=1)

def stabiliser_timestep_4(data, ancilla, leaked_reg, p_leak, pxctrl, pxxctrl, pmot):
    x_type_entangling(data[0], ancilla[1], leaked_reg, d_leak_ind=0, a_leak_ind=1+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=-1)
    x_type_entangling(data[4], ancilla[6], leaked_reg, d_leak_ind=4, a_leak_ind=6+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, s=-1)
    x_type_entangling(data[2], ancilla[3], leaked_reg, d_leak_ind=2, a_leak_ind=3+9, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, s=-1)

def stabiliser_timestep_5(data, ancilla, leaked_reg, p_leak, pxctrl, pyctrl, pxxctrl, pmot):
    z_type_entangling(data[1], ancilla[2], leaked_reg, d_leak_ind=1, a_leak_ind=2+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=1)
    z_type_entangling(data[3], ancilla[5], leaked_reg, d_leak_ind=1, a_leak_ind=2+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=1)
    z_type_entangling(data[7], ancilla[7], leaked_reg, d_leak_ind=7, a_leak_ind=7+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=1)

def stabiliser_timestep_6(data, ancilla, leaked_reg, p_leak, pxctrl, pyctrl, pxxctrl, pmot):
    z_type_entangling(data[2], ancilla[2], leaked_reg, d_leak_ind=2, a_leak_ind=2+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=-1)
    z_type_entangling(data[4], ancilla[5], leaked_reg, d_leak_ind=4, a_leak_ind=5+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=-1)
    z_type_entangling(data[8], ancilla[7], leaked_reg, d_leak_ind=8, a_leak_ind=7+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=-1)

def stabiliser_timestep_7(data, ancilla, leaked_reg, p_leak, pxctrl, pyctrl, pxxctrl, pmot):
    z_type_entangling(data[4], ancilla[2], leaked_reg, d_leak_ind=4, a_leak_ind=2+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=1)
    z_type_entangling(data[6], ancilla[5], leaked_reg, d_leak_ind=6, a_leak_ind=5+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=1)
    z_type_entangling(data[0], ancilla[0], leaked_reg, d_leak_ind=0, a_leak_ind=0+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=1)

def stabiliser_timestep_8(data, ancilla, leaked_reg, p_leak, pxctrl, pyctrl, pxxctrl, pmot):
    z_type_entangling(data[5], ancilla[2], leaked_reg, d_leak_ind=5, a_leak_ind=2+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, v=1, s=-1)
    z_type_entangling(data[7], ancilla[5], leaked_reg, d_leak_ind=7, a_leak_ind=5+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=-1)
    z_type_entangling(data[1], ancilla[0], leaked_reg, d_leak_ind=1, a_leak_ind=0+9, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot, cancel_data_rx=True, v=1, s=-1)

def stabilizer_cycle(data, ancilla, leaked_reg, eng, reset=True, p_leak=0, pxctrl=0, pyctrl=0, pxxctrl=0, pmot=0):
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    stabiliser_timestep_1(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_2(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_3(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_4(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_5(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_6(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_7(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot)
    stabiliser_timestep_8(data, ancilla, leaked_reg, p_leak=p_leak, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl, pmot=pmot)

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
    with open(filename, 'r') as infile:
        correction_table=json.load(infile)
    return correction_table