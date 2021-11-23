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

def x_type_entangling(dataq,ancillaq,leaked_reg,d_leak_ind,a_leak_ind,s=1,p_leak=0,pxctrl=0,pxxctrl=0,pmot=0):
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s*pi/2) | (dataq, ancillaq)
        if random.random() < p_leak:
            leaked_reg[d_leak_ind] = 1 # Insert a leakage error
        if random.random() < p_leak:
            leaked_reg[a_leak_ind] = 1
        if random.random() < pxxctrl:
            insert_2q_ctrl_error(dataq, ancillaq) # Insert control error
        if random.random() < pmot:
            insert_2q_motional_error(dataq,ancillaq) # Insert motional error
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind,a_leak_ind))
    Rx(-s*pi/2) | dataq
    if random.random() < pxctrl:
        insert_1q_x_ctrl_error(dataq)

def z_type_entangling(dataq,ancillaq,leaked_reg,d_leak_ind,a_leak_ind,v=1,s=1,p_leak=0,pyctrl=0,pxctrl=0,pxxctrl=0,pmot=0):
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
            insert_2q_motional_error(dataq,ancillaq)
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind,a_leak_ind))
    Rx(-s * pi / 2) | dataq
    if random.random() < pxctrl:
        insert_1q_x_ctrl_error(dataq)
    Ry(-v * pi / 2) | dataq
    if random.random()<pyctrl:
        insert_1q_y_ctrl_error(dataq)

def stabiliser_timestep_1(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
        x_type_entangling(data[4], ancilla[1], pxctrl=0,pyctrl=0,pmot=0,pd=0)
        x_type_entangling(data[8], ancilla[6], pxctrl=0,pyctrl=0,pmot=0,pd=0)
        x_type_entangling(data[6], ancilla[4], pxctrl=0,pyctrl=0,pmot=0,pd=0)

def stabiliser_timestep_2(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    # TODO:
    # put the qubits in in the right order in the stabiliser timesteps
    x_type_entangling(data[1], ancilla[1], pxctrl=0,pyctrl=0,pmot=0,pd=0)
    x_type_entangling(data[3], ancilla[4], pxctrl=0,pyctrl=0,pmot=0,pd=0)
    x_type_entangling(data[5], ancilla[6], pxctrl=0,pyctrl=0,pmot=0,pd=0)

def stabiliser_timestep_3(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    x_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    x_type_entangling(data[3], ancilla[4], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    x_type_entangling(data[5], ancilla[6], pxctrl=0, pyctrl=0, pmot=0, pd=0)

def stabiliser_timestep_4(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    x_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    x_type_entangling(data[3], ancilla[4], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    x_type_entangling(data[5], ancilla[6], pxctrl=0, pyctrl=0, pmot=0, pd=0)

def stabiliser_timestep_5(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)

def stabiliser_timestep_6(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)

def stabiliser_timestep_7(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)

def stabiliser_timestep_8(data,ancilla,v=1,s=1,pxctrl=0,pyctrl=0,pmot=0,pd=0):
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)
    z_type_entangling(data[1], ancilla[1], pxctrl=0, pyctrl=0, pmot=0, pd=0)

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
    stabiliser_timestep_1(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_2(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_3(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
    stabiliser_timestep_4(data,ancilla,pxctrl,pxxctrl,pyctrl,pmot,pd)
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