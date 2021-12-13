import random
import numpy as np
import json
from math import sqrt, pi
from projectq.ops import All, Measure, X, Y, Z, Rx, Ry, Rxx


def instantiate_error_model(p_set, model):
    """
    :param p_set: list float, probabilities to associate to the different types of errors as per model
    :param model: string, specifies which error model in error_models.json to use
    :return: dict error model probabilities (key string label as given in error model file, value float from p_set)
    """
    error_model_probabilities = {}
    with open('error_models.json', 'r') as infile:
        error_models = json.load(infile)
        error_model = error_models[model]
        try:
            if len(error_model['probabilities']) == len(p_set):
                for i, key in enumerate(error_model['probabilities']):
                    error_model_probabilities[key] = p_set[i]
            else:
                raise ValueError
        except ValueError:
            print('make sure you pass enough probabilities to specify the error model - see error_model.json')
    return error_model, error_model_probabilities


def insert_errors(gate, qubits, error_model, e_model_probs, c_ind=-1, t_ind=-1, d_ind=-1, display=False):
    error_list = error_model['gates'][gate]
    for e in error_list:
        if display:
            print(e[2])
        insert_error(error=e[0], prob=e_model_probs[e[1]], qubits=qubits, c_ind=c_ind, t_ind=t_ind, d_ind=d_ind)


def insert_error(error, prob, qubits, c_ind, t_ind, d_ind):
    if random.random() < prob:
        if error == 'X':
            X | qubits
        if error == 'Y':
            Y | qubits
        if error == 'Z':
            Z | qubits
        if error == 'Zc':
            Z | qubits[0]
        if error == 'Zt':
            Z | qubits[1]
        if error == 'XX':
            X | qubits[0]
            X | qubits[1]
        if error == 'L':
            insert_leakage(d_ind)
        if error == 'Lc':
            insert_leakage(c_ind)
        if error == 'Lt':
            insert_leakage(t_ind)


def insert_leakage(qubit):
    print('not implemented - needs the index of the qubit in question passing down from insert errors')


def x_type_entangling(dataq, ancillaq, leaked_reg, error_model, e_model_probs, d_leak_ind, a_leak_ind,
                      cancel_data_rx=False, s=1):
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s*pi/2) | (dataq, ancillaq)
        insert_errors(gate='Rxx', qubits=[dataq, ancillaq], error_model=error_model, e_model_probs=e_model_probs,
                      c_ind=a_leak_ind, t_ind=d_leak_ind)
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    if not cancel_data_rx:
        Rx(-s*pi/2) | dataq
        insert_errors(gate='Rx', qubits=dataq, error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind)

def z_type_entangling(dataq, ancillaq, leaked_reg, error_model, e_model_probs, d_leak_ind, a_leak_ind,
                      cancel_data_rx=False, s=1, v=1):
    Ry(v*pi/2) | dataq
    insert_errors(gate='Ry', qubits=dataq, error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind)
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s * pi / 2) | (dataq, ancillaq)
        insert_errors(gate='Rxx', qubits=[dataq, ancillaq], error_model=error_model, e_model_probs=e_model_probs,
        c_ind = d_leak_ind, t_ind = a_leak_ind)
    else:
        print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    if not cancel_data_rx:
        Rx(-s * pi / 2) | dataq
        insert_errors(gate='Rxx', qubits=dataq, error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind)
    Ry(-v * pi / 2) | dataq
    insert_errors(gate='Ry', qubits=dataq, error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind)

def stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs):
    x_type_entangling(data[4], ancilla[1], leaked_reg,  error_model, e_model_probs, d_leak_ind=4, a_leak_ind=1+9,
                      cancel_data_rx=True, s=1)
    x_type_entangling(data[8], ancilla[6], leaked_reg, error_model, e_model_probs, d_leak_ind=8, a_leak_ind=6+9,  s=1)
    x_type_entangling(data[6], ancilla[4], leaked_reg, error_model, e_model_probs, d_leak_ind=6, a_leak_ind=4+9, s=1)

def stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs):
    x_type_entangling(data[1], ancilla[1], leaked_reg, error_model, e_model_probs, d_leak_ind=1, a_leak_ind=1+9, s=-1)
    x_type_entangling(data[5], ancilla[6], leaked_reg, error_model, e_model_probs, d_leak_ind=5, a_leak_ind=6+9,
                      cancel_data_rx=True, s=-1)
    x_type_entangling(data[3], ancilla[4], leaked_reg, error_model, e_model_probs, d_leak_ind=3, a_leak_ind=4+9,
                      cancel_data_rx=True, s=-1)

def stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs):
    x_type_entangling(data[3], ancilla[1], leaked_reg, error_model, e_model_probs, d_leak_ind=3, a_leak_ind=1+9,
                      cancel_data_rx=True, s=1)
    x_type_entangling(data[7], ancilla[6], leaked_reg, error_model, e_model_probs, d_leak_ind=7, a_leak_ind=6+9, s=1)
    x_type_entangling(data[5], ancilla[3], leaked_reg, error_model, e_model_probs, d_leak_ind=5, a_leak_ind=3+9,
                      cancel_data_rx=True, s=1)

def stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs):
    x_type_entangling(data[0], ancilla[1], leaked_reg, error_model, e_model_probs, d_leak_ind=0, a_leak_ind=1+9, s=-1)
    x_type_entangling(data[4], ancilla[6], leaked_reg, error_model, e_model_probs, d_leak_ind=4, a_leak_ind=6+9,
                      cancel_data_rx=True, s=-1)
    x_type_entangling(data[2], ancilla[3], leaked_reg, error_model, e_model_probs, d_leak_ind=2, a_leak_ind=3+9, s=-1)

def stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs):
    z_type_entangling(data[1], ancilla[2], leaked_reg, error_model, e_model_probs, d_leak_ind=1, a_leak_ind=2+9,
                      cancel_data_rx=True, v=1, s=1)
    z_type_entangling(data[3], ancilla[5], leaked_reg, error_model, e_model_probs, d_leak_ind=1, a_leak_ind=2+9,
                      v=1, s=1)
    z_type_entangling(data[7], ancilla[7], leaked_reg, error_model, e_model_probs, d_leak_ind=7, a_leak_ind=7+9,
                      cancel_data_rx=True, v=1, s=1)

def stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs):
    z_type_entangling(data[2], ancilla[2], leaked_reg, error_model, e_model_probs, d_leak_ind=2, a_leak_ind=2+9,
                      v=1, s=-1)
    z_type_entangling(data[4], ancilla[5], leaked_reg, error_model, e_model_probs, d_leak_ind=4, a_leak_ind=5+9,
                      cancel_data_rx=True, v=1, s=-1)
    z_type_entangling(data[8], ancilla[7], leaked_reg, error_model, e_model_probs, d_leak_ind=8, a_leak_ind=7+9,
                      v=1, s=-1)

def stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs):
    z_type_entangling(data[4], ancilla[2], leaked_reg, error_model, e_model_probs, d_leak_ind=4, a_leak_ind=2+9,
                      cancel_data_rx=True, v=1, s=1)
    z_type_entangling(data[6], ancilla[5], leaked_reg, error_model, e_model_probs, d_leak_ind=6, a_leak_ind=5+9,
                      v=1, s=1)
    z_type_entangling(data[0], ancilla[0], leaked_reg, error_model, e_model_probs, d_leak_ind=0, a_leak_ind=0+9,
                      v=1, s=1)

def stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs):
    z_type_entangling(data[5], ancilla[2], leaked_reg, error_model, e_model_probs, d_leak_ind=5, a_leak_ind=2+9,
                      v=1, s=-1)
    z_type_entangling(data[7], ancilla[5], leaked_reg, error_model, e_model_probs, d_leak_ind=7, a_leak_ind=5+9,
                      cancel_data_rx=True, v=1, s=-1)
    z_type_entangling(data[1], ancilla[0], leaked_reg, error_model, e_model_probs, d_leak_ind=1, a_leak_ind=0+9,
                      cancel_data_rx=True, v=1, s=-1)

def stabilizer_cycle(data, ancilla, leaked_reg, eng, error_model, e_model_probs, reset=True):
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs)
    stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs)

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