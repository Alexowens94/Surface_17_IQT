import random
import numpy as np
import json
from math import pi
from projectq.ops import All, Measure, X, Y, Z, H, Rx, Ry, Rxx, Rzz


def instantiate_error_model_variable_probs(p_set_list, model):  # think this is unnecesary - other version achieves this
    """
    :param p_set_list: list of list of float, probabilities to associate to the different types of errors as per model
    first index is which type of error, second index is error location (gate index)
    :param model: string, specifies which error model in error_models.json to use
    :return: dict error model probabilities (key string label as given in error model file, value float from p_set)
    """
    error_model_probabilities = {}
    with open('error_models.json', 'r') as infile:
        error_models = json.load(infile)
        error_model = error_models[model]
        try:
            if len(error_model['probabilities']) == len(p_set_list):
                for i, key in enumerate(error_model['probabilities']):
                    error_model_probabilities[key] = p_set_list[
                        i]  # list of probabilities for all locations of error location type i
            else:
                raise ValueError
        except ValueError:
            print('make sure you pass enough probabilities to specify the error model - see error_model.json')
    return error_model, error_model_probabilities


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


def logical_measurement(data, basis, eng, classical_lookup, quiescent, leaked_q_reg):
    if basis == 'X':
        All(H) | data  # change Z -> X basis

    All(Measure) | data
    eng.flush()  # flush all gates (and execute measurements)
    data_meas = [int(q) for q in data]
    # print('leaked reg {}'.format(leaked_q_reg))
    # print('without leak {}'.format(data_meas))
    for i in range(len(data_meas)):
        if leaked_q_reg[i] == 1:
            data_meas[i] = 1

    # print('with leak {}'.format(data_meas))
    # print('d_meas {}'.format(data_meas))
    # classical corrections TODO: implement for X basis prep
    weight = 0
    if basis == 'Z':  # from measuring data qubits in z basis, can infer a z stabiliser measurement,
        # to perform bit flip corrections (if measuring in x basis, can correct phase flips)
        q = [quiescent[0], quiescent[2], quiescent[5], quiescent[7]]
        classical_syndrome = [0, 0, 0, 0]
        classical_syndrome[0] = (data_meas[0] + data_meas[1]) % 2
        classical_syndrome[1] = (data_meas[1] + data_meas[2] + data_meas[4] + data_meas[5]) % 2
        classical_syndrome[2] = (data_meas[3] + data_meas[4] + data_meas[6] + data_meas[7]) % 2
        classical_syndrome[3] = (data_meas[7] + data_meas[8]) % 2
        csyndrome = (np.array(classical_syndrome) - np.array(q)) % 2
        # print(csyndrome)
        weight = return_weight_classical_lookup(csyndrome, classical_lookup)
        # print(weight)
    # ZL1 = (data_meas[0] + data_meas[3] + data_meas[6]) % 2
    # ZL2 = (data_meas[1] + data_meas[4] + data_meas[7]) % 2
    # ZL3 = (data_meas[2] + data_meas[5] + data_meas[8]) % 2
    # print(data_meas)
    logic_measurement = (sum(data_meas) + weight) % 2  # weight adds classical correction
    # if (ZL1+ZL2+ZL3) % 3 != 0:
    #     print('round {}'.format(round))
    #     print('ZLs {} {} {}'.format(ZL1, ZL2, ZL3))
    #     print('l_meas (corrected) {}'.format(logic_measurement))
    #     print('synd {}'.format(csyndrome))
    return logic_measurement


def logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, e_probs, cz_compilation=False):
    '''
    prepare the logical qubit, return the quiescent state
    :param data: list of data qubits
    :param basis: string "X" or "Z", basis to prepare in
    :param state: int 0 or 1, state you wish to prepare interpretting 0(1) in X basis as |+>(|->)
    '''
    if state == 1:
        All(X) | data  # rotate to orthogonal state (all qubits start in 0, Z basis)
    if basis == 'X':
        All(H) | data  # change Z -> X basis

    # Perfect logical state prep (well motivated as per
    # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
    # quiescent prepared perfectly (all errors = 0)
    errorless_dict = {}  # because of the way ive coded probabilities as dictionary with keys known apriori,
    # need to construct the 0 error version at run time to pass as an argument for perfect quiescent state prep
    for key in e_probs.keys():
        errorless_dict[key] = len(e_probs[key]) * [0]
    # quiescent = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, errorless_dict, reset=True))
    if cz_compilation:
        quiescent = np.array(cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, errorless_dict, reset=True))
    else:
        quiescent = np.array(stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, errorless_dict, reset=True))
    return quiescent


def insert_errors(gate, qubits, leaked_reg, error_model, e_model_probs, c_ind=-1, t_ind=-1, d_ind=-1, display=False,
                  rx_ind=0, ry_ind=0, rxx_ind=0, rzz_ind=0):
    num_gates_ms_comp = [12, 18, 24]
    num_gates_cz_comp = [34, 24]
    error_list = error_model['gates'][gate]
    for e in error_list:
        index = None
        if display:
            print(e[2])
        if e[1] == 'p_deph':
            #print('doesnt work right _ split 1q deph and 2q deph up to have 2 different placeholder probs')
            if gate == 'Rx':
                index = rx_ind
            if gate == 'Ry':
                index = num_gates_ms_comp[0] + ry_ind
            if gate == 'Rxx':
                if e[0] == 'Zc':
                    index = num_gates_ms_comp[0]+num_gates_ms_comp[1]+rxx_ind
                    # print('error {} index {}'.format(e[0], index))
                if e[0] == 'Zt':
                    index = num_gates_ms_comp[0] + num_gates_ms_comp[1] + num_gates_ms_comp[2] + rxx_ind
                    # print('error {} index {}'.format(e[0], index))
        if e[1] == 'p_X_ctrl':
            index = rx_ind
        if e[1] == 'p_XX_ctrl' or e[1] == 'p_heat':
            index = rxx_ind
        if e[1] == 'p_Y_ctrl':
            index = ry_ind
        if e[1] == 'p_leak':
            if gate == 'Rx':
                index = rx_ind
            if gate == 'Ry':
                index = num_gates_ms_comp[0] + ry_ind
            if gate == 'Rxx':
                if e[0] == 'Lc':
                    index = num_gates_ms_comp[0] + num_gates_ms_comp[1] + rxx_ind
                    # print('error {} index {}'.format(e[0], index))
                if e[0] == 'Lt':
                    index = num_gates_ms_comp[0] + num_gates_ms_comp[1] + num_gates_ms_comp[2] + rxx_ind
                    # print('error {} index {}'.format(e[0], index))
        #### cz comp error indexing
        if e[1] == 'p_1q_deph':
            index = ry_ind+rx_ind
        if e[1] == 'p_2q_deph':
            index = rzz_ind
        if e[1] == 'p_ctrl_leak':
            index = rzz_ind
        if e[1] == 'p_env_leak':
            if e[0] == 'Lc':
                index = rzz_ind
            if e[0] == 'Lt':
                index = num_gates_cz_comp[1] + rzz_ind
        # print('e1 {}'.format(e[1]))
        insert_error(error=e[0], prob=e_model_probs[e[1]], leaked_reg=leaked_reg, qubits=qubits,
                     c_ind=c_ind, t_ind=t_ind, d_ind=d_ind, index=index)


def insert_error(error, prob, leaked_reg, qubits, c_ind, t_ind, d_ind, index=0):
    # print('{} location index {}'.format(error, index))
    # print('prob {}'.format(prob))
    if random.random() < prob[index]:
        # print(error, index)
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
        if error == 'ZZ':
            Z | qubits[0]
            Z | qubits[1]
        if error == 'L':
            insert_leakage(d_ind, leaked_reg)
        if error == 'LL':
            # print(leaked_reg)
            insert_leakage(c_ind, leaked_reg)
            insert_leakage(t_ind, leaked_reg)
            # print(leaked_reg)
        if error == 'Lc':
            # print(leaked_reg)
            insert_leakage(c_ind, leaked_reg)
            # print(leaked_reg)
        if error == 'Lt':
            # print(leaked_reg)
            insert_leakage(t_ind, leaked_reg)
            # print(leaked_reg)
        if error == 'depol':
            insert_2q_error(qubits[0], qubits[1])


def insert_2q_error(qubit1, qubit2):
    rand = random.random()
    if rand < 1 / 15:
        # I|qubit1
        X | qubit2
    if 1 / 15 < rand < (2 * 1 / 15):
        # I|qubit1
        Y | qubit2
    if (2 * 1 / 15) < rand < (3 * 1 / 15):
        # I|qubit1
        Z | qubit2
    if (3 * 1 / 15) < rand < (4 * 1 / 15):
        X | qubit1
        # I|qubit2
    if (4 * 1 / 15) < rand < (5 * 1 / 15):
        X | qubit1
        X | qubit2
    if (5 * 1 / 15) < rand < (6 * 1 / 15):
        X | qubit1
        Y | qubit2
    if (6 * 1 / 15) < rand < (7 * 1 / 15):
        X | qubit1
        Z | qubit2
    if (7 * 1 / 15) < rand < (8 * 1 / 15):
        Y | qubit1
        # I|qubit2
    if (8 * 1 / 15) < rand < (9 * 1 / 15):
        Y | qubit1
        X | qubit2
    if (9 * 1 / 15) < rand < (10 * 1 / 15):
        Y | qubit1
        Y | qubit2
    if (10 * 1 / 15) < rand < (11 * 1 / 15):
        Y | qubit1
        Z | qubit2
    if (11 * 1 / 15) < rand < (12 * 1 / 15):
        Z | qubit1
        # I|qubit2
    if (12 * 1 / 15) < rand < (13 * 1 / 15):
        Z | qubit1
        X | qubit2
    if (13 * 1 / 15) < rand < (14 * 1 / 15):
        Z | qubit1
        Y | qubit2
    if (14 * 1 / 15) < rand < (15 * 1 / 15):
        Z | qubit1
        Z | qubit2


def insert_leakage(index, leaked_reg):
    leaked_reg[index] = 1
    # print('leakage at {}'.format(index))


def x_type_entangling(dataq, ancillaq, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0,
                      d_leak_ind=0, a_leak_ind=0, cancel_data_rx=False, s=1):
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s * pi / 2) | (dataq, ancillaq)
        insert_errors(gate='Rxx', qubits=[dataq, ancillaq], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs,
                      c_ind=a_leak_ind, t_ind=d_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)

    # else:
    #     print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    rxx_ind += 1
    if not cancel_data_rx:
        Rx(-s * pi / 2) | dataq
        insert_errors(gate='Rx', qubits=dataq, leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
        rx_ind += 1
    return rx_ind, rxx_ind


def z_type_entangling(dataq, ancillaq, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0,
                      d_leak_ind=0, a_leak_ind=0, cancel_data_rx=False, cancel_ry1=False, cancel_ry2=False, s=1, v=1):
    if not cancel_ry1:
        Ry(v * pi / 2) | dataq

        insert_errors(gate='Ry', qubits=dataq, leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
        ry_ind += 1

    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rxx(s * pi / 2) | (dataq, ancillaq)
        insert_errors(gate='Rxx', qubits=[dataq, ancillaq], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs,
                      c_ind=d_leak_ind, t_ind=a_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)

    else:
        print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
    rxx_ind += 1
    if not cancel_data_rx:
        Rx(-s * pi / 2) | dataq
        insert_errors(gate='Rx', qubits=dataq, leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
        rx_ind += 1
    if not cancel_ry2:
        Ry(-v * pi / 2) | dataq
        insert_errors(gate='Ry', qubits=dataq, leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=d_leak_ind,
                      rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
        ry_ind += 1
    return rx_ind, ry_ind, rxx_ind


def cz_entangling(dataq, ancillaq, leaked_reg, error_model, e_model_probs, rzz_ind=0, d_leak_ind=0, a_leak_ind=0):
    # entangling step for use in spin spin gate compilation
    # need to be sandwiched with root ys - see cz_x_stab func
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 0:
        Rzz(pi) | (dataq, ancillaq)
    if leaked_reg[d_leak_ind] == 1 and leaked_reg[a_leak_ind] == 0:
        # print('leak data{} ancilla{}'.format(d_leak_ind, a_leak_ind))
        Z | ancillaq
    if leaked_reg[d_leak_ind] == 0 and leaked_reg[a_leak_ind] == 1:
        Z | dataq
    insert_errors(gate='Rzz', qubits=[dataq, ancillaq], leaked_reg=leaked_reg,
                  error_model=error_model, e_model_probs=e_model_probs,
                  c_ind=d_leak_ind, t_ind=a_leak_ind, rzz_ind=rzz_ind)
    rzz_ind += 1

    return rzz_ind


def stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, xx_inc = x_type_entangling(data[4], ancilla[1], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=4, a_leak_ind=1 + 9, cancel_data_rx=True, s=1)
    # print('stab1 ent1 {} {}'.format(x_inc,xx_inc))
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[8], ancilla[6], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=8, a_leak_ind=6 + 9, s=1)
    # print('stab1 ent2 {} {}'.format(x_inc, xx_inc))
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[6], ancilla[4], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=6, a_leak_ind=4 + 9, s=1)
    # print('stab1 ent3 {} {}'.format(x_inc, xx_inc))
    rx_ind = x_inc
    rxx_ind = xx_inc
    return rx_ind, rxx_ind


def cz_stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[4], ancilla[1], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=4, a_leak_ind=1 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[8], ancilla[6], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=8, a_leak_ind=6 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[6], ancilla[4], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=6, a_leak_ind=4 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[1], ancilla[1], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=1, a_leak_ind=1 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[5], ancilla[6], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=5, a_leak_ind=6 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[3], ancilla[4], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=3, a_leak_ind=4 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[3], ancilla[1], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=3, a_leak_ind=1 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[7], ancilla[6], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=7, a_leak_ind=6 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[5], ancilla[3], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=5, a_leak_ind=3 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[0], ancilla[1], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=0, a_leak_ind=1 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[4], ancilla[6], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=4, a_leak_ind=6 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[2], ancilla[3], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=2, a_leak_ind=3 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[1], ancilla[2], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=1, a_leak_ind=2 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[3], ancilla[5], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=3, a_leak_ind=5 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[7], ancilla[7], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=7, a_leak_ind=7 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[2], ancilla[2], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=2, a_leak_ind=2 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[4], ancilla[5], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=4, a_leak_ind=5 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[8], ancilla[7], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=8, a_leak_ind=7 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[4], ancilla[2], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=4, a_leak_ind=2 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[6], ancilla[5], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=6, a_leak_ind=5 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[0], ancilla[0], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=0, a_leak_ind=0 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def cz_stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=0):
    zz_inc = cz_entangling(data[5], ancilla[2], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=5, a_leak_ind=2 + 9)
    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[7], ancilla[5], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=7, a_leak_ind=5 + 9)

    rzz_ind = zz_inc
    zz_inc = cz_entangling(data[1], ancilla[0], leaked_reg, error_model, e_model_probs, rzz_ind,
                           d_leak_ind=1, a_leak_ind=0 + 9)
    rzz_ind = zz_inc
    return rzz_ind


def stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, xx_inc = x_type_entangling(data[1], ancilla[1], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=1, a_leak_ind=1 + 9, s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[5], ancilla[6], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=5, a_leak_ind=6 + 9, cancel_data_rx=True,
                                      s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[3], ancilla[4], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=3, a_leak_ind=4 + 9, cancel_data_rx=True,
                                      s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    return rx_ind, rxx_ind


def stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, xx_inc = x_type_entangling(data[3], ancilla[1], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=3, a_leak_ind=1 + 9, cancel_data_rx=True, s=1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[7], ancilla[6], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=7, a_leak_ind=6 + 9, s=1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[5], ancilla[3], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=5, a_leak_ind=3 + 9, cancel_data_rx=True, s=1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    return rx_ind, rxx_ind


def stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, xx_inc = x_type_entangling(data[0], ancilla[1], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=0, a_leak_ind=1 + 9, s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[4], ancilla[6], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=4, a_leak_ind=6 + 9, cancel_data_rx=True,
                                      s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = x_type_entangling(data[2], ancilla[3], leaked_reg, error_model, e_model_probs,
                                      rx_ind, ry_ind, rxx_ind, d_leak_ind=2, a_leak_ind=3 + 9, s=-1)
    rx_ind = x_inc
    rxx_ind = xx_inc
    return rx_ind, rxx_ind


def stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, y_inc, xx_inc = z_type_entangling(data[1], ancilla[2], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=1, a_leak_ind=2 + 9, cancel_data_rx=True, cancel_ry2=True, v=1,
                                             s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[3], ancilla[5], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=1, a_leak_ind=2 + 9, v=1, s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[7], ancilla[7], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=7, a_leak_ind=7 + 9, cancel_data_rx=True, cancel_ry2=True, v=1,
                                             s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    return rx_ind, ry_ind, rxx_ind


def stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, y_inc, xx_inc = z_type_entangling(data[2], ancilla[2], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=2, a_leak_ind=2 + 9, v=1, s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[4], ancilla[5], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=4, a_leak_ind=5 + 9, cancel_data_rx=True, cancel_ry2=True, v=1,
                                             s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[8], ancilla[7], leaked_reg, error_model, e_model_probs, rx_ind,
                                             ry_ind, rxx_ind,
                                             d_leak_ind=8, a_leak_ind=7 + 9, v=1, s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    return rx_ind, ry_ind, rxx_ind


def stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, y_inc, xx_inc = z_type_entangling(data[4], ancilla[2], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=4, a_leak_ind=2 + 9, cancel_data_rx=True, cancel_ry1=True, v=1,
                                             s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[6], ancilla[5], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=6, a_leak_ind=5 + 9, v=1, s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[0], ancilla[0], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=0, a_leak_ind=0 + 9, v=1, s=1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    return rx_ind, ry_ind, rxx_ind


def stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs, rx_ind=0, ry_ind=0, rxx_ind=0):
    x_inc, y_inc, xx_inc = z_type_entangling(data[5], ancilla[2], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=5, a_leak_ind=2 + 9, v=1, s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[7], ancilla[5], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=7, a_leak_ind=5 + 9, cancel_data_rx=True, cancel_ry1=True, v=1,
                                             s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = z_type_entangling(data[1], ancilla[0], leaked_reg, error_model, e_model_probs,
                                             rx_ind, ry_ind, rxx_ind,
                                             d_leak_ind=1, a_leak_ind=0 + 9, cancel_data_rx=True, cancel_ry1=True, v=1,
                                             s=-1)
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    return rx_ind, ry_ind, rxx_ind


def stabilizer_cycle(data, ancilla, leaked_reg, eng, error_model, e_model_probs, reset=True):
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    if len(data) != 9:
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
    # print('leaked reg {}'.format(leaked_reg))
    syndrome_t = [int(q) for q in ancilla]
    # print('a without leak {}'.format(syndrome_t))
    # account for measuring a leaked qubit to always give 1
    for i in range(len(syndrome_t)):
        if leaked_reg[9+i] == 1:
            syndrome_t[i] = 1
            leaked_reg[9+i] = 0  # reset leaked ancillas post measurement
    # print('with leak {}'.format(syndrome_t))
    if reset:
        for a in ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
            if int(a) == 1:
                X | a
    return syndrome_t


def cz_stabilizer_cycle(data, ancilla, leaked_reg, eng, error_model, e_model_probs, reset=True):
    """ stabilizer cycle compilation for minimum 1q gates using cz instead of cx
    as native entangling gate (spin-spin coupling gate)"""
    rzz_ind, ry_ind = cz_x_stabilizers(data, ancilla, leaked_reg, error_model, e_model_probs)
    cz_z_stabilizers(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind, ry_ind)
    All(Measure) | ancilla
    eng.flush()
    # print('leaked reg {}'.format(leaked_reg))
    syndrome_t = [int(q) for q in ancilla]
    # print('a without leak {}'.format(syndrome_t))
    # account for measuring a leaked qubit to always give 1
    for i in range(len(syndrome_t)):
        if leaked_reg[9+i] == 1:
            syndrome_t[i] = 1
            leaked_reg[9+i] = 0  # reset leaked ancillas post measurement
    # print('a with leak {}'.format(syndrome_t))
    if reset:
        for a in ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
            if int(a) == 1:
                X | a
    return syndrome_t


def cz_x_stabilizers(data, ancilla, leaked_reg, error_model, e_model_probs):
    All(Ry(-pi / 2)) | ancilla
    ry_ind = 0
    for i in [1, 3, 4, 6]:
        insert_errors(gate='Ry', qubits=ancilla[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=9 + i, ry_ind=ry_ind)
        ry_ind += 1
    All(Ry(-pi / 2)) | data
    for i in range(len(data)):
        insert_errors(gate='Ry', qubits=data[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=i, ry_ind=ry_ind)
        ry_ind += 1
    rzz_ind = cz_stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs)
    rzz_ind = cz_stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    rzz_ind = cz_stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    rzz_ind = cz_stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    All(Ry(pi / 2)) | ancilla
    for i in [1, 3, 4, 6]:
        insert_errors(gate='Ry', qubits=ancilla[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=9 + i, ry_ind=ry_ind)
        ry_ind += 1
    All(Ry(pi / 2)) | data
    for i in range(len(data)):
        insert_errors(gate='Ry', qubits=data[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=i, ry_ind=ry_ind)
        ry_ind += 1
    return rzz_ind, ry_ind


def cz_z_stabilizers(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind, ry_ind):
    All(Ry(-pi / 2)) | ancilla
    for i in [0, 2, 5, 7]:
        insert_errors(gate='Ry', qubits=ancilla[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=9 + i, ry_ind=ry_ind)
        ry_ind += 1
    rzz_ind = cz_stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    rzz_ind = cz_stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    rzz_ind = cz_stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    rzz_ind = cz_stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs, rzz_ind=rzz_ind)
    All(Ry(pi / 2)) | ancilla
    for i in [0, 2, 5, 7]:
        insert_errors(gate='Ry', qubits=ancilla[i], leaked_reg=leaked_reg,
                      error_model=error_model, e_model_probs=e_model_probs, d_ind=9 + i, ry_ind=ry_ind)
        ry_ind += 1
    # print('final indices Rzz {} Ry {}'.format(rzz_ind, ry_ind))
    return rzz_ind


def stabilizer_cycle_error_index(data, ancilla, leaked_reg, eng, error_model, e_model_probs, reset=True):
    '''
    :param data: list of the data qubits
    :param ancilla: list of the ancilla qubits
    :param eng:
    :param reset:
    :return:
    '''
    if len(data) != 9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')

    rx_ind = 0
    ry_ind = 0
    rxx_ind = 0
    x_inc, xx_inc = stabiliser_timestep_1(data, ancilla, leaked_reg, error_model, e_model_probs,
                                          rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("1")
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = stabiliser_timestep_2(data, ancilla, leaked_reg, error_model, e_model_probs,
                                          rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("2")
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = stabiliser_timestep_3(data, ancilla, leaked_reg, error_model, e_model_probs,
                                          rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("3 ryind {}".format(ry_ind))

    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, xx_inc = stabiliser_timestep_4(data, ancilla, leaked_reg, error_model, e_model_probs,
                                          rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("4 ryind {}".format(ry_ind))
    rx_ind = x_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = stabiliser_timestep_5(data, ancilla, leaked_reg, error_model, e_model_probs,
                                                 rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("5")
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = stabiliser_timestep_6(data, ancilla, leaked_reg, error_model, e_model_probs,
                                                 rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("6")
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    x_inc, y_inc, xx_inc = stabiliser_timestep_7(data, ancilla, leaked_reg, error_model, e_model_probs,
                                                 rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("7")
    rx_ind = x_inc
    ry_ind = y_inc
    rxx_ind = xx_inc
    stabiliser_timestep_8(data, ancilla, leaked_reg, error_model, e_model_probs,
                          rx_ind=rx_ind, ry_ind=ry_ind, rxx_ind=rxx_ind)
    # print("8")

    All(Measure) | ancilla
    eng.flush()
    # print('leaked reg {}'.format(leaked_reg))
    syndrome_t = [int(q) for q in ancilla]
    # print('a without leak {}'.format(syndrome_t))
    # account for measuring a leaked qubit to always give 1
    for i in range(len(syndrome_t)):
        if leaked_reg[9 + i] == 1:
            syndrome_t[i] = 1
            leaked_reg[9 + i] = 0  # reset leaked ancillas post measurement
    # print('a with leak {}'.format(syndrome_t))
    if reset:
        for a in ancilla:  # reset the ancillas to 0 at end of stab round (allow for repeat rounds)
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


def lookup(syndrome, table, display=False):
    key = str(syndrome).strip('[,]')
    error_vec = table[key][0]
    if display:
        print('ft syndrome: {}'.format(syndrome))
        print('error vector: {}'.format(error_vec))
    return error_vec


def return_weight_classical_lookup(syndrome, table):
    key = str(syndrome).strip('[,]')
    weight = table[key][1]
    return weight


def apply_correction(error_vec, data):
    for i in range(9):
        if error_vec[i] == 1:
            X | data[i]
        if error_vec[i + 9] == 1:
            Z | data[i]
    return


def load_lookup_table(filename):
    with open(filename, 'r') as infile:
        correction_table = json.load(infile)
    return correction_table
