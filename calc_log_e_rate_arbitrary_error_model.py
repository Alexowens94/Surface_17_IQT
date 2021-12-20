import random
import time
from compiled_surface_code_arbitrary_error_model import *
import numpy as np
from scipy.optimize import curve_fit
from projectq import MainEngine
from projectq.backends import Simulator
import matplotlib.pyplot as plt

def error_test(index, type, data):
    index = index % 9
    if type == 'X':
        X | data[index]
        print('X {}'.format(index))
    if type == 'Y':
        Y | data[index]
        print('Y {}'.format(index))
    if type == 'Z':
        Z | data[index]
        print('Z {}'.format(index))

def ancilla_test(index):
    index = index % 8
    faulty_synd = np.zeros(8, dtype=int)
    faulty_synd[index]=1
    print('a error {}'.format(index))
    return faulty_synd

def calculate_log_e_rate(num_runs, filename, correction_table, e_model, e_probs, bitflip_table, num_bins=20):
    """
    Iterates through a number of runs.
    Each run fails after a number of rounds of error correction.
    Uses these numbers to produce a bar graph: round_num vs. number of runs that failed on this round.
    Fits the bar graph to find the logical error rate.
    Saves the bar graph in file_name.
    """
    round_count_list_X0 = []
    round_count_list_X1 = []
    round_count_list_Z0 = []
    round_count_list_Z1 = []
    syndrome_b_list = []

    round_count_list = []
    leaked_q_reg = 17 * [0]  # a register (initialised as 0s) to track which qubits  leaked, LEAKED == 1 0-8 data,9-16 ancilla

    t_begin = time.perf_counter()
    for _ in range(num_runs):
        # randomly choose basis and initial logical state - keep consistent throughout all rounds within run
        if random.random() < 1:
            basis = 'Z'
        else:
            basis = 'X'
        if random.random() < 1:
            state = 0
        else:
            state = 1
        # print('b = {}, s = {}'.format(basis, state))
        eng = MainEngine(Simulator())
        round_count = 0
        syndrome_b_count = 0
        #max_rounds=1 # TODO: Catch as an exception
        #while round_count < max_rounds:
        while True:
            # Initialise qubit register
            data = eng.allocate_qureg(9)
            ancilla = eng.allocate_qureg(8)
            # prepare logical qubit
            quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, e_probs)
            # print('quiescent {}'.format(quiescent))
            # Error correction cycle
            # error_test(round_count, 'X', data)
            prev_syndrome = np.array(quiescent)
            # print('q {}'.format(prev_syndrome))
            syndrome = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))  #ancilla_test(round_count)#
            # print(syndrome)
            flips_a = (prev_syndrome - syndrome) % 2
            # print('fa {}'.format(flips_a))
            prev_syndrome = syndrome
            if np.all((flips_a == 0)):
                ft_syndrome = flips_a
            else:
                syndrome = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))
                syndrome_b_count += 1
                # print(syndrome)
                flips_b = (prev_syndrome - syndrome) % 2
                # print('fb {}'.format(flips_b))
                # prev_syndrome=syndrome
                ft_syndrome = (flips_a + flips_b) % 2   # this is equivalent to B - quiescent
            # print('ft {}'.format(ft_syndrome))
            error_vec = lookup(ft_syndrome, correction_table)
            apply_correction(error_vec, data)
            round_count += 1

            eng.flush()
            # print('round {}'.format(round_count))

            # #Measure logical qubit
            logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent, round_count)
            # print(logical_meas)
            if logical_meas != state:
                # print('failed round {}'.format(round_count))

                break

        if basis == "X":
            if state == 0:
                round_count_list_X0.append(round_count)
            if state == 1:
                round_count_list_X1.append(round_count)
        if basis == "Z":
            if state == 0:
                round_count_list_Z0.append(round_count)
            if state == 1:
                round_count_list_Z1.append(round_count)

        round_count_list.append(round_count)
        syndrome_b_list.append(syndrome_b_count)
    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    print('total time taken {}'.format(total_time_taken))
    print(round_count_list)

    results = {

        'total_time_taken_{}_runs'.format(num_runs): total_time_taken,
        'probability_set': e_probs,
        'rounds_til_fail_list': round_count_list,
        'syndrome_b_measured_list': syndrome_b_list,
        'rounds_til_fail_list_X0': round_count_list_X0,
        'rounds_til_fail_list_X1': round_count_list_X1,
        'rounds_til_fail_list_Z0': round_count_list_Z0,
        'rounds_til_fail_list_Z1': round_count_list_Z1,
        'error_model': e_model,
    }
    with open("data/" + filename + '.json', 'w') as file:
        json.dump(results, file)
    log_e_rate = plot_log_e_rate_graph(filename, results, num_bins=num_bins, include_b=False)
    return log_e_rate


def plot_log_e_rate_graph(filename, results, num_bins, include_b, save=True):
    data = np.array(results['rounds_til_fail_list'], dtype='float64')
    if include_b:
        bs = np.array(results['syndrome_b_measured_list'], dtype='float64')
        print(bs)
        data += bs
    vals, bins, patches = plt.hist(data, align='left', rwidth=0.5, bins=num_bins)
    popt, pcov = curve_fit(f, bins[:-1], vals)
    plt.ylabel('runs ending after t QEC round')
    plt.xlabel('t')
    plt.plot(bins, f(bins, *popt), 'r-')
    log_e_rate = popt[1]
    title = filename + " log_e_rate = {}, A = {}".format(log_e_rate, popt[0])
    plt.title(title)
    if save:
        plt.savefig("figs/" + filename + ".png")
    else:
        plt.show()
    plt.close()
    return log_e_rate

def plot_log_e_rate_graphs(filename, results, num_bins, save=True):
    log_e_rate = []
    for data in [results['rounds_til_fail_list_X0'],results['rounds_til_fail_list_X1'],
                 results['rounds_til_fail_list_Z0'],results['rounds_til_fail_list_Z1']]:
        vals, bins, patches = plt.hist(data, align='left', rwidth=0.5, bins=num_bins)
        popt, pcov = curve_fit(f, bins[:-1], vals)

        plt.plot(bins, f(bins, *popt))
        log_e_rate.append(popt[1])
    plt.ylabel('runs ending after t QEC round')
    plt.xlabel('t')
    plt.legend(['X0 {}'.format(log_e_rate[0]),
                'X1 {}'.format(log_e_rate[1]),
                'Z0 {}'.format(log_e_rate[2]),
                'Z1 {}'.format(log_e_rate[3])])
    plt.title(filename)

    if save:
        plt.savefig("figs/" + filename + ".png")
    else:
        plt.show()
    plt.close()
    # return log_e_rate

def f(t, A, r):
    return A * np.exp(-r * t)

# filename = "data/classical_corrections_1000_p=0.007"
# with open(filename+'.json', 'r') as file:
#     results = json.load(file)
# plot_log_e_rate_graph(filename, results, 20, include_b=True, save=False)
#
