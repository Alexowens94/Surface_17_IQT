import random
import time
from compiled_surface_code_arbitrary_error_model import *
import numpy as np
from scipy.optimize import curve_fit
from projectq import MainEngine
from projectq.backends import Simulator
import matplotlib.pyplot as plt

def calculate_log_e_rate(num_runs, filename, correction_table, e_model, e_probs, num_bins=20):
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

    round_count_list = []
    leaked_q_reg = 17 * [0]  # a register (initialised as 0s) to track which qubits  leaked, LEAKED == 1 0-8 data,9-16 ancilla

    t_begin = time.perf_counter()
    for _ in range(num_runs):
        # randomly choose basis and initial logical state - keep consistent throughout all rounds within run
        if random.random() < 0.5:
            basis = 'Z'
        else:
            basis = 'X'
        if random.random() < 0.5:
            state = 0
        else:
            state = 1

        eng = MainEngine(Simulator())
        round_count=0
        #max_rounds=1 # TODO: Catch as an exception
        #while round_count < max_rounds:
        while True:
            # Initialise qubit register
            data = eng.allocate_qureg(9)
            ancilla = eng.allocate_qureg(8)
            # prepare logical qubit
            quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model,e_probs)
            # Error correction cycle
            prev_syndrome=np.array(quiescent)

            syndrome=np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))

            flips_a=(prev_syndrome - syndrome) % 2
            prev_syndrome=syndrome
            if np.all((flips_a == 0)):
                ft_syndrome=flips_a
            else:
                syndrome = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))

                flips_b = (prev_syndrome - syndrome) % 2
                # prev_syndrome=syndrome
                ft_syndrome = (flips_a + flips_b) % 2
            error_vec = lookup(ft_syndrome, correction_table)
            apply_correction(error_vec, data)
            eng.flush()
            round_count += 1
            # print('round {}'.format(round_count))

            # #Measure logical qubit
            logical_meas = logical_measurement(data, basis, eng)
            if logical_meas != state:
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
    t_end=time.perf_counter()
    total_time_taken=t_end - t_begin
    print('total time taken {}'.format(total_time_taken))
    print(round_count_list)

    results={

        'total_time_taken_{}_runs'.format(num_runs): total_time_taken,
        'probability_set': e_probs,
        'rounds_til_fail_list': round_count_list,
        'rounds_til_fail_list_X0': round_count_list_X0,
        'rounds_til_fail_list_X1': round_count_list_X1,
        'rounds_til_fail_list_Z0': round_count_list_Z0,
        'rounds_til_fail_list_Z1': round_count_list_Z1,
        'error_model': e_model,
    }
    with open("data/" + filename + '.json', 'w') as file:
        json.dump(results, file)
    log_e_rate = plot_log_e_rate_graph(filename, results, num_bins=num_bins)
    return log_e_rate

def plot_log_e_rate_graph(filename, results, num_bins, save=True):
    data = results['rounds_til_fail_list']
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


def f(t, A, r):
    return A * np.exp(-r * t)

# with open("data\PDD_no_deph_p=0.005.json",'r') as file:
#     results = json.load(file)
# plot_log_e_rate_graph("PDD_no_deph_p=0.005", results,20, save=False)

