import random
import time
import matplotlib.pyplot as plt

import numpy as np

from error import *
from logical_qubit import *

error_subset = {"XCtrlError": 5, "DephasingError": 6}

def calc_error_subset_error_rate(error_subset):
    failed_count = 0
    t_begin = time.perf_counter()
    num_runs = 50
    for _ in range(num_runs):
        leaked_q_reg = 17 * [0]  # tracks which qubits  leak, LEAKED == 1. indices -> 0-8 data,9-16 ancilla

        # randomly choose basis and initial logical state - keep consistent throughout all rounds within run
        if random.random() < 1:
            basis = 'Z'
        else:
            basis = 'X'
        if random.random() < 1:
            state = 0
        else:
            state = 1

        error_list = Error.generate_error_list(error_subset)

        logical_qubit = LogicalQubit()

        # Error correction cycle
        # error_test(round_count, 'X', data)
        prev_syndrome = np.array(logical_qubit.quiescent_state)
        # print('real syndrome meas')
        syndrome = logical_qubit.measure_syndrome([1],error_list=error_list)
        syndrome = np.array(syndrome)
        flips_a = (prev_syndrome - syndrome) % 2
        error_vec = logical_qubit.lookup(flips_a)
        logical_qubit.apply_correction(error_vec)
        logical_qubit.eng.flush()

        # Measure logical qubit
        logical_meas = logical_qubit.measure_qubit()

        if logical_meas != state:
            failed_count += 1
    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    log_e_rate = failed_count / num_runs
    print(log_e_rate)
    print('total_time_taken_{}_runs_{}'.format(num_runs, total_time_taken))
    return log_e_rate

error_subset_list = create_error_subset_list()

def calc_many_error_subset_error_rates(error_subset_list):
    log_e_rate_dict = {}
    for error_subset in error_subset_list:
        print("error subset : {}".format(error_subset))
        log_e_rate = calc_error_subset_error_rate(error_subset)
        log_e_rate_dict[error_subset["DephasingError"]] = log_e_rate
        print("log_e_rate : {}".format(log_e_rate))
    print("log_e_rate_dict : {}".format(log_e_rate_dict))
    return log_e_rate_dict

filename = "importance_sampling_1"
log_e_rate_dict = calc_many_error_subset_error_rates(error_subset_list)

with open("imp_samp/" + filename + '.json', 'w') as file:
    json.dump(log_e_rate_dict, file)

print("error subset list: {}".format(error_subset_list))
plt.bar(log_e_rate_dict.keys(), log_e_rate_dict.values(), width=0.1)
plt.ylabel('logical error rate')
plt.xlabel('error subset')
plt.xticks(rotation='vertical')
plt.show()