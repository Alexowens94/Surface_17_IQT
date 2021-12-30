import random
import time
from compiled_surface_code_arbitrary_error_model import *
import numpy as np
from scipy.optimize import curve_fit
from projectq import MainEngine
from projectq.backends import Simulator
import matplotlib.pyplot as plt
from math import comb

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


def generate_error_locations(num_e, num_e_locations):
    mask = num_e_locations*[0]
    for _ in range(num_e):
        success = False
        while success == False:
            ind = random.randint(0,num_e_locations-1)
            # print('ind {}'.format(ind))
            if mask[ind] == 0:
                mask[ind] = 1
                success = True
    return mask


def calculate_log_e_rate_error_subset(num_runs, filename, correction_table, model, eset, bitflip_table):
    """

    """
    failed_count = 0

    t_begin = time.perf_counter()
    for _ in range(num_runs):
        # in error subset testing over 1 round, so reinitialise leaked_reg within loop
        leaked_q_reg = 17 * [0]  # a register (initialised as 0s) to track which qubits  leaked, LEAKED == 1 0-8 data,9-16 ancilla

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
        # print('eset {}'.format(eset))
        p_set = [generate_error_locations(num_e=eset[0], num_e_locations=12),  # px
                 generate_error_locations(num_e=eset[1], num_e_locations=18),  # py
                 generate_error_locations(num_e=eset[2], num_e_locations=24),  # pxx
                 generate_error_locations(num_e=eset[3], num_e_locations=78),  # dephasing
                 generate_error_locations(num_e=eset[4], num_e_locations=24)  # heating
                 ]
        # print(p_set)
        e_model, model_probs = instantiate_error_model(p_set, model)
        # print(e_model['gates'], model_probs)
        eng = MainEngine(Simulator())

        # Initialise qubit register
        data = eng.allocate_qureg(9)
        ancilla = eng.allocate_qureg(8)
        # prepare logical qubit
        quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, model_probs)
        # Error correction cycle
        # error_test(round_count, 'X', data)
        prev_syndrome = np.array(quiescent)
        # print('real syndrome meas')
        syndrome = np.array(stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, model_probs,
                                                         reset=True))
        flips_a = (prev_syndrome - syndrome) % 2
        error_vec = lookup(flips_a, correction_table)
        apply_correction(error_vec, data)

        eng.flush()
        # print('round {}'.format(round_count))

        # #Measure logical qubit
        logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent)
        # print(logical_meas)
        if logical_meas != state:
            failed_count += 1



    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    log_e_rate = failed_count/num_runs
    print('total_time_taken_{}_runs_{}'.format(num_runs, total_time_taken))
    # results = {
    #
    #     'total_time_taken_{}_runs'.format(num_runs): total_time_taken,
    #     'probability_set': model_probs,
    #     'subset logical error rate': log_e_rate,
    #     'error subset': eset,
    #     'error_model': e_model,
    # }
    # with open("imp_samp/" + filename + '.json', 'w') as file:
    #     json.dump(results, file)
    return log_e_rate


def subset_weight(params):
    weight = 1
    for tuple in params:  # this loop should run over every type of error ya got
        num_gates, prob_error, num_error = tuple[0],tuple[1], tuple[2]
        prob = comb(num_gates, num_error)*prob_error**num_error*(1-prob_error)**(num_gates-num_error)
        weight *= prob
    return weight

def weighted_logical_error_rate(probs, esets, log_e_dict, gates = [12, 18, 24, 78, 24]):
    sum_log_e = 0
    weights = []
    subset_e_rates = []

    for eset in esets:
        # weight = subset_weight([(gates[0], probs[0], eset[0]),
        #                (gates[1], probs[1], eset[1]),
        #                (gates[2], probs[2], eset[2])])
        weight = subset_weight([(gates[i], probs[i], eset[i]) for i in range(len(eset))])
        log_e_subset = weight * log_e_dict[str(eset)]
        sum_log_e += log_e_subset
        weights.append(weight)
        subset_e_rates.append(log_e_subset)
    return sum_log_e, weights, subset_e_rates

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
#
# filename = "data\manual_cancel_ys_after_index_tracking_2000_p=0.007"
# with open(filename+'.json', 'r') as file:
#     results = json.load(file)
# plot_log_e_rate_graph(filename, results, 20, include_b=True, save=False)



with open("imp_samp\\1000_Dress_imp_samp_reset_leakage_each_run.json", 'r') as file:
    log_e_dict = json.load(file)
with open("imp_samp\\significant_subsets_Dress_pl_2e4.json") as file:
    error_subsets = json.load(file)
pts = 50
sum_log_e = np.zeros(pts)
# sum_log_e = np.zeros((pts, pts))
ps = np.geomspace(0.01, 0.00001, pts)
pds = np.linspace(0.001, 0.00001, pts)
pheat = 1e-3
for pd in [2e-6, 2e-5, 2e-4]:
    for i, p in enumerate(ps):
        probs = (p/10, p/10, p, pd, pheat)
        # for j, pd in enumerate(pds):
        #     probs = (p/10, p/10, p, pd, 0.0001)
        sum_log_e[i] = weighted_logical_error_rate(probs=probs, esets=error_subsets, log_e_dict=log_e_dict)[0]
    plt.loglog(ps, sum_log_e, label='logical_e, pd={}'.format(pd))
# X, Y = np.meshgrid(ps, pds)
# plt.contourf(X, Y, sum_log_e, levels = np.linspace(5e-4, 5e-3, 11))
# plt.xscale('log')
# plt.yscale('log')
# plt.colorbar()
plt.title('pdd pheat={}'.format(pheat))
plt.loglog(ps,ps, label='physical error rate')
plt.legend()
plt.show()