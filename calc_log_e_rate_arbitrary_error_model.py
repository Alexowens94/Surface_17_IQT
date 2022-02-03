import time
from compiled_surface_code_arbitrary_error_model import *
import itertools
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


def calculate_log_e_rate(num_runs, filename, correction_table, e_model, e_probs, bitflip_table, num_bins=20,
                         cz_compilation=False):
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
            quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, e_probs, cz_compilation)
            # print('quiescent {}'.format(quiescent))
            # Error correction cycle
            # error_test(round_count, 'X', data)
            prev_syndrome = np.array(quiescent)
            # print('q {}'.format(prev_syndrome))
            if cz_compilation:
                syndrome = np.array(cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))
            else:
                syndrome = np.array(stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))
            flips_a = (prev_syndrome - syndrome) % 2
            # print('fa {}'.format(flips_a))
            prev_syndrome = syndrome
            if np.all((flips_a == 0)):
                ft_syndrome = flips_a
            else:
                if cz_compilation:
                    syndrome = np.array(cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))
                else:
                    syndrome = np.array(stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, e_probs, reset=True))
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
            print('round {}'.format(round_count))

            # #Measure logical qubit
            logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent, leaked_q_reg)
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

def generate_from_distr(probs):
    r = random.random()*sum(probs)
    total=0
    for i,p in enumerate(probs):
        total += p
        if total > r:
            return i

def generate_error_locations_non_uniform(num_e, num_e_locations, probs):
    mask = num_e_locations*[0]
    for _ in range(num_e):
        success = False
        while success == False:
            ind = generate_from_distr(probs)
            # print('ind {}'.format(ind))
            if mask[ind] == 0:
                mask[ind] = 1
                success = True
    return mask

def calculate_log_e_rate_error_subset(num_runs, correction_table, model, eset, locations, bitflip_table,
                                      prob_lists, cz_comp=False):
    """
    """
    failed_count = 0
    t_begin = time.perf_counter()
    for _ in range(num_runs):
        print('run {}'.format(_))
        # in error subset testing over 1 round, so reinitialise leaked_reg within loop
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
        # print('b = {}, s = {}'.format(basis, state))

        # generate random choice of error locations, given num. of each error type (eset) and total possible (locations)
        # print('eset {}'.format(eset))
        if type(prob_lists[0]) == list:

            p_set = [generate_error_locations_non_uniform(eset[i], locations[i], prob_lists[i]) for i in range(len(eset))]
            if _ == 1:
                print('non_uniform')
                print(p_set)
        else:
            # print('uniform')
            p_set = [generate_error_locations(eset[i], locations[i]) for i in range(len(eset))]
        # for p in p_set:
        #     print('error locations {}'.format(len(p)))
        #     print('errors of type {} '.format(sum(p)))
        e_model, model_probs = instantiate_error_model(p_set, model)
        # print(e_model['gates'], model_probs)
        eng = MainEngine(Simulator())

        # Initialise qubit register
        data = eng.allocate_qureg(9)
        ancilla = eng.allocate_qureg(8)

        # prepare logical qubit
        quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, model_probs, cz_comp)

        # Error correction cycle
        # error_test(round_count, 'X', data)
        prev_syndrome = np.array(quiescent)
        # print('real syndrome meas')
        if cz_comp:
            syndrome = cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True)
        else:
            syndrome = stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True)
        syndrome = np.array(syndrome)
        # print('syndrome {}'.format(syndrome))
        flips_a = (prev_syndrome - syndrome) % 2
        error_vec = lookup(flips_a, correction_table)
        # print('error vex {}'.format(error_vec))
        apply_correction(error_vec, data)
        eng.flush()
        # print('round {}'.format(round_count))

        # #Measure logical qubit
        logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent, leaked_q_reg)
        # print(logical_meas)
        if logical_meas != state:
            failed_count += 1
            print('failed in round {}'.format(_))
    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    log_e_rate = failed_count/num_runs
    print(log_e_rate)
    print('total_time_taken_{}_runs_{}'.format(num_runs, total_time_taken))

    return log_e_rate

def s1_calculate_log_e_rate_error_subset(num_runs, correction_table, model, eset, locations, bitflip_table,
                                      prob_lists, cz_comp=False, cooling=False):
    """
    """
    failed_count = 0
    flipsa_not0_count = 0
    t_begin = time.perf_counter()
    for _ in range(num_runs):
        # print('run {}'.format(_))
        # in error subset testing over 1 round, so reinitialise leaked_reg within loop
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
        # print('b = {}, s = {}'.format(basis, state))

        # generate random choice of error locations, given num. of each error type (eset) and total possible (locations)
        p_set = []
        for i in range(len(eset)):
            if len(prob_lists[i]) == 1:
                p_set.append(generate_error_locations(eset[i], locations[i]))
                # if _ == 1:
                    # print('uniform')
            else:
                p_set.append(generate_error_locations_non_uniform(eset[i], locations[i], prob_lists[i]))
                # if _ == 1:
                    # print('non_uniform')
        # print('pset {}'.format(len(p_set)))
        e_model, model_probs = instantiate_error_model(p_set, model)
        # print('e_model_probs {}'.format(model_probs))
        # print(e_model['gates'], model_probs)
        eng = MainEngine(Simulator())

        # Initialise qubit register
        data = eng.allocate_qureg(9)
        ancilla = eng.allocate_qureg(8)

        # prepare logical qubit
        quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, model_probs, cz_comp)

        # Error correction cycle
        # error_test(round_count, 'X', data)
        prev_syndrome = np.array(quiescent)
        # print('prev syndrome {}'.format(prev_syndrome))
        # print('real syndrome meas')
        if cz_comp:
            syndrome = cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True)
        else:
            syndrome = stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True,
                                                    cooling=cooling)
        syndrome = np.array(syndrome)
        # print('syndrome {}'.format(syndrome))
        flips_a = (prev_syndrome - syndrome) % 2
        if np.all((flips_a == 0)):
            error_vec = lookup(flips_a, correction_table)
            # print('error vex {}'.format(error_vec))
            apply_correction(error_vec, data)
            eng.flush()
            # #Measure logical qubit
            logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent, leaked_q_reg)
            # print(logical_meas)
            if logical_meas != state:
                failed_count += 1
                # print('failed in round {}'.format(_))
        else:
            flipsa_not0_count += 1
            All(Measure) | ancilla
            All(Measure) | data
            eng.flush()
    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    try:
        log_e_rate = failed_count/(num_runs-flipsa_not0_count)
    except ZeroDivisionError:
        print('cant assign error rate, 0 runs ending after first stab, placeholder e rate -1')
        log_e_rate = None
    flipsa_not0_rate = flipsa_not0_count/num_runs
    print('counts')
    print('failed {} times'.format(failed_count))
    print('wouldnt stop at s1 {} '.format(flipsa_not0_count))
    print('rates')
    print(log_e_rate)
    print(flipsa_not0_rate)
    print('total_time_taken_{}_runs_{}'.format(num_runs, total_time_taken))

    return log_e_rate, flipsa_not0_rate, failed_count, flipsa_not0_count

def s2_calculate_log_e_rate_error_subset(num_runs, correction_table, model, eset, locations, bitflip_table,
                                      prob_lists, cz_comp=False, cooling=False):
    """
    """
    failed_count = 0
    flipsa_not0_count = 0
    t_begin = time.perf_counter()
    for _ in range(num_runs):
        # print('run {}'.format(_))
        # in error subset testing over 1 round, so reinitialise leaked_reg within loop
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
        # print('b = {}, s = {}'.format(basis, state))

        # generate random choice of error locations, given num. of each error type (eset) and total possible (locations)
        p_set = []
        for i in range(len(eset)):
            if len(prob_lists[i]) == 1:
                p_set.append(generate_error_locations(eset[i], locations[i]))
                # if _ == 1:
                    # print('uniform')
            else:
                p_set.append(generate_error_locations_non_uniform(eset[i], locations[i], prob_lists[i]))
                # if _ == 1:
                    # print('non_uniform')
                #comeback
        e_model, model_probs = instantiate_error_model(p_set, model)
        # print(e_model['gates'], model_probs)
        eng = MainEngine(Simulator())

        # Initialise qubit register
        data = eng.allocate_qureg(9)
        ancilla = eng.allocate_qureg(8)

        # prepare logical qubit
        quiescent = logical_prep(data, basis, state, ancilla, leaked_q_reg, eng, e_model, model_probs, cz_comp)

        # Error correction cycle
        # error_test(round_count, 'X', data)
        prev_syndrome = np.array(quiescent)
        # print('real syndrome meas')
        if cz_comp:  # todo make cz version return indices too
            syndrome = cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True)
        else:
            # syndrome, rxi, ryi, rxxi = stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model,
            #                                                         model_probs, reset=True, return_indices=True)
            syndrome = stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model,
                                                                    model_probs, reset=True, cooling=cooling)

        syndrome = np.array(syndrome)
        # print('syndrome {}'.format(syndrome))
        flips_a = (prev_syndrome - syndrome) % 2
        if not np.all(flips_a == 0):
            flipsa_not0_count += 1
            if cz_comp:
                syndrome = cz_stabilizer_cycle(data, ancilla, leaked_q_reg, eng, e_model, model_probs, reset=True)
            else:
                # print('rx {} ry {} rxx {} indices pre s2'.format(rxi,ryi,rxxi))
                syndrome = stabilizer_cycle_error_index(data, ancilla, leaked_q_reg, eng, e_model, model_probs,
                                                        stab_ind=1, reset=True, cooling=cooling)
            ft_syndrome = np.array(prev_syndrome - syndrome) % 2
            error_vec = lookup(ft_syndrome, correction_table)
            # print('error vex {}'.format(error_vec))
            apply_correction(error_vec, data)
            eng.flush()
            # #Measure logical qubit
            logical_meas = logical_measurement(data, basis, eng, bitflip_table, quiescent, leaked_q_reg)
            # print(logical_meas)
            # print('s2 {}'.format(not np.all(flips_a == 0)))
            if logical_meas != state:
                failed_count += 1
                # print('failed in round {}'.format(_))
                # print('model probs {}'.format(model_probs))
        else:
            All(Measure) | ancilla
            All(Measure) | data
            eng.flush()
    t_end = time.perf_counter()
    total_time_taken = t_end - t_begin
    try:
        log_e_rate = failed_count/(flipsa_not0_count)
    except ZeroDivisionError:
        print('cant assign error rate, 0 runs ending after first stab, placeholder e rate -1')
        log_e_rate = -1
    flipsa_not0_rate = flipsa_not0_count/num_runs
    print('counts')
    print('failed {} times'.format(failed_count))
    print('wouldnt stop at s1 {} '.format(flipsa_not0_count))
    print('rates')
    print(log_e_rate)
    print(flipsa_not0_rate)
    print('total_time_taken_{}_runs_{}'.format(num_runs, total_time_taken))

    return log_e_rate, flipsa_not0_rate, failed_count, flipsa_not0_count


def subset_weight(params):
    weight = 1
    for tuple in params:  # this loop should run over every type of error ya got
        num_gates, prob_error, num_error = tuple[0],tuple[1], tuple[2]
        prob = comb(num_gates, num_error)*prob_error**num_error*(1-prob_error)**(num_gates-num_error)
        weight *= prob
    return weight

def weight_contribution_variable_e_rate(num_gates, prob_error_list, num_error):
    #  difference with subset weight is that prob_error_list has a probability for every error location
    # i.e it should be a list of length num_gates not a single probability as in subset_weight
    sum_prob = 0
    # print('len prob e list {}'.format(len(prob_error_list)))
    for e_locations in itertools.combinations(range(num_gates), num_error):
        prob = 1
        #     print('error locs {}'.format(e_locations))
        no_e_locations = list(range(num_gates))
        # print('len no e {}'.format(len(no_e_locations)))
        for i in e_locations:
            no_e_locations.remove(i)
            prob *= prob_error_list[i]
        for j in no_e_locations:
            prob *= 1 - prob_error_list[j]
        #     print('not error locs {}'.format(no_e_locations))
        #     print(prob)
        sum_prob += prob
    weight = sum_prob
    return weight

def subset_weight_mixed(params):
    weight = 1
    for num_gates, prob_error, num_error in params:
        #if type(prob_error) is int:
        # print(num_gates, prob_error, num_error)
        if len(prob_error) == 1:  # we pass one error_prob if error prob constant (over locations) for this error type
            p = prob_error[0]
            prob = comb(num_gates, num_error) * p ** num_error * (1 - p) ** (num_gates - num_error)
            weight *= prob
        else:
            weight *= weight_contribution_variable_e_rate(num_gates, prob_error, num_error)
    return weight

def subset_weight_variable_e_rate(params):
    ''' calculate the statistical weight of an error subset, allowing for
        different probabilities for the same error type depending on the location
        see https://iopscience.iop.org/article/10.1088/1367-2630/aab341 eq 20-22'''
    weight = 1
    for tuple in params:
        #  difference with subset weight is that prob_error_list has a probability for every error location
        # i.e it should be a list of length num_gates not a single probability as in subset_weight
        num_gates, prob_error_list, num_error = tuple[0], tuple[1], tuple[2]
        sum_prob = 0
        for e_locations in itertools.combinations(range(num_gates), num_error):
            prob = 1
            #     print('error locs {}'.format(e_locations))
            no_e_locations = list(range(num_gates))
            for i in e_locations:
                no_e_locations.remove(i)
                prob *= prob_error_list[i]
            for j in no_e_locations:
                prob *= 1 - prob_error_list[j]
            #     print('not error locs {}'.format(no_e_locations))
            #     print(prob)
            sum_prob += prob
        weight *= sum_prob
    return weight


def weighted_logical_error_rate(probs, esets, log_e_dict, error_locations, variable_e_rate=False):
    sum_log_e = 0
    weights = []
    subset_e_rates = []

    for eset in esets:
        #  recalculate subset weights with probabilities given in argument - to plot log e curves
        if variable_e_rate:
            weight = subset_weight_variable_e_rate([(error_locations[ind], probs[ind], eset[ind]) for ind in range(len(eset))])
        else:
            weight = subset_weight([(error_locations[ind], probs[ind], eset[ind]) for ind in range(len(eset))])
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

def calculate_significant_subsets_incrementally(locations, prob_lists, cutoff_weight=1e-6, display_checked=True):
    results_dict = {}
    error_subset_list = []
    subset_weight_list = []
    already_checked = []
    num_error_types = len(locations)
    esets = [num_error_types*[0]]
    while len(esets) != 0:
        print('esets at start of while loop {}'.format(esets))
        for eset in esets:
            print('checking eset {}'.format(eset))
            esets.remove(eset)  # remove current subset from list to check weights of
            if eset in already_checked:
                print('eset already checked')
                continue  # if already checked weight of this set and added it, skip the weight calc
            params = [(locations[i], prob_lists[i], eset[i]) for i in range(len(eset))]
            eset_weight = subset_weight_mixed(params)  # calculate the weight of the subset: eset
            if eset_weight > cutoff_weight:
                print('eset weight {}'.format(eset_weight))
                print('weight calc parameters {}'.format(params))
                error_subset_list.append(eset)  # if weight is above cutoff, add to significant subsets list
                subset_weight_list.append(eset_weight)
                new_esets = [eset.copy() for j in range(len(eset))]  # make a new eset for each error type
                print('copies {}'.format(new_esets))
                for j in range(len(new_esets)):
                    new_esets[j][j] += 1  # make new error subsets, one for each error type with one additional error
                print('new esets {} '.format(new_esets))
                for new_eset in new_esets:  # add new esets to list of subsets to check weight of
                    esets.append(new_eset)
            already_checked.append(eset)
    if display_checked:
        print('checked')
        print(already_checked)
        print(len(already_checked))
    print('{} sig subsets'.format(len(error_subset_list)))
    results_dict['error_subset_list'] = error_subset_list  # significant subsets
    results_dict['subset_weights'] = subset_weight_list  # statistical weights of the subsets, order matches subset list
    results_dict['cutoff_weight'] = cutoff_weight
    results_dict['locations'] = locations
    results_dict['probabilities'] = prob_lists
    return results_dict



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

def increment(eset_a):
    new_esets = [eset_a.copy() for j in range(len(eset_a))]  # make a new eset for each error type
    for j in range(len(new_esets)):
        new_esets[j][j] += 1  # make new error subsets, one for each error type with one additional error
    return new_esets


def find_significant_subsets_with_more_error(eset_a, error_subsets_list):
    """

    :param eset_a: error subset (list of int)
    :param error_subsets_list: (list of list) list of error subsets
    :return: list of all subsets in error_subsets_list that can be produced by adding errors to eset_a
            returns the list of such subsets, and a list of the difference between eset_a and the list returned
    """

    s2_esets = []
    s2_subtracted = []
    # print('eset a {}'.format(eset_a))
    new_esets = increment(eset_a)
    new_esets.append(eset_a)  # this handles the case you have eset errors in s1, still go to s2, and have 0 in s2
    while len(new_esets) != 0:
        # print('new esets {} '.format(new_esets))
        for new_eset in new_esets:
            # print('checking {}'.format(new_eset))
            new_esets.remove(new_eset)  # remove current item from list to check
            if new_eset in s2_esets:
                # print('already in s2 esets')
                continue
            if new_eset in error_subsets_list:
                subtracted = new_eset.copy()
                for i in range(len(new_eset)):
                    subtracted[i] -= eset_a[i]
                s2_esets.append(new_eset)
                s2_subtracted.append(subtracted)
                next_new_esets = increment(new_eset)
                for ele in next_new_esets:
                    new_esets.append(ele)
    # print('s2 esets {}'.format(s2_esets))
    # print('s2 esets - eset_a {}'.format(s2_subtracted))
    return s2_esets, s2_subtracted


def prob_s2_AND_eset_a_total_errors(eset_a, s2_rate_dict, significant_subsets_s1,
                                    significant_subsets_s1_and_s2, locations1round, prob_lists_s1, prob_lists_s2):
    sump = 0
    ratelist, sw1list, sw2list = [], [], []
    for eset_s1 in significant_subsets_s1:
        s2, s2_subtracted = find_significant_subsets_with_more_error(eset_s1, significant_subsets_s1_and_s2)
        if eset_a in s2:
            ind = s2.index(eset_a)
            eset_subtracted = s2_subtracted[ind]
            print('eset_s1, eset_s2 {} {}'.format(eset_s1, eset_subtracted))
            params1 = [(locations1round[i], prob_lists_s1[i], eset_s1[i]) for i in range(len(eset_s1))]
            params2 = [(locations1round[i], prob_lists_s2[i], eset_subtracted[i]) for i in range(len(eset_subtracted))]
            key_s1 = str(eset_s1)
            sw1 = subset_weight_mixed(params1)
            sw2 = subset_weight_mixed(params2)
            sump += s2_rate_dict[key_s1]*sw1*sw2
            sw1list.append(sw1)
            sw2list.append(sw2)
            ratelist.append(s2_rate_dict[key_s1])
    return sump, ratelist, sw1list, sw2list


def s2sim_calculate_significant_subsets_incrementally(locations, prob_lists_s1, prob_lists_s2, s2_rate_dict, significant_subsets_s1,
                                                      significant_subsets_s1_and_s2,
                                                      cutoff_weight=1e-6, display_checked=True):
    results_dict = {}
    error_subset_list = []
    subset_weight_list = []

    already_checked = []
    num_error_types = len(locations)
    eset1 = num_error_types*[0]
    esets = [eset1.copy() for j in range(len(eset1))] # make a new eset for each error type
    for j in range(len(esets)):
        esets[j][j] += 1  # can't seed with just all 0's as that will stop algo (0 likelihood of seeing s2)
    esets.append(eset1)  # adding all 0 case in (would highlight if a bad compilation can go to s2 with no errors)
    while len(esets) != 0:
        # print('esets at start of while loop {}'.format(esets))
        for eset in esets:
            print('checking eset {}'.format(eset))
            esets.remove(eset)  # remove current subset from list to check weights of
            if eset in already_checked:
                # print('eset already checked')
                continue  # if already checked weight of this set and added it, skip the weight calc
            weight, s2rates, s1w, s2w = prob_s2_AND_eset_a_total_errors(eset, s2_rate_dict, significant_subsets_s1,
                                                                        significant_subsets_s1_and_s2,
                                                                        locations, prob_lists_s1, prob_lists_s2)
            print('s2rate list {} eset1 weights {} eset2 weights {}'.format(s2rates, s1w, s2w))
            if weight > cutoff_weight:
                # print('weight {} eset {}'.format(weight,eset))
                error_subset_list.append(eset)  # if weight is above cutoff, add to significant subsets list
                subset_weight_list.append(weight)
                new_esets = [eset.copy() for j in range(len(eset))]  # make a new eset for each error type
                # print('copies {}'.format(new_esets))
                for j in range(len(new_esets)):
                    new_esets[j][j] += 1  # make new error subsets, one for each error type with one additional error
                # print('new esets {} '.format(new_esets))
                for new_eset in new_esets:  # add new esets to list of subsets to check weight of
                    esets.append(new_eset)
            already_checked.append(eset)
    if display_checked:
        print('checked')
        print(already_checked)
        print(len(already_checked))
    results_dict['error_subset_list'] = error_subset_list  # significant subsets
    results_dict['subset_weights'] = subset_weight_list  # statistical weights of the subsets, order matches subset list
    results_dict['cutoff_weight'] = cutoff_weight
    results_dict['locations'] = locations
    results_dict['error probabilities s1'] = prob_lists_s1
    results_dict['error probabilities s2'] = prob_lists_s2
    # results_dict['s2rates'] = s2rates
    # results_dict['prob k in subset 1'] = s1w
    # results_dict['prob n-k in subset 2'] = s2w    # results_dict['s2rates'] = s2rates
    # results_dict['prob k in subset 1'] = s1w
    # results_dict['prob n-k in subset 2'] = s2w
    return results_dict


def f(t, A, r):
    return A * np.exp(-r * t)
#
# filename = "data\manual_cancel_ys_after_index_tracking_2000_p=0.007"
# with open(filename+'.json', 'r') as file:
#     results = json.load(file)
# plot_log_e_rate_graph(filename, results, 20, include_b=True, save=False)

