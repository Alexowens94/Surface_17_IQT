import time
from compiled_surface_code_pdd_error_model import *
import numpy as np
from scipy.optimize import curve_fit
from projectq import MainEngine
from projectq.backends import Simulator
import matplotlib.pyplot as plt

correction_table=load_lookup_table("C:\\Users\\daisy\\Documents\\Code\\Surface_17_IQT\\correction_table_depolarising.json")

def calculate_log_e_rate(num_runs, filename, correction_table, pxctrl=0, pyctrl=0, pxxctrl=0, pmot=0, pd=0):
    """
    Iterates through a number of runs.
    Each run fails after a number of rounds of error correction.
    Uses these numbers to produce a bar graph: round_num vs. number of runs that failed on this round.
    Fits the bar graph to find the logical error rate.
    Saves the bar graph in file_name.
    """
    result_dict = {}
    round_count_list = []
    t_begin = time.perf_counter()
    for _ in range(num_runs):
        eng = MainEngine(Simulator())
        round_count=0
        #max_rounds=600 # TODO: Catch as an exception
        #while round_count < max_rounds:
        while True:
            t_start=time.perf_counter()

            #Initialise qubit register
            data=eng.allocate_qureg(9)
            ancilla=eng.allocate_qureg(8)

            ##Perfect logical state prep (well motivated as per
            # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
            #All(X) | data

            # quiescent prepared perfectly (all errors = 0)
            quiescent=np.array(stabilizer_cycle(data, ancilla, eng, reset=True,
                                                  pyctrl=0, pxctrl=0, pxxctrl=0, pmot=0, pd=0))

            ## Error correction cycle
            prev_syndrome=np.array(quiescent)

            syndrome=np.array(stabilizer_cycle(data, ancilla, eng, reset=True,
                                                 pyctrl=pyctrl, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, pd=pd))

            flips_a=(prev_syndrome - syndrome) % 2
            prev_syndrome=syndrome
            if np.all((flips_a == 0)):
                ft_syndrome=flips_a
            else:
                syndrome=np.array(stabilizer_cycle(data, ancilla, eng, reset=True,
                                                     pyctrl=pyctrl, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, pd=pd))

                flips_b=(prev_syndrome - syndrome) % 2
                prev_syndrome=syndrome
                ft_syndrome=(flips_a + flips_b) %2
            error_vec=lookup(ft_syndrome, correction_table)
            apply_correction(error_vec, data)
            eng.flush()
            round_count += 1
            #Measure logical qubit
            All(Measure) | data
            eng.flush()  # flush all gates (and execute measurements)
            data_meas=[int(q) for q in data]

            logic_Z_meas=sum(data_meas)%2
            if logic_Z_meas == 1: #incorrect
                print('incorrect logic meas {}'.format(data_meas))
                print('round {}'.format(round_count))
                break
        round_count_list.append(round_count)
    t_end=time.perf_counter()
    total_time_taken=t_end - t_begin
    print('total time taken {}'.format(total_time_taken))
    print(round_count_list)

    results={
        'pyctrl' : pyctrl,
        'pxctrl': pxctrl,
        'pxxctrl': pxxctrl,
        'pmot': pmot,
        'pd' : pd,
        'rounds_til_fail_list': round_count_list,
        'total_time_taken_{}_runs'.format(num_runs): total_time_taken
    }
    with open("data/" + filename + '.json', 'w') as file:
        json.dump(results, file)
    log_e_rate = plot_log_e_rate_graph(filename, results, 10)
    return log_e_rate

def plot_log_e_rate_graph(filename, results, num_bins):
    data = results['rounds_til_fail_list']
    vals, bins, patches = plt.hist(data, bins=num_bins)
    popt, pcov = curve_fit(f, bins[:-1], vals)
    plt.ylabel('runs ending after t QEC round')
    plt.xlabel('t')
    plt.plot(bins, f(bins, *popt), 'r-')
    log_e_rate = popt[1]
    title = filename + " log_e_rate = {}, A = {}".format(log_e_rate, popt[0])
    plt.title(title)
    plt.savefig("figs/" + filename + ".png")
    plt.close()
    return log_e_rate

def f(t, A, r):
    return A * np.exp(-r * t)

#calculate_log_e_rate(100, "this_file5")

