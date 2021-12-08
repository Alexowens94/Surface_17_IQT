import time
from compiled_surface_code_pdd_error_model import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator
eng = MainEngine(Simulator())
ancilla_number=8
data_number=9
result_dict = {}
round_count_list=[]

run_name = 'test_log_e_rate_til_fail_pdd_new_syndrome_rules'

number_of_runs = 200
pxctrl=0.5
pyctrl=0
pxxctrl=0
pmot=0
pd=0
correction_table = load_lookup_table("C:\\Users\\daisy\\Documents\\Code\\Surface_17_IQT\\correction_table_depolarising.json")

t_begin = time.perf_counter()

# Initialise qubit register
data = eng.allocate_qureg(data_number)
ancilla = eng.allocate_qureg(ancilla_number)

##Perfect logical state prep (well motivated as per
# https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
# All(X) | data


for _ in range(number_of_runs):
    # For a run, initialize in quiescent state
    quiescent = np.array(stabilizer_cycle(data, ancilla, eng, reset=True,
                                          pyctrl=0, pxctrl=0, pxxctrl=0, pmot=0, pd=0))
    ft_syndrome = [00000000]  # At this point we believe there are no errors to correct
    prev_syndrome = np.array(quiescent)

    round_count = 0
    # max_rounds = 10
    # while round_count < max_rounds:
    while True:
        #print(round_count)
        t_start = time.perf_counter()


        ## Error correction cycle
        syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset = True,
                                             pyctrl = pyctrl, pxctrl = pxctrl, pxxctrl = pxxctrl, pmot = pmot, pd = pd))

        flips_a = (prev_syndrome - syndrome) % 2
        prev_syndrome = syndrome
        syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset=True,
                                             pyctrl=pyctrl, pxctrl=pxctrl, pxxctrl=pxxctrl, pmot=pmot, pd=pd))
        flips_b = (prev_syndrome - syndrome) % 2

        ft_syndrome = (flips_a + flips_b) % 2 # Daisy convince self of logic of this
        error_vec = lookup(ft_syndrome, correction_table)
        apply_correction(error_vec, data)
        eng.flush()
        round_count += 1

        # Measure logical qubit to check you still have 0, otherwise error correction has failed
        All(Measure) | data
        eng.flush()  # flush all gates (and execute measurements)
        data_meas = [int(q) for q in data] # measurements of the data qubits
        #print("data_meas: {}".format(data_meas))
        logic_Z_meas = sum(data_meas)%2 # logical qubit measurement (Z projection) = sum % 2 of the data qubits
        #print("logic_Z_meas: {}".format(data_meas))
        if logic_Z_meas == 1: #incorrect, the logical qubit was measured as 1 which is incorrect
            print('incorrect logic meas {}'.format(data_meas))
            print('round {}'.format(round_count))
            break

    t_stop = time.perf_counter()
    time_taken = t_stop-t_start

    round_count_list.append(round_count)


t_end = time.perf_counter()
total_time_taken = t_end - t_begin
print('total time taken {}'.format(total_time_taken))
print(round_count_list)

results={
    'pyctrl' : pyctrl,
    'pxctrl': pxctrl,
    'pxxctrl': pxxctrl,
    'pmot': pmot,
    'pd' : pd,
    #'rounds_til_fail_tally': result_dict,
    'rounds_til_fail_list': round_count_list,
    'total_time_taken_{}_runs'.format(number_of_runs): total_time_taken
}


with open(run_name+'.txt','w') as file:
     json.dump(results,file)
