import time
from surface_code import *
from compiled_surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator

eng = MainEngine(Simulator())
ancilla_number=8
data_number=9
result_list = []

run_name = 'test_log_e_rate'
number_of_runs = 100
correction_table = load_lookup_table('correction_table_depolarising.json')

# parameters=[[1,2],[0],[0]]
# param_combinations=list(itertools.product(*parameters))
param_combinations=[1,5,10,15,20]
print('number of parameter combinations to run: {}'.format(len(param_combinations)))
for rounds in param_combinations:
    incorrect_count = 0
    t_start=time.perf_counter()
    for i in range(number_of_runs):

        #Initialise qubit register
        data=eng.allocate_qureg(data_number)
        ancilla = eng.allocate_qureg(ancilla_number)

        ##Perfect logical state prep (well motivated as per
        # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
        #All(X) | data
        quiescent = stabilize_no_noise(eng,ancilla,data)

        ## Error correction cycles
        prev_syndrome=np.array(quiescent)
        for i in range(rounds):
            syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset=True, p1q=0.001, p2q=0.01))
            flips_a = (prev_syndrome - syndrome) % 2
            prev_syndrome = syndrome
            if np.all((flips_a == 0)):
                ft_syndrome = flips_a
            else:
                syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset=True, p1q=0.001, p2q=0.01))
                flips_b = (prev_syndrome - syndrome) % 2
                prev_syndrome = syndrome
                ft_syndrome = (flips_a + flips_b) %2
            error_vec = lookup(ft_syndrome, correction_table)
            apply_correction(error_vec, data)
            eng.flush()

        #Measure logical qubit
        All(Measure) | data
        eng.flush()  # flush all gates (and execute measurements)
        logic_Z_meas = sum([int(q) for q in data])%2
        if logic_Z_meas == 1: #incorrect
            incorrect_count+=1


    t_stop = time.perf_counter()
    time_taken = t_stop-t_start
    result = get_results_log_e(rounds, number_of_runs, incorrect_count, time_taken)
    result_list.append(result)
    print('{} QEC rounds'.format(rounds))
    print('time taken for {} loops: {}'.format(number_of_runs,time_taken))

    print('incorrect count = {}'.format(incorrect_count))

with open(run_name+'.txt','w') as file:
    json.dump(result_list,file)
