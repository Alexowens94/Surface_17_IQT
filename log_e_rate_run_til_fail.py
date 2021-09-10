import time
from surface_code import *
from compiled_surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator

eng = MainEngine(Simulator())
ancilla_number=8
data_number=9
result_dict = {}

run_name = 'test_log_e_rate_til_fail_20_old_order'

number_of_runs = 20
p1q = 0.001
p2q = 0.01
correction_table = load_lookup_table('correction_table_depolarising.json')

t_begin = time.perf_counter()
for _ in range(number_of_runs):
    round_count = 0
    #max_rounds = 2000
    #while round_count < max_rounds:
    while True:

        t_start = time.perf_counter()

        #Initialise qubit register
        data=eng.allocate_qureg(data_number)
        ancilla = eng.allocate_qureg(ancilla_number)

        ##Perfect logical state prep (well motivated as per
        # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
        #All(X) | data
        quiescent = stabilize_no_noise(eng,ancilla,data)

        ## Error correction cycle
        prev_syndrome=np.array(quiescent)

        syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset=True, old_order=True, p1q=p1q, p2q=p2q))
        flips_a = (prev_syndrome - syndrome) % 2
        prev_syndrome = syndrome
        if np.all((flips_a == 0)):
            ft_syndrome = flips_a
        else:
            syndrome = np.array(stabilizer_cycle(data, ancilla, eng, reset=True, old_order=True, p1q=p1q, p2q=p2q))
            flips_b = (prev_syndrome - syndrome) % 2
            prev_syndrome = syndrome
            ft_syndrome = (flips_a + flips_b) %2
        error_vec = lookup(ft_syndrome, correction_table)
        apply_correction(error_vec, data)
        eng.flush()
        round_count += 1

        #Measure logical qubit
        All(Measure) | data
        eng.flush()  # flush all gates (and execute measurements)
        logic_Z_meas = sum([int(q) for q in data])%2
        if logic_Z_meas == 1: #incorrect
            break


    t_stop = time.perf_counter()
    time_taken = t_stop-t_start


    if str(round_count) in result_dict:
        result_dict[str(round_count)] += 1
    else:
        result_dict[str(round_count)] = 1
    #
    # print('{} QEC rounds til failure'.format(round_count))
    # print('time taken {}s'.format(time_taken))
t_end = time.perf_counter()
total_time_taken = t_end - t_begin
print('total time taken {}'.format(total_time_taken))
print(result_dict)

results={
    'p1q_error': p1q,
    'p2q_error': p2q,
    'rounds_til_fail_tally': result_dict,
    'total_time_taken_{}_runs'.format(number_of_runs): total_time_taken
}


with open(run_name+'.txt','w') as file:
     json.dump(results,file)
