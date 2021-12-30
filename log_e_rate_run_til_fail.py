import time
from surface_code import *
from compiled_surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator

eng = MainEngine(Simulator())
ancilla_number = 8
data_number = 9
round_count_list = []

run_name = 'test_log_e_rate_til_fail_p20_ancilla_data_repumped'

number_of_runs = 2
p1q = 0.001
p2q = 0
p_leak = 0.00
bad_hook_order = False  # set true to take stab meas in wrong order, single qubit errors mid cycle can give hook errors

correction_table = load_lookup_table('correction_table_depolarising.json')

t_begin = time.perf_counter()
for _ in range(number_of_runs):
    leaked_q_reg = 17 * [0]  # register (init as 0s) to track which qubits  leaked, LEAKED == 1 0-8 data,9-16 ancilla
    round_count = 0
    #max_rounds = 60
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

        syndrome = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, reset=True,
                                             bad_hook_order=bad_hook_order, p1q=p1q, p2q=p2q, p_leak=p_leak))
        flips_a = (prev_syndrome - syndrome) % 2
        prev_syndrome = syndrome
        if np.all((flips_a == 0)):
            ft_syndrome = flips_a
        else:
            syndrome = np.array(stabilizer_cycle(data, ancilla, leaked_q_reg, eng, reset=True,
                                                 bad_hook_order=bad_hook_order, p1q=p1q, p2q=p2q, p_leak=p_leak))
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
        data_meas = [int(q) for q in data]
        # print('leaked q reg, data meas')
        # print(leaked_q_reg)
        # print(data_meas)
        ### measure leaked states as bright
        for i, q in enumerate(leaked_q_reg[:9]):
            if int(q) == 1:
                data_meas[i] = 1
                #print('leak registered at data meas')
        logic_Z_meas = sum(data_meas)%2
        if logic_Z_meas == 1: #incorrect
            print('incorrect logic meas {}'.format(data_meas))
            break


    t_stop = time.perf_counter()
    time_taken = t_stop-t_start
    round_count_list.append(round_count)

t_end = time.perf_counter()
total_time_taken = t_end - t_begin
print('total time taken {}'.format(total_time_taken))
print(round_count_list)

results={
    'p1q_error': p1q,
    'p2q_error': p2q,
    'p_leak_error': p_leak,
    'rounds_til_fail_list': round_count_list,
    'total_time_taken_{}_runs'.format(number_of_runs): total_time_taken
}


with open(run_name+'.txt','w') as file:
     json.dump(results,file)
