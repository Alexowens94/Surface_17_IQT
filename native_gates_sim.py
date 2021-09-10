import time
import itertools
from surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator  # CommandPrinter,CircuitDrawerMatplotlib#,QEC_Simulator
import json
from projectq.ops import All, CNOT, H, Measure, X, Z, Y, Rx, Ry, Rxx, Deallocate, Command


eng = MainEngine(Simulator(rnd_seed=10))

# config
cycle_number=2
n_1q_errors = 1
n_2q_errors = 0
n_a_meas_errors = 0
number_of_runs = 100
display = False
display_results = True
parameters=[[1,2],[0,1,2],[0,1,2,3,4,5]]

parameters=[[1,2],[0],[0]]
param_combinations=list(itertools.product(*parameters))
print('number of parameter combinations to run: {}'.format(len(param_combinations)))
results_list=[]
for (cycle_number,n_1q_errors,n_2q_errors) in param_combinations:
    print('cycles={},s={},t={}'.format(cycle_number,n_1q_errors,n_2q_errors))
    ancilla_number=8
    data_number=9
    gates_per_func_1q = [32, 24, 24, 32]*cycle_number # not including prep/meas or ancilla meas
    gates_per_func_2q = [6,6,6,6]*cycle_number

    number_stab_steps_1cycle = 4
    number_stab_steps = number_stab_steps_1cycle*cycle_number

    # counters
    a0_count = 0 # number of times the first measured syndrome matches quiescent state
                 # in this case in reality we would end the EC round after 1 cycle,
                 # and consider no correction is necessary
                 # however for importance sampling noise insertion need fix circuit depth beforehand
                 # see 5.1.1 of Colin J Trout et al 2018 New J. Phys. 20 043038
    incorrect_count =0  # number of runs the logical measurement doesn't match intended result
    correct_count =0  # number of runs logical measurement does match intended result
    a0_incorrect_count = 0 # number of runs where logic meas. does NOT match AND should end after 1 EC cycle
    a0_correct_count = 0 # number of runs where logic meas. does match AND should end after 1 EC cycle
    a_not0_incorrect_count = 0 # number of runs where logic meas. does NOT match AND should NOT end after 1 EC cycle
    a_not0_correct_count = 0 # number of runs where logic meas. does match AND should NOT end after 1 EC cycle

    correction_table = load_lookup_table('correction_table_depolarising.json')



    t_start=time.perf_counter()
    for _ in range(number_of_runs):

        #generate random locations for your errors (importance sampling)
        errors_per_func_1q=generate_noise_locations(n_1q_errors,number_stab_steps,gates_per_func_1q)
        errors_per_func_2q=generate_noise_locations(n_2q_errors,number_stab_steps,gates_per_func_2q)
        error_locations_a_meas=generate_noise_locations(n_a_meas_errors, cycle_number*ancilla_number,
                                                        ancilla_number*[1])
        if display:
            print('errors per step 1q {}'.format(errors_per_func_1q))
            print('errors per step 2q {}'.format(errors_per_func_2q))

        #Initialise qubit register
        syndrome=[] # list for storing full syndrome (over all stabilizer cycles)
        data=eng.allocate_qureg(data_number)
        ancilla = eng.allocate_qureg(ancilla_number)

        ##Perfect logical state prep (well motivated as per
        # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
        # see supplementray material
        #All(X) | data
        quiescent = stabilize_no_noise(eng,ancilla,data)

        ## Error correction cycles
        for i in range(cycle_number):
            if i ==0:
                insert_1q_error(data[1])
                insert_1q_error(data[2])
                insert_1q_error(data[0])
            # syndrome.append(stabilize(eng,ancilla,data,errors_per_func_1q,errors_per_func_2q,
            #                           gates_per_func_1q,gates_per_func_2q,cycle=i))

            syndrome.append(stabilize(eng, ancilla, data, errors_per_func_1q, errors_per_func_2q,
                                      gates_per_func_1q, gates_per_func_2q, cycle=i,gateset='native'))

        ### Decode
        a0_increment = decode(quiescent=np.array(quiescent),syndrome=np.array(syndrome),num_cycles=cycle_number,
                     table=correction_table,data=data)
        eng.flush()

        #Measure logical qubit
        All(Measure) | data
        eng.flush()  # flush all gates (and execute measurements)

        #Assess logical error rate for this subset of errors
        a0_count+=a0_increment
        logic_Z_meas = sum([int(q) for q in data])%2
        if logic_Z_meas == 1: #incorrect
            incorrect_count+=1
            if a0_increment == 0: #first syndrome warranted a correction
                a_not0_incorrect_count+=1
            if a0_increment == 1: #should have stopped after 1 syndrome measurement
                a0_incorrect_count+=1
        if logic_Z_meas == 0: #correct
            correct_count+=1
            if a0_increment == 0: #first syndrome warranted a correction
                a_not0_correct_count+=1
            if a0_increment == 1: #should have stopped after 1 syndrome measurement
                a0_correct_count+=1

    t_stop = time.perf_counter()
    print('time taken for {} loops: {}'.format(number_of_runs,t_stop-t_start))
    not_a0_count=number_of_runs-a0_count

    results=get_results(cycle_number, gates_per_func_1q, n_1q_errors, gates_per_func_2q, n_2q_errors,
                    number_of_runs, a0_count, a0_incorrect_count, a_not0_incorrect_count)
    results_list.append(results)



    if display_results:
        for key in results.keys():
            print('{} = {}'.format(key,results[key]))

with open('results/native_gates_no_ancilla_meas_error_runs_{}.json'.format(number_of_runs),'w') as file:
    json.dump(results_list,file)
    # print('number of runs where we shouldnt stop after first syndrome {}'.format(not_a0_count))
    # print('total runs where we stop after first syndrome measurement: {}'.format(a0_count))
    # print('runs stopping after first syndrome with logical error: {}'.format(a0_incorrect_count))
    # print('runs stopping after first syndrome WITHOUT logical error: {}'.format(a0_correct_count))
    # print('runs NOT stopping after first syndrome with logical error: {}'.format(a_not0_incorrect_count))
    # print('runs NOT stopping after first syndrome WITHOUT logical error: {}'.format(a_not0_correct_count))
    # print('failure rate, runs stopping after first syndrome extraction: {}'.format(a0_incorrect_count/a0_count))
    # print('failure rate, runs NOT stopping after 1 syn. extraction: {}'.format(a_not0_incorrect_count /not_a0_count))
    # print('fraction of runs stopping after first syndrome extraction: {}'.format(a0_count/number_of_runs))

