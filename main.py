import time
from surface_code import *
from compiled_surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator, CommandPrinter,CircuitDrawerMatplotlib#,QEC_Simulator

if __name__ == '__main__':
    #drawing_engine=CircuitDrawerMatplotlib()
    #cmd_engine=CommandPrinter()
    eng = MainEngine(Simulator(rnd_seed=10))

    #eng = MainEngine(engine_list=get_engine_list()+[cmd_engine])#Simulator(rnd_seed=10))#for code compilation
    cycle_number=2
    ancilla_number=8
    data_number=9
    n_1q_errors = 1
    n_2q_errors = 0
    number_of_runs = 100
    display = True
    a0_count = 0
    incorrect_count =0
    correct_count =0
    a0_incorrect_count = 0
    a0_correct_count = 0
    a_not0_incorrect_count = 0
    a_not0_correct_count = 0

    correction_table = load_lookup_table('correction_table_depolarising.json')

    gates_per_func_1q = [4,0,0,4]*cycle_number #not including prep/meas or ancilla meas
    gates_per_func_2q = [6,6,6,6]*cycle_number
    number_stab_steps_1cycle = 4
    number_stab_steps = 4*cycle_number

    t_start=time.perf_counter()
    for i in range(number_of_runs):

        #generate random locations for your errors (importance sampling)
        errors_per_func_1q=generate_noise_locations(n_1q_errors,number_stab_steps,gates_per_func_1q)
        errors_per_func_2q=generate_noise_locations(n_2q_errors,number_stab_steps,gates_per_func_2q)

        if display:
            print('errors per step 1q {}'.format(errors_per_func_1q))
            print('errors per step 2q {}'.format(errors_per_func_2q))

        #Initialise qubit register
        syndrome=[] # list for storing full syndrome (over all stabilizer cycles)
        data=eng.allocate_qureg(data_number)
        ancilla = eng.allocate_qureg(ancilla_number)

        ##Perfect logical state prep (well motivated as per
        # https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
        #All(X) | data
        quiescent = stabilize_no_noise(eng,ancilla,data)

        ## Error correction cycles
        for i in range(cycle_number):
            # if i ==1:
            #     X|data[7] # inject a bitflip error
            syndrome.append(stabilizer_cycle(data, ancilla, eng, reset=True, p1q=0.001, p2q=0.01)) # compiled for ion traps)

        ### Decode
        #print(syndrome)
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




        #drawing_engine.draw()
        #plt.show()

        # print('full syndrome {}'.format(syndrome))
        # print('data qubits {}'.format([int(q) for q in data]))
        # print('logical Z meas result {}'.format(logic_Z_meas))

    t_stop = time.perf_counter()
    print('time taken for {} loops: {}'.format(number_of_runs,t_stop-t_start))
    not_a0_count=number_of_runs-a0_count
    print('number of runs where we shouldnt stop after first syndrome {}'.format(not_a0_count))
    print('total runs where we stop after first syndrome measurement: {}'.format(a0_count))
    print('runs stopping after first syndrome with logical error: {}'.format(a0_incorrect_count))
    print('runs stopping after first syndrome WITHOUT logical error: {}'.format(a0_correct_count))
    print('runs NOT stopping after first syndrome with logical error: {}'.format(a_not0_incorrect_count))
    print('runs NOT stopping after first syndrome WITHOUT logical error: {}'.format(a_not0_correct_count))
    print('failure rate, runs stopping after first syndrome extraction: {}'.format(a0_incorrect_count/a0_count))
    print('failure rate, runs NOT stopping after 1 syn. extraction: {}'.format(a_not0_incorrect_count /not_a0_count))
    print('fraction of runs stopping after first syndrome extraction: {}'.format(a0_count/number_of_runs))

