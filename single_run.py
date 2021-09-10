from math import pi
from surface_code import *
from compiled_surface_code import *
import numpy as np
from projectq import MainEngine
from projectq.backends import Simulator, CircuitDrawer,ResourceCounter, CommandPrinter,CircuitDrawerMatplotlib
from projectq.setups.default import get_engine_list
from projectq.ops import All, CNOT, H, Measure, X, Z, Y, Rz, Rx, Ry, Rxx, Deallocate, Command
from projectq.setups import restrictedgateset, trapped_ion_decomposer

g1q = (Rx,Ry,Rz)
g2q = (Rxx(pi/4),Rxx(-pi/4))

# engines=restrictedgateset.get_engine_list(one_qubit_gates=g1q,
#                                                two_qubit_gates=(Rxx,),
#                                           compiler_chooser=trapped_ion_decomposer.chooser_Ry_reducer)
# drawing_engine=CircuitDrawer()#Matplotlib()
resource_counter = ResourceCounter()
#cmd_engine=CommandPrinter()
eng = MainEngine(Simulator())
eng = MainEngine(engine_list=get_engine_list()+[resource_counter])#Simulator(rnd_seed=10))#for code compilation

# config
cycle_number=1
n_1q_errors = 0
n_2q_errors = 0
n_a_meas_errors = 0
display = False
display_results = True


ancilla_number=8
data_number=9
gates_per_func_1q = [32, 24, 24, 32]*cycle_number # not including prep/meas or ancilla meas
gates_per_func_2q = [6,6,6,6]*cycle_number

number_stab_steps_1cycle = 8
number_stab_steps = number_stab_steps_1cycle*cycle_number

correction_table = load_lookup_table('correction_table_depolarising.json')
#generate random locations for your errors (importance sampling)
#errors_per_func_1q=generate_noise_locations(n_1q_errors,number_stab_steps,gates_per_func_1q)
#errors_per_func_2q=generate_noise_locations(n_2q_errors,number_stab_steps,gates_per_func_2q)
#error_locations_a_meas=generate_noise_locations(n_a_meas_errors, cycle_number*ancilla_number,
#                                                ancilla_number*[1])
# if display:
#     print('errors per step 1q {}'.format(errors_per_func_1q))
#     print('errors per step 2q {}'.format(errors_per_func_2q))

#Initialise qubit register
syndrome=[] # list for storing full syndrome (over all stabilizer cycles)
data = eng.allocate_qureg(data_number)
ancilla = eng.allocate_qureg(ancilla_number)

##Perfect logical state prep (well motivated as per
# https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf)
# see supplementray material
#All(H) | data #x basis prep
# quiescent = stabilize_no_noise(eng,ancilla,data)
#print(quiescent)
#X | data[0]
# Z | data[3]
# Y | data[6]
## Error correction cycles
for i in range(cycle_number):
    # if i ==1:
    #     X|data[7] # inject a bitflip error
    syndrome.append(stabilizer_cycle(data, ancilla, eng, reset=True, p1q=0, p2q=0))  # compiled for ion traps
    print(resource_counter)
    #syndrome.append(stabilize_no_noise(eng,ancilla,data,reset=False))
    eng.flush()

### Decode
# a0_increment = decode(quiescent=np.array(quiescent),syndrome=np.array(syndrome),num_cycles=cycle_number,
#              table=correction_table,data=data)
#eng.flush()

#Measure logical qubit
#All(H) | data #x basis meas
All(Measure) | data
eng.flush()  # flush all gates (and execute measurements)
#print(drawing_engine.get_latex())

#print(resource_counter)

#drawing_engine.draw()
#plt.show()
print('full syndrome {}'.format(syndrome))
print('data qubits {}'.format([int(q) for q in data]))
print('logical Z meas result {}'.format(sum([int(q) for q in data])%2))