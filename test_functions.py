# old tests - I leave it here incase it's useful in rewriting a test suite for the current version

# from .surface_code import *
# import pytest
# import numpy as np
# from projectq import MainEngine
# from projectq.backends import Simulator, CommandPrinter,CircuitDrawerMatplotlib#,QEC_Simulator
# from projectq.ops import All, CNOT, H, Measure, X, Z, Y, Rx, Ry, Rxx, Deallocate, Command
#
#
#
# @pytest.mark.parametrize(
#     "input0,input1,input2,expected",
#     [(10, 5, 5*[2], 10), (5, 10, 10*[2],5), (5, 5, 5*[1], 5)]
# )
# def test_partition_total_errors(input0,input1,input2, expected):
#     part=generate_noise_locations(input0,input1,input2)
#     print(part)
#     assert sum(generate_noise_locations(input0,input1,input2)) == expected
#
# @pytest.mark.parametrize(
#     "input0,input1,input2,expected",
#     [(10, 5, 5*[2], 5), (5, 10, 10*[2], 10), (10, 11, 11*[1], 11)])
# def test_partition_parts(input0,input1,input2, expected):
#     assert len(generate_noise_locations(input0,input1,input2)) == expected
#
# @pytest.mark.parametrize(
#     "input0,input1,input2,expected",
#     [(10, 5, 5*[2], 3), (5, 10, 10*[2], 3), (10, 15, 15*[1], 2)])
# def test_partition_max(input0,input1,input2, expected):
#     a_part=np.array(generate_noise_locations(input0,input1,input2))
#     print(a_part)
#     np.testing.assert_array_less(a_part,expected)
#
# def test_0error_syndrome():
#     '''run some cycles, see the syndrome doesnt change in absence of error'''
#     eng = MainEngine(Simulator())
#     cycle_number = 10
#     ancilla_number = 8
#     data_number = 9
#     number_stab_steps_1cycle = 4
#     number_stab_steps = number_stab_steps_1cycle * cycle_number
#     gates_per_func_1q = [4, 0, 0, 4] * cycle_number
#     gates_per_func_2q = [6, 6, 6, 6] * cycle_number
#
#     n_1q_errors = 0  # don't insert any noise
#     n_2q_errors = 0
#
#     errors_per_func_1q = generate_noise_locations(n_1q_errors, number_stab_steps, gates_per_func_1q)
#     errors_per_func_2q = generate_noise_locations(n_2q_errors, number_stab_steps, gates_per_func_2q)
#
#     data = eng.allocate_qureg(data_number)
#     ancilla = eng.allocate_qureg(ancilla_number)
#     syndrome = []
#     for i in range(cycle_number):
#         syndrome.append(stabilize(eng, ancilla, data, errors_per_func_1q, errors_per_func_2q,
#                                   gates_per_func_1q, gates_per_func_2q, cycle=i,gateset='circuit'))
#     All(Measure) | data
#     eng.flush()
#
#
#     original_syndrome = [ syndrome[0] for i in range(cycle_number-1)] #create copies of first syndrome
#     syndromes = np.array(syndrome[1:])  # subsequently measured syndromes
#     np.testing.assert_equal(syndromes,original_syndrome)
# def test_logical0():
#     eng = MainEngine(Simulator())
#     cycle_number = 3
#     ancilla_number = 8
#     data_number = 9
#     number_stab_steps_1cycle = 4
#     number_stab_steps = number_stab_steps_1cycle * cycle_number
#     gates_per_func_1q = [4, 0, 0, 4] * cycle_number
#     gates_per_func_2q = [6, 6, 6, 6] * cycle_number
#
#     n_1q_errors = 0 #don't insert any noise
#     n_2q_errors = 0
#
#     errors_per_func_1q = generate_noise_locations(n_1q_errors, number_stab_steps, gates_per_func_1q)
#     errors_per_func_2q = generate_noise_locations(n_2q_errors, number_stab_steps, gates_per_func_2q)
#
#     data = eng.allocate_qureg(data_number)
#     ancilla = eng.allocate_qureg(ancilla_number)
#     syndrome=[]
#     for i in range(cycle_number):
#         syndrome.append(stabilize(eng, ancilla, data, errors_per_func_1q, errors_per_func_2q,
#                                   gates_per_func_1q, gates_per_func_2q, cycle=i,gateset='circuit'))
#     All(Measure) | data
#     eng.flush()
#
#     logical_state = sum([int(q) for q in data])  %2
#     assert logical_state == 0
#
# def test_logical1():
#     eng = MainEngine(Simulator())
#     cycle_number = 3
#     ancilla_number = 8
#     data_number = 9
#     number_stab_steps_1cycle = 4
#     number_stab_steps = number_stab_steps_1cycle * cycle_number
#     gates_per_func_1q = [4, 0, 0, 4] * cycle_number
#     gates_per_func_2q = [6, 6, 6, 6] * cycle_number
#
#     n_1q_errors = 0 #don't insert any noise
#     n_2q_errors = 0
#
#     errors_per_func_1q = generate_noise_locations(n_1q_errors, number_stab_steps, gates_per_func_1q)
#     errors_per_func_2q = generate_noise_locations(n_2q_errors, number_stab_steps, gates_per_func_2q)
#
#     data = eng.allocate_qureg(data_number)
#     ancilla = eng.allocate_qureg(ancilla_number)
#     syndrome=[]
#
#     for i in range(cycle_number):
#         syndrome.append(stabilize(eng, ancilla, data, errors_per_func_1q, errors_per_func_2q,
#                                   gates_per_func_1q, gates_per_func_2q, cycle=i,gateset='circuit'))
#     X | data[0]  #perform a logical X gate
#     X | data[1]
#     X | data[2]
#     All(Measure) | data
#     eng.flush()
#
#     logical_state = sum([int(q) for q in data])  %2
#     print(sum([int(q) for q in data]))
#     assert logical_state == 1
#
# def test_hadamard_decomp():
#     eng = MainEngine()
#     qubits = eng.allocate_qureg(8)
#     X | qubits[0]
#     X | qubits[1]
#     X | qubits[2]
#     X | qubits[3]
#     for q in qubits:
#         native_noisy_hadamard(q, error_locations_1q=[0,0])
#     All(H) | qubits
#     All(Measure) | qubits
#     eng.flush()
#     result = [int(q) for q in qubits]
#     assert result == [1, 1, 1, 1, 0, 0, 0, 0]
# def test_native_noisy_ancilla_to_x():
#     eng = MainEngine()
#     qbits=eng.allocate_qureg(8)
#     X | qbits[1]
#     X | qbits[3]
#     native_noisy_ancilla_to_x_basis(qbits,[0]*8)
#
#     H | qbits[1]
#     H | qbits[3]
#     H | qbits[4]
#     H | qbits[6]
#
#     All(Measure) | qbits
#     eng.flush()
#
#     assert [int(q) for q in qbits] == [0,1,0,1,0,0,0,0]
#
# @pytest.mark.parametrize(
#     "s,v,expected",
#     [(1, 1, [0,0,0,1,1,0,1,1]), (1, -1, [0,0,0,1,1,0,1,1]),
#      (-1, 1, [0,0,0,1,1,0,1,1]), (-1, -1, [0,0,0,1,1,0,1,1])]
# )
# def test_native_CNOT(s,v,expected):
#     eng = MainEngine()
#     qbits=eng.allocate_qureg(8)
#     #I | qbits[0]
#     # I | qbits[1]
#
#     # I | qbits[2]
#     X | qbits[3]
#
#     X | qbits[4]
#     # I | qbits[5]
#
#     X | qbits[6]
#     X | qbits[7]
#
#     nat_CNOT(qbits[0], qbits[1],s=s,v=v)
#     CNOT | (qbits[0], qbits[1])
#     nat_CNOT(qbits[2], qbits[3],s=s,v=v)
#     CNOT | (qbits[2], qbits[3])
#     nat_CNOT(qbits[4], qbits[5],s=s,v=v)
#     CNOT | (qbits[4], qbits[5])
#     nat_CNOT(qbits[6], qbits[7],s=s,v=v)
#     CNOT | (qbits[6], qbits[7])
#
#     All(Measure) | qbits
#     eng.flush()
#     assert [int(q) for q in qbits] == expected
#
# @pytest.mark.parametrize(
#     "s,v,expected",
#     [(1, 1, [0,0,0,1,1,0,1,1]), (1, -1, [0,0,0,1,1,0,1,1]),
#      (-1, 1, [0,0,0,1,1,0,1,1]), (-1, -1, [0,0,0,1,1,0,1,1])]
# )
# def test_native_noisy_CNOT(s,v,expected):
#     eng = MainEngine()
#     qbits=eng.allocate_qureg(8)
#     #I | qbits[0]
#     # I | qbits[1]
#
#     # I | qbits[2]
#     X | qbits[3]
#
#     X | qbits[4]
#     # I | qbits[5]
#
#     X | qbits[6]
#     X | qbits[7]
#
#     native_noisy_CNOT(qbits[0], qbits[1],s=s,v=v)
#     CNOT | (qbits[0], qbits[1])
#     native_noisy_CNOT(qbits[2], qbits[3],s=s,v=v)
#     CNOT | (qbits[2], qbits[3])
#     native_noisy_CNOT(qbits[4], qbits[5],s=s,v=v)
#     CNOT | (qbits[4], qbits[5])
#     native_noisy_CNOT(qbits[6], qbits[7],s=s,v=v)
#     CNOT | (qbits[6], qbits[7])
#
#     All(Measure) | qbits
#     eng.flush()
#     assert [int(q) for q in qbits] == expected
#
# def test_native_entangle1():
#     eng = MainEngine()
#     data=eng.allocate_qureg(9)
#     ancilla=eng.allocate_qureg(8)
#     d_indices = [1,3,5,2,4,8]
#     a_indices = [1,4,6,2,5,7]
#     native_noisy_entangle(data,d_indices,ancilla,a_indices,error_locations_1q=24*[0],
#                           error_locations_2q=6*[0])
#     #entangle
#     CNOT | (ancilla[1],data[1])
#     CNOT | (ancilla[4],data[3])
#     CNOT | (ancilla[6],data[5])
#
#
#     CNOT | (data[2],ancilla[2])
#     CNOT | (data[4],ancilla[5])
#     CNOT | (data[8],ancilla[7])
#
#     All(Measure) | data
#     All(Measure) | ancilla
#
#     assert [int(d) for d in data]+[int(a) for a in ancilla] == 17*[0]
#
# # def test_insert_1q_error():
# #     eng = MainEngine()
# #     qubit = eng.allocate_qubit()
# #     insert_1q_error(qubit)
#
