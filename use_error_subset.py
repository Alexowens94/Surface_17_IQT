from compiled_surface_code_error_subsets import *
from generate_error_subset import *
from projectq import MainEngine
from projectq.backends import Simulator

# Create qubits
eng = MainEngine(Simulator())
data = eng.allocate_qureg(9)
ancilla = eng.allocate_qureg(8)

# Generate an error subset (a list of errors and their corresponding random locations)
# error_number_dictionary = {"XCtrlError" : 5, "YCtrlError" : 6, "XXCtrl" : 1}
errors_to_generate = {"XCtrlError" : 5, "DephasingError" : 6}

# Create an error subset
error_subset = generate_error_subset(errors_to_generate)
# x_error_1 = XCtrlError(location=[0,7,2,1])
for error in error_subset:
    print("name: {}, location: {}".format(error.name, error.location))
# error_subset_2 = [x_error_1]
# print(error_subset_1[0].location)
# print(error_subset_2[0].location)

# Create syndrome with a randomly placed error
syndrome0 = StabiliserCycle(location=[0])
syndrome1 = StabiliserCycle(location=[1])

syndrome0.run(data, ancilla, eng, error_subset=error_subset)
syndrome1.run(data, ancilla, eng, error_subset=error_subset)

All(Measure) | data
All(Measure) | ancilla
