"""

Here we keep the error model logic.

"""

# Pulsed Dynamical Decoupling Error Model

# Gate                 |       Associated errors        |   Effect of error     |   Probability of error
#----------------------------------------------------------------------------------------------------------
# Rx                   |       Controlled-x             |   X                   |   pxctrl
#                      |       Dephasing                |   Z                   |   pd
# Ry                   |       Controlled-y             |   Y                   |   pyctrl
#                      |       Dephasing                |   Z                   |   pd
# Rxx                  |       Rotation error           |   X on both qubits    |   Pxx      (XX error is projected into two X errors)
#                      |       Motional error           |   X on both qubits    |   Ph       (XX error is projected into two X errors)
#                      |       Dephasing error          |   Z on both qubits    |   Pd


# Dressed-state error model

# Gate                 |       Associated errors        |   Effect of error     |   Probability of error
#----------------------------------------------------------------------------------------------------------
# Rx                   |       Controlled-x             |   X                   |   pxctrl
#                      |       Leakage                  |   L                   |   pl
# Ry                   |       Controlled-y             |   Y                   |   pyctrl
#                      |       Leakage                  |   L                   |   pl
# Rxx                  |       Rotation error           |   X on both qubits    |   Pxx      (XX error is projected into two X errors)
#                      |       Motional error           |   X on both qubits    |   Ph       (XX error is projected into two X errors)
#                      |       Leakage error            |   L on both qubits    |   Pl
#                      |                                |   Removes all future MS gates


# * Use of capital P to indicate probability of 2 qubit error
