import random
import matplotlib.pyplot as plt
import numpy as np
from surface_code import *
import pypartitions
import json
from math import sqrt
from projectq.ops import All, CNOT, H, Measure, X, Z, Y, Rx, Ry, Rxx, Deallocate, Command, SqrtX, MatrixGate
from projectq import MainEngine
from projectq.backends import Simulator, CommandPrinter,CircuitDrawerMatplotlib#,QEC_Simulator
from projectq.backends import CommandPrinter
from projectq import MainEngine
from projectq.meta import Loop
from projectq.setups.default import get_engine_list
cmd_engine = CommandPrinter()
drawing_engine=CircuitDrawerMatplotlib()
pi2Y = MatrixGate([[1 / sqrt(2), -1 / sqrt(2)], [1 / sqrt(2), 1 / sqrt(2)]])
minpi2Y = MatrixGate([[1 / sqrt(2), 1 / sqrt(2)], [-1 / sqrt(2), 1 / sqrt(2)]])
pi2X = MatrixGate([[1 / sqrt(2), -1j / sqrt(2)], [-1j / sqrt(2), 1 / sqrt(2)]])
minpi2X = MatrixGate([[1 / sqrt(2), 1j / sqrt(2)], [1j / sqrt(2), 1 / sqrt(2)]])
pi4XX = MatrixGate([[1 / sqrt(2), 0, 0, -1j / sqrt(2)],
                    [0, 1 / sqrt(2), -1j / sqrt(2), 0],
                    [0, -1j / sqrt(2), 1 / sqrt(2), 0],
                    [-1j / sqrt(2), 0, 0, 1 / sqrt(2)]
                    ])
minpi4XX = MatrixGate([[1 / sqrt(2), 0, 0, 1j / sqrt(2)],
                       [0, 1 / sqrt(2), 1j / sqrt(2), 0],
                       [0, 1j / sqrt(2), 1 / sqrt(2), 0],
                       [1j / sqrt(2), 0, 0, 1 / sqrt(2)]
                       ])

pi4ms = np.matrix([[1 / sqrt(2), 0, 0, -1j / sqrt(2)], #pi4xx
                    [0, 1 / sqrt(2), -1j / sqrt(2), 0],
                    [0, -1j / sqrt(2), 1 / sqrt(2), 0],
                    [-1j / sqrt(2), 0, 0, 1 / sqrt(2)]
                    ])
negpi4ms = np.matrix([[1 / sqrt(2), 0, 0, 1j / sqrt(2)],
                       [0, 1 / sqrt(2), 1j / sqrt(2), 0],
                       [0, 1j / sqrt(2), 1 / sqrt(2), 0],
                       [1j / sqrt(2), 0, 0, 1 / sqrt(2)]
                       ])
m3 = np.matrix([[1,0,0,0],
                [0,-1,0,0],
                [0,0,1,0],
                [0,0,0,-1]])
ix = np.matrix([[0,1,0,0],
               [1,0,0,0],
               [0,0,0,1],
               [0,0,1,0]])

iypi2 = np.matrix([[1/sqrt(2),-1/sqrt(2),0,0],
               [1/sqrt(2),1/sqrt(2),0,0],
               [0,0,1/sqrt(2),-1/sqrt(2)],
               [0,0,1/sqrt(2),1/sqrt(2)]])
ipi2x=np.matrix([[1/sqrt(2),-1j/sqrt(2),0,0],
                 [-1j/sqrt(2),1j/sqrt(2),0,0],
                 [0,0,1/sqrt(2),-1j/sqrt(2)],
                 [0,0,-1j/sqrt(2),1/sqrt(2)]])
pi2xi=np.matrix([[1/sqrt(2),0,-1j/sqrt(2),0],
                 [0,1/sqrt(2),0,-1j/sqrt(2)],
                 [-1j/sqrt(2),0,0,1/sqrt(2)],
                 [0,0,-1j/sqrt(2),1/sqrt(2)]])
pi2y=np.matrix([[1/sqrt(2),-1/sqrt(2)],
               [1/sqrt(2),1/sqrt(2)]])
pi2x=np.matrix([[1/sqrt(2),-1j/sqrt(2)],
               [-1j/sqrt(2),1/sqrt(2)]])

negpi2y=np.matrix([[1/sqrt(2),1/sqrt(2)],
               [-1/sqrt(2),1/sqrt(2)]])
negpi2x=np.matrix([[1/sqrt(2),1j/sqrt(2)],
               [1j/sqrt(2),1/sqrt(2)]])

A = np.matmul(ipi2x,pi4ms)
B = np.matmul(pi2xi,pi4ms)
print('Rxpi/2 becomes  {}'.format(np.matmul(negpi4ms,A)))

print('Rxpi/2 on left becomes  {}'.format(np.matmul(negpi4ms,B)))



#
# product1 = np.matmul(pi4ms,ipi2x)
# print(product1)
# print(product1)
# product2 = np.matmul(ipi2x,pi4ms)
# print('com {}'.format(product1-product2))
# error1 = np.matmul(product1,pi2x)
# error2 = np.matmul(product1,negpi2x)
# print('error1 {}'.format(error1))
# print('error2 {}'.format(error2))
# print('product3 {}'.format(product3))
# print('comm12 {}'.format(comm))
# print('error {}'.format(error))

eng = MainEngine(engine_list=get_engine_list()+[drawing_engine])
qbits = eng.allocate_qureg(8)

# I | qbits[0]
# I | qbits[1]

#I | qbits[2]
X | qbits[3]

X | qbits[4]
# I | qbits[5]

X | qbits[6]
X | qbits[7]
display = True

with Loop(eng,10):
    nat_CNOT(qbits[0], qbits[1])
    CNOT | (qbits[0], qbits[1])
    nat_CNOT(qbits[2], qbits[3])
    CNOT | (qbits[2], qbits[3])
    nat_CNOT(qbits[4], qbits[5])
    CNOT | (qbits[4], qbits[5])
    nat_CNOT(qbits[6], qbits[7])
    CNOT | (qbits[6], qbits[7])

All(Measure) | qbits
eng.flush()
# drawing_engine.draw()
# plt.show()
#print([int(q) for q in qbits])  # == [0, 0, 0, 1, 1, 0, 1, 1]

    #print([int(qb) for qb in q])