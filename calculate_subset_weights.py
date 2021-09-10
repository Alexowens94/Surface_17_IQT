import itertools
from math import comb

cycles=2
num_1q_gates=30#8*cycles
num_2q_gates=24*cycles
p1q = 1e-4
p2q = 0.005

lists= [[0,1,2,3],[0,1,2,3,4,5,6]]
error_subsets=list(itertools.product(*lists))
print(error_subsets)
for ele in error_subsets:

    ele_prob1 = comb(num_1q_gates, ele[0])*p1q**ele[0]*(1-p1q)**(num_1q_gates-ele[0])
    ele_prob2 = comb(num_2q_gates, ele[1])*p2q**ele[1]*(1-p2q)**(num_2q_gates-ele[1])
    ele_prob = ele_prob1 * ele_prob2
    if ele_prob>1e-6:
        print(ele)
        print(ele_prob)
