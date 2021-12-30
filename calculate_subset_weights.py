import itertools
from math import comb
import json
import matplotlib.pyplot as plt
from calc_log_e_rate_arbitrary_error_model import weighted_logical_error_rate


filename = "significant_subsets_Dress_pl_2e4"
num_x_gates = 12
num_y_gates = 18
num_2q_gates = 24
num_dephasing_locations = 78  # 24x2 from 2 qubit gates, 30 from single qubit gates

p1q = 0.0005
p2q = 0.005
pdep = 0.0002
pheat = 0.001

error_subset_list = []
lists = [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4]]
error_subsets=list(itertools.product(*lists))
for ele in error_subsets:

    ele_prob1 = comb(num_x_gates, ele[0])*p1q**ele[0]*(1-p1q)**(num_x_gates-ele[0])
    ele_prob2 = comb(num_y_gates, ele[1]) * p1q ** ele[1] * (1 - p1q) ** (num_y_gates - ele[1])
    ele_prob3 = comb(num_2q_gates, ele[2])*p2q**ele[2]*(1-p2q)**(num_2q_gates-ele[2])
    ele_prob4 = comb(num_dephasing_locations, ele[3]) * pdep ** ele[3] * (1 - pdep) ** (num_dephasing_locations - ele[3])
    ele_prob5 = comb(num_2q_gates, ele[4]) * pheat ** ele[4] * (1 - pheat) ** (num_2q_gates - ele[4])
    ele_prob = ele_prob1 * ele_prob2 * ele_prob3 * ele_prob4 * ele_prob5
    if ele_prob>1e-6:
        error_subset_list.append(ele)
        # print(ele)
        # print(ele_prob)
# print(len(error_subset_list))
# with open("imp_samp/" + filename + '.json', 'w') as file:
#     json.dump(error_subset_list, file)

with open("imp_samp\\1000_Dress_imp_samp.json", 'r') as file:
    log_e_dict = json.load(file)

with open("imp_samp\\significant_subsets_Dress_pl_2e4.json") as file:
    error_subsets = json.load(file)

probs = (p1q, p1q, p2q, pdep, pheat)
pl, ws, pl_subsets = weighted_logical_error_rate(probs=probs, esets=error_subsets, log_e_dict=log_e_dict)

plt.bar(x=[str(e) for e in error_subsets],height=pl_subsets, width=0.1)
plt.ylabel('logical error rate')
plt.xlabel('error subset')
plt.xticks(rotation='vertical')
plt.show()