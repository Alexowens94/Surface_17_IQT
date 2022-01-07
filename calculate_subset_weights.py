import itertools
import json
import matplotlib.pyplot as plt
from calc_log_e_rate_arbitrary_error_model import subset_weight, weighted_logical_error_rate


filename = "significant_subsets_pdd"
cutoff_weight = 1e-6  # consider error subsets of statistical weights greater than this significant
types = ['Rx_ctrl]', 'Ry_ctrl', 'Rxx_ctrl', 'leak', 'heat']  # pdd model ms comp
         #['Ry_ctrl', 'Ry_dephase', 'CZ_2q_dephase', 'CZ_2q_leak', '1Q_leak_sigmaZ_noise'] #spin spin cz comp
locations = [34, 34, 24, 24, 48]#[12, 18, 24, 78, 24]#
probs = [0.0001, 0.0005, 0.001, 0.001, 0.0005]

error_subset_list = []
subset_weight_list = []
results_dict = {}
results_dict['error_type'] = types  # human readable description of the errors
results_dict['error_locations'] = locations  # number of places a given error can occur
results_dict['error_probability'] = probs  # probabilities used to calculate subset weights

lists = [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4]]
error_subsets = list(itertools.product(*lists))


for eset in error_subsets:
    params = [(locations[i], probs[i], eset[i]) for i in range(len(eset))]
    eset_weight = subset_weight(params)
    if eset_weight > cutoff_weight:
        print(eset_weight)
        print(params)
        error_subset_list.append(eset)
        subset_weight_list.append(eset_weight)
print(len(error_subset_list))
results_dict['error_subset_list'] = error_subset_list  # significant subsets
results_dict['subset_weights'] = subset_weight_list  # statistical weights of the subsets, order matches subset list
results_dict['cutoff_weight'] = cutoff_weight

with open("imp_samp/" + filename + '.json', 'w') as file:
    json.dump(results_dict, file)


# ######
# with open("imp_samp\\ImpSamp_z0_leakageredux_1000test_significant_subsets_cz.json", 'r') as file:
#     log_e_dict = json.load(file)
#
# with open("imp_samp\\test_significant_subsets_cz.json") as file:
#     data = json.load(file)
#     error_subsets = data['error_subset_list']
#
# probs = [0.0005, 0.0005, 0.005, 0.005, 0.0002]
# print(probs)
# print(error_subsets[0])
# print(locations)
# pl, ws, pl_subsets = weighted_logical_error_rate(probs, error_subsets, log_e_dict, locations)
#
# plt.bar(x=[str(e) for e in error_subsets],height=pl_subsets, width=0.1)
# plt.ylabel('logical error rate')
# plt.xlabel('error subset')
# plt.xticks(rotation='vertical')
# plt.show()