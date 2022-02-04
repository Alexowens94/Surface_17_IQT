from calc_log_e_rate_arbitrary_error_model import *
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *
import time

correction_table = load_lookup_table("correction_table_depolarising.json")
bitflip_table = load_lookup_table("correction_table_classical_bitflip.json")

# run 1 round of error correction, inserting a specified number of each type of error, in random locations



runs = 1000
model = 'Dress_cooling'
filename = 'pl_{}'.format(model)
# cz_comp = False


subsets_file_s1 = "incremental_subsets_s1_{}".format(model) # sim would end at s1
subsets_file_s1_and_s2 = "incremental_subsets_s1ands2_{}".format(model) # 2 rounds locations, no account of if end after 1 or 2 rounds
subsets_file_s2_only = 'incremental_subsets_s2_{}'.format(model)
s1is = "ImpSamp_1000_syndrome_process_s1subsets_{}".format(model)
s2is = "ImpSamp_1000_syndrome_process_s2subsets_{}".format(model)
# s1is = "ImpSamp_1000_syndrome_process_s1subsets_{}_nbarheat".format(model)
# s2is = "ImpSamp_1000_syndrome_process_s2subsets_{}_nbarheat".format(model)
# subsets_file_s1 = "incremental_subsets_s1_{}_nbarheat".format(model) # sim would end at s1
# subsets_file_s1_and_s2 = "incremental_subsets_s1ands2_{}_nbarheat".format(model) # 2 rounds locations, no account of if end after 1 or 2 rounds
# subsets_file_s2_only = 'incremental_subsets_s2_{}_nbarheat'.format(model)


with open("imp_samp\\"+subsets_file_s1+".json") as file:
    data = json.load(file)
    s1_significant_subsets = data['error_subset_list']
    s1_weights = data['subset_weights']
    locations_1round = data['locations']
    # prob_lists_s1 = data['probabilities']
with open("imp_samp\\" + subsets_file_s1_and_s2 + ".json") as file:
    data = json.load(file)
    s1_and_s2_significant_subsets = data['error_subset_list']
with open("imp_samp\\" + subsets_file_s2_only + ".json") as file:
    data = json.load(file)
    s2_significant_subsets = data['error_subset_list']
    s2_weights = data['subset_weights']
with open("imp_samp\\"+s1is+".json") as file:
    s1data = json.load(file)
    s1_log_e_rates = s1data["s1_log_e_rate_dict"]
    s2_rate_dict = s1data['s2_rate_dict']
with open("imp_samp\\" + s2is + ".json") as file:
    s2data = json.load(file)
    s2_log_e_rates = s2data["s2_log_e_rate_dict"]
    s2_counts = s2data['s2count']
    s2_fails = s2data['failcount']



def pl_from_file():
    print('pl from file')
    p = 0
    for i, eset in enumerate(s1_significant_subsets):
        print(i, eset)
        erate = s1_log_e_rates[str(eset)]   # prob error GIVEN eset errors and stopping after 1 syndrome
        if erate == -1:
            continue
        weight = s1_weights[i]  # statistical weighting of eset p(eset AND stop after 1 syndrome)
        # print('S1: erate {} weight {}'.format(erate, weight))
        p += erate * weight
    for i, eset in enumerate(s2_significant_subsets):
        print(i, eset)
        erate = s2_fails[str(eset)]/s2_counts[str(eset)]  # prob error GIVEN eset errors and stopping after 2 syndrome
        weight = s2_weights[i]  # statistical weighting of eset p(eset AND stop after 2 syndrome)
        # print('S2: erate {} weight {}'.format(erate, weight))
        p += erate * weight
    return p


def pl(prob_lists_s1, prob_lists_s2):
    print('pl from input probs')

    p = 0
    for i, eset in enumerate(s1_significant_subsets):
        print(i, eset)
        erate = s1_log_e_rates[str(eset)]   # prob error GIVEN eset errors and stopping after 1 syndrome
        if erate == None:  # never stops after one syndrome so can't assign an error rate
            continue
        params = [(locations_1round[i], prob_lists_s1[i], eset[i]) for i in range(len(eset))]
        weight = subset_weight_mixed(params)  # statistical weighting of eset p(eset AND stop after 1 syndrome)
        print('S1: erate {} weight {}'.format(erate, weight))
        p += erate * weight
    for i, eset in enumerate(s2_significant_subsets):
        print(i, eset)
        erate = s2_fails[str(eset)]/s2_counts[str(eset)]  # prob error GIVEN eset errors and stopping after 2 syndrome
        weight = prob_s2_AND_eset_a_total_errors(eset, s2_rate_dict, s1_significant_subsets,
                                                 s1_and_s2_significant_subsets, locations_1round,
                                                 prob_lists_s1, prob_lists_s2)[0]
        # statistical weighting of eset p(eset AND stop after 2 syndrome)
        print('S2: erate {} weight {}'.format(erate, weight))
        p += erate * weight
    return p

pdep = 0.2e-4
pheat = 1e-4
pidle = 1e-4
pts = 15
pls = np.zeros(pts)
p_sum = np.zeros(pts)
# n2bar = [ 1,2.57,4.7,7.4, 10.7, 14.57,19,24]
ps = np.geomspace(1e-2, 1e-6, pts)
t1 = time.perf_counter()
num_steps = 8
gates_per_timestep = 3
for i, p in enumerate(ps):
    pxctrl, pyctrl, pctrl = p/10, p/10, p
    prob_list = [pxctrl, pyctrl, pctrl, pdep, pheat, pidle]
    pxxctrl_heating_s1 = []
    pxxctrl_heating_s2 = []
    for j in range(num_steps):
        pxxctrl_heating_s1 += gates_per_timestep * [(j + 1) * prob_list[2]]
        pxxctrl_heating_s2 += gates_per_timestep * [(j + 1 + num_steps) * prob_list[2]]
    # print('pxxctrl_heating_s1 {}'.format(pxxctrl_heating_s1))
    # print(len(pxxctrl_heating_s1))
    # print('pxxctrl_heating_s2 {}'.format(pxxctrl_heating_s2))
    # print(len(pxxctrl_heating_s2))
    # print('prob_list {}'.format(prob_list))
    prob_lists_s1 = [[prob_list[0]], [prob_list[1]], [prob_list[2]], [prob_list[3]], [prob_list[4]], [prob_list[5]]]
    prob_lists_s2 = [[prob_list[0]], [prob_list[1]], [prob_list[2]], [prob_list[3]], [prob_list[4]], [prob_list[5]]]

    # prob_lists_s1 = [[prob_list[0]], [prob_list[1]], pxxctrl_heating_s1, [prob_list[3]], [prob_list[4]]]
    print('prob_lists_s1 {}'.format(prob_lists_s1))
    # prob_lists_s2 = [[prob_list[0]], [prob_list[1]], pxxctrl_heating_s2, [prob_list[3]], [prob_list[4]]]
    pls[i] = pl(prob_lists_s1, prob_lists_s2)
    p_sum[i] = pheat+pdep+pctrl+pidle
t2 = time.perf_counter()
results = {'p_physical_ctrl': list(ps), 'p_logical': list(pls), 'p_physical_sum': list(p_sum)}
print(results)
print('time taken {}'.format(t2-t1))
with open(filename+'.json', 'w') as file:
    json.dump(results, file)

# ps = [0.01, 0.0089, 0.0078, 0.0067, 0.0056, 0.0045, 0.0034000000000000002, 0.0023, 0.0011999999999999997, 0.0001]
# pls = [0.4219390011185434, 0.43029362846802954, 0.4212656397690426, 0.3917726208626029, 0.3409808513179788, 0.27152632411810496, 0.19037054516999474, 0.10881872489535545, 0.04136016281352608, 0.0033538007645720967]
# p_sum = [0.012, 0.0109, 0.0098, 0.0087, 0.0076, 0.0065, 0.0054, 0.0043, 0.0031999999999999997, 0.0021]

plt.loglog(ps, ps, label='physical')
plt.loglog(ps, pls, label='logical')
plt.loglog(ps, p_sum, label='sum')
plt.ylabel('Error rate')
plt.xlabel('MS control error rate')
plt.title('p_heat {} p_dephase {}'.format(pheat, pdep))
plt.legend()
plt.show()