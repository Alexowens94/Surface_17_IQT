import time
from calc_log_e_rate_arbitrary_error_model import *
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *

correction_table = load_lookup_table("correction_table_depolarising.json")
bitflip_table = load_lookup_table("correction_table_classical_bitflip.json")
runs = 1000
model = 'PDD_cooling'
filename_s1_impsamp = "ImpSamp_1000_syndrome_process_s1subsets_{}".format(model)
filename_s2_impsamp = "ImpSamp_1000_syndrome_process_s2subsets_{}".format(model)
subsets_file_s1 = "incremental_subsets_s1_{}".format(model) # sim would end at s1
subsets_file_s1_and_s2 = "incremental_subsets_s1ands2_{}".format(model) # 2 rounds locations, no account of if end after 1 or 2 rounds
subsets_file_s2_only = 'incremental_subsets_s2_{}'.format(model)

cz_comp = False
cooling = True
types = ['Rx_ctrl]', 'Ry_ctrl', 'Rxx_ctrl', 'dephase', 'heat', 'Idle']
locations = [12, 18, 24, 78, 24, 119]  # locations errors can occur in pdd model in round 1
locations_both_rounds = [2*12, 2*18, 2*24, 2*78, 2*24, 2*119] # locations for errors to occur, 2 rounds
prob_list = [1e-4, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3]
# since cooling in between timesteps, don't model gate errors changing throughout circuit
prob_lists = [[prob_list[0]],
              [prob_list[1]],
              [prob_list[2]],
              [prob_list[3]],
              [prob_list[4]],
              [prob_list[5]]]
prob_lists_both_rounds = [[prob_list[0]],
              [prob_list[1]],
              [prob_list[2]],
              [prob_list[3]],
              [prob_list[4]],
              [prob_list[5]]]

# Generate significant subsets (if already generated can comment out below)
#########################################
# generate significant subsets for runs with 1 round of stabilizer measurement
t_begin = time.perf_counter()
results_dict = calculate_significant_subsets_incrementally(locations, prob_lists, cutoff_weight=1e-6)
t_end = time.perf_counter()
print('time taken to calc s1 significant subsets {}'.format(t_end-t_begin))
with open("imp_samp\\"+subsets_file_s1+".json", 'w') as file:
    json.dump(results_dict, file)
# generate significant subsets for runs with 2 rounds of stabilizer measurement
# (not accounting for how likely the error subset is to result in stopping after round 1)
t_begin = time.perf_counter()
results_dict = calculate_significant_subsets_incrementally(locations_both_rounds, prob_lists_both_rounds,
                                                           cutoff_weight=1e-6)
t_end = time.perf_counter()
print('time taken to calc s1 and s2 significant subsets {}'.format(t_end-t_begin))
with open("imp_samp\\"+subsets_file_s1_and_s2+".json", 'w') as file:
    json.dump(results_dict, file)
#########################################
# Importance Sampling
# run 1 round of error correction, inserting a specified number of each type of error, in random locations
s1_log_e_rate_dict = {}  # values = e rates given 'key' errors in syndrome 1
s2_rate_dict = {}  # values equal rate a second syndrome is measured, given 'key' errors in syndrome 1
results_list = []
fail_count_dict = {}
s2_count_dict={}
esets = []
with open("imp_samp\\"+subsets_file_s1+".json") as file:
    data = json.load(file)
    s1_significant_subsets = data['error_subset_list']
    locations = data['locations']
    prob_lists = data['probabilities']
with open("imp_samp\\"+subsets_file_s1_and_s2+".json") as file:
    data = json.load(file)
    significant_subsets_s1_and_s2 = data['error_subset_list']
    locations_both_rounds = data['locations']
    prob_lists_both_rounds = data['probabilities']
# importance sampling simulate s1 circuit
################################################################
for i, eset in enumerate(s1_significant_subsets):
    print('sim number {}'.format(i))
    print("calc log e rate of {} in S1".format(eset))
    eset_log_e_rate, eset_s2_rate, fails, s2count = s1_calculate_log_e_rate_error_subset(runs, correction_table, model, eset, locations,
                                                        bitflip_table, prob_lists, cz_comp, cooling=cooling)
    key = str(eset)
    s1_log_e_rate_dict[key] = eset_log_e_rate
    s2_rate_dict[key] = eset_s2_rate
    fail_count_dict[key] = fails
    s2_count_dict[key] = s2count
    esets.append(key)
s1_results = {'s1_log_e_rate_dict': s1_log_e_rate_dict, 's2_rate_dict': s2_rate_dict,
              'failcount': fail_count_dict, 's2count': s2_count_dict}
with open("imp_samp\\"+filename_s1_impsamp+".json", 'w') as file:
    json.dump(s1_results, file)
################################################################################
# load s1 importance sampling results
with open("imp_samp\\"+filename_s1_impsamp+".json", 'r') as file:
    results = json.load(file)
    s2_rate_dict = results['s2_rate_dict']
# calculate significant subsets for s2 simulations from the s2 rates the impsamp for s1 spits out
t1 = time.perf_counter()
print('starting s2 sig sets calc')
s2res = s2sim_calculate_significant_subsets_incrementally(locations, prob_lists, s2_rate_dict, s1_significant_subsets,
                                                          significant_subsets_s1_and_s2)
t2 = time.perf_counter()
print('time taken to calc s2 sig subsets{}'.format(t2-t1))
with open("imp_samp\\"+subsets_file_s2_only+".json", 'w') as file:
    json.dump(s2res, file)
######################################################################################
# load s2 only significant subsets (calculated from s1 imp samp data)
with open("imp_samp\\"+subsets_file_s2_only+".json", 'r') as file:
    s2result = json.load(file)
s2_significant_subsets = s2result['error_subset_list']
print('s2 sig subsets')
print(s2_significant_subsets)
print(len(s2_significant_subsets))
# s2 importance sampling
s2_s2_rate_dict = {}
s2_log_e_rate_dict = {}
fail_count_dict = {}
s2_count_dict = {}
for i, eset in enumerate(s2_significant_subsets):
    print('sim number {}'.format(i))
    print("calc log e rate of {} for 2 syndrome measurements (S2)".format(eset))
    #replace with loading loc from subset file
    eset_log_e_rate, eset_s2_rate, fails, s2count = s2_calculate_log_e_rate_error_subset(runs, correction_table, model,
                                                                                         eset, locations_both_rounds,
                                                                                         bitflip_table, prob_lists_both_rounds,
                                                                                         cz_comp, cooling=cooling)
    key = str(eset)
    s2_log_e_rate_dict[key] = eset_log_e_rate
    s2_s2_rate_dict[key] = eset_s2_rate
    fail_count_dict[key] = fails
    s2_count_dict[key] = s2count
    esets.append(key)
s2_results = {'s2_log_e_rate_dict': s2_log_e_rate_dict, 's2_rate_dict': s2_rate_dict, 'failcount': fail_count_dict,
              's2count': s2_count_dict}

with open("imp_samp\\"+filename_s2_impsamp+".json", 'w') as file:
    json.dump(s2_results, file)