import calc_log_e_rate_arbitrary_error_model
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *


correction_table = load_lookup_table("correction_table_depolarising.json")
bitflip_table = load_lookup_table("correction_table_classical_bitflip.json")

# run 1 round of error correction, inserting a specified number of each type of error, in random locations

runs = 1000
filename = "1000_Dress_imp_samp_reset_leakage_each_run"

log_e_rate_dict = {}
log_e_rate_list = []
results_list = []
esets = []

model = 'Dress_model'
with open("imp_samp\\significant_subsets_Dress_pl_2e4.json") as file:
    error_subsets = json.load(file)

for eset in error_subsets:

    eset_log_e_rate = calc_log_e_rate_arbitrary_error_model.calculate_log_e_rate_error_subset(runs, filename, correction_table,
                                                                                      model, eset, bitflip_table)
    key = str(eset)
    log_e_rate_dict[key] = eset_log_e_rate
    log_e_rate_list.append(eset_log_e_rate)
    esets.append(key)


with open("imp_samp/" + filename + '.json', 'w') as file:
    json.dump(log_e_rate_dict, file)
plt.bar(esets, log_e_rate_list, width=0.1)
plt.ylabel('logical error rate')
plt.xlabel('error subset')
plt.xticks(rotation='vertical')
plt.show()