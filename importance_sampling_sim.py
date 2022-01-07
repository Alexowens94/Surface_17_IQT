from calc_log_e_rate_arbitrary_error_model import calculate_log_e_rate_error_subset
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *


correction_table = load_lookup_table("correction_table_depolarising.json")
bitflip_table = load_lookup_table("correction_table_classical_bitflip.json")

# run 1 round of error correction, inserting a specified number of each type of error, in random locations
runs = 1000
subsets_file = "significant_subsets_dress"
model = 'Dress_model'
cz_comp = False

filename = "ImpSamp_leakageredux_{}".format(runs)+subsets_file
print(filename)
log_e_rate_dict = {}
log_e_rate_list = []
results_list = []
esets = []

with open("imp_samp\\"+subsets_file+".json") as file:
    data = json.load(file)
    error_subsets = data['error_subset_list']
    locations = data['error_locations']


for eset in error_subsets:
    print(eset)
    eset_log_e_rate = calculate_log_e_rate_error_subset(runs, correction_table, model, eset, locations,
                                                        bitflip_table, cz_comp)
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
