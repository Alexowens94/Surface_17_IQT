import calc_log_e_rate_arbitrary_error_model
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *


correction_table = load_lookup_table("correction_table_depolarising.json")
bitflip_table = load_lookup_table("correction_table_classical_bitflip.json")

ps = [0.01, 0.007, 0.005]
# ps = [0]
runs = 2000
run_name = "test_cz_comp"

log_e_rate_list = []
for p in ps:
    p_set = [12*[p/10],  # px
             18*[p/10],  # py
             24*[p],  # pxx
             30*[0],  # dephasing
             24*[0]]  # heating
    e_model, model_probs = instantiate_error_model(p_set, 'PDD_model')
    print(e_model['gates'], model_probs)
    filename = run_name+"_p={}".format(p)
    log_e_rate = calc_log_e_rate_arbitrary_error_model.calculate_log_e_rate(runs, filename, correction_table, e_model,
                                                                            model_probs, bitflip_table,
                                                                            cz_compilation=False)
    log_e_rate_list.append(log_e_rate)

#log_e_rate = calculate_log_e_rate.calculate_log_e_rate(100, "this_file3", 0.7)

plt.scatter(ps, log_e_rate_list)
plt.ylabel('logical error rate')
plt.xlabel('p')
plt.show()