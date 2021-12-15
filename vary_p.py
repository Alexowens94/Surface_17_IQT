import calc_log_e_rate_arbitrary_error_model
import matplotlib.pyplot as plt
from compiled_surface_code_arbitrary_error_model import *


correction_table=load_lookup_table("correction_table_depolarising.json")

ps = [0.005, 0.007, 0.01, 0.02]


log_e_rate_list = []
for p in ps:
    p_set = [p/10, p/10, p, 0, 0]
    e_model, model_probs = instantiate_error_model(p_set, 'PDD_model')
    print(e_model['gates'], model_probs)
    filename = "PDD_no_deph_rand_init_p={}".format(p)
    log_e_rate = calc_log_e_rate_arbitrary_error_model.calculate_log_e_rate(100, filename, correction_table, e_model,
                                                                            model_probs)
    log_e_rate_list.append(log_e_rate)

#log_e_rate = calculate_log_e_rate.calculate_log_e_rate(100, "this_file3", 0.7)

plt.scatter(ps, log_e_rate_list)
plt.ylabel('logical error rate')
plt.xlabel('p')
plt.show()