import calculate_log_e_rate
import matplotlib.pyplot as plt
from compiled_surface_code_pdd_error_model import *


correction_table=load_lookup_table("C:\\Users\\daisy\\Documents\\Code\\Surface_17_IQT\\correction_table_depolarising.json")

pxctrls = [5e-4, 1e-3]
pyctrls = pxctrls
pxxctrls = 10*np.array(pxctrls)

log_e_rate_list = []
for p in range(len(pxctrls)):
    pxctrl = pxctrls[p]
    pyctrl = pyctrls[p]
    pxxctrl = pxxctrls[p]
    filename = "pxctrl = {}".format(pxctrls[p])
    log_e_rate = calculate_log_e_rate.calculate_log_e_rate(100, filename, correction_table, pxctrl=pxctrl, pyctrl=pyctrl, pxxctrl=pxxctrl)
    log_e_rate_list.append(log_e_rate)

#log_e_rate = calculate_log_e_rate.calculate_log_e_rate(100, "this_file3", 0.7)

plt.scatter(pxctrls, log_e_rate_list)
plt.ylabel('logical error rate')
plt.xlabel('pxctrl')
plt.show()