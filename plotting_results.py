import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.optimize import curve_fit

import numpy as np

with open('test_log_e_rate_til_fail_2000_old_order.txt', 'r') as file:
    results=json.load(file)
tally = results["rounds_til_fail_tally"]
num_bins = 25
bins = num_bins*[0]
max=250
bin_size = max/num_bins


for i in range(num_bins):
    for key in tally.keys():
        if (int(key) <= (i+1)*bin_size) and (i*bin_size < int(key)):
            bins[i] += int(tally[key])
a_bins = np.array(bins)
print(sum(a_bins))

def f(t, A, r):
    return A*np.exp(-r*t)
xdata = np.linspace(0,max,num_bins)
popt, pcov = curve_fit(f, xdata, bins)
print(popt)
plt.ylabel('runs ending after t QEC round')
plt.xlabel('t')
plt.plot(xdata, f(xdata, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f')
plt.plot(xdata, bins)
plt.title('old_order pL={:2f}%'.format(popt[1]*100))
plt.show()

