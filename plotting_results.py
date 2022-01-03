import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit
import numpy as np

filename = 'blah'
with open(filename, 'r') as file:
    results=json.load(file)
#tally = results["rounds_til_fail_tally"]
data = results["rounds_til_fail_list"]

num_bins = 10

def f(t, A, r):
    return A*np.exp(-r*t)
vals,bins,patches = plt.hist(data,bins=num_bins)
popt, pcov = curve_fit(f, bins[:-1], vals)
print(popt)
plt.ylabel('runs ending after t QEC round')
plt.xlabel('t')
plt.plot(bins, f(bins, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f')
plt.title('pL={:2f}%'.format(popt[1]*100))
plt.show()

