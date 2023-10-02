import glob
import matplotlib.pyplot as plt
import numpy as np

# Read in txt files from this folder. 
path = "./data/*.txt"
files = glob.glob(path)

# Dictionary to store the y data values (time values). 
# Key is number of the method in the C++ code (1, 2, 3). 
# Values: List of two lists.
# First list contains first time results, i.e. code run for the first time.
# Second list contains second time results, i.e. code run for the second time. 
d = {
    0: [[]], # Serial code. 
    1: [[], []],
    2: [[], []],
    3: [[], []],
}

for f in range(len(files)):
    # Get the file name.
    fname = files[f]
    # Get the computation version. 
    vers = int(fname[3])
    # Check for first or second run. 
    if len(fname) > 13:
        first = 1 # second run. 
    else:
        first = 0 # first run. 

    _file = open(fname, 'r')
    for line in _file:
        if line[0:4] == 'time':
            # Get the time. 
            t = float(line[6:])
            d[vers][first].append(t)
            # Break loop as only this one line is of relevance. 
            break

# Plot and save results.

x = [1, 2, 4, 8] # Tested thread values (x-axis values).
labels = ['No arrays', 'Arrays (no padding)', 'Arrays (padding)']
colors = ['orange', 'blue', 'darkviolet']

for i in range(len(d)):
    if i == 0:
        # Plot reference line of serial code. 
        plt.axhline(d[i][0][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.text(2, 2.45, 'Serial program', fontsize=8)
    else:
        plt.plot(x, d[i][0], marker='.', markersize=9, color=colors[i-1], label=labels[i-1])
        plt.plot(x, d[i][1], marker='.', markersize=9, color=colors[i-1], linestyle='--')

plt.xlabel('No. Threads [1]')
plt.ylabel('Time [s]')
plt.xlim(1, 8)
plt.ylim(0, 3)
plt.xticks(np.arange(1, 8.5, 1))
plt.yticks(np.arange(0, 3.5, 0.5))
plt.legend(loc='upper right', frameon=False)

p = './results/Exercise_2_Question_1.png'
plt.savefig(p, bbox_inches='tight')

plt.close()
