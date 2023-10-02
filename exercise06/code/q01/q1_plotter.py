"""
Plots the timing results of the cpp code. 
"""
import matplotlib.pyplot as plt
import numpy as np

# Data for plots. 
grids = [1535, 2048, 2560, 3072, 3584, 4096, 4608, 5120, 5632, 6144]
# Wrote down the results from the terminal output. 
time = [0.057, 0.026, 0.100, 0.042, 0.151, 0.057, 0.218, 0.075, 0.294, 0.107,
        0.421, 0.224, 0.489, 0.165, 0.610, 0.202, 0.738, 0.237, 0.897, 0.444]
serialTime = [time[i] for i in range(len(time)) if i % 2 == 0]
parallelTime = [time[i] for i in range(len(time)) if i % 2 != 0]
su = []
for i in range(len(serialTime)):
    su.append(serialTime[i]/parallelTime[i])

# Plotting. 
# Time. 
# Plot concentration over time. 
plt.plot(grids, serialTime, marker='o', label='serial')
plt.plot(grids, parallelTime, marker='o', label='parallel')

plt.xlabel('Grid size NxN [1]')
plt.ylabel('Time [s]')
plt.xlim(1536, 6144)
plt.xticks(np.arange(1536, 6145, 512))
plt.ylim(0, 1)
plt.yticks(np.arange(0, 1.1, 0.2))
plt.legend(loc='upper left', frameon=False)
plt.title('I/O File write performance (8 MPI ranks)')

path = './results/concentration_over_time.png' 
plt.savefig(path, bbox_inches='tight')
plt.show()
plt.close()

# Speed-up. 
plt.plot(grids, su, marker='o', label='SU')

plt.xlabel('Grid size NxN [1]')
plt.ylabel('Time [s]')
plt.xlim(1536, 6144)
plt.xticks(np.arange(1536, 6145, 512))
plt.ylim(0, 4)
plt.yticks(np.arange(0, 4.1, 1))
plt.legend(frameon=False)
plt.title('I/O File write performance SU (8 MPI ranks)')

path = './results/serial_parallel_su.png' 
plt.savefig(path, bbox_inches='tight')
plt.show()
plt.close()
