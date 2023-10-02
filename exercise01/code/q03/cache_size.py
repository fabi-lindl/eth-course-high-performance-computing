import csv
import matplotlib.pyplot as plt
import numpy as np

# Read in data.
n, size, mode0, mode1, mode2 = [], [], [], [], []
cnt = 0
with open('./data/test_plots.csv', newline='') as csvFile:
    f = csv.reader(csvFile, delimiter=";")
    for row in f:
        if cnt != 0:
            n.append(int(row[0]))
            size.append(float(row[1].replace(',', '.')))
            mode0.append(float(row[2].replace(',', '.')))
            mode1.append(float(row[3].replace(',', '.')))
            mode2.append(float(row[4].replace(',', '.')))
        cnt+=1

# Plot and save results. 

plt.scatter(n, mode0, color="white", edgecolors="blue", label="Variant 1")
plt.scatter(n, mode1, color="white", edgecolors="red", label="Variant 2")
plt.scatter(n, mode2, color="white", edgecolors="orange", label="Variant 3")

plt.xlim([100, 1000000])
plt.ylim([0, 0.4])
plt.xscale("log")
plt.yticks(np.arange(0, 0.41, 0.05))
plt.xlabel("Array size [No. elements]")
plt.ylabel("Operations per second [10^9]")
plt.legend(loc="lower left", frameon=False)

# Vertical line for the size display, to inidcate temporal cache locality.  
plt.axvline(7500, color="grey", linestyle="--", linewidth=0.6, zorder=0) # L1
plt.axvline(61291, color="grey", linestyle="--", linewidth=0.6, zorder=0) # L2 

# Indicate the different cache sections. 
plt.text(1000, 0.36, "L1", fontsize=10) # L1 section
plt.text(20000, 0.36, "L2", fontsize=10) # L2 section
plt.text(200000, 0.36, "L3", fontsize=10) # L3 section
plt.text(8300, 0.04, "32 KB", fontsize=8)
plt.text(66000, 0.04, "256 KB", fontsize=8)

path = "./results/Cache_size.png"
plt.savefig(path, bbox_inches="tight")

plt.show()
plt.close()
