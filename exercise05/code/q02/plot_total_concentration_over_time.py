import matplotlib.pyplot as plt
import numpy as np

f = open("./data/concentration_over_time.txt", "r")

# Read lines and store time values and concentration values
# in two separate lists.
t, c = [], []
cnt = 0
for line in f:
    # Format data and save it to a list.
    if cnt % 2 == 0:
        x = line.split(' ')
        t.append(float(x[0])) # Store the time value.
        c.append(float(x[1])) # Store the concentration value.
    cnt+=1
print(len(c))

plt.plot(t, c)

plt.xlabel('Time [s]')
plt.ylabel('Total Concentration [1]')
plt.xlim(0, 1)
plt.xticks(np.arange(0, 1.2, .2))
plt.ylim(0, 1.2)
plt.yticks(np.arange(0, 1.3, 0.2))

path = './results/concentration_over_time.png'
plt.savefig(path, bbox_inches='tight')
plt.show()

plt.close()
