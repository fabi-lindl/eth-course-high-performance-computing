import matplotlib.pyplot as plt
import numpy as np

def autolabel(rects):
    """
    Attach a text label above each bar in *rects*, displaying its height.
    """
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

# Data was stored in lsf.oxxx files.
data = [
    ( 512, (0.0003, 0.0005)),
    (1024, (0.0013, 0.0036)),
    (2048, (0.0054, 0.0226)),
    (4096, (0.0213, 0.1147)),
    (8192, (0.1034, 0.5933)),
]
labels = [i[0] for i in data]

# The 3 smallest matrices.
x = np.arange(3)
width = 0.45

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, [data[i][1][0] for i in range(3)],
                width, label='Row major')
rects2 = ax.bar(x + width/2, [data[i][1][1] for i in range(3)],
                width, label='Column major')

ax.set_xlabel('Matrix dimension (NxN)')
ax.set_ylabel('Time [s]')
ax.set_xticks(x)
ax.set_xticklabels(labels)
plt.ylim(0, 0.025)
ax.legend(frameon=False)

autolabel(rects1)
autolabel(rects2)

fig.tight_layout()
plt.show()

# All matrices.
x = np.arange(len(data))
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, [i[1][0] for i in data], width, label='Row major')
rects2 = ax.bar(x + width/2, [i[1][1] for i in data], width, label='Col major')

ax.set_xlabel('Matrix dimension (NxN)')
ax.set_ylabel('Time [s]')
ax.set_xticks(x)
ax.set_xticklabels(labels)
plt.ylim(0, 0.7)
ax.legend(loc='upper left', frameon=False)

autolabel(rects1)
autolabel(rects2)

fig.tight_layout()
plt.show()
