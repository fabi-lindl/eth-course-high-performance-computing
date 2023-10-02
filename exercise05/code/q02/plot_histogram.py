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
    (1, 1760),
    (2, 1352),
    (3, 1152),
    (4, 1056),
    (5,  908),
    (6,  868),
    (7,  796),
    (8,  736),
    (9,  716),
    (10, 656)
]
labels = [i[0] for i in data]

# The 3 smallest matrices. 
# Plot and save results. 
x = np.arange(10)
width = 0.45

fig, ax = plt.subplots()
rects1 = ax.bar(x, [data[i][1] for i in range(10)], width)

ax.set_xlabel('Concentration distribution [1]')
ax.set_ylabel('Grid elements [1]')
ax.set_xticks(x)
ax.set_xticklabels(labels)
plt.ylim(0, 2000)
ax.legend(frameon=False)

autolabel(rects1)

fig.tight_layout()

path = './results/histogram_results.png'
plt.savefig(path, bbox_inches='tight')

plt.show()
