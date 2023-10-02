import matplotlib.pyplot as plt
import numpy as np

f = open("./data/transpose_times.txt", "r")

# Read lines and remove all space characters. 
d = {}
cnt = 0
for line in f:
    # Format data and save it to a list. 
    tmpList = list(filter(None, line.split(' ')))
    # Remove the last element from the list ('\n'). 
    tmpList.pop()
    if cnt != 0:
        d[tmpList[0]] = [float(tmpList[1]), [float(x) for x in tmpList[2:]]]
    cnt += 1

# List of block values. 
blocks = [2, 4, 8, 16, 32, 64, 128]

# Plot and save results. 
cnt = 256
for i in range(6):
    print(cnt)
    if cnt == 256:
        plt.scatter(blocks, d['256'][1], color='white', edgecolors='blue')
        plt.axhline(d['256'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.0004, 0.0009)
        plt.yticks(np.arange(0.0004, 0.00091, 0.0001))
        plt.text(5, 0.00078, 'Non optimized transposition', fontsize=8)
    elif cnt == 512:
        plt.scatter(blocks, d['512'][1], color='white', edgecolors='red')
        plt.axhline(d['512'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.002, 0.004)
        plt.yticks(np.arange(0.002, 0.0041, 0.0005))
        plt.text(25, 0.00315, 'Non optimized transposition', fontsize=8)
    elif cnt == 1024:
        plt.scatter(blocks, d['1024'][1], color='white', edgecolors='orange')
        plt.axhline(d['1024'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.012, 0.017)
        plt.yticks(np.arange(0.012, 0.0171, 0.001))
        plt.text(80, 0.014, 'Non optimized transposition', fontsize=8)
    elif cnt == 2048:
        plt.scatter(blocks, d['2048'][1], color='white', edgecolors='darkviolet')
        plt.axhline(d['2048'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.04, 0.08)
        plt.yticks(np.arange(0.04, 0.081, 0.01))
        plt.text(5, 0.0713, 'Non optimized transposition', fontsize=8)
    elif cnt == 4096:
        plt.scatter(blocks, d['4096'][1], color='white', edgecolors='olivedrab')
        plt.axhline(d['4096'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.2, 0.32)
        plt.yticks(np.arange(0.2, 0.321, 0.04))
        plt.text(5, 0.308, 'Non optimized transposition', fontsize=8)
    else:
        plt.scatter(blocks, d['8192'][1], color='white', edgecolors='saddlebrown')
        plt.axhline(d['8192'][0], color='grey', linestyle='--', linewidth=0.6, zorder=0)
        plt.ylim(0.8, 1.4)
        plt.yticks(np.arange(0.8, 1.41, 0.1))
        plt.text(5, 1.31, 'Non optimized transposition', fontsize=8)

    plt.xlabel('Blocksize [Bytes]')
    plt.ylabel('Time [s]')
    plt.xlim(0, 140)
    plt.xticks(np.arange(0, 141, 20))

    path = './results/Matrix_transpose_%s.png' % (cnt) 
    plt.savefig(path, bbox_inches='tight')
    plt.close()

    cnt *= 2

    # plt.show()
