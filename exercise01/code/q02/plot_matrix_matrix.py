import matplotlib.pyplot as plt
import numpy as np

f = open("./data/matrix_matrix_times.txt", "r")

# Read lines and remove all space characters. 
d = {}
cnt = 0
for line in f:
    # Format data and save it to a list. 
    tmpList = list(filter(None, line.split(' ')))
    # Remove last element ('\n'). 
    if cnt != 0:
        tmpList.pop()
        x = [float(k) for k in tmpList[2:]]
        d[tmpList[0]] = [float(tmpList[1]), x]
    else:
        d[tmpList[0]] = [tmpList[1], tmpList[2:]]
    cnt += 1

# Adapt the block value labelling. 
nums = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
formatBr, formatBc = [], []
for i in d['N'][1]:
    # Keep only the numbers in the labels. 
    formattedString = ''
    for char in i:
        if char in nums: 
            formattedString += char
    # Store formatted value into corresponding list. 
    if i[0:2] == 'br':
        formatBr.append(int(formattedString))
    else:
        formatBc.append(int(formattedString))
# Store created lists into data dict. 
d['N'][1] = [formatBr, formatBc]

# Split block row and block column data into seperate lists. 
numBr = len(formatBr)
for i in d:
    if i != 'N':
        blockRows = d[i][1][0:numBr]
        blockCols = d[i][1][numBr:]
        d[i][1] = [blockRows, blockCols]

# Plot and save results. 
cnt = 256
for i in range(4):
    if cnt == 256:
        plt.scatter(d['N'][1][0], d['256'][1][0], color='blue',
                    label='B (row major)')
        plt.scatter(d['N'][1][0], d['256'][1][1], color='white',
                    edgecolors='blue', label="B (column major)")
        plt.axhline(d['256'][0], color='grey', linestyle='--',
                    linewidth=0.6, zorder=0)
        plt.ylim(0.01, 0.04)
        plt.yticks(np.arange(0.01, 0.041, 0.005))
        plt.text(5, 0.0285, 'Non optimized multiply', fontsize=8)
    elif cnt == 512:
        plt.scatter(d['N'][1][0], d['512'][1][0], color='red',
                    label='B (row major)')
        plt.scatter(d['N'][1][0], d['512'][1][1], color='white',
                    edgecolors='red', label='B (column major)')
        plt.axhline(d['512'][0], color='grey', linestyle='--',
                    linewidth=0.6, zorder=0)
        plt.ylim(0.1, 0.4)
        plt.yticks(np.arange(0.1, 0.4, 0.05))
        plt.text(5, 0.25, 'Non optimized multiply', fontsize=8)
    elif cnt == 1024:
        plt.scatter(d['N'][1][0], d['1024'][1][0], color='orange',
                    label='B (row major)')
        plt.scatter(d['N'][1][0], d['1024'][1][1], color='white',
                    edgecolors='orange', label='B (column major)')
        plt.axhline(d['1024'][0], color='grey', linestyle='--',
                    linewidth=0.6, zorder=0)
        plt.ylim(1.0, 4.0)
        plt.yticks(np.arange(1.0, 4.0, 0.5))
        plt.text(5, 3.3, 'Non optimized multiply', fontsize=8)
    else:
        plt.scatter(d['N'][1][0], d['2048'][1][0], color='darkviolet',
                    label='B (row major')
        plt.scatter(d['N'][1][0], d['2048'][1][1], color='white',
                    edgecolors='darkviolet', label='B (column major)')
        plt.axhline(d['2048'][0], color='grey', linestyle='--',
                    linewidth=0.6, zorder=0)
        plt.ylim(0, 60)
        plt.yticks(np.arange(0, 60, 10))
        plt.text(5, 47, 'Non optimized multiply', fontsize=8)

    plt.xlabel('Blocksize [Bytes]')
    plt.ylabel('Time [s]')
    plt.xlim(0, 140)
    plt.xticks(np.arange(0, 141, 20))
    plt.legend(loc='upper right', frameon=False)

    path = './results/Matrix_matrix_multiply_%s.png' % (cnt) 
    plt.savefig(path, bbox_inches='tight')
    plt.close()

    cnt *= 2

    # plt.show()
