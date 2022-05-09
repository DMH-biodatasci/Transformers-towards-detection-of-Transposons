import numpy as np

#Arrays have been stored as 1D-arrays
frequencies = np.loadtxt('dataset/frequencies.csv', delimiter=',').astype(int) 
mcls = np.loadtxt('dataset/mcls.csv', delimiter=',').astype(str)
labels = np.loadtxt('dataset/labels.csv', delimiter=',').astype(int)

#Bring 1D-arrays into the correct shape:
frequencies = frequencies.reshape(807447, 150, 1)
mcls = mcls.reshape(807447, 150, 1)
labels = labels.reshape(807447, 150, 1)

print(frequencies.shape, mcls.shape, labels.shape)
