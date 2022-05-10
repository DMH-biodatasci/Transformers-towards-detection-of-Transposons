# Boundary detection of composite transposable elements:

<p align="center">
    <img src="https://github.com/DMH-dutte/Detection_of_composite_transposable_elements/blob/main/preview/frequency_landscape.png" width="400" />
</p>

## Tasks:

1. Find model architectures that would fit our approach
2. Develop an approach using the gene frequencies -> 1D-data
3. Develop an approach using mcl class labels -> multidimensional Embedding


## Getting started:


### Padding: 0.0, 0 or '0.0' in frequencies or mcls means that this chunk has been padded

### Chunk size: 150
### Overlap: 15

```python
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
```



## Description:
...


