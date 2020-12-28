import numpy as np


a = np.array([1,2,np.nan])

b = np.array([1,2,np.inf])

print(a)

print(b)

print(np.isnan(a).any())

print(np.isinf(b).any())