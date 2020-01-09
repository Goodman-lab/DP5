import numpy as np

from scipy.stats import gaussian_kde, norm

from matplotlib import pyplot as plt

def l(x,p,a):

    return (a)/(1+(x-p)**2)


x = np.linspace(0,100,1000)

y = np.zeros(len(x))

y+= l(x,20,0.5) + l(x,21.7,1)

plt.plot(x,y)

plt.plot(x[2:],np.diff(y,2)/np.max(np.diff(y,2)))

plt.show()