import numpy as np
import pickle
from matplotlib import pyplot as plt

cor = pickle.load(open("/Users/Maidenhair/Desktop/c_w_kde_mean_s_0.025.p","rb"))

incor = pickle.load(open("/Users/Maidenhair/Desktop/i_w_kde_mean_s_0.025.p","rb"))

x = np.linspace(0,1,100)

plt.plot(x,cor.pdf(x))

plt.plot(x,incor.pdf(x))

plt.show()