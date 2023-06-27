import numpy as np
import matplotlib.pyplot as plt
plt.style.use("seddon_TUD")

a = np.loadtxt("182x.txt")
b = np.loadtxt("182y.txt")

c = np.loadtxt("102x.txt")
d = np.loadtxt("102y.txt")

e = np.loadtxt("52x.txt")
f = np.loadtxt("52y.txt")
plt.plot(a, b, label = "182")
plt.plot(c, d, label = "102")
plt.plot(e, f, label = "52")
plt.legend()
plt.show()
