#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

u = np.genfromtxt("a.csv", delimiter=',', names=True)

theta = u['theta']
grad = u['grad']
exact = u['exact']

plt.plot(theta, grad)
plt.plot(theta, exact, 'k--')
plt.savefig("a.pdf")


