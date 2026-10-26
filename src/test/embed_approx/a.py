#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

u = np.genfromtxt("a.csv", delimiter=",", names=True)

theta = u["theta"]
grad = u["grad"]
exact = u["exact"]

plt.plot(theta, grad)
plt.plot(theta, exact, "k--")
plt.savefig("a.pdf")
