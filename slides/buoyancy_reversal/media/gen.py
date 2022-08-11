#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 1, 3)

yy = [
    (1 - x, 'reversal0.svg', 'C0'),
    (np.maximum(1 - x * 4, -2 + x * 2), 'reversal1.svg', 'C1'),
]

for y, path, c in yy:
    fig, ax = plt.subplots(figsize=(1.5,1.5))
    ax.set_xticks([0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_ylim(-1, 1)
    ax.set_xlabel('$\chi$', labelpad=-5)
    ax.set_ylabel('$B$', labelpad=0)
    ax.plot(x, y, c=c, lw=2, clip_on=False)
    print(path)
    fig.savefig(path)
    plt.close(fig)
