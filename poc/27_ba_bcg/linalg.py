#!/usr/bin/env python

import numpy as np


def Solve(a, b):
    print(a, b)
    x = np.linalg.solve(a, b)
    print(x)

Solve(
        np.array([
            [1, 1, 1],
            [0.5, 1, 2],
            [0.25, 1, 4],
            ])
        ,
        np.array([
            0,
            1,
            2,
            ])
        )

# du/dx(0.5) = k0 * u(0.5) + k1 * u(1) + k2 * u(2)
Solve(
        np.array([
            [1, 1, 1],
            [0.5, 1, 2],
            [0.25, 1, 4],
            ])
        ,
        np.array([
            0,
            1,
            1,
            ])
        )
