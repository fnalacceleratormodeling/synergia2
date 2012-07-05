#!/usr/bin/env python
import numpy as np
from matplotlib  import pyplot

def vertex():
    (u0, u1) = (2.2355, 0.0)
    (v0, v1) = (2.811, 2.811)
    (r0, r1) = (0.5755, 2.250 + 0.0435)
    (ang1, ang2) = (78.43, 66.86)
    d2r = np.pi / 180.0

    x = np.zeros(20,'d')
    y = np.zeros(20,'d')

    for i in range (0, 4):
        j = i
        rot_ang = (i + 1) * ang1 / 4
        x[j] = u0 + r0 * np.cos(rot_ang * d2r)
        y[j] = u1 + r0 * np.sin(rot_ang * d2r)

    for i in range(0, 10):
        j = i + 4
        rot_ang = (270 - (90 - ang1)) - (i + 1) * ang2 / 10
        x[j] = v0 + r1 * np.cos(rot_ang * d2r)
        y[j] = v1 + r1 * np.sin(rot_ang * d2r)

    for i in range (0, 3):
        j = i + 14
        rot_ang = (90- ang1) + (i + 1) * ang1 / 4
        x[j] = u1 + r0 * np.cos(rot_ang * d2r)
        y[j] = u0 + r0 * np.sin(rot_ang * d2r)
    return x * 0.0254, y * 0.0254

if __name__ == "__main__":
    (x, y) = vertex()
    pyplot.plot(x, y)
    pyplot.show()
