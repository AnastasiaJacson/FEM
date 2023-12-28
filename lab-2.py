import matplotlib.pyplot as plt
import numpy as np

import triangle as tr

A = dict(vertices=[[0, 0], [1, 0], [1, 1], [0, 1], [2, 0.5]])
#B = tr.triangulate(A, 'a0.2')
B = tr.triangulate(A, 'qa0.2c')
print(B)
for i, v in enumerate(B['vertices']):
    print(f'{i}: {v}')

tr.compare(plt, A, B)
plt.show()