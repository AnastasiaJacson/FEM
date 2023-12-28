import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import math
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LightSource
import matplotlib

import triangle as tr

a11 = 1
a22 = 1
d = 1
#з граничної умови наприклад:
sigma = 1e5
beta = 1e-5
psi = 0
fe = 3 * np.ones(3).reshape(-1, 1)

 
a = [
    lambda t: t[1][0]*t[2][1] - t[2][0]*t[1][1],
    lambda t: t[2][0]*t[0][1] - t[0][0]*t[2][1],
    lambda t: t[0][0]*t[1][1] - t[1][0]*t[0][1]
]

b = [
    lambda t: t[1][1] - t[2][1],
    lambda t: t[2][1] - t[0][1],
    lambda t: t[0][1] - t[1][1]
]

c = [
    lambda t: t[2][0] - t[1][0],
    lambda t: t[0][0] - t[2][0],
    lambda t: t[1][0] - t[0][0]
]

phi = [
    lambda t: 1/delta(t) * (a[0](t) + b[0](t)*t[0][0] + c[0](t)*t[0][1]),
    lambda t: 1/delta(t) * (a[1](t) + b[1](t)*t[1][0] + c[1](t)*t[1][1]),
    lambda t: 1/delta(t) * (a[2](t) + b[2](t)*t[2][0] + c[2](t)*t[2][1])
] 

def delta(t):
    area = 0.5 * (t[0][0] * (t[1][1] - t[2][1]) + t[1][0] * (t[2][1] - t[0][1]) + t[2][0]
                  * (t[0][1] - t[1][1]))
    return 2*area

def dist(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

#функція для отримання бареоцентричних координат для заданої точки і задоно трикутника
def bareocentr(triangle, point):
    A = np.vstack([triangle.transpose(), np.array([1, 1, 1])])
    B = np.vstack([point.reshape(-1, 1), np.array([1])])
    bc_point = np.linalg.solve(A, B)

    return bc_point

def find_triangle(triangles, point):
    for i, triangle in enumerate(triangles):
        bc_point = bareocentr(triangle, point)
        if not np.any(bc_point < 0):
            # index of the triangle
            return i, bc_point
        
        




A = dict(vertices=[[0, 0], [1, 0], [1, 1], [0, 1]])
# B = tr.triangulate(A, 'c')
#B = tr.triangulate(A, 'qa0.3c')
B = tr.triangulate(A, 'qa0.1c')

segments = B['segments']

N = len(B['vertices'])
# Глобальна матриця 
Glo = np.zeros((N,N))  

# Права частина  
b_glo = np.zeros((N,1))

fig, (trax, plax) = plt.subplots(2, 1, sharex=True)

for i, p in enumerate(B['vertices']):
    trax.text(p[0], p[1] + .01, i, path_effects=[pe.withStroke(linewidth=2, foreground="white")])

for i, t in enumerate(B['triangles']):
    points = np.array([B['vertices'][j] for j in t])
    midpoint = points.mean(axis=0)
    trax.text(midpoint[0], midpoint[1], i, color='red')

tr.plot(trax, **B)

for t in B['triangles']:
    # t - масив індексів точок, t_points - масив точок трикутника
    # t = B['triangles'][0]
    t_points = list(map(lambda x: B['vertices'][x], t))


    #переробила по людськи гамму
    t_pairs = [[t[0], t[1]], [t[1], t[2]], [t[0], t[2]]]
    t_idx_on_segments = []
    for pair in t_pairs:
        for segment in segments:
            if (pair[0] == segment[0] and pair[1] == segment[1]) or (pair[0] == segment[1] and pair[1] == segment[0]):
                t_idx_on_segments.append(pair)

    if len(t_idx_on_segments) > 0:
        Gamma = sum([dist(B['vertices'][segment[0]], B['vertices'][segment[1]]) for segment in t_idx_on_segments])
    else:
        Gamma = 0


    bi, bj, bm = b[0](t_points), b[1](t_points), b[2](t_points)
    ci, cj, cm = c[0](t_points), c[1](t_points), c[2](t_points)

    Ke = 1 / (delta(t_points)**2) * np.array([
        [a11*bi**2 + a22*ci**2, a11*bi*bj + a22*ci*cj, a11*bi*bm + a22*ci*cm],
        [a11*bi*bj + a22*ci*cj, a11*bj**2 + a22*cj**2, a11*bj*bm + a22*cj*cm],
        [a11*bi*bm + a22*ci*cm, a11*bj*bm + a22*cj*cm, a11*bm**2 + a22*cm**2]
    ])

    Me = d * (delta(t_points) / 24)* np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2]
    ])

    Re = (sigma / beta)* (Gamma / 6) * np.array([
        [2, 1],
        [1, 2]
    ])

    Pe = (psi / beta) * (Gamma / 2)* np.array([
        [1],
        [1]
    ])

    Qe = (np.array(Me) / d) @ fe



    # Індекси вузлів 
    idx = t
    # Формування глобальної матриці 
    for i in range(3):
        for j in range(3):
            Glo[idx[i], idx[j]] += Ke[i, j] + Me[i, j]

    for segment in t_idx_on_segments:
        k, p = segment
        Glo[k, k] += Re[0, 0]
        Glo[k, p] += Re[0, 1]
        Glo[p, p] += Re[1, 1]
        Glo[p, k] += Re[1, 0]




    # накопичення правої частини
    for i in range(3):
        b_glo[idx[i]] += Qe[i]
    
    for segment in t_idx_on_segments:
        b_glo[segment[0]] += Pe.flatten()[0]
        b_glo[segment[1]] += Pe.flatten()[1]

mult = np.ones(Glo.shape)
mult[0,0] = 1e+3
Glo = Glo * mult
#Glo = ((np.diag(np.ones(N)) * (10**4 - 1)) + 1) * Glo

# Розв'язання СЛАР:  
u = np.linalg.solve(Glo, b_glo)


y = .3
trax.plot([0, 1], [y, y], 'g')

triangles = [np.array([B['vertices'][point_idx] for point_idx in triangle]) for triangle in B['triangles']]

points_x = np.linspace(0, 1, 100)
interpolation_results = []

for point_x in points_x:
    i, bc_point = find_triangle(triangles, np.array([point_x, y]))
    t = B['triangles'][i]
    u_in_triangle = u[[t[0], t[1], t[2]]]
    interpolation_result = bc_point.reshape(1, -1) @ u_in_triangle
    interpolation_results.append(interpolation_result.flatten()[0])

plax.plot(points_x, interpolation_results)
plt.show()

x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
z = np.zeros((100, 100))

for i, x_i in enumerate(x):
    for j, y_j in enumerate(y):
        i, bc_point = find_triangle(triangles, np.array([x_i, y_j]))
        t = B['triangles'][i]
        u_in_triangle = u[[t[0], t[1], t[2]]]
        interpolation_result = bc_point.reshape(1, -1) @ u_in_triangle
        z[i,j] = interpolation_result.flatten()[0]

x, y = np.meshgrid(x, y)

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

ls = LightSource(270, 45)
# To use a custom hillshading mode, override the built-in shading and pass
# in the rgb colors of the shaded surface calculated from "shade".
rgb = ls.shade(z, cmap=matplotlib.cm.gist_earth, vert_exag=0.1, blend_mode='soft')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=False, shade=False)

plt.show()

print("Глобальна матриця:")
print(pd.DataFrame(Glo))

print("\nВектор правої частини:")
print(b_glo)

print("\Р-к:")
print(u)

print(u.sum())


