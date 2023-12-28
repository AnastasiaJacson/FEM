import matplotlib.pyplot as plt
import numpy as np
import math

import triangle as tr

a11 = 1
a22 = 1
d = 1
#з граничної умови наприклад:
sigma = 1e5
beta = 1e-5
psi = 0
fe = 3


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


A = dict(vertices=[[0, 0], [1, 0], [1, 1], [0, 1]])
B = tr.triangulate(A, 'qa0.1c')
print(B['segments'])
#були масиві B['segments'] пари індексі, зробили пари точок
segments = list(map(lambda s: (B['vertices'][s[0]], B['vertices'][s[1]]), B['segments']))
# зробили з пар плоский масив бо по ньому зручніше шукати 
segments = np.array([point for segment in segments for point in segment])
# видаляємо дублі
segments = np.unique(segments, axis=0)
# елементи унікальні по осі рядків


tr.compare(plt, A, B)


for t in B['triangles'][:1]:
# t - масив індексів точок, t_points - масив точок трикутника
# t = B['triangles'][0]
    t_points = list(map(lambda x: B['vertices'][x], t))
    # print(t_points)

    #фільтруємо точки і вибираємо лише ті що є на сегменті
    t_points_on_segments = list(filter(lambda point: any(np.equal(segments, point).all(1)) , t_points))
    #якщо точток на сегменті менеше ніж одна то межа = 0
    if len(t_points_on_segments) < 2:
        Gamma = 0
    else:
        # в іншому випадку повертає гамма відстань між цтми двома точками
        Gamma = dist(t_points_on_segments[0], t_points_on_segments[1])
    
    # Gamma = sum(map(lambda p: dist([p[0]], B['vertices'][p[1]]), B['segments']))

    bi, bj, bm = b[0](t_points), b[1](t_points), b[2](t_points)
    ci, cj, cm = c[0](t_points), c[1](t_points), c[2](t_points)

    Ke = 1 / (delta(t_points)**2) * np.array([
        [a11*bi**2 + a22*ci**2, a11*bi*bj + a22*ci*cj, a11*bi*bm + a22*ci*cm],
        [a11*bi*bj + a22*ci*cj, a11*bj**2 + a22*cj**2, a11*bj*bm + a22*cj*cm],
        [a11*bi*bm + a22*ci*cm, a11*bj*bm + a22*cj*cm, a11*bm**2 + a22*cm**2]
    ])
    print(Ke)

    Me = d * (delta(t_points) / 24)* np.array([
        [2, 1, 1],
        [1, 2, 2],
        [1, 1, 2]
    ])
    print(Me)

    Re = (sigma / beta)* (Gamma / 6) * np.array([
        [2, 1],
        [1, 2]
    ])
    print(Re)

    Pe = (psi / beta) * (Gamma / 2)* np.array([
        [1],
        [1]
    ])
    print(Pe)

    Qe = fe* np.array(Me) / d
    print (Qe)


#проходимось по кожній точці трикутника t і малюємо її червоним 

# for t_point in t_points:
#     plt.plot(t_point[0], t_point[1], 'r.')
plt.show()