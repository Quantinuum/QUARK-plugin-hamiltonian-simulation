import numpy as np
from .square_lattice import createCouplingsSquare


def createScheduling(
    partition_start: list, partition_end: list, n_trot: int, dt: float
) -> list:
    """returns the list of times used for the adiabatic evolution"""

    def varphi(x):
        return (1 + np.tanh(np.tan((x - 0.5) * np.pi))) / 2

    T = []
    n = len(partition_start)

    if partition_end.count(0) == n - 1:
        decrease_dt = False
    else:
        decrease_dt = True

    for t in range(n_trot):
        a = [0] * n
        for j in range(n):
            a[j] = (
                partition_start[j]
                + varphi((t + 0.5) / n_trot) * (partition_end[j] - partition_start[j])
            ) * dt
            if decrease_dt:
                a[j] *= (1 - (t + 0.5) / n_trot) ** 0.5
        T.append(a)

    return T


def createCircuitSquareIsing(
    lx: int, ly: int, h: float, periodic: bool, n_trot: int, dt: float
) -> tuple[list, list]:
    """returns the list of all the circuits to run and the number of shots to perform"""
    L = createCouplingsSquare(lx, ly, periodic)  # couplings of square lattice
    V = L[0]
    E = L[1]
    T = createScheduling([1, 0], [h, 1], n_trot, dt)  # angles of gates

    U = []
    for j in V:
        U.append(["h", j])
        U.append(["z", j])

    for t in T:  # time evolution of ising
        for e in E:
            U.append(["rzz", t[1], e[0], e[1]])
        for j in V:
            U.append(["rx", t[0], j])

    Uzz = list(U)
    Ux = list(U)

    for j in V:
        Ux.append(["h", j])
        Ux.append(["measure", j, j])
        Uzz.append(["measure", j, j])

    resU = [Uzz, Ux]  # whether we measure ZZ or X

    shot_proportions = [
        abs(h) / (1 + abs(h)),
        1 / (1 + abs(h)),
    ]  # proportion of number of shots for every circuit

    for j in V:  # appends a X on site j at the beginning of every circuit
        resU.append([["x", j]] + Uzz)
        resU.append([["x", j]] + Ux)
        shot_proportions = shot_proportions + [
            abs(h) / (1 + abs(h)) / len(V),
            1 / (1 + abs(h)) / len(V),
        ]

    U_mirror = []  # mirror circuit
    for j in V:
        U_mirror.append(["h", j])
        U_mirror.append(["z", j])

    for t in T[: len(T) // 2]:
        for e in E:
            U_mirror.append(["rzz", t[1], e[0], e[1]])
        for j in V:
            U_mirror.append(["rx", t[0], j])
    for t in T[: len(T) // 2][::-1]:
        for j in V:
            U_mirror.append(["rx", -t[0], j])
        for e in E:
            U_mirror.append(["rzz", -t[1], e[0], e[1]])

    for j in V:
        U_mirror.append(["h", j])
        U_mirror.append(["measure", j, j])

    resU.append(U_mirror)
    shot_proportions.append(1)

    for j in V:  # mirror circuit with X on site j at the beginning
        resU.append([["x", j]] + U_mirror)
        shot_proportions.append(1 / len(V))

    return resU, shot_proportions
