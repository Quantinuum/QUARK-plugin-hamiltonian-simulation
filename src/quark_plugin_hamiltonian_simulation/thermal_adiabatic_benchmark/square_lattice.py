def coordinates(x, y, lx, ly):
    return (x % lx) + (y % ly) * lx


def createCouplingsSquare(lx, ly, periodic):
    l = lx * ly
    E = []

    for k in range(ly):
        for j in range(0, lx, 2):
            if j + 1 < lx or periodic:
                E.append([coordinates(j, k, lx, ly), coordinates(j + 1, k, lx, ly)])
    for k in range(0, ly, 2):
        for j in range(lx):
            if k + 1 < ly or periodic:
                E.append([coordinates(j, k, lx, ly), coordinates(j, k + 1, lx, ly)])
    for k in range(ly):
        for j in range(1, lx, 2):
            if j + 1 < lx or periodic:
                E.append([coordinates(j, k, lx, ly), coordinates(j + 1, k, lx, ly)])
    for k in range(1, ly, 2):
        for j in range(lx):
            if k + 1 < ly or periodic:
                E.append([coordinates(j, k, lx, ly), coordinates(j, k + 1, lx, ly)])

    return [[j for j in range(l)], E]
