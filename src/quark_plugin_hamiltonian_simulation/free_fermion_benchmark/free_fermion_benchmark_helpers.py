import logging

logger = logging.getLogger()


def coordinates(x: int, y: int, lx: int, ly: int, periodic: bool) -> int:
    """Coordinate of site (x, y)"""
    if not periodic:
        if x < 0 or x > lx - 1 or y < 0 or y > ly - 1:
            return -1  # no lattice site
    return (x % lx) + (y % ly) * lx


def coordinates_ancilla(x: int, y: int, lxp: int, lyp: int, periodic: bool) -> int:
    """Coordinate of ancilla (x, y)"""
    if periodic:
        return (x % lxp) + (y % lyp) * lxp
    else:
        if x < 0 or y < 0:
            return -1  # no ancilla
        elif y > lyp - 2:
            return -1
        elif y % 2 == 0 and x > lxp - 2:
            return -1
        elif y % 2 == 1 and x > lxp - 1:
            return -1
        else:
            return x + (y // 2) * (2 * lxp - 1) + (y % 2) * (lxp - 1)


def toric_code_preparation(
    u: list, lx: int, ly: int, periodic: bool, shift: int
) -> None:
    """toric code ground state preparation"""

    lxp = lx // 2
    lyp = ly

    for k in range(lyp // 2 - 1):  # initial Hadamards
        for j in range(lxp):
            f1 = coordinates_ancilla(j, lyp - 2 * k - 3, lxp, lyp, periodic)
            if f1 >= 0:
                u.append(["h", shift + f1])

    for k in range(lyp // 2 - 1):  # sequence of CXs
        for shifting in [[0, 2], [-1, 1], [0, 1]]:
            for j in range(lxp):
                f1 = coordinates_ancilla(j, lyp - 2 * k - 3, lxp, lyp, periodic)
                f2 = coordinates_ancilla(
                    j + shifting[0], lyp - 2 * k - 3 + shifting[1], lxp, lyp, periodic
                )
                if f1 >= 0 and f2 >= 0:
                    u.append(["cx", shift + f1, shift + f2])

    for j in range(lxp - 1):  # CXs on the last row
        f1 = coordinates_ancilla(lxp - j - 2, 0, lxp, lyp, periodic)
        f2 = coordinates_ancilla(lxp - j - 2 + 1, 1, lxp, lyp, periodic)
        f3 = coordinates_ancilla(lxp - j - 2 + 1, -1, lxp, lyp, periodic)
        f4 = coordinates_ancilla(lxp - j - 2 + 1, 0, lxp, lyp, periodic)

        if f1 >= 0:
            u.append(["h", shift + f1])
        if f1 >= 0 and f2 >= 0:
            u.append(["cx", shift + f1, shift + f2])
        if f1 >= 0 and f3 >= 0:
            u.append(["cx", shift + f1, shift + f3])
        if f1 >= 0 and f4 >= 0:
            u.append(["cx", shift + f1, shift + f4])

    for k in range(
        lyp
    ):  # change of basis to make all the stabilizers into identical XXYY strings
        for j in range(lxp):
            a = coordinates_ancilla(j, k, lxp, lyp, periodic)
            if a >= 0:
                if k % 2 == 1:
                    u.append(["sdg", shift + a])
                    u.append(["h", shift + a])
                if k % 2 == 0:
                    u.append(["s", shift + a])
                    u.append(["h", shift + a])
                    u.append(["s", shift + a])


def inverse_toric_code_preparation(
    u: list, lx: int, ly: int, periodic: bool, shift: int
) -> None:
    """inverse toric code ground state preparation"""

    lxp = lx // 2
    lyp = ly

    for k in range(lyp)[::-1]:
        for j in range(lxp)[::-1]:
            a = coordinates_ancilla(j, k, lxp, lyp, periodic)
            if a >= 0:
                if k % 2 == 1:
                    u.append(["h", shift + a])
                    u.append(["s", shift + a])
                if k % 2 == 0:
                    u.append(["sdg", shift + a])
                    u.append(["h", shift + a])
                    u.append(["sdg", shift + a])

    for j in range(lxp - 1)[::-1]:
        f1 = coordinates_ancilla(lxp - j - 2, 0, lxp, lyp, periodic)
        f2 = coordinates_ancilla(lxp - j - 2 + 1, 1, lxp, lyp, periodic)
        f3 = coordinates_ancilla(lxp - j - 2 + 1, -1, lxp, lyp, periodic)
        f4 = coordinates_ancilla(lxp - j - 2 + 1, 0, lxp, lyp, periodic)

        if f1 >= 0 and f4 >= 0:
            u.append(["cx", shift + f1, shift + f4])
        if f1 >= 0 and f3 >= 0:
            u.append(["cx", shift + f1, shift + f3])
        if f1 >= 0 and f2 >= 0:
            u.append(["cx", shift + f1, shift + f2])
        if f1 >= 0:
            u.append(["h", shift + f1])

    for k in range(lyp // 2 - 1)[::-1]:
        for shifting in [[0, 2], [-1, 1], [0, 1]][::-1]:
            for j in range(lxp)[::-1]:
                f1 = coordinates_ancilla(j, lyp - 2 * k - 3, lxp, lyp, periodic)
                f2 = coordinates_ancilla(
                    j + shifting[0], lyp - 2 * k - 3 + shifting[1], lxp, lyp, periodic
                )
                if f1 >= 0 and f2 >= 0:
                    u.append(["cx", shift + f1, shift + f2])

    for k in range(lyp // 2 - 1)[::-1]:
        for j in range(lxp)[::-1]:
            f1 = coordinates_ancilla(j, lyp - 2 * k - 3, lxp, lyp, periodic)
            if f1 >= 0:
                u.append(["h", shift + f1])


def plaquette_trotter_step(
    u: list,
    theta: float,
    x0: int,
    x1: int,
    x2: int,
    x3: int,
    a: int,
    shift: int,
    shift_ancilla: int,
) -> None:
    """Implements one plaquette hopping with coordinates x0,x1,x2,x3 and ancilla a"""
    if a >= 0 and x0 >= 0:
        u.append(["cx", shift + x0, shift_ancilla + a])
    if a >= 0 and x1 >= 0:
        u.append(["cx", shift + x1, shift_ancilla + a])
    if x0 >= 0 and x2 >= 0:
        u.append(["rxx", -theta / 2, shift + x0, shift + x2])
    if x0 >= 0 and x2 >= 0:
        u.append(["ryy", -theta / 2, shift + x0, shift + x2])
    if a >= 0 and x0 >= 0:
        u.append(["cx", shift + x0, shift_ancilla + a])
    if x1 >= 0 and x3 >= 0:
        u.append(["rxx", theta / 2, shift + x1, shift + x3])
    if x1 >= 0 and x3 >= 0:
        u.append(["ryy", theta / 2, shift + x1, shift + x3])
    if a >= 0 and x1 >= 0:
        u.append(["cx", shift + x1, shift_ancilla + a])

    if a >= 0:
        u.append(["sdg", shift_ancilla + a])
    if a >= 0 and x0 >= 0:
        u.append(["cx", shift + x0, shift_ancilla + a])
    if a >= 0 and x2 >= 0:
        u.append(["cx", shift + x2, shift_ancilla + a])
    if x0 >= 0 and x1 >= 0:
        u.append(["rxx", theta / 2, shift + x0, shift + x1])
    if x0 >= 0 and x1 >= 0:
        u.append(["ryy", theta / 2, shift + x0, shift + x1])
    if a >= 0 and x0 >= 0:
        u.append(["cx", shift + x0, shift_ancilla + a])
    if x2 >= 0 and x3 >= 0:
        u.append(["rxx", theta / 2, shift + x2, shift + x3])
    if x2 >= 0 and x3 >= 0:
        u.append(["ryy", theta / 2, shift + x2, shift + x3])
    if a >= 0 and x2 >= 0:
        u.append(["cx", shift + x2, shift_ancilla + a])
    if a >= 0:
        u.append(["s", shift_ancilla + a])


def trotter_step(
    u: list, dt: float, lx: int, ly: int, periodic: bool, two_spin_species: bool
) -> None:
    """implements one Trotter step on the system"""
    if two_spin_species:
        spin_species = 2
    else:
        spin_species = 1

    shift_ancilla = spin_species * lx * ly

    for stage_species in range(spin_species):  # on which spin species we act
        for stage in [0, 1]:  # whether 'P' or 'Q' plaquettes
            for k in range(
                -stage, ly // 2
            ):  # loops over plaquettes. Negative values are for the non-periodic case where edge sites don't form a complete face
                for j in range(-1 + stage, lx // 2):
                    if j >= 0 or (not periodic):
                        if k >= 0 or (not periodic):
                            a = coordinates_ancilla(
                                j, 2 * k + stage, lx // 2, ly, periodic
                            )
                            x0 = coordinates(
                                2 * j + 1 - stage, 2 * k + stage, lx, ly, periodic
                            )
                            x1 = coordinates(
                                2 * j + 2 - stage, 2 * k + stage, lx, ly, periodic
                            )
                            x2 = coordinates(
                                2 * j + 1 - stage, 2 * k + 1 + stage, lx, ly, periodic
                            )
                            x3 = coordinates(
                                2 * j + 2 - stage, 2 * k + 1 + stage, lx, ly, periodic
                            )

                            plaquette_trotter_step(
                                u,
                                dt,
                                x0,
                                x1,
                                x2,
                                x3,
                                a,
                                ((stage_species + stage) % 2)
                                * (spin_species - 1)
                                * lx
                                * ly,
                                shift_ancilla,
                            )
        if spin_species == 2:
            for j in range(lx * ly):
                u.append(["cz", j, j + lx * ly])


def droplet(lx: int, ly: int, periodic: bool) -> list:
    """initialize state with a 'droplet' of fermions in the center of the lattice"""
    toApply = []
    a = 1
    for j in range(ly):
        for i in range(a):
            toApply.append(coordinates(lx // 2 - 1 - i, j, lx, ly, periodic))
            toApply.append(coordinates(lx // 2 + i, j, lx, ly, periodic))
        if j < ly // 2 - 1 and j < lx // 2 - 1:
            a += 1
        else:
            a -= 1
    return toApply


def adjacent_ancilla(j0: int, j1: int, lx: int, ly: int, periodic: bool) -> int:
    """returns the ancilla adjacent to edge j0-j1"""
    x0 = j0 % lx
    x1 = j1 % lx
    y0 = j0 // lx
    y1 = j1 // lx

    if x1 == x0 + 1:
        a = coordinates_ancilla(
            x0 // 2, y0 - ((x0 + y0 + 1) % 2), lx // 2, ly, periodic
        )
    elif y1 == y0 + 1:
        a = coordinates_ancilla((x0 - 1 + (y0 % 2)) // 2, y0, lx // 2, ly, periodic)
    elif periodic and x1 == x0 - lx + 1:
        a = coordinates_ancilla(
            x0 // 2, y0 - ((x0 + y0 + 1) % 2), lx // 2, ly, periodic
        )
    elif periodic and y1 == y0 - ly + 1:
        a = coordinates_ancilla((x0 - 1 + (y0 % 2)) // 2, y0, lx // 2, ly, periodic)
    else:
        raise ValueError("vertices non-adjacent")

    return a


def change(s: str, j: int, new: str) -> str:
    """simple routine to change a character in a string"""
    return s[:j] + new + s[(j + 1) :]


def initialize_product_state(
    u: list, initial: list, lx: int, ly: int, periodic: bool, two_spin_species: bool
) -> None:
    """initialize the product state by putting fermions on the sites in initial"""
    if two_spin_species:  # if there are two spin species, the initialization is trivial and does not require acting on the ancillas as it commutes with the stabilizers
        for j in initial:
            u.append(["x", j])
            u.append(["x", j + lx * ly])
    else:  # otherwise, one must act on the stabilizers appropriately
        initial_editable = list(initial)

        while len(initial_editable) > 0:
            x0 = (
                initial_editable.pop()
            )  # take one vertex on which we must put a fermion
            lowest_distance = lx + ly
            lowest_index = -1
            for b in range(len(initial_editable)):  # look for the closest other fermion
                j = initial_editable[b]
                distance = abs(x0 // lx - j // lx) + abs((x0 % lx) - (j % lx))
                if distance < lowest_distance:
                    lowest_distance = distance
                    lowest_index = b

            x1 = initial_editable.pop(lowest_index)

            u.append(["x", x0])
            u.append(["x", x1])

            # putting two fermions creates two defects on the toric code, that must be brought together to be suppressed

            if x0 % lx > x1 % lx:
                x0 = x0 + x1
                x1 = x0 - x1
                x0 = x0 - x1
            while x0 % lx != x1 % lx:  # move one defect horizontally
                a = adjacent_ancilla(x0, x0 + 1, lx, ly, periodic)
                if a >= 0:
                    u.append(["y", lx * ly + a])
                x0 = x0 + 1

            if x0 // lx > x1 // lx:
                x0 = x0 + x1
                x1 = x0 - x1
                x0 = x0 - x1
            while x0 // lx != x1 // lx:  # move one defect vertically
                a = adjacent_ancilla(x0, x0 + lx, lx, ly, periodic)
                if a >= 0:
                    u.append(["x", lx * ly + a])
                x0 = x0 + lx


def create_circuit(
    lx: int, ly: int, dt: float, periodic: bool, two_spin_species: bool, n_trot: int
) -> list:
    """returns an abstract circuit in the form of a list of operations"""
    logger.info(f"Creating simulation circuit for {n_trot} Trotter steps")
    u: list[list[object]] = []

    if two_spin_species:
        shift_ancilla = 2 * lx * ly
    else:
        shift_ancilla = lx * ly

    toric_code_preparation(u, lx, ly, periodic, shift_ancilla)

    initial = droplet(lx, ly, periodic)

    initialize_product_state(u, initial, lx, ly, periodic, two_spin_species)

    for t in range(n_trot):
        trotter_step(u, dt, lx, ly, periodic, two_spin_species)

    u.append(["barrier"])

    inverse_toric_code_preparation(u, lx, ly, periodic, shift_ancilla)

    if two_spin_species:
        N = 2 * lx * ly
    else:
        N = lx * ly

    if periodic:
        N += lx * ly // 2
    else:
        N += (lx // 2 - 1) * ly // 2 + lx // 2 * (ly // 2 - 1)

    for j in range(N):
        u.append(["measure", j, j])

    return u
