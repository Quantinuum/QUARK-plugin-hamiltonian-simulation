from .free_fermion_benchmark_helpers import (
    adjacent_ancilla,
    toric_code_preparation,
    change,
)


def stabilizers(lx: int, ly: int, periodic: bool, two_spin_species: bool) -> list:
    """returns the list of stabilizers under the form of Pauli strings"""
    S = []

    if two_spin_species:
        N = 2 * lx * ly
        shift_ancilla = 2 * lx * ly
    else:
        N = lx * ly
        shift_ancilla = lx * ly

    if periodic:
        N += lx * ly // 2
    else:
        N += (lx // 2 - 1) * ly // 2 + lx // 2 * (ly // 2 - 1)

    for y in range(ly):  # loops over all odd faces of the square lattice
        for x in range(lx // 2):
            if periodic or (y < ly - 1 and (x + y % 2) < lx // 2):
                # starts creating new stabilizer
                s = "I" * N
                x0 = (2 * x + y % 2) % lx + y * lx  # the four sites of the face
                x1 = (2 * x + y % 2 + 1) % lx + y * lx
                x2 = (2 * x + y % 2) % lx + ((y + 1) % ly) * lx
                x3 = (2 * x + y % 2 + 1) % lx + ((y + 1) % ly) * lx
                s = change(s, x0, "Z")
                s = change(s, x1, "Z")
                s = change(s, x2, "Z")
                s = change(s, x3, "Z")
                if (
                    two_spin_species
                ):  # if two species are present, couples to other other species as well
                    s = change(s, x0 + lx * ly, "Z")
                    s = change(s, x1 + lx * ly, "Z")
                    s = change(s, x2 + lx * ly, "Z")
                    s = change(s, x3 + lx * ly, "Z")
                a = adjacent_ancilla(
                    x0, x1, lx, ly, periodic
                )  # adds ancilla adjacent to x0-x1, etc
                if a >= 0:
                    if s[a + shift_ancilla] == "I":
                        s = change(s, a + shift_ancilla, "Y")
                    else:
                        s = change(s, a + shift_ancilla, "I")
                a = adjacent_ancilla(x2, x3, lx, ly, periodic)
                if a >= 0:
                    if s[a + shift_ancilla] == "I":
                        s = change(s, a + shift_ancilla, "Y")
                    else:
                        s = change(s, a + shift_ancilla, "I")
                a = adjacent_ancilla(x0, x2, lx, ly, periodic)
                if a >= 0:
                    if s[a + shift_ancilla] == "I":
                        s = change(s, a + shift_ancilla, "X")
                    else:
                        s = change(s, a + shift_ancilla, "I")
                a = adjacent_ancilla(x1, x3, lx, ly, periodic)
                if a >= 0:
                    if s[a + shift_ancilla] == "I":
                        s = change(s, a + shift_ancilla, "X")
                    else:
                        s = change(s, a + shift_ancilla, "I")
                S.append(s)
    return S


def evolve(string: str, operation: list) -> tuple[str, float]:
    """evolve a Pauli string through a list of Clifford operations"""
    name = operation[0]
    if name == "x":
        j = operation[1]
        result = str(string)
        if string[j] == "I" or string[j] == "X":
            coeff = 1
        else:
            coeff = -1
    elif name == "y":
        j = operation[1]
        result = str(string)
        if string[j] == "I" or string[j] == "Y":
            coeff = 1
        else:
            coeff = -1
    elif name == "z":
        j = operation[1]
        result = str(string)
        if string[j] == "I" or string[j] == "Z":
            coeff = 1
        else:
            coeff = -1
    elif name == "s":
        j = operation[1]
        if string[j] == "I" or string[j] == "Z":
            result = str(string)
            coeff = 1
        elif string[j] == "X":
            result = change(string, j, "Y")
            coeff = -1
        elif string[j] == "Y":
            result = change(string, j, "X")
            coeff = 1
    elif name == "sdg":
        j = operation[1]
        if string[j] == "I" or string[j] == "Z":
            result = str(string)
            coeff = 1
        elif string[j] == "X":
            result = change(string, j, "Y")
            coeff = 1
        elif string[j] == "Y":
            result = change(string, j, "X")
            coeff = -1
    elif name == "h":
        j = operation[1]
        if string[j] == "I":
            result = str(string)
            coeff = 1
        elif string[j] == "X":
            result = change(string, j, "Z")
            coeff = 1
        elif string[j] == "Y":
            result = change(string, j, "Y")
            coeff = -1
        elif string[j] == "Z":
            result = change(string, j, "X")
            coeff = 1
    elif name == "cx":
        j = operation[1]
        k = operation[2]
        result = str(string)
        coeff = 1
        if not (
            (string[j] == "I" or string[j] == "Z")
            and (string[k] == "I" or string[k] == "X")
        ):
            if string[j] == "I":
                if string[k] == "Y":
                    result = change(string, j, "Z")
                    coeff = 1
                elif string[k] == "Z":
                    result = change(string, j, "Z")
                    coeff = 1
            elif string[j] == "X":
                if string[k] == "I":
                    result = change(string, k, "X")
                    coeff = 1
                elif string[k] == "X":
                    result = change(string, k, "I")
                    coeff = 1
                elif string[k] == "Y":
                    result = change(change(string, j, "Y"), k, "Z")
                    coeff = 1
                elif string[k] == "Z":
                    result = change(change(string, j, "Y"), k, "Y")
                    coeff = -1
            elif string[j] == "Y":
                if string[k] == "I":
                    result = change(string, k, "X")
                    coeff = 1
                elif string[k] == "X":
                    result = change(string, k, "I")
                    coeff = 1
                elif string[k] == "Y":
                    result = change(change(string, j, "X"), k, "Z")
                    coeff = -1
                elif string[k] == "Z":
                    result = change(change(string, j, "X"), k, "Y")
                    coeff = 1
            elif string[j] == "Z":
                if string[k] == "Y":
                    result = change(string, j, "I")
                    coeff = 1
                elif string[k] == "Z":
                    result = change(string, j, "I")
                    coeff = 1
    else:
        raise ValueError("non-Clifford gate encountered")
    return result, coeff


def stabilizers_after_toric_code(
    lx: int, ly: int, periodic: bool, two_spin_species: bool
) -> list:
    """computes the list of stabilizers after applying the toric code state preparation"""
    S = stabilizers(
        lx, ly, periodic, two_spin_species
    )  # list of stabilizers before toric code state preparation
    U = []
    if two_spin_species:
        shift = 2 * lx * ly
    else:
        shift = lx * ly

    toric_code_preparation(
        U, lx, ly, periodic, shift
    )  # U is the circuit that prepares the toric code

    Safter = []

    for s in S:
        safter = str(s)
        coeff = 1
        for g in U[::-1]:
            safter, multiplier = evolve(
                safter, g
            )  # computes the Pauli string obtained after toric code preparation for every stabilizer
            coeff *= multiplier
        Safter.append([safter, coeff])

    return Safter
