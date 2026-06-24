import numpy as np
from .square_lattice import create_couplings_square


def extract_simulation_results_square_ising(
    lx: int, ly: int, periodic: bool, h: float, counts_per_circuit: list[dict[str, int]]
) -> list[float]:
    """computes expectation value of energy or X for all circuits"""
    L = create_couplings_square(lx, ly, periodic)
    V = L[0]
    E = L[1]

    def meanZ(counts: dict[str, int]) -> tuple[float, float]:
        res: float = 0
        var: float = 0
        n_shots = np.sum(list(counts.values()))
        for s in counts:
            a: float = 0
            for j in V:
                if s[j] == "1":
                    a += -1
                else:
                    a += 1
            res += a * counts[s]
            var += a**2 * counts[s]
        res = res / n_shots
        var = np.sqrt(var / n_shots - res**2) / np.sqrt(n_shots - 1)
        return res, var

    def meanZZ(counts: dict[str, int]) -> tuple[float, float]:
        res: float = 0
        var: float = 0
        n_shots = np.sum(list(counts.values()))
        for s in counts:
            a: float = 0
            for e in E:
                if s[e[0]] == s[e[1]]:
                    a += 1
                else:
                    a += -1
            res += a * counts[s]
            var += a**2 * counts[s]
        res = res / n_shots
        var = np.sqrt(var / n_shots - res**2) / np.sqrt(n_shots - 1)
        return res, var

    temp, var_temp = meanZZ(counts_per_circuit[0])
    energy = temp
    var_energy = var_temp**2

    temp, var_temp = meanZ(counts_per_circuit[1])
    energy += h * temp
    var_energy += h**2 * var_temp**2
    var_energy = np.sqrt(var_energy)

    energy_excited: float = 0
    var_energy_excited: float = 0
    for j in V:
        temp, var_temp = meanZZ(counts_per_circuit[2 * j])
        energy_excited += temp / len(V)
        var_energy_excited += var_temp**2 / len(V)

        temp, var_temp = meanZ(counts_per_circuit[2 * j + 1])
        energy_excited += h * temp / len(V)
        var_energy_excited += h**2 * var_temp**2 / len(V)

    var_energy_excited = np.sqrt(var_energy_excited)

    mirror, var_mirror = meanZ(counts_per_circuit[2 * len(V) + 2])

    mirror_excited: float = 0
    var_mirror_excited: float = 0
    for j in V:
        temp, var_temp = meanZ(counts_per_circuit[2 * len(V) + 2 + j])

        mirror_excited += temp / len(V)
        var_mirror_excited += var_temp**2 / len(V)
    var_mirror_excited = np.sqrt(var_mirror_excited)

    return [
        energy / len(V),
        var_energy / len(V),
        energy_excited / len(V),
        var_energy_excited / len(V),
        mirror / len(V),
        var_mirror / len(V),
        mirror_excited / len(V),
        var_mirror_excited / len(V),
    ]


def thermodynamics(
    results: list[float],
) -> tuple[float, float, float, float, float, float, float, float]:
    energy = results[0]
    var_energy = results[1]

    x = (1 + results[4]) / 2 + 0.0000000001
    entropy = -x * np.log(x) - (1 - x) * np.log(1 - x)  # entropy function
    var_entropy = np.log((1 - x) / x) * results[5] / 2  # variance of entropy measured

    def f(e: float, ep: float, m: float, mp: float) -> float:
        return (
            2 * (e - ep) / (m - mp) / np.log((1 - m) / (1 + m))
        )  # function computing the temperature

    values: list[float] = []
    for n in range(10000):
        xi1 = np.random.normal()
        xi2 = np.random.normal()
        xi3 = np.random.normal()
        xi4 = np.random.normal()
        value = f(
            results[0] + results[1] * xi1,
            results[2] + results[3] * xi2,
            min(max(results[4] + results[5] * xi3, -0.99999999), 0.99999999),
            min(max(results[6] + results[7] * xi4, -0.99999999), 0.99999999),
        )
        values.append(max(value, 0))

    temperature = float(np.mean(values))
    var_temperature = float(np.sqrt(np.var(values)))

    free_energy = float(energy - temperature * entropy)  # free energy of the system
    var_free_energy = float(
        np.sqrt(
            var_energy**2
            + temperature**2 * var_entropy**2
            + var_temperature**2 * entropy**2
        )
    )

    return (
        energy,
        var_energy,
        float(entropy),
        float(var_entropy),
        temperature,
        var_temperature,
        free_energy,
        var_free_energy,
    )
