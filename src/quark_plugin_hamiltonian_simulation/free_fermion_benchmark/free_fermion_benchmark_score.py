# Copyright Quantinuum
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import numpy as np
import numpy.typing as npt
import logging
from scipy.stats import chi2
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from .free_fermion_benchmark_helpers import droplet, coordinates
from .free_fermion_benchmark_stabilizers import stabilizers_after_toric_code

logger = logging.getLogger()

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]


class FreeFermionSolver:
    def __init__(self, j, k, s, lx, ly, bound_hor, bound_vert, n):
        self.j = j
        self.k = k
        self.s = s
        self.lx = lx
        self.ly = ly
        self.n = n
        self.sig = 1
        if self.j < 0 or self.k < 0:
            self.sig *= 0
        if j > k:
            abs_jk = abs(j - k)
            if bound_hor == 1 and abs_jk == lx - 1:
                self.sig *= -1
            if bound_vert == 1 and abs_jk >= lx:
                self.sig *= -1
            if bound_hor == -1 and abs_jk == lx - 1:
                self.sig *= 0
            if bound_vert == -1 and abs_jk >= lx:
                self.sig *= 0

    def dc(
        self, c: ComplexArray, d: ComplexArray
    ) -> ComplexArray:  # derivative of evolution of c=<c_i^\dagger c_j>
        deriv: ComplexArray = np.zeros((self.n, self.n), dtype=np.complex128)
        deriv[self.j, :] += 1j * (c[self.k, :] - self.s * d[self.k, :]) * self.sig
        deriv[self.k, :] += 1j * (c[self.j, :] + self.s * d[self.j, :]) * self.sig
        deriv[:, self.j] += (
            1j * (-c[:, self.k] + self.s * np.conj(d[self.k, :])) * self.sig
        )
        deriv[:, self.k] += (
            1j * (-c[:, self.j] - self.s * np.conj(d[self.j, :])) * self.sig
        )
        return deriv

    def dd(
        self, c: ComplexArray, d: ComplexArray
    ) -> ComplexArray:  # derivative of evolution of D=<c_i c_j>
        deriv: ComplexArray = np.zeros((self.n, self.n), dtype=np.complex128)
        deriv[self.j, :] += 1j * (-d[self.k, :] + self.s * c[self.k, :]) * self.sig
        deriv[self.k, :] += 1j * (-d[self.j, :] - self.s * c[self.j, :]) * self.sig
        deriv[:, self.j] += 1j * (-d[:, self.k] - self.s * c[self.k, :]) * self.sig
        deriv[self.k, self.j] += 1j * self.s * self.sig
        deriv[:, self.k] += 1j * (-d[:, self.j] + self.s * c[self.j, :]) * self.sig
        deriv[self.j, self.k] += -1j * self.s * self.sig
        return deriv

    def diff(
        self, t: float, cvec: ComplexArray
    ) -> ComplexArray:  # function to call for the differential equation
        """Function to pass to solve_ivp"""
        call: ComplexArray = cvec.reshape((2 * self.n, self.n))
        return (
            np.concatenate(
                (
                    self.dc(call[: self.n], call[self.n :]),
                    self.dd(call[: self.n], call[self.n :]),
                )
            )
        ).reshape(2 * self.n**2)


def exact_values_and_variance(
    n_trot: int, dt: float, lx: int, ly: int, periodic: bool, two_spin_species: bool
) -> FloatArray:
    """computes the exact output of the benchmark"""
    l_tot = lx * ly
    n = 2 * l_tot
    # index 0: number of steps; index 1: expectation value of imbalance; index
    # 2: expectation value of square of imbalance
    res: FloatArray = np.zeros((n_trot + 1, 3))

    if periodic:
        if two_spin_species:
            boundary_list = [
                [0, 0, 0],
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 0, 1],
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 1],
            ]
        else:
            boundary_list = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
    else:
        if two_spin_species:
            boundary_list = [[-1, -1, 0], [-1, -1, 1]]
        else:
            boundary_list = [[-1, -1, 0]]

    initial = droplet(lx, ly, periodic)

    for bb in boundary_list:  # loops over the boundary conditions and spin species
        boundary_vert = bb[0]
        boundary_hor = bb[1]
        spin_species = bb[2]

        c: ComplexArray = np.zeros((n, n), dtype=np.complex128)
        d: ComplexArray = np.zeros((n, n), dtype=np.complex128)

        for j in range(l_tot):  # initialize in the product state
            if j in initial:
                c[j, j] = 1

        res[0, 1] += -1
        res[0, 2] += 1

        order = [[0, 0, 0, 1], [1, 0, 1, 1], [0, 0, 1, 0], [0, 1, 1, 1]]
        f: list = [1 / l_tot if j in initial else -1 / l_tot for j in range(l_tot)]
        for t in range(n_trot):  # loop over Trotter steps
            for decal in [
                spin_species,
                1 - spin_species,
            ]:  # loop over the 2 plaquette types
                for k in range(-1, ly // 2):  # loop over vertical coordinate
                    for j in range(-1, lx // 2):  # loop over horizontal coordinate
                        if j >= 0 or (not periodic):
                            if k >= 0 or (not periodic):
                                for (
                                    o
                                ) in order:  # loop over the 4 edges of the plaquette
                                    jcur = coordinates(
                                        2 * j + decal + o[0] + 1,
                                        2 * k + decal + o[1],
                                        lx,
                                        ly,
                                        periodic,
                                    )
                                    kcur = coordinates(
                                        2 * j + decal + o[2] + 1,
                                        2 * k + decal + o[3],
                                        lx,
                                        ly,
                                        periodic,
                                    )
                                    solver = FreeFermionSolver(
                                        jcur,
                                        kcur,
                                        0,
                                        lx,
                                        ly,
                                        boundary_hor,
                                        boundary_vert,
                                        n,
                                    )
                                    cc: ComplexArray = (
                                        solve_ivp(
                                            solver.diff,
                                            [0, dt],
                                            np.concatenate((c, d)).reshape(2 * n**2),
                                            atol=1e-9,
                                            rtol=1e-9,
                                        ).y
                                    )[:, -1].reshape((2 * n, n))
                                    c = cc[:n]
                                    d = cc[n:]
            a: float = np.sum([f[j] * (1 - 2 * c[j, j]) for j in range(l_tot)])
            var: float = 4 * np.sum(
                [
                    f[i] * f[j] * c[i, i] * c[j, j]
                    for i in range(l_tot)
                    for j in range(l_tot)
                ]
            )
            var += -4 * np.sum(
                [
                    f[i] * f[j] * c[i, j] * c[j, i]
                    for i in range(l_tot)
                    for j in range(l_tot)
                ]
            )
            var += 4 * np.sum([f[i] ** 2 * c[i, i] for i in range(l_tot)])
            var += 4 * np.sum(
                [
                    f[i] * f[j] * abs(d[i, j]) ** 2
                    for i in range(l_tot)
                    for j in range(l_tot)
                ]
            )
            var += 2 * np.sum(f) * np.real(a)
            res[t + 1, 0] += t + 1
            res[t + 1, 1] += np.real(a)
            res[t + 1, 2] += np.real(var) - np.real(a) ** 2

    res = res / len(boundary_list)
    res[:, 2] = np.sqrt(res[:, 2])  # standard deviation per shot
    return res


def value_shot(s: str, l_tot: int, initial: list[int]) -> float:
    a: float = 0
    for j in range(l_tot):
        if j in initial:
            if s[-1 - j] == "1":
                a += -1 / l_tot
            else:
                a += 1 / l_tot
        else:
            if s[-1 - j] == "1":
                a += 1 / l_tot
            else:
                a += -1 / l_tot
    return a


def extract_simulation_results(
    dt: float,
    lx: int,
    ly: int,
    counts_per_circuit: list[dict[str, int]],
    periodic: bool,
    two_spin_species: bool,
) -> tuple[
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
    list[list[float]],
]:
    """Returns the simulation results.

    For every time step returns the time, expectation value and standard deviation as a tuple for that step.
    """
    l_tot = lx * ly

    initial = droplet(lx, ly, True)

    Safter = stabilizers_after_toric_code(
        lx, ly, periodic, two_spin_species
    )  # generates the Pauli strings of the stabilizers after toric code preparation. These must be only Z strings.

    results = []
    results_post = []
    stabilizers = []
    for n, counts in enumerate(counts_per_circuit):
        res: float = 0
        var: float = 0
        n_shots = np.sum(list(counts.values()))
        for s in counts:  # computes the expectation value of the observable
            a = value_shot(s, l_tot, initial)
            res += a * counts[s]
            var += a**2 * counts[s]
        res = res / n_shots
        var = var / n_shots

        results.append(
            (
                dt * n,
                res,
                np.sqrt(max(var - res**2, 0.000000000001)) / np.sqrt(n_shots - 1),
            )
        )

        stabilizer_count = [dt * n] + [0] * (
            len(Safter) + 1
        )  # now computes the value of the stabilizers for each shot
        res_post: float = 0
        var_post: float = 0
        n_shots_post = 0
        for s in counts:
            wrong_stabi = 0
            for a in range(
                len(Safter)
            ):  # computes the value of all stabilizers in shot s, and count the number of wrong stabilizers
                coeff = Safter[a][1]
                for j in range(len(s)):
                    if Safter[a][0][j] == "Z" and s[-1 - j] == "1":
                        coeff *= -1
                if coeff == -1:
                    wrong_stabi += 1
            stabilizer_count[1 + wrong_stabi] += counts[s] / n_shots

            if (
                wrong_stabi == 0
            ):  # computes the expectation value of the observable when post-selecting onto 0 wrong stabilizers
                a = value_shot(s, l_tot, initial)
                res_post += a * counts[s]
                var_post += a**2 * counts[s]
                n_shots_post += counts[s]

        if n_shots_post > 1:
            res_post = res_post / n_shots_post
            var_post = var_post / n_shots_post
            results_post.append(
                (
                    dt * n,
                    res_post,
                    np.sqrt(max(var_post - res_post**2, 0.000000000001))
                    / np.sqrt(n_shots_post - 1),
                )
            )

        stabilizers.append(stabilizer_count)

    return results, results_post, stabilizers


def computes_score_values(
    delta: FloatArray, std_exp: FloatArray, std: FloatArray, l_tot: int
) -> tuple[int, int, int]:
    """Computes score values.

    Returns the score in terms of
     1.) Number of gates
     2.) Number of shots
     3.) Number of Trotter steps
    """
    n: int = len(delta)
    delta_corrected = np.zeros(n)
    for j in range(n):
        if std[j] > 0:
            delta_corrected[j] = max(abs(delta[j]), std_exp[j]) / std[j]
            # redefines delta as the maximum between the experimental standard deviation and the measured value,
            # normalized by the theoretical standard deviation
    rewards: float = float(delta_corrected[0]) ** 2
    opt: int = 0
    for j in range(1, n):  # looks for the time point opt with maximal reward
        temp: float = float(delta_corrected[j]) ** 2 / (j + 1)
        if temp > rewards:
            rewards = temp
            opt = j

    def ff(x):
        return chi2.cdf(delta_corrected[opt] ** 2 * x, df=n - 1) - 0.997

    x: float = float(fsolve(ff, n / delta_corrected[opt] ** 2)[0])
    # looks for x such that chi2.cdf(delta[opt]**2*x*L,df=1)=0.997

    return (
        6 * int(np.floor(x) + 1) * (opt + 1) * l_tot,
        int(np.floor(x) + 1),
        int(np.floor(x) + 1) * (opt + 1),
    )


def create_and_show_plot(
    n_trot: int,
    dt: float,
    simulation_results,
    simulation_results_post,
    exact_results,
    score_gates: int,
) -> None:
    plt.plot(
        np.array(list(range(n_trot + 1))) * dt,
        exact_results[:, 1],
        color="black",
        label="exact",
    )
    plt.errorbar(
        simulation_results[:, 0],
        simulation_results[:, 1],
        yerr=simulation_results[:, 2],
        label="simulated",
    )
    plt.errorbar(
        simulation_results_post[:, 0],
        simulation_results_post[:, 1],
        yerr=simulation_results_post[:, 2],
        label="post-selected",
    )
    plt.title("SCORE = " + str(score_gates) + " gates")
    plt.xlabel("Time")
    plt.ylabel("Imbalance")
    plt.legend()
    plt.show()
    plt.close()


def create_and_show_histogram(stabilizers_results):
    for s in range(len(stabilizers_results)):
        stabilizers = stabilizers_results[s]
        plt.stairs(
            edges=[j - 0.5 for j in range(len(stabilizers))],
            values=stabilizers[1:],
            fill=True,
            color="C" + str(s),
        )
        plt.xlabel("Wrong stabilizers")
        plt.ylabel("Probability")
        plt.title("t=" + str(stabilizers[0]))
        plt.show()
        plt.close()
