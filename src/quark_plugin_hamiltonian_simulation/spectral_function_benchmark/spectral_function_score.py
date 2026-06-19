import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp

FloatArray = npt.NDArray[np.float64]


def extract_simulation_results(
    l_tot: int, counts_per_circuit: list[dict[str, int]]
) -> FloatArray:
    """compute A(k,omega) from the measurement outcomes"""

    results: FloatArray = np.zeros((l_tot, l_tot, 4))

    for omega in range(l_tot):
        for filling in range(2):
            NS = 0
            counts = counts_per_circuit[2 * omega + filling]
            for s in counts:
                NS += counts[s]
                for k in range(l_tot):
                    if filling == 0 and s[-1 - k - l_tot] == "1":
                        results[omega, k, 0] += counts[s]
                    if filling == 1 and s[-1 - k - l_tot] == "0":
                        results[omega, k, 2] += counts[s]
            results[omega, :, 2 * filling] = results[omega, :, 2 * filling] / NS
            results[omega, :, 2 * filling + 1] = np.sqrt(
                results[omega, :, 2 * filling] - results[omega, :, 2 * filling] ** 2
            ) / np.sqrt(NS - 1)

    results2: FloatArray = np.zeros((l_tot, l_tot, 2))

    results2[:, :, 0] = results[:, :, 0] + results[:, :, 2]
    results2[:, :, 1] = np.sqrt(results[:, :, 1] ** 2 + results[:, :, 3] ** 2)

    return results2


def exact_values(
    L: int, dt: float, Ntrot: int, epsilon: float, domega: float, n_omega: int
) -> FloatArray:
    """compute the exact values obtained for a noiseless implementation of the circuit"""

    def dC(C):
        deriv = np.zeros((N, N)) * 1j

        if mode == 0:
            parity = 0
        if mode == 2:
            parity = 1

        if mode == 0 or mode == 2:
            for j in range(L):
                for k in range(L):
                    if j % 2 == parity and j + 1 < L:
                        deriv[j, k] += -1j * (C[(j + 1) % L, k])
                    if j % 2 == 1 - parity and j - 1 > -1:
                        deriv[j, k] += -1j * (C[(j - 1) % L, k])
                    if k % 2 == parity and k + 1 < L:
                        deriv[j, k] += -1j * (-C[j, (k + 1) % L])
                    if k % 2 == 1 - parity and k - 1 > -1:
                        deriv[j, k] += -1j * (-C[j, (k - 1) % L])
                    if k % 2 == parity and k + 1 < L:
                        deriv[j + L, k] += -1j * (-C[j + L, (k + 1) % L])
                    if k % 2 == 1 - parity and k - 1 > -1:
                        deriv[j + L, k] += -1j * (-C[j + L, (k - 1) % L])
                    if j % 2 == parity and j + 1 < L:
                        deriv[j, k + L] += -1j * (C[(j + 1) % L, k + L])
                    if j % 2 == 1 - parity and j - 1 > -1:
                        deriv[j, k + L] += -1j * (C[(j - 1) % L, k + L])

        if mode == 1:
            for j in range(L):
                for k in range(L):
                    if j + 1 == L:
                        deriv[j, k] += -1j * (C[(j + 1) % L, k])
                    if j - 1 == -1:
                        deriv[j, k] += -1j * (C[(j - 1) % L, k])
                    if k + 1 == L:
                        deriv[j, k] += -1j * (-C[j, (k + 1) % L])
                    if k - 1 == -1:
                        deriv[j, k] += -1j * (-C[j, (k - 1) % L])
                    if k + 1 == L:
                        deriv[j + L, k] += -1j * (-C[j + L, (k + 1) % L])
                    if k - 1 == -1:
                        deriv[j + L, k] += -1j * (-C[j + L, (k - 1) % L])
                    if j + 1 == L:
                        deriv[j, k + L] += -1j * (C[(j + 1) % L, k + L])
                    if j - 1 == -1:
                        deriv[j, k + L] += -1j * (C[(j - 1) % L, k + L])

        if mode == 3:
            for j in range(L):
                for k in range(L):
                    deriv[j, k] += -1j * (
                        epsilon / 2 * C[L + j, k] - epsilon / 2 * C[j, k + L]
                    )
                    deriv[j + L, k] += -1j * (
                        epsilon / 2 * C[j, k] - epsilon / 2 * C[j + L, k + L]
                    )
                    deriv[j, k + L] += -1j * (
                        epsilon / 2 * C[j + L, k + L] - epsilon / 2 * C[j, k]
                    )
                    deriv[j + L, k + L] += -1j * (
                        epsilon / 2 * C[j, k + L] - epsilon / 2 * C[j + L, k]
                    )

        if mode == 4:
            for j in range(L):
                for k in range(L):
                    deriv[j + L, k] += -1j * (omega * C[j + L, k])
                    deriv[j, k + L] += -1j * (-omega * C[j, k + L])

        return deriv

    def diff(t, Cvec):
        Call = Cvec.reshape((N, N))
        return (dC(Call)).reshape(N**2)

    def modes(C):
        e = np.zeros(L) * 1j
        for k in range(L):
            for j in range(L):
                for ell in range(L):
                    e[k] += (
                        C[L + j, L + (j + ell) % L]
                        * np.exp(-1j * 2 * np.pi * (k + 0.0) / L * ell)
                        / L
                    )

        return np.real(e)

    N = 2 * L

    toPlot = []
    for omega0 in range(n_omega):
        omega = (omega0 - (n_omega - 1) / 2) * domega

        C = np.zeros((N, N)) * 1j
        K = [j for j in range(L)]

        for j in range(L):
            for ell in range(L):
                for k in K:
                    C[j, (j + ell) % L] += (
                        np.exp(1j * 2 * np.pi * (k + 0.0) / L * ell) / L
                    )

        for nn in range(Ntrot):
            for mm in [0, 1, 2, 3, 4]:
                mode = mm
                Cc = (
                    solve_ivp(diff, [0, dt], C.reshape(N**2), atol=1e-8, rtol=1e-8).y
                )[:, -1].reshape((N, N))
                C = Cc
        toPlot.append(modes(Cc))
    return np.array(toPlot, dtype=float)


def compute_score(
    result: FloatArray, exact: FloatArray, domega: float, n_omega: int, l_tot: int
) -> tuple[float, FloatArray, float, float]:
    fidelity = np.sum(np.array(result) * np.array(exact)) / np.sqrt(
        np.sum(np.array(result) ** 2) * np.sum(np.array(exact) ** 2)
    )

    dispersion: list[list[float]] = []
    omega_array = np.array(
        [domega * (omega - (n_omega - 1) / 2) for omega in range(n_omega)]
    )
    for k in range(l_tot):
        omega = np.sum(result[:, k] ** 4 * omega_array) / np.sum(result[:, k] ** 4)
        omega2 = np.sum(np.array(exact)[:, k] ** 4 * omega_array) / np.sum(
            np.array(exact)[:, k] ** 4
        )
        dispersion.append([k + 0.5, omega, omega2])
    dispersion_array: FloatArray = np.array(dispersion, dtype=float)
    dispersion_fidelity = np.sum(
        dispersion_array[:, 1] * dispersion_array[:, 2]
    ) / np.sqrt(
        np.sum(dispersion_array[:, 1] ** 2) * np.sum(dispersion_array[:, 2] ** 2)
    )

    exact_dispersion = 2 * np.cos(
        2 * np.pi * np.array([k for k in range(l_tot)]) / l_tot
    )
    dispersion_fidelity_exact = np.sum(
        dispersion_array[:, 1] * exact_dispersion
    ) / np.sqrt(np.sum(dispersion_array[:, 1] ** 2) * np.sum(exact_dispersion**2))

    return (
        float(fidelity),
        dispersion_array,
        float(dispersion_fidelity),
        float(dispersion_fidelity_exact),
    )
