from dataclasses import dataclass, field
from typing import override, Dict

from quark.core import Core, Result, Data
from quark.interface_types import InterfaceType, Other

import matplotlib.pyplot as plt
import numpy as np
from quark_plugin_quantinuum.interfaces.backend_result import (  # type: ignore[import-untyped]
    BackendResult,
)
from .spectral_function_helpers import create_sequence_circuit
from .spectral_function_score import (
    extract_simulation_results,
    exact_values,
    compute_score,
)
from ..abstract_classes.abstract_circuit import (
    AbstractCircuits,
)
import logging

logger = logging.getLogger()


@dataclass
class SpectralFunction(Core):
    l_tot: int = 8

    dt: float = 0.5
    n_trot: int = 10
    n_shots: int = 200

    mode: str = "log2N"

    create_plot: bool = True

    epsilon: float = np.pi * 0.8 / (n_trot * dt)

    metrics: Dict[str, float | int] = field(init=False, default_factory=dict)

    @override
    def preprocess(self, data: InterfaceType) -> Result:
        if not (
            self.mode == "log2N"
            or self.mode == "CZ_pyramid"
            or self.mode == "CX_ladder"
        ):
            raise ValueError(
                "mode must be either log2N (default), CZ_pyramid or CX_ladder"
            )

        U_list = create_sequence_circuit(
            self.epsilon,
            6 / (self.l_tot - 1),
            self.l_tot,
            self.dt,
            self.n_trot,
            self.mode,
            self.l_tot,
        )

        return Data(
            Other(
                AbstractCircuits(
                    [
                        {
                            "circuit": U_list[n],
                            "number_of_qubits": self.l_tot * 2,
                            "shot_proportions": 1,
                            "n_shots": self.n_shots,
                        }
                        for n in range(len(U_list))
                    ],
                    self.benchmark_tag(),
                )
            )
        )

    @override
    def postprocess(self, input_data: Other[BackendResult]) -> Result:
        backend_result = input_data.data

        counts_per_circuit = backend_result.counts
        results = extract_simulation_results(self.l_tot, counts_per_circuit)

        exact = exact_values(
            self.l_tot,
            self.dt,
            self.n_trot,
            self.epsilon,
            6 / (self.l_tot - 1),
            self.l_tot,
        )

        fidelity, dispersion, dispersion_fidelity, dispersion_fidelity_exact = (
            compute_score(
                results[:, :, 0], exact, 6 / (self.l_tot - 1), self.l_tot, self.l_tot
            )
        )

        logger.info(f"Benchmark score (global fidelity): {fidelity}")
        logger.info(
            f"Benchmark score (noiseless dispersion fidelity): {dispersion_fidelity}"
        )
        logger.info(
            f"Benchmark score (exact dispersion fidelity): {dispersion_fidelity_exact}"
        )
        self.metrics.update(
            {
                "application_score_fidelity": fidelity,
                "application_score_noiseless_dispersion": dispersion_fidelity,
                "application_score_exact_dispersion": dispersion_fidelity_exact,
            }
        )

        if self.create_plot:
            # plt.figure()

            f, axarr = plt.subplots(1, 2)

            axarr[0].imshow(
                results[:, :, 0],
                aspect=self.l_tot / 6,
                extent=[0, self.l_tot, -3, 3],
                origin="lower",
                vmin=0,
                vmax=1,
            )
            axarr[0].plot(dispersion[:, 0], dispersion[:, 1], color="red")
            axarr[0].set_title("results")
            axarr[0].set_xlabel("k")
            axarr[0].set_ylabel("$\omega$")

            im1 = axarr[1].imshow(
                exact,
                aspect=self.l_tot / 6,
                extent=[0, self.l_tot, -3, 3],
                origin="lower",
                vmin=0,
                vmax=1,
            )
            axarr[1].plot(dispersion[:, 0], dispersion[:, 2], color="red")
            axarr[1].set_title("exact")
            axarr[1].set_xlabel("k")
            f.colorbar(im1, ax=f.get_axes(), shrink=0.5)

            plt.show()

        return Data(
            Other(
                {
                    "application_score_fidelity": fidelity,
                    "application_score_noiseless_dispersion": dispersion_fidelity,
                    "application_score_exact_dispersion": dispersion_fidelity_exact,
                }
            )
        )

    def get_metrics(self) -> dict:
        return self.metrics

    def benchmark_tag(self) -> str:
        return f"spectral_function_{self.l_tot}_{self.n_trot}_{self.dt}_{self.epsilon}_{self.mode}"
