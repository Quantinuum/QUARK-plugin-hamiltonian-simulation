from quark.core import Core, Result, Data
from quark.interface_types import InterfaceType, Other
from dataclasses import dataclass, field
from typing import override, Dict
from quark_plugin_quantinuum.interfaces.backend_result import (  # type: ignore[import-untyped]
    BackendResult,
)
from .thermal_adiabatic_helpers import createCircuitSquareIsing
from .thermal_adiabatic_score import (
    extract_simulation_results_square_ising,
    thermodynamics,
)
from quark_plugin_hamiltonian_simulation.abstract_classes.abstract_circuit import (
    AbstractCircuits,
)
import logging

logger = logging.getLogger()


@dataclass
class ThermalAdiabaticSquareIsing(Core):
    lx: int = 2
    ly: int = 2
    h: float = 1
    periodic: bool = False

    dt: float = 0.2
    n_trot: int = 2
    n_shots: int = 200

    metrics: Dict[str, float | int] = field(init=False, default_factory=dict)

    @override
    def preprocess(self, data: InterfaceType) -> Result:
        U_list, shot_proportions = createCircuitSquareIsing(
            self.lx, self.ly, self.h, self.periodic, self.n_trot, self.dt
        )

        return Data(
            Other(
                AbstractCircuits(
                    [
                        {
                            "circuit": U_list[n],
                            "number_of_qubits": self.lx * self.ly,
                            "shot_proportions": shot_proportions[n],
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
        results = extract_simulation_results_square_ising(
            self.lx, self.ly, self.periodic, self.h, counts_per_circuit
        )

        (
            energy,
            var_energy,
            entropy,
            var_entropy,
            temperature,
            var_temperature,
            free_energy,
            var_free_energy,
        ) = thermodynamics(results)

        logger.info(f"Benchmark score (energy): {energy} +/- {var_energy}")
        logger.info(f"Benchmark score (entropy): {entropy} +/- {var_entropy}")
        logger.info(
            f"Benchmark score (temperature): {temperature} +/- {var_temperature}"
        )
        logger.info(
            f"Benchmark score (free energy): {free_energy} +/- {var_free_energy}"
        )
        self.metrics.update(
            {
                "application_score_energy": energy,
                "application_score_energy_var": var_energy,
                "application_score_entropy": entropy,
                "application_score_entropy_var": var_entropy,
                "application_score_temperature": temperature,
                "application_score_temperature_var": var_temperature,
                "application_score_free_energy": free_energy,
                "application_score_free_energy_var": var_free_energy,
            }
        )

        return Data(
            Other(
                {
                    "energy": energy,
                    "energy_std": var_energy,
                    "entropy": entropy,
                    "entropy_std": var_entropy,
                    "temperature": temperature,
                    "temperature_std": var_temperature,
                    "free_energy": free_energy,
                    "free_energy_std": var_free_energy,
                }
            )
        )

    def get_metrics(self) -> dict:
        return self.metrics

    def benchmark_tag(self) -> str:
        return f"square_ising_{self.lx}_{self.ly}_{self.n_trot}_{self.dt}"
