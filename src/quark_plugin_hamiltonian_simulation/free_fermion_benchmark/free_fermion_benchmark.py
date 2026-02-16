from dataclasses import dataclass
from typing import override

from quark.core import Core, Result,Data
from quark.interface_types import InterfaceType,Other
from .free_fermion_benchmark_helpers import create_circuit
from .free_fermion_benchmark_score import (
    exact_values_and_variance,
    computes_score_values,
    extract_simulation_results,
    create_and_show_plot,
    create_and_show_histogram
)
import matplotlib.pyplot as plt
import numpy as np
import logging
logger = logging.getLogger()
from dataclasses import dataclass, field
from typing import override, Dict

from quark_plugin_hamiltonian_simulation.backends.backend_result import BackendResult
from quark_plugin_hamiltonian_simulation.abstract_classes.abstract_circuit import AbstractCircuit

@dataclass
class FreeFermionBenchmark(Core):

    lx: int = 2
    ly: int = 2
    periodic: bool=False
    two_spin_species: bool=False
    create_plot: bool = False
    n_trot: int=2
    dt: float=0.5
    metrics: Dict[str, float | int] = field(init=False, default_factory=dict)

    @override
    def preprocess(self, data: InterfaceType) -> Result:
        
        if(self.two_spin_species):
            N=2*self.lx*self.ly
        else:
            N=self.lx*self.ly

        if(self.periodic):
            N+=self.lx*self.ly//2
        else:
            N+=(self.lx//2-1)*self.ly//2+self.lx//2*(self.ly//2-1)

        u=[{'number_of_qubits':N,'circuit':create_circuit(self.lx,self.ly,self.dt,self.periodic,self.two_spin_species,n),'shot_proportions':1} for n in range(self.n_trot+1)]

        return Data(Other(AbstractCircuit(u)))

    @override
    def postprocess(self, input_data: Other[BackendResult]) -> Result:
        backend_result = input_data.data

        counts_per_circuit = backend_result.counts

        simulation_results_list,simulation_results_post_list, stabilizers = extract_simulation_results(self.dt, self.lx, self.ly, counts_per_circuit,self.periodic,self.two_spin_species)
        
        simulation_results=np.array(simulation_results_list)
        simulation_results_post=np.array(simulation_results_post_list)

        exact_results = np.real(np.array(exact_values_and_variance(self.n_trot, self.dt, self.lx, self.ly
                                                                 ,self.periodic,self.two_spin_species)))
        
        score_gate, score_shot, score_runtime = computes_score_values(exact_results[:, 1] - simulation_results[:, 1],
                                                                      simulation_results[:, 2],
                                                                      exact_results[:, 2], self.lx * self.ly)
        logger.info(f"Benchmark score (number of gates): {score_gate}")
        logger.info(f"Benchmark score (number of shots): {score_shot}")
        logger.info(f"Benchmark score (number of trotter steps): {score_runtime}")
        self.metrics.update({
            "application_score_value": score_gate,
            "application_score_value_gates": score_gate,
            "application_score_value_shots": score_shot,
            "application_score_value_trotter_steps": score_runtime,
            "application_score_unit": "N_gates",
            "application_score_type": "int"
        })
        if self.create_plot:
            create_and_show_plot(
                self.n_trot,
                self.dt,
                simulation_results,
                simulation_results_post,
                exact_results,
                score_gate,
            )
            create_and_show_histogram(stabilizers)
        return Data(Other(score_gate))

    def get_metrics(self) -> dict:
        return self.metrics
