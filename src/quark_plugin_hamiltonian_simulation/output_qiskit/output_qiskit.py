from dataclasses import dataclass
from typing import override

from quark.core import Core, Result,Data
from quark.interface_types import InterfaceType,Other
from .output_qiskit_helpers import create_circuit_qiskit
from quark_plugin_hamiltonian_simulation.backends.backend_input import BackendInput
from quark_plugin_hamiltonian_simulation.abstract_classes.abstract_circuit import AbstractCircuit

@dataclass
class OutputQiskit(Core):
    """
    This is an example module following the recommended structure for a quark module.

    A module must have a preprocess and postprocess method, as required by the Core abstract base class.
    A module's interface is defined by the type of data parameter those methods receive and return, dictating which other modules it can be connected to.
    Types defining interfaces should be chosen form QUARKs predefined set of types to ensure compatibility with other modules. TODO: insert link
    """

    @override
    def preprocess(self, data: Other[AbstractCircuit]) -> Result:
        U=[create_circuit_qiskit(circuit['circuit'],circuit['number_of_qubits']) for circuit in data.data.circuits]
        shot_proportions=[circuit['shot_proportions'] for circuit in data.data.circuits]

        return Data(Other(BackendInput(circuits=U,shot_proportions=shot_proportions)))

    @override
    def postprocess(self, data: InterfaceType):
        return Data(data)
