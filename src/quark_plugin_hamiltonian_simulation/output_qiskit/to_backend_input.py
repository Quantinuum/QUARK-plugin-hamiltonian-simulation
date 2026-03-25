from dataclasses import dataclass
from math import ceil
from typing import override
import logging

from qiskit import QuantumCircuit
from quark.core import Core, Result, Data
from quark.interface_types import InterfaceType, Other
from quark_plugin_hamiltonian_simulation.abstract_classes.abstract_circuit import (
    AbstractCircuits,
)
from quark_plugin_quantinuum.interfaces.backend_input_pytket import BackendInputPytket
from quark_plugin_quantinuum.interfaces.backend_input_qiskit import BackendInputQiskit
from pytket.extensions.qiskit import qiskit_to_tk

logger = logging.getLogger()


@dataclass
class ToBackendInput(Core):
    """
    Represents a dataclass used to convert circuit data to a specified backend input
    type with preprocessing and postprocessing capabilities.

    The main purpose of this class is to handle the preparation of circuit data for
    backend processing and to convert the output into a compatible format.

    Attributes:
        circuit_type: Specifies the type of the circuit representation
            to be used for backend processing, either "pytket" or "qiskit". Defaults to "pytket".

    Methods:
        preprocess(data: Other[AbstractCircuit]) -> Result:
            Takes circuit data and processes it into a Quantinuum backend-compatible
            format. This includes creating the required circuit representation and
            extracting shot proportions.

        postprocess(data: InterfaceType):
            Processes the backend output data into a consistent format for use.
    """

    circuit_type: str = "pytket"

    @override
    def preprocess(self, input_data: Other[AbstractCircuits]) -> Result:
        abstract_circuits = input_data.data.circuits
        benchmark_tag = input_data.data.benchmark_name

        circuits = [
            create_circuit_qiskit(circuit["circuit"], circuit["number_of_qubits"])
            for circuit in abstract_circuits
        ]

        shots_per_circuit = [
            ceil(circuit["n_shots"] * circuit["shot_proportions"])
            for circuit in abstract_circuits
        ]

        match self.circuit_type:
            case "pytket":
                circuits_pytket = [qiskit_to_tk(circ) for circ in circuits]
                return Data(
                    Other(
                        BackendInputPytket(
                            circuits_pytket, shots_per_circuit, benchmark_tag
                        )
                    )
                )
            case "qiskit":
                return Data(
                    Other(
                        BackendInputQiskit(circuits, shots_per_circuit, benchmark_tag)
                    )
                )
            case _:
                raise ValueError(
                    f"Invalid output circuit type: {self.circuit_type} only qiskit and pytket are supported"
                )

    @override
    def postprocess(self, data: InterfaceType):
        return Data(data)


def create_circuit_qiskit(U, N):
    V = QuantumCircuit(N, N)

    for u in U:
        if u[0] == "x":
            V.x(u[1])
        elif u[0] == "y":
            V.y(u[1])
        elif u[0] == "z":
            V.z(u[1])
        elif u[0] == "h":
            V.h(u[1])
        elif u[0] == "s":
            V.s(u[1])
        elif u[0] == "sdg":
            V.sdg(u[1])
        elif u[0] == "rx":
            if u[1] != 0:
                V.rx(2 * u[1], u[2])
        elif u[0] == "ry":
            if u[1] != 0:
                V.ry(2 * u[1], u[2])
        elif u[0] == "rz":
            if u[1] != 0:
                V.rz(2 * u[1], u[2])
        elif u[0] == "cx":
            V.cx(u[1], u[2])
        elif u[0] == "cy":
            V.cy(u[1], u[2])
        elif u[0] == "cz":
            V.cz(u[1], u[2])
        elif u[0] == "rxx":
            if u[1] != 0:
                V.rxx(2 * u[1], u[2], u[3])
        elif u[0] == "ryy":
            if u[1] != 0:
                V.ryy(2 * u[1], u[2], u[3])
        elif u[0] == "rzz":
            if u[1] != 0:
                V.rzz(2 * u[1], u[2], u[3])
        elif u[0] == "measure":
            V.measure(u[1], u[2])
        elif u[0] == "barrier":
            V.barrier()

    return V
