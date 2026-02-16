from dataclasses import dataclass
from typing import override

from quark.core import Core, Result
from quark.interface_types import InterfaceType


@dataclass
class ExampleModule(Core):
    """
    This is an example module following the recommended structure for a quark module.

    A module must have a preprocess and postprocess method, as required by the Core abstract base class.
    A module's interface is defined by the type of data parameter those methods receive and return, dictating which other modules it can be connected to.
    Types defining interfaces should be chosen form QUARKs predefined set of types to ensure compatibility with other modules. TODO: insert link
    """

    @override
    def preprocess(self, data: InterfaceType) -> Result:
        raise NotImplementedError

    @override
    def postprocess(self, data: InterfaceType) -> Result:
        raise NotImplementedError
