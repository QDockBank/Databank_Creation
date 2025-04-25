# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0.
# You may obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

"""An interface for sampling problems."""
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Union

from qiskit_algorithms.minimum_eigensolvers import MinimumEigensolverResult

from qiskit.quantum_info import SparsePauliOp, Pauli

if TYPE_CHECKING:
    from .protein_folding_result import ProteinFoldingResult


class SamplingProblem(ABC):
    """An interface for sampling problems."""

    @abstractmethod
    def qubit_op(self) -> Union[Pauli, SparsePauliOp]:
        """Returns a qubit operator that represents a Hamiltonian encoding the sampling problem."""

    @abstractmethod
    def interpret(self, raw_result: MinimumEigensolverResult) -> "ProteinFoldingResult":
        """Interprets fig of an optimization."""
