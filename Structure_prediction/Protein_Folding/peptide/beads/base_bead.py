"""An abstract class defining a bead of a peptide."""

from abc import ABC, abstractmethod
from typing import Tuple, Union, Optional
from qiskit.quantum_info import SparsePauliOp
from ..pauli_ops_builder import _build_full_identity
from ...residue_validator import _validate_residue_symbol

class BaseBead(ABC):
    """An abstract class defining a bead of a peptide."""

    def __init__(
        self,
        chain_type: str,
        main_index: int,
        residue_type: Optional[str],
        vector_qubits: Tuple[SparsePauliOp, SparsePauliOp],
    ):
        """
        Args:
            chain_type: Type of the chain, either "main_chain" or "side_chain".
            main_index: Index of the bead on the main chain in a peptide.
            residue_type: A character representing the type of a residue for the bead. An empty
                          string in case of non-existing side bead.
            vector_qubits: A tuple of two Pauli operators that encodes the turn following from a
                         given bead index.
        """
        self.chain_type = chain_type
        self.main_index = main_index
        self._residue_type = residue_type
        _validate_residue_symbol(residue_type)
        self._turn_qubits = vector_qubits

        if self._residue_type and self.turn_qubits is not None:
            self._full_id = _build_full_identity(vector_qubits[0].num_qubits)

    @property
    def turn_qubits(self) -> Tuple[SparsePauliOp, SparsePauliOp]:
        """Returns the list of two qubits that encode the turn following from the bead."""
        return self._turn_qubits

    @property
    def residue_type(self) -> Optional[str]:
        """Returns a residue type."""
        return self._residue_type

    @property
    def indicator_functions(
        self,
    ) -> Union[None, Tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]]:
        """
        Returns all turn indicator functions for the bead.

        Returns:
            A tuple of all turn indicator functions for the bead.
        """
        if self.turn_qubits is None:
            return None
        return tuple(self.get_turn_indicator_function(i) for i in range(4))

    def get_turn_indicator_function(self, turn_index: int) -> SparsePauliOp:
        """Returns the turn indicator function for the specified turn index."""
        turn_vectors = self._get_turn_vectors()
        turn_vector = turn_vectors[turn_index]
        return self._build_turn_indicator_fun(turn_vector)

    @abstractmethod
    def _get_turn_vectors(self):
        """Returns a dictionary of turn vectors. Must be implemented in subclasses."""
        pass

    @abstractmethod
    def _build_turn_indicator_fun(self, turn_vector) -> SparsePauliOp:
        """Builds the turn indicator function based on the turn vector. Must be implemented in subclasses."""
        pass



