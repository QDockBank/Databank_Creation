"""A class defining a main bead of a peptide."""

from typing import Tuple
from qiskit.quantum_info import SparsePauliOp
from .base_bead import BaseBead
from ..chains.side_chain import SideChain


class MainBead(BaseBead):
    """A class defining a main bead of a peptide."""

    def __init__(
        self,
        main_index: int,
        residue_type: str,
        vector_qubits: Tuple[SparsePauliOp, SparsePauliOp],
        side_chain: SideChain,
    ):
        """
        Args:
            main_index: Index of the bead on the main chain in a peptide.
            residue_type: A character representing the type of a residue for the bead.
            vector_qubits: A tuple of two SparsePauliOp operators that encodes the turn following from a
                         given bead index.
            side_chain: An object representing a side chain attached to this main bead.
        """
        super().__init__(
            "main_chain",
            main_index,
            residue_type,
            vector_qubits,
        )
        self._side_chain = side_chain

    # def __str__(self):
    #     return self.chain_type + "_" + str(self.main_index)
    #
    # def __hash__(self):
    #     return hash(str(self))
    #
    # def __eq__(self, other):
    #     if not isinstance(other, MainBead):
    #         return False
    #     return (
    #         self.main_index == other.main_index and self.chain_type == other.chain_type
    #     )

    def _get_turn_vectors(self):
        return {
            0: [0, 0],
            1: [0, 1],
            2: [1, 0],
            3: [1, 1],
        }

    def _build_turn_indicator_fun(self, turn_vector) -> SparsePauliOp:
        operator = self._full_id
        for i, qubit in enumerate(self._turn_qubits):
            if turn_vector[i] == 0:
                operator = operator @ (self._full_id - qubit)
            else:
                operator = operator @ qubit
        return (self._full_id ^ operator).simplify()

    @property
    def side_chain(self) -> SideChain:
        """Returns the side chain attached to this main bead."""
        return self._side_chain

