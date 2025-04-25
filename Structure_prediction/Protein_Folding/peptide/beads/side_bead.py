"""A class defining a side bead of a peptide."""

from typing import Tuple, Optional
from qiskit.quantum_info import SparsePauliOp
from .base_bead import BaseBead


class SideBead(BaseBead):
    """A class defining a side bead of a peptide."""

    def __init__(
        self,
        main_index: int,
        side_index: int,
        residue_type: Optional[str],
        vector_qubits: Tuple[SparsePauliOp, SparsePauliOp],
    ):
        """
        Args:
            main_index: Index of the bead on the main chain in a peptide to which the side
                        chain of this side bead is attached.
            side_index: Index of the bead on the related side chain in a peptide.
            residue_type: A character representing the type of a residue for the bead. Empty
                          string if a side bead does not exist.
            vector_qubits: A tuple of two Pauli operators that encodes the turn following from a given
                         bead index.
        """
        super().__init__(
            "side_chain",
            main_index,
            residue_type,
            vector_qubits,
        )
        self.side_index = side_index

    # def __str__(self):
    #     return (
    #         self.chain_type
    #         + "_"
    #         + str(self.side_index)
    #         + "_main_chain_ind_"
    #         + str(self.main_index)
    #     )
    #
    # def __hash__(self):
    #     return hash(str(self))
    #
    # def __eq__(self, other):
    #     if not isinstance(other, SideBead):
    #         return False
    #     return (
    #         self.main_index == other.main_index
    #         and self.side_index == other.side_index
    #         and self.chain_type == other.chain_type
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
        return (operator ^ self._full_id).simplify()

