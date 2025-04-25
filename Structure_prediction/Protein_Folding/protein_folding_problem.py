# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0.
# You may obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

"""Defines a protein folding problem that can be passed to algorithms."""
from __future__ import annotations
from typing import TYPE_CHECKING, List, Union
from qiskit.quantum_info import Pauli, SparsePauliOp
from .interactions.interaction import Interaction
from .penalty_parameters import PenaltyParameters
from .peptide.Peptide import Peptide
from .qubit_op_builder import QubitOpBuilder
from .qubit_utils import qubit_number_reducer
from .sampling_problem import SamplingProblem


if TYPE_CHECKING:
    from .protein_folding_result import ProteinFoldingResult


class ProteinFoldingProblem(SamplingProblem):
    """Defines a protein folding problem that can be passed to algorithms. Example initialization:

    .. code-block:: python

        penalty_terms = PenaltyParameters(15, 15, 15)
        main_chain_residue_seq = "SAASSASAAG"
        side_chain_residue_sequences = ["", "", "A", "A", "A", "A", "A", "A", "S", ""]
        peptide = Peptide(main_chain_residue_seq, side_chain_residue_sequences)
        mj_interaction = MiyazawaJerniganInteraction()
        protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
        qubit_op = protein_folding_problem.qubit_op()
    """

    def __init__(
        self,
        peptide: Peptide,
        interaction: Interaction,
        penalty_parameters: PenaltyParameters,
    ):
        """
        Args:
            peptide: A peptide object that defines the protein subject to the folding problem.
            interaction: A type of interaction between the beads of the peptide.
            penalty_parameters: Parameters that define the strength of constraints enforcing in
                                the problem.
        """
        self._peptide = peptide
        self._interaction = interaction
        self._penalty_parameters = penalty_parameters
        self._pair_energies = interaction.calculate_energy_matrix(
            peptide.get_main_chain.main_chain_residue_sequence
        )
        self._qubit_op_builder = QubitOpBuilder(
            self._peptide, self._pair_energies, self._penalty_parameters
        )
        self._unused_qubits: List[int] = []

    def qubit_op(self) -> Union[SparsePauliOp, Pauli]:
        """
        Builds a qubit operator for the Hamiltonian encoding a protein folding problem. The
        number of qubits needed for optimization is optimized (compressed), if possible.
        To obtain the full qubit operator for a Hamiltonian, use the method `qubit_op_full`.

        Returns:
            A qubit operator for the Hamiltonian encoding a protein folding problem on an
            optimized number of qubits.
        """
        qubit_operator, unused_qubits = qubit_number_reducer.remove_unused_qubits(
            self._qubit_op_full()
        )
        self._unused_qubits = unused_qubits
        return qubit_operator

    def _qubit_op_full(self) -> Union[Pauli, SparsePauliOp]:
        """
        Builds a full qubit operator for the Hamiltonian encoding a protein folding problem. Full
        means that the number of qubits needed for optimization is not optimized and may be
        larger that necessary. To ensure the optimal number of qubits, use the method `qubit_op`.

        Returns:
            A qubit operator for the Hamiltonian encoding a protein folding problem.
        """
        qubit_operator = self._qubit_op_builder.build_qubit_op()
        return qubit_operator

    def interpret(self, binary_probs: dict) -> "ProteinFoldingResult":
        """
        Interprets a binary string probability distribution, in the context of this problem,
        and returns a ProteinFoldingResult. The returned class can plot the protein and generate a
        .xyz file with the coordinates of each of its atoms.

        Args:
            binary_probs: A dictionary where keys are binary strings (e.g., "0101") representing
                          states, and values are probabilities associated with each state.

        Returns:
            A :class:`~Protein_Folding.ProteinFoldingResult`
            instance that contains the protein folding result_10_11.
        """
        # pylint: disable=import-outside-toplevel
        from .protein_folding_result import ProteinFoldingResult

        # Find the binary sequence with the highest probability
        best_turn_sequence = max(binary_probs, key=binary_probs.get)

        # Construct and return the ProteinFoldingResult
        return ProteinFoldingResult(
            unused_qubits=self.unused_qubits,
            peptide=self.peptide,
            vector_sequence=best_turn_sequence,
        )

    @property
    def unused_qubits(self) -> List[int]:
        """Returns the list of indices for qubits in the original problem formulation that were
        removed during compression."""
        return self._unused_qubits

    @property
    def peptide(self) -> Peptide:
        """Returns the peptide defining the protein subject to the folding problem."""
        return self._peptide
