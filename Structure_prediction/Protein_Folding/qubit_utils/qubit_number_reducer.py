# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0.
# You may obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

"""Removes qubit registers that are not relevant for the problem."""
from typing import Union, List, Dict, Tuple

import numpy as np
from qiskit.quantum_info import PauliList, SparsePauliOp, Pauli, Operator


def remove_unused_qubits(
    total_hamiltonian: Union[SparsePauliOp, Pauli]
) -> Tuple[Union[SparsePauliOp, Pauli], List[int]]:
    """
    Removes those qubits from a total Hamiltonian that are equal to an identity operator across
    all terms, i.e. they are irrelevant for the problem. It makes the number of qubits required
    for encoding the problem smaller or equal.

    Args:
        total_hamiltonian: A full Hamiltonian for the protein folding problem.

    Returns:
        Tuple consisting of the total_hamiltonian compressed to an equivalent Hamiltonian and
        indices of qubits in the original Hamiltonian that were unused as optimization variables.
    """
    unused_qubits = _find_unused_qubits(total_hamiltonian)
    num_qubits = total_hamiltonian.num_qubits


    if isinstance(total_hamiltonian, SparsePauliOp):
        return (
            _compress_pauli_sum_op(num_qubits, total_hamiltonian, unused_qubits),
            unused_qubits,
        )
    # return None, None

def _compress_pauli_sum_op(
        num_qubits: int,
        total_hamiltonian: Union[SparsePauliOp, Pauli],
        unused_qubits: List[int],
) -> Union[SparsePauliOp, Pauli]:
    """
    Compress the given Pauli or SparsePauliOp Hamiltonian by removing unused qubits.

    Args:
        num_qubits: Total number of qubits in the system.
        total_hamiltonian: The Hamiltonian as a Pauli or SparsePauliOp.
        unused_qubits: List of qubit indices that are unused and should be removed.

    Returns:
        Compressed SparsePauliOp or Pauli.
    """
    new_paulis = []
    new_coeffs = []

    if isinstance(total_hamiltonian, SparsePauliOp):
        # Use .to_list() to get the terms and their coefficients
        terms = total_hamiltonian.to_list()
        for pauli_str, coeff in terms:
            pauli_obj = Pauli(pauli_str)  # Create a Pauli object from the string

            # Extract Z and X as boolean arrays from Pauli object
            table_z = pauli_obj.z
            table_x = pauli_obj.x

            # Remove the unused qubits
            new_table_z, new_table_x = _calc_reduced_pauli_tables(
                num_qubits, table_x, table_z, unused_qubits
            )

            # Construct new Pauli from the reduced Z and X tables
            new_pauli = Pauli((new_table_z, new_table_x))
            new_paulis.append(new_pauli)
            new_coeffs.append(coeff)

    # Create a new PauliList and return a compressed SparsePauliOp
    new_pauli_list = PauliList(new_paulis)
    total_hamiltonian_compressed = SparsePauliOp(new_pauli_list, coeffs=new_coeffs).simplify()

    return total_hamiltonian_compressed


# def _calc_reduced_pauli_tables(
#     num_qubits: int, table_x, table_z, unused_qubits
# ) -> Tuple[List[bool], List[bool]]:
#     new_table_z = []
#     new_table_x = []
#     for ind in range(num_qubits):
#         if ind not in unused_qubits:
#             new_table_z.append(table_z[ind])
#             new_table_x.append(table_x[ind])
#
#     return new_table_z, new_table_x

def _calc_reduced_pauli_tables(num_qubits, table_x, table_z, unused_qubits):
    """
    Helper function to reduce the Pauli tables by removing the unused qubits.
    """
    # Remove the unused qubits from the Pauli tables
    new_table_z = np.delete(table_z, unused_qubits)
    new_table_x = np.delete(table_x, unused_qubits)

    return new_table_z, new_table_x


def _find_unused_qubits(total_hamiltonian: Union[SparsePauliOp, Pauli]) -> List[int]:
    used_map: Dict[int, bool] = {}
    unused = []
    num_qubits = total_hamiltonian.num_qubits

    if isinstance(total_hamiltonian, SparsePauliOp):
        for term in total_hamiltonian:
            table_z = term.paulis.z[0]
            _update_used_map(num_qubits, table_z, used_map)

    else:
        raise ValueError("Invalid input PDBbind_data for Pauli.")

    for ind in range(num_qubits):
        if ind not in used_map:
            unused.append(ind)

    return unused


def _update_used_map(num_qubits: int, table_z, used_map: Dict[int, bool]):
    for ind in range(num_qubits):
        if table_z[ind]:
            used_map[ind] = True
