"""Builds Pauli operators of a given size."""
from typing import Set

from qiskit.quantum_info import SparsePauliOp, Pauli

def _build_full_identity(num_qubits: int) -> SparsePauliOp:
    """
    Builds a full identity operator of a given size.

    Args:
        num_qubits: number of qubits on which a full identity operator will be created.

    Returns:
        A full identity operator of a given size.
    """
    I = SparsePauliOp('I')

    full_identity = I
    for _ in range(1, num_qubits):
        full_identity = I ^ full_identity

    return full_identity


def _build_pauli_z_op(num_qubits: int, pauli_z_indices: Set[int]) -> SparsePauliOp:
    """
    Builds a Pauli operator of a given size with Pauli Z operators on indicated positions and
    identity operators on other positions.

    Args:
        num_qubits: number of qubits on which a Pauli operator will be created.
        pauli_z_indices: a set of indices in a Pauli operator on which a Pauli Z operator shall
                        appear.

    Returns:
        A Pauli operator of a given size with Pauli Z operators on indicated positions and
        identity operators on other positions.
    """
    I = SparsePauliOp('I')
    Z = SparsePauliOp('Z')

    if 0 in pauli_z_indices:
        operator = Z
    else:
        operator = I
    for i in range(1, num_qubits):
        if i in pauli_z_indices:
            operator = Z ^ operator
        else:
            operator = I ^ operator

    return operator


def _build_full_identity_Pauli(num_qubits: int) -> Pauli:

    I = Pauli('I')

    full_identity = I
    for _ in range(1, num_qubits):
        full_identity = I ^ full_identity

    return full_identity


def _build_pauli_z_Pauli(num_qubits: int, pauli_z_indices: Set[int]) -> Pauli:

    I = Pauli('I')
    Z = Pauli('Z')

    if 0 in pauli_z_indices:
        operator = Z
    else:
        operator = I
    for i in range(1, num_qubits):
        if i in pauli_z_indices:
            operator = Z ^ operator
        else:
            operator = I ^ operator

    return operator
