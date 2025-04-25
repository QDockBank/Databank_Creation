# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0.
# You may obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

"""Changes certain qubits to fixed values."""
from typing import Union

import numpy as np
from qiskit.quantum_info import SparsePauliOp, Pauli
from qiskit.quantum_info import Operator, PauliList


def _fix_qubits(
        operator: Union[int, SparsePauliOp, Pauli, Operator],
        has_side_chain_second_bead: bool = False,
) -> Union[int, SparsePauliOp, Pauli, Operator]:
    """
    Assigns predefined values for turn qubits on positions 0, 1, 2, 3, 5 in the main chain
    without the loss of generality. Qubits on these positions are considered fixed and not subject
    to optimization.

    Args:
        operator: an operator whose qubits shall be fixed.

    Returns:
        An operator with relevant qubits changed to fixed values.
    """
    # If the type is int, return
    if isinstance(operator, int):
        return operator

    new_tables = []

    # Handle SparsePauliOp
    if isinstance(operator, SparsePauliOp):
        table_z = np.copy(operator.paulis.z)  # Get Z
        table_x = np.copy(operator.paulis.x)  # Get X

        # binary vals
        _preset_binary_vals(table_z, has_side_chain_second_bead)

        # Change Z and X to String
        for i in range(table_z.shape[0]):
            pauli_str = ''.join([
                'Z' if z and not x else 'X' if not z and x else 'Y' if z and x else 'I'
                for z, x in zip(table_z[i], table_x[i])
            ])
            new_tables.append(Pauli(pauli_str))
        new_coeffs = operator.coeffs

    # Handel Pauli
    elif isinstance(operator, Pauli):
        table_z = np.copy(operator.z)  # Get Z
        table_x = np.copy(operator.x)  # Get X
        _preset_binary_vals(table_z, has_side_chain_second_bead)

        # Create a new Pauli
        return Pauli((table_z, table_x))

    # Else
    else:
        raise ValueError("Unsupported operator type")

    # Buliding a new SparsePauliOp
    new_pauli_table = PauliList(new_tables)
    operator_updated = SparsePauliOp(new_pauli_table, coeffs=new_coeffs)

    # Simplify Operator
    operator_updated = operator_updated.simplify()

    return operator_updated

def _calc_updated_coeffs(
    hamiltonian: Union[SparsePauliOp, Pauli], table_z, has_side_chain_second_bead: bool
) -> np.ndarray:
    coeffs = np.copy(hamiltonian.coeffs[0])
    if len(table_z) > 1 and table_z[1] == np.bool_(True):
        coeffs = -1 * coeffs
    if (
        not has_side_chain_second_bead
        and len(table_z) > 6
        and table_z[5] == np.bool_(True)
    ):
        coeffs = -1 * coeffs
    return coeffs


def _preset_binary_vals(table_z, has_side_chain_second_bead: bool):
    main_beads_indices = [0, 1, 2, 3]
    if not has_side_chain_second_bead:
        main_beads_indices.append(5)
    for index in main_beads_indices:
        _preset_single_binary_val(table_z, index)


def _preset_single_binary_val(table_z, index: int):
    try:
        table_z[index] = np.bool_(False)
    except IndexError:
        pass
