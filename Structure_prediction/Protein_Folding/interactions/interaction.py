# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0.
# You may obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

"""An abstract class defining an interaction between beads of a peptide."""
from abc import ABC, abstractmethod

import numpy as np

# pylint: disable=too-few-public-methods


class Interaction(ABC):
    """An abstract class defining an interaction between beads of a peptide."""

    @abstractmethod
    def calculate_energy_matrix(self, residue_sequence: str) -> np.ndarray:
        """
        Calculates an energy matrix for a particular interaction type.

        Args:
            residue_sequence: A string that contains characters defining residues for
                            a chain of proteins.

        Returns:
            Numpy array of pair energies for amino acids.
        """
