# --*-- conding:utf-8 --*--
# @Time : 11/13/24 12:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : measure.py

from qiskit import QuantumCircuit
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from typing import Dict

class StateCalculator:
    def __init__(self,service ,min_qubit_num, ansatz: QuantumCircuit):
        """
        Initialize the ProbabilityDistributionCalculator with a predefined ansatz circuit.
        """
        self.service = service
        self.ansatz = ansatz
        # self.backend = self._select_backend(min_quits=min_qubit_num)
        self.min_qubits = min_qubit_num

    # def _select_backend(self, min_quits):
    #     """
    #             Selects the least busy IBM Quantum backend with enough qubits and returns it.
    #             This ensures that the quantum job is processed faster by using a backend with fewer queued jobs.
    #
    #             Returns:
    #             - backend: The IBM Quantum backend with the least busy queue.
    #     """
    #     backend = self.service.least_busy(simulator=False, operational=True, min_num_qubits=min_quits)
    #     return backend

    def get_probability_distribution(self, optimized_params) -> Dict:
        """
        Calculate the probability distribution by assigning optimized parameters to the ansatz,
        executing the circuit, and returning the measurement fig.

        Parameters:
        optimized_params: The parameters to assign to the ansatz circuit.

        Returns:
        Dict: A dictionary with binary strings as keys and their probabilities as values.
        """
        circuit = self.ansatz.assign_parameters(optimized_params)
        circuit.measure_all()

        backend = self.service.least_busy(simulator=False, operational=True, min_num_qubits=self.min_qubits)

        pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
        isa_circuit = pm.run(circuit)

        sampler = Sampler(backend,options={"default_shots": 100000})
        job = sampler.run([isa_circuit])
        result = job.result()

        threshold = 0 #

        pub_result = result[0]

        counts = pub_result.data.meas.get_counts()

        # total_shots = sum(counts.values())
        measure_result = {key: value for key, value in counts.items() if value > threshold}

        return measure_result
