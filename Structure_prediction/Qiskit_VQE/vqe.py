# --*-- conding:utf-8 --*--
# @Time : 10/23/24 PM1:19
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : vqe.py

import numpy as np
from qiskit.circuit.library import EfficientSU2
from scipy.optimize import minimize
from qiskit_ibm_runtime import Session
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.primitives import Sampler
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager


class VQE:
    """
        Variational Quantum Eigensolver (VQE) class for performing quantum simulations
        using Qiskit and the IBM Quantum runtime environment. It defines the quantum
        ansatz, sets up the optimization process, and computes the minimum eigenvalue
        of a given Hamiltonian using classical-quantum hybrid optimization.
    """
    def __init__(self, service, hamiltonian, optimization_level=3, shots=200, min_qubit_num=100, maxiter=20):
        """
                Initializes the VQE class with the necessary quantum service, backend, and
                Hamiltonian information.

                Parameters:
                - service: QiskitRuntimeService object for interacting with IBM Quantum services.
                - instance: The specific IBM Quantum instance.
                - token: Authentication token to access IBM Quantum services.
                - hamiltonian: Pauli terms defining the Hamiltonian.
                - optimization_level: Integer representing the optimization level for transpiling circuits (default: 3).
                - shots: Number of shots (repeated measurements) to be performed per circuit execution (default: 1000).
        """
        self.service = service
        self.shots = shots
        self.backend = self._select_backend(min_quits=min_qubit_num)
        self.hamiltonian = hamiltonian
        self.ansatz = EfficientSU2(self.hamiltonian.num_qubits)
        self.optimization_level = optimization_level
        self.cost_history_dict = {"prev_vector": None, "iters": 0, "cost_history": []}
        self.energy_list = []
        self.maxiter = maxiter

    def _select_backend(self, min_quits):
        """
                Selects the least busy IBM Quantum backend with enough qubits and returns it.
                This ensures that the quantum job is processed faster by using a backend with fewer queued jobs.

                Returns:
                - backend: The IBM Quantum backend with the least busy queue.
        """
        backend = self.service.least_busy(simulator=False, operational=True, min_num_qubits=min_quits)
        return backend

    def _generate_pass_manager(self):
        """
               Generates a preset pass manager to optimize the quantum circuit for the selected backend.
               The pass manager optimizes the circuit by reducing gate depth and improving efficiency.

               Returns:
               - pm: The preset pass manager object.
        """
        target = self.backend.target
        pm = generate_preset_pass_manager(target=target, optimization_level=self.optimization_level)
        return pm

    def cost_func(self, params, ansatz_isa, hamiltonian_isa, estimator):
        """
                The cost function for the VQE optimization. This function estimates the energy
                for a given set of parameters by executing the ansatz circuit and measuring
                the energy expectation value of the Hamiltonian.

                Parameters:
                - params: Array of ansatz parameters to be optimized.
                - ansatz_isa: Quantum circuit representing the ansatz after pass manager optimizations.
                - hamiltonian_isa: The Hamiltonian applied with the ansatz layout.
                - estimator: The EstimatorV2 object used to estimate the energy from the quantum circuit.

                Returns:
                - energy: The estimated energy value for the given parameters.
        """
        pub = (ansatz_isa, [hamiltonian_isa], [params])
        result = estimator.run(pubs=[pub]).result()
        energy = result[0].data.evs[0]

        self.energy_list.append(energy)

        self.cost_history_dict["iters"] += 1
        self.cost_history_dict["prev_vector"] = params
        self.cost_history_dict["cost_history"].append(energy)
        print(f"Iters. done: {self.cost_history_dict['iters']} [Current cost: {energy}]")

        return energy

    def get_probability_distribution(self, optimized_params) -> 'Dict':

        circuit = self.ansatz.assign_parameters(optimized_params)
        circuit.measure_all()

        sampler = Sampler()

        job_result = sampler.run([circuit]).result()

        # print(f'job result_10_11:{job_result}')

        data = [q.binary_probabilities() for q in job_result.quasi_dists]

        measure_result = data[0]

        return measure_result


    def run_vqe(self):
        """
                Executes the VQE algorithm. This method:
                1. Generates the optimized quantum circuit using the pass manager.
                2. Prepares the Hamiltonian for computation.
                3. Initializes the optimization process using COBYLA.
                4. Returns the result_10_11 of the optimization (minimum eigenvalue of the Hamiltonian).

                Returns:
                - res: The result_10_11 of the classical optimizer containing the optimized parameters and minimum energy value.
        """
        pm = self._generate_pass_manager()
        ansatz_isa = pm.run(self.ansatz)
        hamiltonian_isa = self.hamiltonian.apply_layout(layout=ansatz_isa.layout)

        x0 = np.random.random(self.ansatz.num_parameters)

        with Session(backend=self.backend) as session:
            estimator = Estimator(mode=session)
            estimator.options.default_shots = self.shots

            res = minimize(self.cost_func, x0, args=(ansatz_isa, hamiltonian_isa, estimator), method="cobyla", options={'maxiter': self.maxiter}) #type

        # prob_distribution = self.get_probability_distribution(res.x)

        return self.energy_list, res.x, self.ansatz