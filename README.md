# Submission
This directory contains the submission files for the [2023 Quantum Computing for Drug Discovery Challenge at ICCAD](https://qccontest.github.io/QC-Contest/index.html). The goal of this contest is to calculate the ground state energy of hydroxyl cation (·OH) by using the variational quantum eigensolver (VQE) algorithm on noisy intermediate-scale quantum (NISQ) devices. To run a VQE on NISQ devices sucessfully, we have to (1) reduce the problem size, (2) design an efficient ansatz, and (3) mitigate the hardware noises. Thus, we applied the following methods to achieve these goals:

* **Problem size reduction**
    * **ZZ symmetry tapering**: This method exploits the symmetry of the Hamiltonian to eliminate redundant qubits. It reduces the number of qubits from 12 to 8.
    * **Single-X-Removal tapering**: This is a novel method that removes the non-symmetric terms in the Hamiltonian, which restores the symmetry and allows further tapering. It reduces the number of qubits from 8 to 6.
    * **Pauli grouping**: This method groups the qubit-wise commuting Pauli terms in the Hamiltonian to reduce the number of measurements.
* **Efficient ansatz circuit**
    * We design an ansatz circuit that simulates the excitation of electrons using Givens rotations, which are particle-conserving transformations. The ansatz only contains the excitation gates that correspond to the important terms in the Hamiltonian, which reduces the depth and the number of parameters of the circuit.
* **Mitigating quantum errors**
    * **Z-basis exchanging**: This method exchanges the Z-basis of the Hamiltonian, which reduces the readout error rate by exploiting the asymmetry of the error from 1 to 0 and from 0 to 1.
    * **Error mitigations**: We apply the readout error mitigation and the reference state error mitigation [1] to correct the output of the quantum device using calibration matrices and reference states, respectively.

By applying these methods, we obtain following results:
|Backend|Accuracy (%)|Duration|shots|
|------|---|---|---|
|cairo|99.9412|2784|1799984|
|montreal|99.9612|5056|1799984|
|kolkata|99.9196|3904|1799984|

<!--
## Directory structure
```
|-- algorithm_seeds
|   `-- requiredseeds.txt
|-- hamiltonian
|   `-- OHhamiltonian.txt
|-- noise_model
|   |-- fakecairo.pkl
|   |-- fakekolkata.pkl
|   `-- fakemontreal.pkl
|-- qasm
|   |-- cairo_QASM
|   |-- cairo_ref_QASM
|   |-- kolkata_QASM
|   |-- kolkata_ref_QASM
|   |-- montreal_QASM
|   `-- montreal_ref_QASM
|-- Qubit_taper.py
|-- How_to_use.ipynb
|-- evaluation.ipynb
|-- result.json
`-- README.md
```
-->

## Requirements
```
conda create -n iccad python=3.11
conda activate iccad
pip install qiskit, qiskit_aer, qiskit_nature, networkx, matplotlib, pylatexenc
pip install --prefer-binary pyscf
```

## How to run
### Optimization
Follow the notebook `How_to_use.ipynb` to run the proposed hamiltonian tapering method and the optimization. This contains the detailed explanation of our methods including the Single-X-Removal tapering, Z-basis exchanging, Pauli grouping, efficient chemistry-inspired ansatz and more.

### Evaluation
Run the notebook `evaluation.ipynb` to evaluate the performance. The result will be saved in `result.json`. The format of `result.json` is shown below. Note that the `metric` is the objective of the contest, but calculated for each seed.
```json
[
    {
        "backend": "cairo",
        "seed": 20,
        "expval": -78.70870653629866,
        "metric": 99.9410990978365,
        "duration": 2784,
        "shot": 1799984
    },
    ...
]
```
The result of the evaluation (the last section of `evaluation.ipynb`) is shown below. 
```
Backend: cairo
  Avg. expval: -78.70880023201202
  Metric: 99.94122505469606 (%)
Backend: montreal
  Avg. expval: -78.72368420157821
  Metric: 99.96123384695724 (%)
Backend: kolkata
  Avg. expval: -78.69275194853438
  Metric: 99.91965105396831 (%)
```

## Team: Quantum Korea
- Youngseok Lee ([pop756](https://github.com/pop756))
- Youngdu Kim ([ydkim0119](https://github.com/ydkim0119))
- Dohun Kim ([yh08037](https://github.com/yh08037))

## References
[1] Lolur, Phalgun, et al. "Reference-State Error Mitigation: A Strategy for High Accuracy Quantum Computation of Chemistry." *Journal of Chemical Theory and Computation* 19.3 (2023): 783-789.
