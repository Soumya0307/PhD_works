import cirq
import numpy as np
from openfermion.circuits import trotter
from openfermion import MolecularData, get_sparse_operator, jordan_wigner
from openfermionpyscf import run_pyscf

# small mol
geo = [('H', (0.0, 0.0, 0.0)), ('H', (0.0, 0.0, 0.74))]
molecule = MolecularData(geo, 'sto-3g', 1, 0)
molecule = run_pyscf(molecule, run_scf=True, run_fci=True)
H_int = molecule.get_molecular_hamiltonian(occupied_indices=[], active_indices=[0,1])

system_qubits = [cirq.LineQubit(i) for i in range(4)]
custom_algorithm = trotter.LowRankTrotterAlgorithm(final_rank=2)

try:
    ops = trotter.simulate_trotter(
            system_qubits, H_int,
            time=0.1,
            omit_final_swaps=True,
            algorithm=custom_algorithm)
    print("One step without n_steps arg works, length:", len(list(ops)))
except Exception as e:
    print("Error 1:", e)

# look if trotter_steps is an arg
import inspect
sig = inspect.signature(trotter.simulate_trotter)
print("simulate_trotter signature:", sig)
