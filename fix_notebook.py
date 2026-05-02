import json
import logging

try:
    with open("Our_circuit_and_resource_est.ipynb", "r") as f:
        nb = json.load(f)

    # The failing cell is likely the last or second to last.
    # We will find the cell that contains `geometries.items():`
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            if "for name, atom_str in geometries.items():" in source:
                target_cell = cell
                break
    else:
        raise Exception("Could not find the target cell")

    new_source = []
    # Let's iterate through lines and insert/replace what we need
    lines = target_cell["source"]
    for i, line in enumerate(lines):
        if "steps_H, errors_H = trotter_error_vs_steps_H(" in line:
            # Insert the qubit definition
            new_source.extend([
                "    n_system = int(np.log2(H_mat.shape[0]))\n",
                "    ancilla = cirq.LineQubit(0)\n",
                "    system_qubits = [cirq.LineQubit(i) for i in range(1, n_system + 1)]\n",
                "    # Create PauliSum for build_random_W_stack\n",
                "    H_pauli = csr_to_paulisum(H_sparse, system_qubits)\n",
                "\n"
            ])
            new_source.append(line)
        elif "H_pauli=fermion_ham" in line:
            new_source.append(line.replace("fermion_ham", "H_pauli"))
        else:
            new_source.append(line)

    target_cell["source"] = new_source

    with open("Our_circuit_and_resource_est.ipynb", "w") as f:
        json.dump(nb, f, indent=2)

    print("Notebook successfully fixed!")

except Exception as e:
    logging.exception("Failed to fix notebook")
