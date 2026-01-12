import numpy as np
import pennylane as qml


char2qml = {
    "I": qml.I,
    "X": qml.X,
    "Y": qml.Y,
    "Z": qml.Z,
}
paulis = ("I", "X", "Y", "Z")


def get_rand_statevector(nqubits, seed=None):
    rng = np.random.default_rng(seed=seed)
    real = rng.random(2**nqubits)
    imag = rng.random(2**nqubits)
    psi = real + 1j * imag
    return psi / np.linalg.norm(psi)


def get_rand_pauliword(nqubits, paulis=paulis, seed=None):
    rng = np.random.default_rng(seed=seed)
    pauli_list = rng.choice(paulis, size=nqubits)
    return "".join(pauli_list)


def get_rand_hamiltonian(nqubits, num_terms, paulis=paulis, seed=None):
    rng = np.random.default_rng(seed=seed)
    coeffs = rng.random(num_terms)
    observables = []
    for _ in range(num_terms):
        word = get_rand_pauliword(nqubits=nqubits, paulis=paulis, seed=seed)
        obs = [char2qml[p](i) for i, p in enumerate(word)]
        observables.append(qml.prod(*obs))
    ham = sum(c * obs for c, obs in zip(coeffs, observables))
    return qml.simplify(ham)
