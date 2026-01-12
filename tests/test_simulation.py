import pytest
import numpy as np
import scipy as sp
import pennylane as qml
from quantlop import Hamiltonian, evolve
from quantlop import get_rand_hamiltonian


@pytest.mark.parametrize("nqubits", range(1, 9))
def test_scipy(nqubits):
    psi = np.zeros(2**nqubits, dtype=complex)
    psi[0] = 1
    op = get_rand_hamiltonian(nqubits=nqubits, num_terms=nqubits * 5)

    mat = op.matrix(range(nqubits))
    psi_scipy = sp.linalg.expm(-1j * mat) @ psi
    ham = Hamiltonian.from_pennylane(op, nqubits=nqubits)
    psi_linop = evolve(ham, psi)
    assert np.allclose(psi_scipy, psi_linop)


@pytest.mark.parametrize("nqubits", range(1, 11))
def test_pennylane(nqubits):
    psi = np.zeros(2**nqubits, dtype=complex)
    psi[0] = 1
    op = get_rand_hamiltonian(nqubits=nqubits, num_terms=nqubits * 5)

    @qml.qnode(qml.device("default.qubit"))
    def circuit():
        qml.evolve(op)
        return qml.state()

    psi_pennylane = circuit()
    ham = Hamiltonian.from_pennylane(op, nqubits=nqubits)
    psi_linop = evolve(ham, psi)
    assert np.allclose(psi_pennylane, psi_linop)
