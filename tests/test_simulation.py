import pytest
import numpy as np
import scipy as sp
import pennylane as qml
from quantlop import Hamiltonian, evolve
from quantlop import get_rand_hamiltonian


@pytest.mark.parametrize("num_qubits", range(1, 9))
@pytest.mark.parametrize("method", ("higham", "krylov"))
def test_scipy(num_qubits, method):
    psi = np.zeros(2**num_qubits, dtype=complex)
    psi[0] = 1
    op = get_rand_hamiltonian(nqubits=num_qubits, num_terms=num_qubits * 5)
    mat = op.matrix(range(num_qubits))
    psi_scipy = sp.linalg.expm(-1j * mat) @ psi
    ham = Hamiltonian.from_pennylane(op, num_qubits=num_qubits)
    psi_linop = evolve(ham, psi, method=method)
    assert np.allclose(psi_scipy, psi_linop)


@pytest.mark.parametrize("num_qubits", range(1, 11))
@pytest.mark.parametrize("method", ("higham", "krylov"))
def test_pennylane(num_qubits, method):
    psi = np.zeros(2**num_qubits, dtype=complex)
    psi[0] = 1
    op = get_rand_hamiltonian(nqubits=num_qubits, num_terms=num_qubits * 5)

    @qml.qnode(qml.device("default.qubit"))
    def circuit():
        qml.evolve(op)
        return qml.state()

    psi_pennylane = circuit()
    ham = Hamiltonian.from_pennylane(op, num_qubits=num_qubits)
    psi_linop = evolve(ham, psi, method=method)
    assert np.allclose(psi_pennylane, psi_linop)
