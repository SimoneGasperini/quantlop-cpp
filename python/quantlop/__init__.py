from .quantlop_cpp import PauliWord as PauliWord_cpp
from .quantlop_cpp import Hamiltonian as Hamiltonian_cpp
from .quantlop_cpp import evolve as evolve_cpp
from .utils import get_rand_hamiltonian


class Hamiltonian(Hamiltonian_cpp):
    @classmethod
    def from_pennylane(cls, operator_pl, num_qubits):
        pws = []
        for pauliword_pl, coeff in operator_pl.pauli_rep.items():
            string = "".join(pauliword_pl.get(i, "I") for i in range(num_qubits))
            pws.append(PauliWord_cpp(coeff=coeff, string=string))
        return cls(pauli_words=pws)


def evolve(ham, psi, coeff=1, method="higham"):
    return evolve_cpp(ham, psi, coeff, method)
