from ._core import PauliWord as _PauliWord
from ._core import Hamiltonian as _Hamiltonian
from ._core import evolve as _evolve
from quantlop.utils import get_rand_hamiltonian


class Hamiltonian(_Hamiltonian):
    @classmethod
    def from_pennylane(cls, operator_pl, nqubits):
        pws = []
        for pauliword_pl, coeff in operator_pl.pauli_rep.items():
            string = "".join(pauliword_pl.get(i, "I") for i in range(nqubits))
            pws.append(_PauliWord(coeff=coeff, string=string))
        return cls(pauli_words=pws)


def evolve(ham, psi, coeff=1):
    return _evolve(ham, psi, coeff)
