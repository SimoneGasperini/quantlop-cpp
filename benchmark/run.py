import time
import json
import numpy as np
import pennylane as qml
from tqdm import trange
from memory_profiler import memory_usage
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from quantlop import Hamiltonian, evolve
from quantlop.utils import get_rand_hamiltonian


@qml.qnode(qml.device("default.qubit"))
def circuit(op):
    qml.evolve(op)
    return qml.state()


def pennylane_simulation(op):
    # same as scipy_dense
    return circuit(op)


def scipy_dense_simulation(nq, op, psi):
    dense = op.matrix(wire_order=range(nq))
    return expm(-1j * dense) @ psi


def scipy_sparse_simulation(nq, op, psi):
    sparse = op.sparse_matrix(wire_order=range(nq))
    return expm_multiply(-1j * sparse, psi, traceA=0)


def quantlop_simulation(nq, op, psi):
    linop = Hamiltonian.from_pennylane(op, nqubits=nq)
    return evolve(linop, psi)


def runtime_and_memory(func, *args, reps, interval=0.0005):
    runtime = []
    memory = []
    for _ in trange(reps, ncols=80):
        t1 = time.perf_counter()
        mem, result = memory_usage(
            (func, (args), {}),
            interval=interval,
            retval=True,
            include_children=True,
            multiprocess=True,
            max_iterations=1,
        )
        t2 = time.perf_counter()
        runtime.append(t2 - t1)
        memory.append(max(mem) - mem[0])
    return runtime, memory, result


def run_benchmark(num_qubits, time_fname, mem_fname, reps):
    runtime = {"SciPy dense": {}, "SciPy sparse": {}, "QuantLop": {}}
    memory = {"SciPy dense": {}, "SciPy sparse": {}, "QuantLop": {}}
    for nq in num_qubits:
        print(f"\nRunning simulation for {nq} qubits:")
        op = get_rand_hamiltonian(nqubits=nq, num_terms=5 * nq)
        psi = np.zeros(2**nq, dtype=complex)
        psi[0] = 1
        if nq < 15:
            time, mem, result1 = runtime_and_memory(
                scipy_dense_simulation, nq, op, psi, reps=reps
            )
            runtime["SciPy dense"][nq] = time
            memory["SciPy dense"][nq] = mem
        time, mem, result2 = runtime_and_memory(
            scipy_sparse_simulation, nq, op, psi, reps=reps
        )
        runtime["SciPy sparse"][nq] = time
        memory["SciPy sparse"][nq] = mem
        time, mem, result3 = runtime_and_memory(
            quantlop_simulation, nq, op, psi, reps=reps
        )
        runtime["QuantLop"][nq] = time
        memory["QuantLop"][nq] = mem
        if nq < 15:
            assert np.allclose(result1, result2)
        assert np.allclose(result2, result3)
        with open(time_fname, "w") as file:
            json.dump(runtime, file, indent=4)
        with open(mem_fname, "w") as file:
            json.dump(memory, file, indent=4)


if __name__ == "__main__":
    run_benchmark(
        num_qubits=range(1, 23),
        time_fname="runtime.json",  # sec
        mem_fname="memory.json",  # MB
        reps=7,
    )
