import json
import numpy as np
import pylab as plt


with open("runtime.json", "r") as file:
    data = json.load(file)

lab1, lab2, lab3 = data.keys()
col1, col2, col3 = "tab:blue", "tab:orange", "tab:green"
alpha = 0.4
yticks = [0.01, 0.1, 1, 10, 100]
qubits = sorted(int(i) for i in data[lab3])
runtime1 = np.array([data[lab1][str(q)] for q in qubits if str(q) in data[lab1]])
runtime2 = np.array([data[lab2][str(q)] for q in qubits if str(q) in data[lab2]])
runtime3 = np.array([data[lab3][str(q)] for q in qubits])

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(
    qubits[: len(runtime1)],
    np.mean(runtime1, axis=1),
    marker="s",
    label=lab1,
    color=col1,
)
for k in range(runtime1.shape[1]):
    ax.plot(qubits[: len(runtime1)], runtime1[:, k], color=col1, alpha=alpha)
ax.plot(
    qubits[: len(runtime2)],
    np.mean(runtime2, axis=1),
    marker="s",
    label=lab2,
    color=col2,
)
for k in range(runtime2.shape[1]):
    ax.plot(qubits[: len(runtime2)], runtime2[:, k], color=col2, alpha=alpha)
ax.plot(qubits, np.mean(runtime3, axis=1), marker="s", label=lab3, color=col3)
for k in range(runtime3.shape[1]):
    ax.plot(qubits, runtime3[:, k], color=col3, alpha=alpha)
ax.set_xlabel("Qubits", fontsize=18)
ax.set_xticks(qubits)
ax.set_xticklabels(qubits, fontsize=14)
ax.set_ylabel("Runtime [s]", fontsize=18)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks, fontsize=14)
ax.set_yscale("log")
ax.legend(fontsize=18)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid()
fig.savefig("runtime.pdf")
plt.show()


with open("memory.json", "r") as file:
    data = json.load(file)

yticks = [0.1, 1, 10, 100, 1000]
memory1 = np.array([data[lab1][str(q)] for q in qubits if str(q) in data[lab1]])
memory2 = np.array([data[lab2][str(q)] for q in qubits if str(q) in data[lab2]])
memory3 = np.array([data[lab3][str(q)] for q in qubits])

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(
    qubits[: len(memory1)], np.mean(memory1, axis=1), marker="s", label=lab1, color=col1
)
for k in range(memory1.shape[1]):
    ax.plot(qubits[: len(memory1)], memory1[:, k], color=col1, alpha=alpha)
ax.plot(qubits, np.mean(memory2[: len(memory2)], axis=1), marker="s", label=lab2, color=col2)
for k in range(memory2.shape[1]):
    ax.plot(qubits[: len(memory2)], memory2[:, k], color=col2, alpha=alpha)
ax.plot(qubits, np.mean(memory3, axis=1), marker="s", label=lab3, color=col3)
for k in range(memory3.shape[1]):
    ax.plot(qubits, memory3[:, k], color=col3, alpha=alpha)
ax.set_xlabel("Qubits", fontsize=18)
ax.set_xticks(qubits)
ax.set_xticklabels(qubits, fontsize=14)
ax.set_ylabel("Memory [MB]", fontsize=18)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks, fontsize=14)
ax.set_yscale("log")
ax.legend(fontsize=18)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid()
fig.savefig("memory.pdf")
plt.show()
