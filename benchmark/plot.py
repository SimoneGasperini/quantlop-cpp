import json
import numpy as np
import pylab as plt


def plot_series(ax, qubits, values, label, color, alpha):
    means = np.mean(values, axis=1)
    ax.plot(qubits, means, marker="s", label=label, color=color)
    for k in range(values.shape[1]):
        ax.plot(qubits, values[:, k], color=color, alpha=alpha)


with open("runtime.json", "r") as file:
    data = json.load(file)

labels = list(data.keys())
colors = ("tab:blue", "tab:orange", "tab:green")
alpha = 0.4
yticks = [0.01, 0.1, 1, 10, 100]
qubits = sorted({int(q) for label in labels for q in data[label]})

fig, ax = plt.subplots(figsize=(9, 6))
for label, color in zip(labels, colors):
    label_qubits = [q for q in qubits if str(q) in data[label]]
    values = np.array([data[label][str(q)] for q in label_qubits])
    plot_series(ax, label_qubits, values, label, color, alpha)
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

fig, ax = plt.subplots(figsize=(9, 6))
for label, color in zip(labels, colors):
    label_qubits = [q for q in qubits if str(q) in data[label]]
    values = np.array([data[label][str(q)] for q in label_qubits])
    plot_series(ax, label_qubits, values, label, color, alpha)
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
