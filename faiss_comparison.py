import matplotlib.pyplot as plt

BASELINE_RUNTIME = 94.39
BASELINE_N_CLUSTERS = 1057

# Data from benchmarking runs (time in seconds, number of clusters)
data = {
    "flat k-NN (k=2)": (12.13, 78712),
    "flat k-NN (k=4)": (12.13, 20258),
    "flat k-NN (k=16)": (12.23, 1750),
    "flat k-NN (k=64)": (12.54, 1185),
    "flat k-NN (k=256)": (14.42, 1087),
    "flat k-NN (k=1024)": (21.88, 1064),
    "flat ranged": (17.86, 1057),
}

fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xscale("log")
ax.set_yscale("log")
ax.grid(True, which="both", linestyle="--", linewidth=0.5)

for index_type, (time, n_clusters) in data.items():
    ax.scatter(n_clusters, time, marker="o")
    # place label slightly to the right of the point
    ax.text(n_clusters * 1.03, time * 1.03, index_type, va="center", fontsize=8)

ax.axhline(BASELINE_RUNTIME, color="gray", linestyle="--")
ax.axvline(BASELINE_N_CLUSTERS, color="gray", linestyle="--",)
ax.text(
    BASELINE_N_CLUSTERS * 1.03,
    BASELINE_RUNTIME / 1.03,
    "baseline (full distance matrix with numpy)",
    va="center",
    fontsize=8,
)

ax.set_xlabel("time (s)")
ax.set_ylabel("n_clusters")

fig.tight_layout()

plt.show()
