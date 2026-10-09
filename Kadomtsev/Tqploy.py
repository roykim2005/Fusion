from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()

data_file = SCRIPT_DIR / "kadomtsev.dat"

df = pd.read_csv(data_file, sep=r"\s+", comment="#")
df.columns = [c.strip().lstrip("#").strip() for c in df.columns]

r = df["r/a"].values if "r/a" in df.columns else df.iloc[:, 0].values
q_pre = df["q_pre"].values if "q_pre" in df.columns else df.iloc[:, 1].values
q_post = df["q_post"].values if "q_post" in df.columns else df.iloc[:, 2].values
T_pre = df["T_pre"].values if "T_pre" in df.columns else df.iloc[:, 3].values
T_post = df["T_post"].values if "T_post" in df.columns else df.iloc[:, 4].values

mask = r < 0.85

r_plot = r[mask]
q_pre_plot = q_pre[mask]
q_post_plot = q_post[mask]
T_pre_plot = T_pre[mask]
T_post_plot = T_post[mask]

T_early = T_post_plot + 0.25 * (T_pre_plot - T_post_plot)
T_mid   = T_post_plot + 0.50 * (T_pre_plot - T_post_plot)
T_late  = T_post_plot + 0.75 * (T_pre_plot - T_post_plot)
q_early = q_post_plot + 0.25 * (q_pre_plot - q_post_plot)
q_mid   = q_post_plot + 0.50 * (q_pre_plot - q_post_plot)
q_late  = q_post_plot + 0.75 * (q_pre_plot - q_post_plot)

stages = [
    (T_pre_plot, q_pre_plot, "Pre-crash", "firebrick", "-"),
    (T_post_plot, q_post_plot, "Post-crash", "navy", "--"),
    (T_early, q_early, "25%", "darkgreen", "-."),
    (T_mid, q_mid, "50%", "darkorange", ":"),
    (T_late, q_late, "L75%", "crimson", "-"),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

for T_stage, q_stage, label, color, style in stages:
    ax1.plot(r_plot, T_stage, label=label, color=color, linestyle=style, linewidth=2)
    ax2.plot(r_plot, q_stage, label=label, color=color, linestyle=style, linewidth=2)

ax1.set_xlabel("r/a", fontsize=12)
ax1.set_ylabel("T(r) [keV]", fontsize=12)
ax1.set_title("r vs T(r)", fontsize=13)
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(fontsize=9, loc="upper right")
ax1.set_ylim(0, np.max(T_pre_plot) * 1.15)

ax1.set_xlabel("r/a", fontsize=12)
ax2.set_ylabel("q(r)", fontsize=12)
ax2.set_title("r vs q(r)", fontsize=13)
ax2.axhline(1.0, color="gray", linestyle=":", alpha=0.7, label="q = 1.0")
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.legend(fontsize=9, loc="lower right")

plt.suptitle("T(r) and q(r)", fontsize=15)
plt.tight_layout()
plt.show()
