from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()

data_file = SCRIPT_DIR / "kadomtsev.dat"
print("Reading:", data_file)

df = pd.read_csv(data_file, sep=r"\s+", comment="#", header=None)

r = df.iloc[:, 0].values
q_pre = df.iloc[:, 1].values
q_post = df.iloc[:, 2].values
T_pre = df.iloc[:, 3].values
T_post = df.iloc[:, 4].values
changed = (np.abs(q_pre - q_post) > 1e-6) | (np.abs(T_pre - T_post) > 1e-6)
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
    (T_late, q_late, "75%", "crimson", "-"),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

for T_stage, q_stage, label, color, style in stages:
    ax1.plot(r_plot, T_stage, label=label, color=color, linestyle=style, linewidth=2)
    ax2.plot(r_plot, q_stage, label=label, color=color, linestyle=style, linewidth=2)

ax1.set_xlabel("r/a", fontsize=12)
ax1.set_ylabel("T(r)/T_c [keV]", fontsize=12)
ax1.set_title("r vs T(r)/T_c", fontsize=13)
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(fontsize=9, loc="upper right")
ax1.set_ylim(0, np.max(T_pre_plot) * 1.15)

ax1.set_xlabel("r/a", fontsize=12)
ax2.set_ylabel("q(r)", fontsize=12)
ax2.set_title("r vs q(r)", fontsize=13)
ax2.axhline(1.0, color="gray", linestyle=":", alpha=0.7, label="q = 1.0")
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.legend(fontsize=9, loc="lower right")


'''
r0_file = r[changed].max()
ax1.set_xlim(0, 1.5 * r0_file)
ax2.set_xlim(0, 1.5 * r0_file)
q_zoom = q_pre_plot[r_plot < 1.5 * r0_file]
q_zoom_post = q_post_plot[r_plot < 1.5 * r0_file]
lo = min(q_zoom.min(), q_zoom_post.min(), 1.0)
hi = max(q_zoom.max(), q_zoom_post.max(), 1.0)
pad = 0.1 * (hi - lo) if hi > lo else 0.05
ax2.set_ylim(lo - pad, hi + pad)'''

plt.suptitle("T(r)/T_c and q(r)", fontsize=15)
plt.tight_layout()
plt.show()


'''10/16 - resistivity timescale too large due to high temp, try to do like 5 kev 2 kev. normalize. write up what i did so far.'''
