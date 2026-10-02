from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()
FILE_PATH = SCRIPT_DIR / "time_history.dat"

data = pd.read_csv(FILE_PATH, sep=r"\s+", comment="#", names=["t", "T0", "q0"])
steady_state = data[data["t"] * 1000 > 20]
fig, ax1 = plt.subplots(figsize=(10, 5))
ax1.set_xlabel("Time [ms]")
ax1.set_ylabel("Central Temperature T(0) [keV]")
ax1.plot(steady_state["t"] * 1000, steady_state["T0"], linewidth=1.5, color = 'red')
ax1.tick_params(axis="y")
t_min, t_max = steady_state["T0"].min() - 5, steady_state["T0"].max() + 5
ax1.set_ylim(t_min, t_max)
ax2 = ax1.twinx()
ax2.set_ylabel("Safety Factor q(0)")
ax2.plot(steady_state["t"] * 1000, steady_state["q0"], linestyle="--", alpha=0.7)
ax2.tick_params(axis="y")

plt.title("Kadomtsev sawtooth crash")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
