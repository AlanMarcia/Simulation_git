import sys
import csv
from collections import defaultdict

import matplotlib.pyplot as plt


# Uso: python plot_proton_trajectories.py [csv_path]
# Se non passi nulla, usa proton_trajectories.csv nella cartella corrente.
csv_path = sys.argv[1] if len(sys.argv) > 1 else "proton_trajectories.csv"

# Leggi il CSV e raggruppa per proton_id
trajectories = defaultdict(list)
with open(csv_path, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    required = {"proton_id", "pos_x_m", "pos_y_m"}
    if not required.issubset(reader.fieldnames or {}):
        missing = required - set(reader.fieldnames or {})
        raise SystemExit(f"Colonne mancanti nel CSV: {', '.join(sorted(missing))}")
    for row in reader:
        pid = int(row["proton_id"])
        x = float(row["pos_x_m"])
        y = float(row["pos_y_m"])
        trajectories[pid].append((x, y))

if not trajectories:
    raise SystemExit("Nessuna traiettoria da plottare.")

# Plot semplice delle traiettorie
plt.figure(figsize=(8, 6))
for pid in sorted(trajectories.keys()):
    pts = trajectories[pid]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    plt.plot(xs, ys, linewidth=1.0, label=f"proton {pid}")

plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Proton trajectories")
plt.grid(True, linestyle="--", alpha=0.4)
plt.axis("equal")
if len(trajectories) <= 20:
    plt.legend()

plt.tight_layout()
plt.show()
