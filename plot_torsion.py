import matplotlib.pyplot as plt
import pandas as pd

# Carga el archivo generado por torsion_histogram.c
data = pd.read_csv("torsion_phi.dat", header=None, names=["Phi_deg"])

# Histograma simple de ϕ en grados
plt.hist(data["Phi_deg"], bins=36, range=(-180, 180), edgecolor="black")

plt.xlabel("Dihedral angle φ (degrees)")
plt.ylabel("Frequency")
plt.title("Distribution of torsion angle φ")
plt.tight_layout()
plt.show()
