import pandas as pd
import matplotlib.pyplot as plt

# === Load MSD CSV ===
filename = "MSD.csv"  

# Skip the first header line starting with "#"
df = pd.read_csv(filename, comment="#")

# === Plot MSD vs time_lag ===
plt.figure(figsize=(7,5))
plt.plot(df["time_lag"], df["MSD"], marker="o", linestyle="-")

plt.xlabel("Time lag")
plt.ylabel("Mean Squared Displacement (MSD)")
plt.xscale('log')
plt.yscale('log')
plt.title("MSD vs Time Lag")
plt.grid(True)

plt.tight_layout()
plt.show()
