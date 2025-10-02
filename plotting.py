import matplotlib.pyplot as plt
import pandas as pd

### PLOTTING THE VELOCITIES IN A HISTOGRAM ###
# Load data into pandas
data = pd.read_csv("speeds.dat", header=None, names=["Value"])

# Plot histogram
plt.hist(data["Value"], bins=20, edgecolor="black")
plt.xlabel("Speed [m/s]")
plt.ylabel("Frequency")
plt.title("Speed of the particiles from a MD simulation")
plt.tight_layout()
plt.show()


### PLOTTING OF THE DIAGNOSTICS ###

df2 = pd.read_csv("diagnostics.csv")

fig, axes = plt.subplots(2, 1, figsize=(10,8), sharex=True)

# 1) Energies
axes[0].plot(df2["time"], df2["kinetic_energy"], label="Kinetic Energy", color="red")
axes[0].plot(df2["time"], df2["potential_energy"], label="Potential Energy", color="blue")
axes[0].plot(df2["time"], df2["total_energy"], label="Total Energy", color="green")

axes[0].set_ylabel("Energy (kJ/mol)")
axes[0].set_title("Energies vs Time")
axes[0].legend()
axes[0].grid(True)

# 2) Temperature
axes[1].plot(df2["time"], df2["temperature"], label="Temperature", color="orange")

axes[1].set_xlabel("Time (ps)")
axes[1].set_ylabel("Temperature (K)")
axes[1].set_title("Temperature vs Time")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()