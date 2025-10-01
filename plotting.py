import matplotlib.pyplot as plt
import pandas as pd

# Load data into pandas
data = pd.read_csv("speeds.dat", header=None, names=["Value"])

# Plot histogram
plt.hist(data["Value"], bins=20, edgecolor="black")
plt.xlabel("Speed [m/s]")
plt.ylabel("Frequency")
plt.title("Speed of the particiles from a MD simulation")
plt.tight_layout()
plt.show()