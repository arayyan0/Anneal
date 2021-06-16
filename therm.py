import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt('therm.out', dtype=float)

plt.plot(file[1:,0],file[1:,1])
plt.title(f'Temperature = {file[0,0]}')

plt.xlabel("Monte Carlo sweep")
plt.ylabel(f"Energy")

plt.show()
