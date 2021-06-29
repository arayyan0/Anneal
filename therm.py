import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt('therm.out', dtype=float)

end = file.shape[0]

eavg, e2avg = file[end-1,1:2+1]
# print()

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)

ax1.plot(file[0:end-1,0],file[0:end-1,1])
ax1.axhline(eavg,color='grey')

ax2.plot(file[0:end-1,0],file[0:end-1,2])
ax2.axhline(e2avg,color='grey')
cv = (e2avg - eavg**2)/(36**2)/(0.05)**2
plt.title(f'cv={cv}')
plt.show()
