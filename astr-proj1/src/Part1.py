import numpy as np
import matplotlib.pyplot as plt
import math

t = np.arange(0, 10**9,10**8)

plt.plot(t,3.5*np.log10(t),'r')

plt.title('Difference in Abs Magnitude over Time')
plt.xlabel('Time')
plt.ylabel('M(t)-M(0)')
plt.show()

