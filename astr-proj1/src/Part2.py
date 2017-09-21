import numpy as np
import matplotlib.pyplot as plt
import math

t = np.arange(0, 10**9,10**8)

plt.plot(t,np.log(((t**-.4) +2)/(t**-.4))-.65,'b')

plt.title('g-r value over time')
plt.xlabel('Time')
plt.ylabel('g-r')
plt.show()