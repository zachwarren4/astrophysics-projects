import numpy as np
import matplotlib.pyplot as plt

t = np.array([10**7,10**8,10**9,2*(10**9),5*(10**9),10**10])

A=np.log(((t**-.4) +2)/(t**-.4))-.65
B=3.5*np.log10(t)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(A,B,'g')
for xy in zip(A,B):
    ax.annotate('10 Myr', xy=(6.49,24.5), textcoords = 'offset points')
    ax.annotate('100 Myr', xy=(7.41,28), textcoords = 'offset points')
    ax.annotate('1 Gyr', xy=(7.69,29), textcoords = 'offset points')
    ax.annotate('2 Gyr', xy=(8.33,31.5), textcoords = 'offset points')
    ax.annotate('5 Gyr', xy=(8.61,32.55), textcoords = 'offset points')
    ax.annotate('10 Gyr', xy=(9.25,35), textcoords = 'offset points')

plt.title('g-r vs Abs Magnitude Difference')
plt.xlabel('g-r')
plt.ylabel('Abs Magnitude Difference')
plt.show()