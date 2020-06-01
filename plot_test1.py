import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 4)

ax[0].plot(range(10), 'r') #row=0, col=0
ax[1].plot(range(10), 'b') #row=1, col=0
ax[2].plot(range(10), 'g') #row=0, col=1
ax[3].plot(range(10), 'k') #row=1, col=1
plt.show()
