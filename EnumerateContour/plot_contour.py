"""
    plot countour
"""
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from numpy.core.shape_base import block
import pandas as pd

plt.rcParams["font.family"] = "Times New Roman"
data_1 = pd.read_csv("DemandData_50_005.txt")
data_2 = pd.read_csv("DemandData_50_010.txt")

x1 = np.unique(data_1["b52"]).tolist()
y1 = np.unique(data_1["b73"]).tolist()
x2 = np.unique(data_2["b52"]).tolist()
y2 = np.unique(data_2["b73"]).tolist()

Z1 = []
Z2 = []
for yv in y1:
    z1 = []
    for xv in x1:
        for r in range(0, data_1.shape[0]):
            if data_1["b73"][r] == yv and data_1["b52"][r] == xv:
                z1.append(data_1["sum"][r])
    Z1.append(z1)

for yv in y2:
    z2 = []
    for xv in x2:
        for r in range(0, data_2.shape[0]):
            if data_2["b73"][r] == yv and data_2["b52"][r] == xv:
                z2.append(data_2["sum"][r])
    Z2.append(z2)

f, ax = plt.subplots(ncols=2, nrows=1, figsize=(9,4))
X1, Y1 = np.meshgrid(x1, y1)
X2, Y2 = np.meshgrid(x2, y2)
lvs = [0,150,300,450,600,750,900,1050,1300]
CS = ax[0].contour(X1, Y1, Z1,levels=lvs)
ax[0].clabel(CS, inline=True, fmt='%4.1f',fontsize=8,use_clabeltext=True)
CS = ax[1].contour(X2, Y2, Z2,levels=lvs)
ax[1].clabel(CS, inline=True, fmt='%4.1f',fontsize=8,use_clabeltext=True)

ax[0].set_title("(a) Scale Paraemter = 0.05",fontsize=10,fontname='Times New Roman')
ax[1].set_title("(b) Scale Parameter = 0.10",fontsize=10,fontname='Times New Roman')

ax[0].set_xlabel("Capacity of Link 5-2",fontsize=10,fontname="Times New Roman")
ax[1].set_xlabel("Capacity of Link 5-2",fontsize=10,fontname="Times New Roman")

ax[0].set_ylabel("Capacity of Link 7-3",fontsize = 10, fontname = 'Times New Roman')
ax[1].set_ylabel("Capacity of Link 7-3",fontsize = 10, fontname = 'Times New Roman')

plt.tight_layout()
plt.savefig("Contour.png",bbox_inches='tight',dpi=600)
plt.show(block = False)
plt.pause(2)
plt.close()
plt.show()
