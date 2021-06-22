
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker
import pandas as pd
import math

data = pd.read_csv("Results.txt",sep="\t")
y1 = [v for v in data["0.3"] if not (math.isnan(v))]
x1 = np.arange(0,len(y1),1)
y2 = [v for v in data["0.6"] if not (math.isnan(v))]
x2 = np.arange(0,len(y2),1)
y3 = [v for v in data["0.9"] if not (math.isnan(v))]
x3 = np.arange(0,len(y3),1)
y4 = [v for v in data["10"] if not (math.isnan(v))]
x4 = np.arange(0,len(y4),1)
y5 = [v for v in data["30"] if not (math.isnan(v))]
x5 = np.arange(0,len(y5),1)
y6 = [v for v in data["50"] if not (math.isnan(v))]
x6 = np.arange(0,len(y6),1)
x=[x1,x2,x3,x4,x5,x6]
y=[y1,y2,y3,y4,y5,y6]
f, ax = plt.subplots(ncols=2, nrows=1, figsize=(9,4))
f1=sns.lineplot(x[0],y[0],label=r'$\rho_{r}=0.3$', ax=ax[0])
sns.lineplot(x[1],y[1],label=r'$\rho_{r}=0.6$', ax=ax[0])
sns.lineplot(x[2],y[2],label=r'$\rho_{r}=0.9$', ax=ax[0])
f2=sns.lineplot(x[3],y[3],label=r'$o_{r,min}^{+,k}=o_{r,max}^{+,k}=10$',ax=ax[1])
sns.lineplot(x[4],y[4],label=r'$o_{r,min}^{+,k}=o_{r,max}^{+,k}=30$',ax=ax[1])
sns.lineplot(x[5],y[5],label=r'$o_{r,min}^{+,k}=o_{r,max}^{+,k}=50$',ax=ax[1])
ax[0].legend(prop={"family":"Times New Roman","size":8})
ax[1].legend(prop={"family":"Times New Roman","size":8})
ax[0].set_title("(a) Adaptive Step Size",fontsize=10,fontname='Times New Roman')
ax[1].set_title("(b) Constant Step Size",fontsize=10,fontname='Times New Roman')
ax[0].set_ylabel("Convergence Measurement",fontsize = 10, fontname = 'Times New Roman')
ax[1].set_ylabel("Convergence Measurement",fontsize = 10, fontname = 'Times New Roman')
ax[0].set_xlabel("No. of Iterations",fontsize=10,fontname="Times New Roman")
ax[1].set_xlabel("No. of Iterations",fontsize=10,fontname="Times New Roman")
ax[0].set_ylim([0,0.3])
ax[1].set_ylim([0,0.7])
for i in [0,1]:
    xmajorFormatter = plt.FormatStrFormatter('%.0f')
    ymajorFormatter = plt.FormatStrFormatter('%.2f')
    xtick = ax[i].get_xticks()
    ytick = ax[i].get_yticks()
    ax[i].set_xticklabels(xtick, fontsize=8,fontname='Times New Roman')
    ax[i].set_xticklabels(xtick, fontsize=8,fontname='Times New Roman')
    ax[i].set_yticklabels(ytick, fontsize=8,fontname='Times New Roman')
    ax[i].yaxis.set_major_formatter(ymajorFormatter)
    ax[i].xaxis.set_major_formatter(xmajorFormatter)
plt.tight_layout()
plt.savefig("2fig.png",bbox_inches='tight',dpi=600)
plt.show(block = False)
plt.pause(2)
plt.close()
