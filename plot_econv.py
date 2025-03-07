import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

colors=sns.color_palette('tab20b')

stoints = np.arange(2,6.1,1)
stoenergies = [-2.702157146211008,-2.8077839566141956,-2.8357259896578073,-2.843769817330659,-2.8462920947813597]
stotimes = [0.10, 0.11, 0.12,  0.17, 0.25]

hs=[]
energies=[]
times=[]

with open('helium_conv.txt', 'r') as f:
    for line in f:
        h,time,e = line.split('\t')
        energies.append(float(e))
        times.append(float(time))
        hs.append(float(h))

fig, axs = plt.subplots(1,2)

axs[0].plot(stoints,np.array(stoenergies)+2.8617, label='STO-nG', color=colors[5], linestyle=':',marker='o')
axs[0].plot(stoints,np.array(energies)+2.8617, label='DI$(\Delta)$', color=colors[1], linestyle=':', marker='X')
axs[0].axhline(y=0, color='k', ls='--', linewidth = 0.5, label="$E_{\\text{HF}}=-2.861(7)$")
axs[0].legend(fontsize='small')
axs[0].set_xlabel("n for STO-nG or $-\\log (\\Delta)$ for DI")
axs[0].set_ylabel("$E-E_{\\text{HF}}$ ($E_{\\text{h}}$)")

axs[1].plot(stoints,stotimes, label='STO-nG', color=colors[5], linestyle=':',marker='o' )
axs[1].plot(stoints,times, label='DI$(\Delta)$', color=colors[1], linestyle=':', marker='X')
axs[1].set_yscale('log')
axs[1].legend(fontsize='small')
axs[1].set_xticks(stoints)
axs[1].set_xlabel("n for STO-nG or $-\\log (\\Delta)$ for DI")
axs[1].set_ylabel('runtime (seconds)')

plt.subplots_adjust(wspace=0.35)
plt.subplots_adjust(right=0.99)
plt.subplots_adjust(left=0.1)
#plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(top = 0.99)
plt.savefig('econv.pdf')
