from MDAnalysis import Universe
from MDContentAnalysis import MDContentAnalysis as MDCA

#Save folder for all files
saveFolderPath = "data"
exeFilepath="~/Desktop/WORK/Programs/hole2/exe/hole"

#Load MD-system into Universe object
Uni = Universe("TestTraj.pdb", "TestTraj.xtc")
protein=Uni.select_atoms('protein and not (name H*)')

#Initiation analysis
analysis = MDCA.TunnelAnalysis(Uni, 0, -1, 1, saveFolderPath, temp_in_memory=True)

start=[['A',112],['A',230],['A',343]] #S3
end=[['B',586],['B',266],['B',261]] #BindingSight

#Finding the tunnel
tunnels, tunnelFilters = analysis.tunnelAnalysis(protein=protein, 
                                                 start=start,
                                                 end=end,
                                                 hole2Filepath=exeFilepath,
                                                 random_seed=42,
                                                 saveTempfiles=True)

#Analysing water content inside the tunnel
waterU = Uni.select_atoms('resname SOL')
waterIds, waterNum = analysis.contentAnalysis(waterU, ['OW', 'HW1', 'HW2'], 'water', False)

waterOxygenIds = []
for ids in waterIds:
    waterOxygenIds.append(Uni.atoms[ids].select_atoms('name OW').ix) #Get only the oxygen's ids

waterOxygenDistribution  = analysis.atomDistribution(waterOxygenIds)



#Writing .txt files
import numpy as np
with open(f'waterNumPrFrame.txt','w+') as f:
    for i,j in zip(range(analysis.trajLen),analysis.contentNums['water']):
        f.write(f'{i}, {j}\n')
f.close()

with open(f'waterOxygenDistribution.txt','w+') as f:
    for i in waterOxygenDistribution:
        f.write(f'{i}\n')
f.close()

hist, bin_edges = np.histogram(waterOxygenDistribution, bins=500, range=[-10,50])
with open(f'waterOxygenHistogram.txt','w+') as f:
    for i, j in zip(bin_edges[0:-1],hist/analysis.trajLen):
        f.write(f'{i}, {j}\n')
f.close()

#How to read text files
histR = []
bin_edgesR = []
with open('waterOxygenHistogram.txt','r') as f:
    lines=f.read()
    listli=lines.split('\n')
    for i in listli[0:-1]:
        j,k = i.split(',')
        bin_edgesR.append(float(j))
        histR.append(float(k))
f.close()



#Ploting the water density distribution along the tunnel
from matplotlib import pyplot as plt
fig = plt.figure(layout="constrained")
gs = plt.GridSpec(1, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])

ax1.hist(waterOxygenDistribution, range=(-10,50), bins=500, density=True)
ax1.set_xlabel('Length along tunnel direction [Å]', fontsize=16)
ax1.set_ylabel('PDF', fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=12)
plt.margins(0)
ax1.grid()

fig.suptitle('Water Distribution Along Tunnel', fontsize=20)
plt.show()