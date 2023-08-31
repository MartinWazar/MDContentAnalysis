from MDAnalysis import Universe
from MDContentAnalysis import MDContentAnalysis as MDCA

#Save folder for all files
saveFolderPath = "data"
exeFilepath="~/Desktop/WORK/Programs/hole2/exe/hole"

#Load MD-system into Universe object
Uni = Universe("TestTraj.pdb", "TestTraj.xtc")
protein=Uni.select_atoms('protein and not (name H*)')

#Initiation analysis
analysis = MDCA.TunnelAnalysis(Uni, 0, -1, 1, saveFolderPath)

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
waterDistrebution  = analysis.atomPDF(waterIds)



#Ploting the water density distribution along the tunnel
from matplotlib import pyplot as plt
fig = plt.figure(layout="constrained")
gs = plt.GridSpec(1, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])

ax1.hist(waterDistrebution, range=(-10,50), bins=500, density=True)
ax1.set_xlabel('Length along tunnel direction [Å]', fontsize=16)
ax1.set_ylabel('PDF', fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=12)
plt.margins(0)
ax1.grid()

fig.suptitle('Water Distribution Along Tunnel', fontsize=20)
plt.show()