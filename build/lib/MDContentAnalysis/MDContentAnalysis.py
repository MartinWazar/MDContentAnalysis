# %%
from MDAnalysis import Universe, Writer, Merge
import numpy as np
import math
import os

#Function that returns the mean possition of the specified residues
def residueMeanPos(protein, residueList:list, usechain:bool = True) -> list:
    """
    Returns the mean of the residues center of geomitry.
    If usechain = True: residueList have the form [['ChainID1',resid1],['ChainID2',resid2],...]
    If usechain = False: residueList have the form [resid1,resid2,...]
    """
    sumResidueCoord = [0,0,0]
    for resseq in residueList:
        if usechain:
            selection = f"chainID {resseq[0]} and resid {resseq[1]}"
            sumResidueCoord += protein.select_atoms(selection).center_of_geometry()
        else:
            selection = f"resid {resseq}"
            sumResidueCoord += protein.select_atoms(selection).center_of_geometry()

    meanResidueCoord = sumResidueCoord/len(residueList)
    return meanResidueCoord

#Projection of u on v.
proj_len = lambda u,v : (np.dot(u, v)/np.linalg.norm(v))

#pojects the atoms in a frame on to the vector and with tunnel start at 0
def projectAtoms(atoms, vec:list, start:list):
    startTile = np.tile(start,[len(atoms),1])
    atomVecPos = proj_len(atoms.positions-startTile,vec)
    return atomVecPos

#Class collecting all the methods
class TunnelAnalysis:
    def __init__(self,
                 Uni,
                 startFrame:int=0, 
                 endFrame:int=-1, 
                 interval:int=1, 
                 saveFolderPath:str=None):
        #Seting up a TunnelAnalysis object
        self.Uni = Uni

        self.startFrame = startFrame
        self.endFrame = endFrame
        self.interval = interval

        self.lastframe = Uni.trajectory[startFrame:endFrame:interval][-1].frame
        self.frameOrder = int(math.ceil(math.log10(self.lastframe)))
        self.trajLen = len(Uni.trajectory[startFrame:endFrame:interval])

        self.frameNum = 0
        self.tunnelIdList = [None]*self.trajLen
        self.rawTunnelList = [None]*self.trajLen
        self.tunnelList = [None]*self.trajLen
        self.posSave = [None]*self.trajLen
        self.tunnelVect = [None]*self.trajLen

        self.contentIdLists = {}
        self.contentIdDicts = {}
        self.contentNums = {}

        if saveFolderPath == None:
            self.saveFolderPath = "save"
        else:
            self.saveFolderPath = saveFolderPath

        if not os.path.exists(self.saveFolderPath):
            os.mkdir(self.saveFolderPath)

        self.tempFolderPath = f"{saveFolderPath}/tempfiles"
        if not os.path.exists(self.tempFolderPath):
            os.mkdir(self.tempFolderPath)

    def filterChannel(self,
                      centerPoint, 
                      startpos, 
                      index:int, 
                      findChannel:bool=False):
        """
        Filters channel file from HOLE2.0.
        This is done to get a tunnel from the channel HOLE outputs.
        Hole writes a lot of atoms with residu number -888
        this just signifies the atoms it is outside the protein.

        There is made one asumption: The tunnel is within a sphere formed by the starting 
        and ending point of the tunnel and where the spheres center is halv way inbetween.
        If findChannel = True: This asumption is not made and the channel is returned

        Returns the atoms of the tunnel (tunnelStructure) and the
        indexes of the atoms that define the tunnel (tunnelIdList).
        """
        #Load rawTunnel structure from HOLE2.0 into Universe object
        tunnelPath = f"{self.tempFolderPath}/tunnel{self.frameNum:0{self.frameOrder}d}.pdb"
        self.rawTunnelList[index] = Universe(tunnelPath) 

        #Filter channel removing spheres outside the protein
        tunnelStructure = self.rawTunnelList[index].select_atoms("not resnum -888")

        #Filter channel into tunnel
        if not findChannel:
            centerPointStr = ' '.join(map(str, centerPoint))
            radius = np.linalg.norm(centerPoint-startpos)
            tunnelStructure = tunnelStructure.select_atoms(f"point {centerPointStr} {radius}")
        tunnelStructure.atoms.write(tunnelPath)
        tunnelIdList = tunnelStructure.ix
        return tunnelIdList, tunnelStructure

    def tunnelAnalysis(self,
                       protein,
                       start:list, 
                       end:list=None,
                       hole2Filepath:str=None,
                       mole2Filepath:str=None, 
                       xmlFilepath:str="MOLE2.xml",
                       random_seed:int=None, 
                       saveTempfiles:bool=False,
                       usechain:bool=True,
                       findChannel:bool=False,
                       tunnelVect=[None,None,None]):
        """
        Runs HOLE2.0 or MOLE2.5 on trajectory dependent on which filepath is specified.
        Then combines all the tunnels in one trajectory pdb file.
        Only one of the filepaths for either HOLE or MOLE can be specified.

        protein must be an atomgroupe contaning the protein structure to be analysed.

        start/end should be a list of residues that defines the start/end point of the tunnel.
        The usechain should be set depending on topology have chainID's or not.
        If usechain = True: start/end expects the form [['ChainID1',resid1],['ChainID2',resid2],...]
        If usechain = False: start/end expects  the form [resid1,resid2,...]

        startFrame = 0: start form first frame of trajectory.
        endFrame = -1: end on last frame of trajectory.
        interval = 1: use every frame, 2 would mean every other frame of trajectory.

        If saveFolderPath is not specified a default foulder named 'save' 
        will be created. xmlFilepath is to specify the name and saving 
        place of the xml file created when running MOLE.

        hole2Filepath: abselute path to the hole.exe file
        mole2Filepath: relative or abselute path to the hole.exe file

        random_seed: Seed for the pseudo-random functions in HOLE.

        If saveTempfiles = True: temperary tunnel files will NOT be deleted.
        If saveTempfiles = False: temperary tunnel files will be deleted.
        saveTempfiles = False, may give [WinError 32] on windows (MOLE doesn't 
        close it's files properly).

        If findChannel = True: Then a channel is found instead of a tunnel, 
        this option is only for HOLE2.0, MOLE2.5 will fail using this option.
        When findChannel = True a tunnelVect must be specified but end dose 
        not need to be.

        The minimum input for finding a tunnel is (protein, start, end) and 
        hole2Filepath or mole2Filepath. The minimum input for finding a channel  
        is (protein, start, tunnelVect, findChannel=True) and hole2Filepath.

        TunnelTraj.pdb contains the tunnel/channel for each frame. There is some 
        fake atoms in the file. The fake atoms is added to keep the number of 
        atoms constant for all the frames. This is to ensure the file can be read 
        by MDanalysis again and so it can be displayed with vmd. All the fake 
        atoms are added at the coordinates [0 0 0] and have 'name', 'resname', 
        and 'segid' set to 'QZW'.

        Returnes a list of universes contaning the tunnel/channel found for each 
        frame in the trajectory and a list of indices that in the case of using 
        HOLE2.0 must be used as a filter to get the decired tunnel/channel output. 
        Like this: rawTunnelList[i].atoms[tunnelIdList[i]] where i is the frame.
        """
        #Set Variables to use in function
        Uni = self.Uni
        lastframe = self.lastframe
        trajLen = self.trajLen
        saveFolderPath = self.saveFolderPath
        tempFolderPath = self.tempFolderPath
        startFrame = self.startFrame
        endFrame = self.endFrame
        interval = self.interval
        
        proteinTempFile = f"{tempFolderPath}/ProteinTemp.pdb"

        #Error handeling when finding channels.
        if findChannel:
            if tunnelVect[0]==None:
                print('Please specify tunnelVect when looking for channels (sorry for the confusing name)')
                return None, None
            elif (hole2Filepath==None) and (mole2Filepath!=None):
                print('HOLE2.0 should be used when finding channels insted of MOLE2.5')
                return None, None

        #Error handeling for not specifing path for HOLE or MOLE
        if (hole2Filepath==None) and (mole2Filepath==None):
            print('No filepath for HOLE or MOLE specified!')
            return None, None
        
        #Error handeling for trying to run both HOLE and MOLE
        elif(hole2Filepath!=None) and (mole2Filepath!=None):
            print('Only one filepath for either HOLE or MOLE can be specified!')
            return None, None

        #HOLE2.0 Preperation
        elif hole2Filepath!=None:
            print('Using HOLE2.0')
            from MDAnalysis.analysis import hole2
            fakeName = 'SPH'

        #MOLE2.5 Preperation
        elif mole2Filepath!=None:
            print('Using MOLE2.5')
            from MDContentAnalysis import xmlMOLE
            WorkingDirectory = f"{saveFolderPath}/MOLE2_output"
            if not os.path.exists(WorkingDirectory):
                os.makedirs(WorkingDirectory)

            xmlMOLE.generateXML(xmlFilepath, proteinTempFile, "./"+WorkingDirectory+"/", start, end)
            pathDir = f"{WorkingDirectory}/pdb/profile/path_1.pdb"
            fakeName = 'TUN'

        #Looping over trajectory
        Uni.trajectory.rewind()
        for ts in Uni.trajectory[startFrame:endFrame:interval]:
            i = ts.frame
            self.frameNum = i
            index = int((i-startFrame)/interval)

            #Writing temperary protine file
            protein.write(proteinTempFile)

            startPos = residueMeanPos(protein, start, usechain)
            if not findChannel:
                endPos = residueMeanPos(protein, end, usechain)
                centerPoint = (startPos+endPos)/2
                self.posSave[index] = [centerPoint,startPos,endPos]
                tunnelVect = endPos-startPos
            else:
                self.posSave[index][1] = startPos
                centerPoint = startPos
            self.tunnelVect[index] = tunnelVect

            tunnelDir = f"{tempFolderPath}/tunnel{i:0{self.frameOrder}d}.pdb"

            #Mole Analasys
            if mole2Filepath != None:
                os.system(f'cmd /c "{mole2Filepath} {xmlFilepath}"')
                
                if os.path.exists(pathDir):
                    os.replace(pathDir,tunnelDir)
                
                self.rawTunnelList[index] = Universe(tunnelDir)
                self.tunnelList[index] = self.rawTunnelList[index].atoms
                self.tunnelIdList[index] = self.tunnelList[index].ix

            #Hole Analasys
            if hole2Filepath != None:
                hole2.hole(proteinTempFile, 
                        executable=hole2Filepath, 
                        cpoint=centerPoint, 
                        cvect=tunnelVect, 
                        random_seed=random_seed)

                if os.path.exists("hole.sph"):
                    os.replace("hole.sph",tunnelDir)

                self.tunnelIdList[index], self.tunnelList[index] = self.filterChannel(centerPoint, 
                                                                                      startPos, 
                                                                                      index, 
                                                                                      findChannel)
                
                os.remove("hole.out")
                os.remove("simple2.rad")

            #Remove temporary files if specified
            if saveTempfiles != True:
                os.remove(tunnelDir)
                os.remove(proteinTempFile)
            else:
                os.replace(proteinTempFile, f"{tempFolderPath}/ProteinTemp{i:0{self.frameOrder}d}.pdb")

            if i != lastframe:
                print(f'Frame: {i}')
            else:
                print(f'Frame: {i}')
                print(f'Tunnel calculations Done!')

        #Finding the maximum amount of spheres in the trajectory
        maxSPH = 0
        for i in range(trajLen):
            numSPH = len(self.tunnelList[i])
            if maxSPH < numSPH:
                maxSPH = numSPH

        #Making a universe with fake atoms
        tempU = Universe.empty(maxSPH,
                            n_residues=1,
                            atom_resindex=[0]*maxSPH,
                            residue_segindex=[0],
                            trajectory=True) # necessary for adding coordinates
        
        tempU.add_TopologyAttr('name', ['x']*maxSPH)
        tempU.add_TopologyAttr('resname', [fakeName])
        tempU.add_TopologyAttr('segid', ['QZW'])
        tempU.add_TopologyAttr('tempfactor', [0]*maxSPH)

        #Writing a tunnel trajectory file with a constant amount of spheres/atoms for each frame.
        Uni.trajectory.rewind()
        with Writer("TunnelTraj.pdb", multiframe=True) as pdb:
            for ts in Uni.trajectory[startFrame:endFrame:interval]:
                i = ts.frame
                index = int((i-startFrame)/interval)

                tunnel = self.tunnelList[index]
                if len(tunnel) == maxSPH:
                    pdb.write(tunnel)
                else:
                    test2 = Merge(tunnel, tempU.atoms[len(tunnel)-1:-1])
                    pdb.write(test2.atoms)

                if i != lastframe:
                    print(f'Frame: {i}')
                else:
                    print(f'Frame: {i}')
                    print('Trajectory File Done!')
        os.replace("TunnelTraj.pdb",f"{saveFolderPath}/TunnelTraj.pdb")

        return self.rawTunnelList, self.tunnelIdList

    def loadTunnel(self,
                   tunnelPath:str,
                   protein,
                   start:list, 
                   end:list,
                   usechain:bool=True):
        
        tunnelU = Universe(tunnelPath)
        tunnel = tunnelU.select_atoms('not (segid QZW)', updating=True)
        for ts in tunnelU.trajectory:
            i = ts.frame
            self.tunnelList[i] = tunnel

        Uni = self.Uni
        for ts in Uni.trajectory[startFrame:endFrame:interval]:
            i = ts.frame
            self.frameNum = i
            index = int((i-startFrame)/interval)
            startPos = residueMeanPos(protein, start, usechain)
            endPos = residueMeanPos(protein, end, usechain)
            centerPoint = (startPos+endPos)/2
            self.posSave[index] = [centerPoint,startPos,endPos]
            self.tunnelVect[index] = endPos-startPos
        print('Tunnel loaded')
        return 0

    def closeContent(self,
                     Atoms, 
                     centerPoint, 
                     startPos,
                     moleculeSelector:list=['OW','HW1','HW2'],
                     moleculeName:str='water',
                     saveTempfiles:bool=False):
        """
        Finds molecules that overlab with a sphere with center in centerPoint
        and radius |centerPoint-startPos|. A molecule overlabs with the sphere
        if any of the atoms are inside the sphere. This is to remove most of the 
        molecules in the system to speed up following calculations.

        Returns an atomgroupe with the molecules that overlab with the sphere 
        and a dictionary to lookup the id of atom in Atoms from the id of atoms 
        writen in the temperary file.
        """
        moleculeLen = len(moleculeSelector)

        centerPointStr =  ' '.join(map(str, centerPoint))
        radius = np.linalg.norm(centerPoint-startPos)
        tempCloseContent = Atoms.select_atoms(f"point {centerPointStr} {radius}")
        
        PartsList = [[None]]*moleculeLen
        tempList = []
        for i, atomName in zip(range(moleculeLen),moleculeSelector):
            PartsList[i] = tempCloseContent.select_atoms(f"name {atomName}").ix
            for id in PartsList[i]:
                ids = [None]*moleculeLen
                for j in range(moleculeLen):
                    ids[j] = id+j-i
                tempList.extend(ids)

        idlist = sorted(set(tempList))
        closeContent = Atoms.universe.atoms[idlist]
        if saveTempfiles:
            if len(idlist) != 0:
                closeContent.write(f"{self.tempFolderPath}/close{moleculeName}{self.frameNum:0{self.frameOrder}d}.pdb")

        idDict = dict(zip(range(len(idlist)), idlist))
        return idDict, closeContent

    def contentInTunnel(self, 
                        closeContent, 
                        tunnelStructure,
                        moleculeSelector:list=['OW','HW1','HW2'],
                        moleculeName:str='water',
                        saveTempfiles:bool=False):
        """
        Finds the moleclules which atleast one atom overlabing with the tunnel spheres.

        Returns the a list with the id's of the atoms in the molecules that overlab with 
        the tunnel. It also returns the number of molecules that overlab with the tunnel.
        """
        moleculeLen = len(moleculeSelector)

        #Find the over overlabing atoms and add their indexes to a list.
        idList = []
        ids = np.empty(1)
        for sph in tunnelStructure:
            coord = sph.position
            radius = sph.tempfactor
            coordStr = ' '.join(map(str, coord))
            ids = closeContent.select_atoms(f"point {coordStr} {radius}").ix
            if ids.size != 0:
                if ids[0] != -1:
                    idList.extend(ids)

        #Remove repeats and add the indexes of the other atoms for each molecule
        idList = sorted(set(idList))

        tempList = []
        for id in idList:
            for i in range(moleculeLen):
                if closeContent.universe.atoms[id].name == moleculeSelector[i]:
                    ids = [None]*moleculeLen
                    for j in range(moleculeLen):
                        ids[j] = id+j-i
                    tempList.extend(ids)
                
        #Remove repeated indexes
        idList = sorted(set(tempList))

        #Get the number of molecules and write a file
        contentNum = len(idList)/moleculeLen
        if saveTempfiles:
            if len(idList) != 0:
                contentInTunnelPath = f"{self.tempFolderPath}/{moleculeName}InTunnel{self.frameNum:0{self.frameOrder}d}.pdb"
                closeContent.universe.atoms[idList].write(contentInTunnelPath)
        #print(f'Number of molecules: {contentNum}')

        return idList, contentNum

    def contentAnalysis(self,
                        contentU,
                        moleculeSelector:list=['OW', 'HW1', 'HW2'],
                        moleculeName:str='water',
                        saveTempfiles:bool=False):
        """
        Finds molecules that overlab with the tunnel found in tunnelAnalysis() for 
        each frame. Do NOT work for channels.
        If findChannel was specified to True in tunnelAnalysis() this function will crash.
        startFrame, endFrame, and interval must be the same as used in tunnelAnalysis().

        IF saveFolderPath is not specified a default foulder named 'save' will be created.
        saveFolderPath don't have to be the same as used in tunnelAnalysis().

        If saveTempfiles = True: temperary tunnel files will NOT be deleted.
        If saveTempfiles = False: temperary tunnel files will be deleted.

        {moleculeName}InTunnelTraj.pdb is written in the save folder and contains the atoms 
        of the molecules which ovarlab with the tunne in each frame. Note that 
        the number of molecules is not guaranteed to be the same in each frame. 
        {moleculeName}InTunnelTraj.pdb can therfor not be loaded using MDanalysis.

        {moleculeName}TunnelTraj.pdb saved in the same folder on the other hand have the same 
        amount of atoms in each frame. For all the frames all the molecules that 
        overlabe with the tunnel at som point in the trajectory is writing to 
        {moleculeName}TunnelTraj.pdb. 

        Returns a list of lists with the id's of the atoms in the molecules that 
        overlab with the tunnel. One list of id's for each frame. It also returns a 
        list of numbers that decribes the number of molecules overlabing with 
        the tunnel in each frame.
        """

        self.contentIdLists[moleculeName] = [None]*self.trajLen
        self.contentIdDicts[moleculeName] = [None]*self.trajLen
        self.contentNums[moleculeName] = [None]*self.trajLen

        Uni = self.Uni
        lastframe = self.lastframe
        saveFolderPath = self.saveFolderPath
        startFrame = self.startFrame
        endFrame = self.endFrame
        interval = self.interval

        #Content analysis
        TrackIdfile = []
        Uni.trajectory.rewind()
        with Writer(f"{moleculeName}InTunnelTraj.pdb", multiframe=True) as pdb:
            for ts in Uni.trajectory[startFrame:endFrame:interval]:
                i = ts.frame
                index = int((i-startFrame)/interval)
                self.frameNum = i
                
                #Load data from tunnelAnalysis()
                centerPoint = self.posSave[index][0]
                startPos = self.posSave[index][1]
                tunnel = self.tunnelList[index]

                #Find molecules close to the tunnel
                idDict, closeContent = self.closeContent(contentU, 
                                                  centerPoint, 
                                                  startPos,
                                                  moleculeSelector,
                                                  moleculeName,  
                                                  saveTempfiles)
                
                #Find molecules overlabing with the tunnel
                idList, self.contentNums[moleculeName][index] = self.contentInTunnel(closeContent, 
                                                                                     tunnel,
                                                                                     moleculeSelector,
                                                                                     moleculeName,
                                                                                     saveTempfiles)

                #Save data to global variables and 
                self.contentIdDicts[moleculeName][index] = idDict
                self.contentIdLists[moleculeName][index] = idList
                if len(idList) != 0:
                    pdb.write(Uni.atoms[idList])
                else:
                    pdb.write(contentU[-3:-1])

                #Track which molecules is in the tunnel at som point
                TrackIdfile.extend(idList)

                if i != lastframe:
                    if index%10 == 0:
                        print(f'Frame: {i}')
                else:
                    print(f'Frame: {i}')
                    print('Content analysis Done!')
        os.replace(f"{moleculeName}InTunnelTraj.pdb",f"{saveFolderPath}/{moleculeName}InTunnelTraj.pdb")


        # Save trajectory file with
        totContent = 0
        for num in self.contentNums[moleculeName]:
            totContent += num

        if not (totContent < 1):
            SortTrackIdfile = sorted(set(TrackIdfile))
            Uni.trajectory.rewind()
            with Writer(f"{moleculeName}TunnelTraj.pdb", multiframe=True) as pdb:
                for ts in Uni.trajectory[startFrame:endFrame:interval]:
                    i = ts.frame
                    index = int((i-startFrame)/interval)
                    pdb.write(Uni.atoms[SortTrackIdfile])
                    if i != lastframe:
                        if index%100 == 0:
                            print(f'Frame: {i}')
                    else:
                        print(f'Frame: {i}')
                        print('Trajectory File Done!')
            os.replace(f"{moleculeName}TunnelTraj.pdb",f"{saveFolderPath}/{moleculeName}TunnelTraj.pdb")
        else:
            print('No content in channel at any time')

        return self.contentIdLists[moleculeName], self.contentNums[moleculeName]

    def atomPDF(self,
                IdList:list):
        """
        Projects the position of atoms onto the general directions of the protein tunnel.
        Uni: The Universe contaning the MD simulation
        IdList: A list of lists. Where each list contains the id's for the atoms to be 
        projected for each frame.

        startFrame, endFrame, and interval must fit the dimentions of the IdList.

        Returnes a list contaning how far the atoms is along the general tunnel direction.
        """
        
        Uni = self.Uni
        lastframe = self.lastframe
        startFrame = self.startFrame
        endFrame = self.endFrame
        interval = self.interval

        atomTotal = 0
        for Ids in IdList:
            atomTotal += len(Ids)

        atomDensList = np.zeros(atomTotal)

        atomNum = 0
        atomCount = 0
        for ts in Uni.trajectory[startFrame:endFrame:interval]:
            i = ts.frame
            index = int((i-startFrame)/interval)

            atomU = Uni.atoms[IdList[index]]
            startPos = self.posSave[index][1]
            tunnelVec = self.tunnelVect[index]

            atomCount += atomNum
            atomNum = len(IdList[index])
            atomDensList[atomCount:atomCount+atomNum] = projectAtoms(atomU,tunnelVec,startPos)

            if i != lastframe:
                if index%100 == 0:
                    print(f'Frame: {i}')
            else:
                print(f'Frame: {i}')
                print('Atom Coordinate Transformations Calculated!')
        return atomDensList


# %% Running Tunnel analysis (example)
if __name__ == '__main__':
    #References to analsys program
    hole2Filepath = "~/Desktop/WORK/Programs/hole2/exe/hole"    #(when using Hole 2.0 to analyse)
    #mole2Filepath = "Programs\Mole2\mole2.exe"                 #(when using Mole 2.5 to analyse)

    #Save folder for all files
    saveFolderPath = "data/analasys/FilReliansTest/Test2_Full/xtc"

    #Load MD-system
    SysFilepath = "data/short_md_center.xtc"
    TopFilepath = "data/Component_top/short_md_c_sys.pdb"
    Uni = Universe(TopFilepath, SysFilepath) #Load protein structure into Universe object

    #Split system in protein and molecules
    proteinU = Uni.select_atoms('protein')# and not (name H*)')

    #Tunnel starting and ending point
    start = [['A',112],['A',230],['A',343]] #S3
    end = [['B',586],['B',266],['B',261]] #BindingSight

    #Variables for selecting which frames in the trajectory to process
    startFrame = 0      #0 is first frame
    endFrame = -1       #-1 is last frame
    interval = 1        #1 is every frame

    #Initiation analysis
    analysis = TunnelAnalysis(Uni,
                              startFrame,
                              endFrame,
                              interval,
                              saveFolderPath)

    #Finding the tunnel
    tunnels, tunnelFilters = analysis.tunnelAnalysis(protein=proteinU, 
                                                        start=start, 
                                                        end=end, 
                                                        hole2Filepath=hole2Filepath,
                                                        saveTempfiles=True)
    
    #Analysing water content inside the tunnel
    waterU = Uni.select_atoms('resname SOL')
    waterIds, waterNum = analysis.contentAnalysis(waterU, ['OW', 'HW1', 'HW2'], 'water', False)
    waterDistrebution  = analysis.atomPDF(waterIds)

    #Analysing the potassium content inside the tunnel
    potassiumU = Uni.select_atoms('resname POT')
    KIds, KNum = analysis.contentAnalysis(potassiumU, ['POT'], 'potassiumT', False)
    potasiumDistrebution  = analysis.atomPDF(KIds)

    #Ploting number of watermolecules inside the tunnel along the trajectory
    from matplotlib import pyplot as plt
    fig = plt.figure(layout="constrained")
    gs = plt.GridSpec(1, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])

    ax1.plot(np.arange(len(waterNum)), waterNum)
    ax1.set_xlabel('Length along tunnel direction [Å]', fontsize=16)
    ax1.set_ylabel('PDF', fontsize=16)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    plt.margins(0)
    ax1.grid()

    fig.suptitle('Water Distribution Along Tunnel', fontsize=20)
    plt.show()

    #Ploting the water density distribution alon the tunnel
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

    #The analysis oblect can 
    import pydataman as pd
    pd.save(f"analysisXTC", analysis)


