# %%
import xml.etree.ElementTree as ET

def generateXML(filename:str, 
                proteinFile:str, 
                WorkingDirectory:str, 
                Start:list, 
                End:list,
                ProbeRadius:float=3,
                InteriorThreshold:float=1.25):
    """
    Generates a XML file for MOLE. It is made to finde tunnel using (Path).
    -----------------------------------##----------------------------------
    filename: The name and path given to the XML file
    proteinFile: File name and path of the protein structure (must be PDB format)
    WorkingDirectory: Name given to the directory where MOLEs temperary and output files will be saved

    Start: List of residues defining the start of the tunnel
    End: List of residues defining the end of the tunnel
    Start and End parameters must have the form [['ChainID1',resid1],['ChainID2',resid2],...]
    
    ProbeRadius: The value for ProbeRadius in MOLE given in Å
    InteriorThreshold: The value for InteriorThreshold in MOLE given in Å
    InteriorThreshold must be smaller than ProbeRadius for the algorithm to work right.
    """

    #This is the parent (root) tag onto 
    #which other tags are created
    data = ET.Element('Tunnels')
    
    #Adding a subtags inside root tag
    element1 = ET.SubElement(data, 'Input')
    element2 = ET.SubElement(data, 'WorkingDirectory')
    element3 = ET.SubElement(data, 'Params')
    element4 = ET.SubElement(data, 'Export')
    element5 = ET.SubElement(data, 'Paths')
    
    #Adding text to elements
    element1.text = proteinFile
    element2.text = WorkingDirectory

    #Adding subtags inside:
    #`Params`
    s_elem2 = ET.SubElement(element3, 'Cavity')
    #`Export`
    s_elem3 = ET.SubElement(element4, 'Formats')
    s_elem4 = ET.SubElement(element4, 'Types')
    #`Paths`
    s_elem5 = ET.SubElement(element5, 'Path')
    #`Path`
    ss_elem1 = ET.SubElement(s_elem5, 'Start')
    ss_elem2 = ET.SubElement(s_elem5, 'End')
    #`Start`
    #A tag is added for each element in Start
    sss_elem1 = []
    for i in range(len(Start)):
        sss_elem1.append(ET.SubElement(ss_elem1, 'Residue'))
        sss_elem1[i].set('Chain', Start[i][0])
        sss_elem1[i].set('SequenceNumber', f'{Start[i][1]}')
    #`End`
    #A tag is added for each element in End
    sss_elem2 = []
    for i in range(len(End)):
        sss_elem2.append(ET.SubElement(ss_elem2, 'Residue'))
        sss_elem2[i].set('Chain', End[i][0])
        sss_elem2[i].set('SequenceNumber', f'{End[i][1]}')

    #Adding attributes to the tags: 
    #'Cavity'
    s_elem2.set('ProbeRadius', str(ProbeRadius))
    s_elem2.set('InteriorThreshold', str(InteriorThreshold))
    #'Formats'
    s_elem3.set('PDBProfile', '1')
    #'Types'
    s_elem4.set('Tunnels', '0')
    
    #Converting the xml data to byte object,
    #for allowing flushing data to file stream
    b_xml = ET.tostring(data)
    
    #Opening a file under the name filename,
    #with operation mode `wb` (write + binary)
    with open(filename, "wb") as f:
        f.write(b_xml)
    return 0

def updateXML(filename:str, Input:str):
    tree = ET.ElementTree(file=filename)
    root = tree.getroot()

    for oldInput in root.iter("Input"):
        oldInput.text = Input

    b_xml = ET.tostring(root)
    
    # Opening a file, with operation mode
    # `wb` (write + binary)
    with open(filename, "wb") as f:
        f.write(b_xml)
    return 0

if __name__ == '__main__':
    #1TQN test of generateXML
    filename = "TestXML.xml"
    protein = "./Programs/Mole2/tests/1TQN.pdb.gz"
    WorkingDirectory = "./test_output_xml/"
    Start = [["A",308],["A",309]]
    generateXML(filename, protein, WorkingDirectory, Start, [])

# %%
