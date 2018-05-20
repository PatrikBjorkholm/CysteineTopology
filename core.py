
def NumberOfMembraneRegions(topology):
    numberTMs = 0;

    for i in range(1,len(topology)):
        if(topology[i] == 'M'):
            if (topology[i-1] != 'M'):
                numberTMs = numberTMs + 1

    if(topology[0] == 'M'):
        numberTMs = numberTMs + 1
    #print(str(numberTMs) + "\t" + topology.strip('\n'))

    return numberTMs


def OrganizeCSS_palm(rawData):
    Ids = []
    Pos = []

    for i in range(1,len(rawData)):
        Ids.append(rawData[i].split('|')[1])
        Pos.append(rawData[i].split('\t')[1])

    return Ids,Pos



#Function to Retrieve information from the file
def Retrieve_Data_From_File(Path):
    with open(Path,'r' ) as ins:
        array = []
        for line in ins:
            line.encode('utf-8')
            array.append(line)
    return array


def Function_to_SP(rawSP):
    sakf = ""
def FunctionToGetSP_pos(rawData, scampiName, ScampiTop,SeqName,Sequence,Ids,Pos):
    UniprotID = []
    UniprotSP_pos = []

    InnerSequence = []
    InnerSequenceName = []

    for i in range(1,len(rawData)):
        rawData[i] = rawData[i].replace("              ", " ")
        rawData[i] = rawData[i].replace("             ", " ")
        rawData[i] = rawData[i].replace("            ", " ")
        rawData[i] = rawData[i].replace("           ", " ")
        rawData[i] = rawData[i].replace("          ", " ")
        rawData[i] = rawData[i].replace("         ", " ")
        rawData[i] = rawData[i].replace("        ", " ")
        rawData[i] = rawData[i].replace("       ", " ")
        rawData[i] = rawData[i].replace("      ", " ")
        rawData[i] = rawData[i].replace("     ", " ")
        rawData[i] = rawData[i].replace("    ", " ")
        rawData[i] = rawData[i].replace("   ", " ")
        rawData[i] = rawData[i].replace("  ", " ")

        temp = rawData[i].split(" ")
        nameT = temp[0].split("_")[1]


#        print(nameT + "\t" + temp[2] + "\t" + temp[9])

# Get Proper AA index
        ProtID = ""
        protIndex = -1
        ScId = ""
        ScIndex = -1
        for k in range(len(SeqName)):
            TempName = SeqName[k].split('|')[1]
            TempSc = SCname[k].split('|')[1]
            if(nameT == (TempName)):
                ProtID = TempName
                protIndex = k
            if(nameT == TempSc):
                ScId = TempSc
                ScIndex = k
                #print(TempSc + "\t" + nameT)

        #print(ProtID + "\t" + ScId)
        SP_start = 0
        if(temp[9] == 'Y'):
            SP_start = int(temp[2])

        CSS_palm_counter = 0
        CSS_palm_counter_All = 0




        tempSequence = ""
        seq = Sequence[protIndex]
        top = ScampiTop[ScIndex]

        for k in range(len(Ids)):
            if(Ids[k] == nameT):
                CSS_palm_counter_All = CSS_palm_counter_All + 1
                if(top[int(Pos[k])] == 'I'):
                    CSS_palm_counter = CSS_palm_counter+1


        CysCount = 0
        for l in range(SP_start, len(top)):
            if(top[l] == 'I'):

                tempSequence = tempSequence + seq[l]
                if(seq[l] == 'C'):
                    CysCount = CysCount + 1
        if(len(tempSequence)>0):
            # create function to see if
            numTMs = NumberOfMembraneRegions(top)
            print(ScId + "\t" + str(CysCount) + "\t" + temp[9] + "\t" + str(CSS_palm_counter) + "\t" + str(CSS_palm_counter_All) + "\t" + str(numTMs))
            #print(">"+ScId)
            #print(tempSequence)

# Create an integrated approach by creating proteins of inner mebrane domains.

IgnoreID = "Q8WXI7"
AminoAcidData = Retrieve_Data_From_File("Data/E_golgi.fasta")
ScampiData = Retrieve_Data_From_File("Data/E_golgi_scampi.txt")
SPdata = Retrieve_Data_From_File("Data/E_golgi_SP.txt")
CSS_data_raw = Retrieve_Data_From_File("Data/E_golgi_css.txt")




# here I get Amino Acid Count
#-----------------------------

AAcount = 0
for i in range(len(AminoAcidData)):
    TempString = AminoAcidData[i]
    if(TempString[0] == '>'):
        AAcount = AAcount + 1

print(AAcount)
AAname = []
AAseq = []

for i in range(AAcount):
    AAname.append("")
    AAseq.append("")
AAcount = -1

for i in range(len(AminoAcidData)):
    TempString = AminoAcidData[i]
    if(TempString[0] == '>'):
        AAcount = AAcount + 1
        AAname[AAcount] = TempString
    else:
        AAseq[AAcount] = AAseq[AAcount] + TempString.strip('\n')

SCcount = 0
for i in range(len(ScampiData)):
    TempString = ScampiData[i]
    if(TempString[0] == '>'):
        SCcount = SCcount + 1
print(SCcount)

SCname = []
SCseq = []

for i in range(SCcount):
    SCname.append("")
    SCseq.append("")
SCcount = -1
for i in range(len(ScampiData)):
    TempString = ScampiData[i]
    if(TempString[0] == '>'):
        SCcount = SCcount + 1
        SCname[SCcount] = TempString
    else:
        SCseq[SCcount] = SCseq[SCcount] + TempString

# Sequences Imported and organized


for i in range(len(AAname)):
    ProtName = AAname[i].split('|')[1]

    SCindex = -1
    for k in range(len(SCname)):
        TempName = SCname[k].split('|')[1]
        if(TempName == ProtName):
            SCindex = k
    ProtSeq = AAseq[i].strip('\n')
    Topology = SCseq[SCindex].strip('\n')
    innerCysteine = 0
    for k in range(len(ProtSeq)):
        if(Topology[k] == 'I'):
            if(ProtSeq[k] == 'C'):
                innerCysteine = innerCysteine + 1

   # print(ProtName +"\t"+ str(innerCysteine))
  #  print(ProtSeq)
  #  print(Topology)

#Function_to_SP(SPdata)
Ids, Pos = OrganizeCSS_palm(CSS_data_raw)

FunctionToGetSP_pos(SPdata ,SCname,SCseq,AAname,AAseq, Ids, Pos)