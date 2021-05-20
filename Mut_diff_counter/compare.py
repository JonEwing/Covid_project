import pandas as pd
lista = ["Founderlist.fasta"]
import math

#DNA to amino Acids
def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for x in range(len(seq)):
        if x%3 == 0:
            try:
                codon = seq[x] + seq[x+1] + seq[x+2]
                if codon in table:
                    protein += table[codon]
                else:
                    protein += "-"
            except:
                protein += "-"
    return protein

#Read data
df = pd.read_csv("SARS-CoV-2 Reference Sequence_All Regions.csv")
genename = df["name"]
genlocs = df["start"]
genloce = df["end"]

for x in lista:
    fp = open(x, 'r')

    holdarray = fp.readlines() 

    holdstr = ""
    holdarr = []
    for x in range(len(holdarray)):
        if x == 0:
            name = holdarray[x][1:len(holdarray[x])-1]
            if name == 'e':
                name = holdarray[x][1:len(holdarray[x])-4]
        elif ">" in holdarray[x]:
            holdarr.append([name,holdstr])
            name = holdarray[x][1:len(holdarray[x])-1]
            holdstr = ""
        else:
            holdstr += holdarray[x][:len(holdarray[x])-1]
    holdarr.append([name,holdstr])

#Find Mutation Diffrences
totalarray = []
for x in range(len(holdarr)):
    if x != 0:
        partalarray = []
        for y in range(len(holdarr[x][1])):
            if holdarr[0][1][y] != holdarr[x][1][y]:
                tempstr = ""
                for z in range(len(genlocs)):
                    if y + 77 >= genlocs[z] and y + 78 <= genloce[z]:
                        tempstr += genename[z] + ", "       
                if tempstr:
                    partalarray.append(tempstr + "at " + str(y + 78) + ": " + holdarr[0][1][y] + "->" + holdarr[x][1][y])
                else:
                    partalarray.append(str(y + 77) + ": " + holdarr[0][1][y] + "->" + holdarr[x][1][y])
        
        partalarray = [holdarr[x][0]] + partalarray
        totalarray.append(partalarray)

df = pd.DataFrame(list(totalarray))
df.to_csv('mutations.csv', index=False)

#Remove bases until ORF1
replaceh = []
for x in holdarr:
    replaceh.append([x[0],x[1][x[1].find('ATGGAGAGC'):]])

#Convert DNA to amino acids
proarray = []
for x in replaceh:
    totalpro = ""
    for y in range(len(genename)):
        x[1] = x[1].replace("\n", "") 
        x[1] = x[1].replace("\r", "")
        if "ORF" in genename[y]:
            if "ab" in genename[y]:
                pro = translate(x[1][genlocs[y] - 266 : 13469 - 266])
                totalpro += pro
                pro = translate(x[1][13468 - 266 : genloce[y] - 266])
                totalpro += pro
            else:
                pro = translate(x[1][genlocs[y] - 266 : genloce[y] - 266])
                totalpro += pro
    proarray.append([x[0],totalpro])

#Read protein data
df = pd.DataFrame(list(proarray))
df.to_csv('protein.csv', index=False)

#Find diffrences in protein samples
totalarray = []
for x in range(len(proarray)):
    if x != 0:
        partalarray = []
        for y in range(len(proarray[x][1])):
            if proarray[0][1][y] != proarray[x][1][y]:
                tempstr = ""
                for z in range(len(genlocs)):
                    if y + 88 >= genlocs[z] / 3 and y + 88 <= genloce[z] / 3:
                        tempstr += genename[z] + ", "         
                if tempstr:
                    partalarray.append(tempstr + "at " + str(y + 88) + ": " + proarray[0][1][y] + "->" + proarray[x][1][y])
        
        partalarray = [proarray[x][0]] + partalarray
        totalarray.append(partalarray)

#Output
df = pd.DataFrame(list(totalarray))
df.to_csv('protein_diff.csv', index=False)