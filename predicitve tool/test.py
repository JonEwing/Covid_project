import pandas as pd
from numpy import genfromtxt
import shutil

#Put Sample1, Sample2, Sample3, Sample4, ...
samples = [940927, 516773, 475645, 475607]
#Put mvx, mvy, ...
name = ["mv16", "mv5"]

#Open data
while samples:
    mutarray = []
    with open("../data/AlignedFilteredUS.fdi") as fp:
        line = fp.readline()
        holdarray = []
        previous = ""
        while line:
            linesplit = line.split(';')
            if linesplit[0] == "LINK_TAXON1":
                lista = [linesplit[1].replace(" ", "").strip(), linesplit[3].replace(" ", "").strip(), int(linesplit[9])]
                for x in range(int(linesplit[9])):
                    if x == 0:
                        counter = 11
                        lista += [linesplit[counter]]
                    else:
                        counter += 6
                        lista += [linesplit[counter]]
                mutarray.append(lista)
            line = fp.readline()
    fp.close()
    pd.DataFrame(mutarray).to_csv("mutations.csv", index = False)

#Read Sequence data
    seqarray = []
    with open("../data/AlignedFilteredUS.nex") as fp:
        while True:
            line = fp.readline()
            line = line.split(" ")
            try:
                seqarray.append([line[0], line[1]])
            except: 
                break

    fp.close()
    pd.DataFrame(seqarray).to_csv("sequences.csv", index = False)

    for x in mutarray:
        if str(x[0]) in str(samples[0]) and str(x[1]) in str(name[0]):
            muts1 = x
        if str(x[0]) in str(samples[1]) and str(x[1]) in str(name[0]):
            muts2 = x

    for x in seqarray:
        if str(samples[0]) in str(x[0]):
            bases1 = x[1]
        if  str(samples[1]) in str(x[0]):
            bases2 = x[1]

    matchingmuts = False
    matchval = 0
    for x in range(len(muts1)):
        if x > 2:
            for y in range(len(muts2)):
                if y > 2:
                    if muts1[x] == muts2[y]:
                        matchingmuts = True
                        matchval = muts1[x]

#Output
    if matchingmuts == True:
        print("Matching Mutation at", matchval, "Pick a new sequence")
    else:
        f = open("out.txt", "a")
        f.write("Using sequences "+ str(samples[0]) + " and " + str(samples[1]) + " for " + name[0] + "\n")
        if muts1[2] <= muts2[2]:
            mvseq = bases1
            for x in muts1[3:]:
                f.write("Sequence Change at "+ x + mvseq[int(x) - 1]+ " -> "+ bases2[int(x) - 1] + "\n")
                mvseq = mvseq[:int(x) - 1] + bases2[int(x) - 1] + mvseq[int(x):]
        else:
            mvseq = bases2
            for x in muts2[3:]:
                f.write("Sequence Change at "+ x + mvseq[int(x) - 1]+ " -> "+ bases1[int(x) - 1] + "\n")
                mvseq = mvseq[:int(x) - 1] + bases1[int(x) - 1] + mvseq[int(x):]

    f.write("\n" + name[0] + " " + mvseq)
    f.write("\n#########################################################################################\n\n")
    f.close()

    f = open("../data/AlignedFilteredUS.nex", "a")
    f.write(name[0] + " " + mvseq)

    samples.pop(0)
    samples.pop(0)
    name.pop(0)