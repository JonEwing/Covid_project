fp = open('data/5000out.fasta', 'r')
Nfp = open('update_file.nex', 'w')

holdarray = fp.readlines() 

holdstr = ""
holdarr = []
for x in holdarray:
    if ">Reference" in x:
        name = x[1:7]
    elif ">" in x:
        holdarr.append([name,holdstr])
        name = x[9:len(x)-1]
        holdstr = ""
    else:
        holdstr += x[:len(x)-1]
holdarr.append([name,holdstr])

Nfp.write("#NEXUS\n\n")
Nfp.write("begin data;\n")
ntax = len(holdarr)
nchar = len(holdarr[0][1])
Nfp.write("dimensions ntax=" + str(ntax) + " nchar=" + str(nchar) + ";\n")
Nfp.write("format datatype=dna gap=-;\n")
Nfp.write("matrix\n")

for x in holdarr:
    Nfp.write(x[0] + " " + x[1] + "\n")

Nfp.write(";\n")
Nfp.write("end;\n")
