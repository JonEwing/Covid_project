import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from networkx.drawing.nx_agraph import graphviz_layout
import community as community_louvain

############################################################################## adj list
filepath = "../data/AlignedFilteredUS.fdi"
fout = open("adj.adjlist", "w")
fout2 = open("directed.adjlist", "w")
holdlist = []
dirarray = []

#Read FDI file
dataarray = []
with open(filepath) as fp:
    line = fp.readline()
    holdarray = []
    previous = ""
    while line:
        linesplit = line.split(';')
        if linesplit[0] == "LINK_TAXON1":
            dataarray.append([linesplit[1].replace(" ", "").strip(), linesplit[3].replace(" ", "").strip(), int(linesplit[9])])
            if previous == linesplit[1].replace(" ", ""):
                holdarray.append(linesplit[3].replace(" ", ""))

            else:
                if holdarray:
                    fout.write(" ".join(holdarray)+ "\n")

                    holdlist.append([holdarray[0], len(holdarray) - 1])

                    for x in range(len(holdarray)):
                        if x != 0:
                            fout2.write(holdarray[0] + " " + holdarray[x] + "\n")
                            dirarray.append([holdarray[0],holdarray[x]])
                holdarray = []
                holdarray.append(linesplit[1].replace(" ", ""))
                holdarray.append(linesplit[3].replace(" ", ""))
                previous = linesplit[1].replace(" ", "")
        
        line = fp.readline()

df = pd.DataFrame(holdlist, columns=["Name", "link Amount"])
df.to_csv('edge_count.csv', index=False)

fout.close()
fp.close()

pd.DataFrame(dataarray, columns=['Node1', 'Node2', 'Mutation ammount']).to_csv("mutation_number.csv", index = False)

########################################################################################## Graph Setup

df = pd.read_csv("../data/FinalUSMeta.csv")

names = df["Accession ID"].to_list()
symptoms = df["Symptoms"].to_list()

G = nx.Graph()
for x in dataarray:
    G.add_nodes_from([x[0],x[1]])
    G.add_edge(x[0], x[1], m = x[2])

size_map = []
for node in G:
    match = False
    for x in holdlist:
        if x[0] == node:
            match = True
            size_map.append((x[1]/3) * 100)
    if match == False:
        size_map.append(100)

posa = graphviz_layout(G)

##################################################################################### States

location = df["Location"].to_list()
states = []
for x in location:
    temp = x.split(' / ')[2]
    states.append(temp)

color_map = []
for node in G:
    match = False

    if node == "REFERE":
        color_map.append('Yellow')
        match = True
    elif node == "ENGLAN":
        color_map.append('Gold')
        match = True
    elif node == "SOUTHA":
        color_map.append('khaki')
        match = True
    elif node == "BRAZIL":
        color_map.append('darkkhaki')
        match = True

    for x in range(len(names)):
        if node in names[x]:
            if "California" in states[x]:
                    color_map.append("Red")
            elif "Montana" in states[x]:
                    color_map.append("lawngreen")
            elif "New York" in states[x]:
                    color_map.append("Purple")
            elif "Puerto Rico" in states[x]:
                    color_map.append("Pink")
            elif "South Carolina" in states[x]:
                    color_map.append("Blue")
            elif "Louisiana" in states[x]:
                    color_map.append("Cyan")
            elif "Texas" in states[x]:
                    color_map.append("Brown")
            else:
                    color_map.append("darkGreen")
            match = True

    if match == False:
        color_map.append('black')   

plt.figure(figsize=(64,36))
nx.draw(G, pos = posa, node_color=color_map, node_size = size_map, font_weight='bold')
nx.draw_networkx_edge_labels(G, pos = posa)
plt.savefig("output/location.png")
plt.clf()
nx.draw(G, pos = posa, with_labels = True, font_size = 7, node_color="white", edge_color = "White")
plt.savefig("output/names.png")
plt.clf()


########################################################################################### Dates

datedash = df["Collection date"].to_list()
dates = []
for x in datedash:
    temps = x.split('-')
    y = int(temps[0])
    m = int(temps[1])
    if y == 2021:
        m += 12
    dates.append(m)


color_map = []
for node in G:
    match = False

    if node == "REFERE":
        color_map.append('Yellow')
        match = True
    elif node == "ENGLAN":
        color_map.append('Gold')
        match = True
    elif node == "SOUTHA":
        color_map.append('khaki')
        match = True
    elif node == "BRAZIL":
        color_map.append('darkkhaki')
        match = True

    for x in range(len(names)):
        if node in names[x]:
            if 2 == dates[x] or 3 == dates[x]:
                    color_map.append("purple")
            elif 4 == dates[x]:
                    color_map.append("lawngreen")
            elif 5 == dates[x] or 6 == dates[x]:
                    color_map.append("lightcoral")
            elif 7 == dates[x]:
                    color_map.append("Pink")
            elif 8 == dates[x] or 9 == dates[x]:
                    color_map.append("blue")
            elif 10 == dates[x] or 11 == dates[x]:
                    color_map.append("cyan")
            elif 12 == dates[x]:
                    color_map.append("fuchsia")
            elif 13 == dates[x] or 14 == dates[x]:
                    color_map.append("darkred")
            else:
                    color_map.append("darkgreen")
            match = True

    if match == False:
        color_map.append('black')   

plt.figure(figsize=(64,36))
nx.draw(G, pos = posa, node_color=color_map, node_size = size_map, font_weight='bold')
plt.savefig("output/dates.png")
plt.clf()

########################################################################################### Symptamitc

symps = df["Symptoms"].to_list()

color_map = []
for node in G:
    match = False

    if node == "REFERE":
        color_map.append('Yellow')
        match = True
    elif node == "ENGLAN":
        color_map.append('Gold')
        match = True
    elif node == "SOUTHA":
        color_map.append('khaki')
        match = True
    elif node == "BRAZIL":
        color_map.append('darkkhaki')
        match = True

    for x in range(len(names)):
        if node in names[x]:
            if symps[x] == "Symptomatic Moderate to Severe Symptoms":
                    color_map.append("Red")
            elif symps[x] == "Asymptomatic or Mild Symptoms":
                    color_map.append("Blue")
            else:
                    color_map.append("green")
            match = True

    if match == False:
        color_map.append('black')   

plt.figure(figsize=(64,36))
nx.draw(G, pos = posa, node_color=color_map, node_size = size_map, font_weight='bold')
plt.savefig("output/symptom.png")
plt.clf()

########################################################################################### Outcome

outs = df["Outcome"].to_list()

color_map = []
for node in G:
    match = False

    if node == "REFERE":
        color_map.append('Yellow')
        match = True
    elif node == "ENGLAN":
        color_map.append('Gold')
        match = True
    elif node == "SOUTHA":
        color_map.append('khaki')
        match = True
    elif node == "BRAZIL":
        color_map.append('darkkhaki')
        match = True

    for x in range(len(names)):
        if node in names[x]:
            if outs[x] == "Alive":
                    color_map.append("Blue")
            elif outs[x] == "Deceased":
                    color_map.append("purple")
            else:
                    color_map.append("green")
            match = True

    if match == False:
        color_map.append('black')   

plt.figure(figsize=(64,36))
nx.draw(G, pos = posa, node_color=color_map, node_size = size_map, font_weight='bold')
plt.savefig("output/outcome.png")
plt.clf()

########################################################################################### Community Detection

partition = community_louvain.best_partition(G, resolution = 2.725)

values = list(partition.values()) 
keys = list(partition.keys()) 

percentile_list = pd.DataFrame({'Sample Name': list(keys), 'Community': list(values)})
percentile_list.to_csv("community_val_15.csv", index = False)

cmap = cm.get_cmap('tab20', max(partition.values()) + 1)
nx.draw_networkx_nodes(G, posa, partition.keys(), node_size=size_map, cmap=cmap, node_color=values)
nx.draw_networkx_edges(G, posa, alpha=0.5)
plt.savefig("output/community_15.png")
plt.clf()

################################################################### Clustering

comm = pd.read_csv("communities.csv")
comnames = comm["Sample Name"].to_list()
comunity = comm["Community"].to_list()

for ccounter in range(15):
    G = nx.Graph()
    for x in dataarray:
        node1 = False
        node2 = False
        for y in range(len(comnames)):
            if x[0] == comnames[y] and comunity[y] == ccounter:
                node1 = True
            if x[1] == comnames[y] and comunity[y] == ccounter:
                node2 = True
        if node1 == True and node2 == True:
            G.add_nodes_from([x[0],x[1]])
            G.add_edge(x[0], x[1])

    size_map = []
    for node in G:
        match = False
        for x in holdlist:
            if x[0] == node:
                match = True
                size_map.append((x[1]/3) * 100)
        if match == False:
            size_map.append(100)

############################################################################ ids

    plt.figure(figsize=(16,9))
    nx.draw(G, pos = posa, with_labels = True, font_size = 10, node_color="white", edge_color = "lightgray")
    plt.savefig("clustering/id/cluster" + str(ccounter + 1) + ".png")
    plt.clf()

############################################################################ Symps

    color_map = []
    for node in G:
        match = False
        if node == "REFERE":
            color_map.append('Yellow')
            match = True
        elif node == "ENGLAN":
            color_map.append('Gold')
            match = True
        elif node == "SOUTHA":
            color_map.append('khaki')
            match = True
        elif node == "BRAZIL":
            color_map.append('darkkhaki')
            match = True

        for x in range(len(names)):
            if node in names[x]:
                if symps[x] == "Symptomatic Moderate to Severe Symptoms":
                        color_map.append("Red")
                elif symps[x] == "Asymptomatic or Mild Symptoms":
                        color_map.append("Blue")
                else:
                        color_map.append("green")
                match = True

        if match == False:
            color_map.append('black')

    plt.figure(figsize=(16,9))
    nx.draw(G, pos = posa, node_color=color_map, node_size=size_map, font_weight='bold')
    plt.savefig("clustering/symptoms/cluster" + str(ccounter + 1) + ".png")
    plt.clf()

############################################################################ outcomes

    color_map = []
    for node in G:
        match = False

        if node == "REFERE":
            color_map.append('Yellow')
            match = True
        elif node == "ENGLAN":
            color_map.append('Gold')
            match = True
        elif node == "SOUTHA":
            color_map.append('khaki')
            match = True
        elif node == "BRAZIL":
            color_map.append('darkkhaki')
            match = True

        for x in range(len(names)):
            if node in names[x]:
                if outs[x] == "Alive":
                        color_map.append("Blue")
                elif outs[x] == "Deceased":
                        color_map.append("purple")
                else:
                        color_map.append("green")
                match = True

        if match == False:
            color_map.append('black')

    plt.figure(figsize=(16,9))
    nx.draw(G, pos = posa, node_color=color_map, node_size=size_map, font_weight='bold')
    plt.savefig("clustering/outcomes/cluster" + str(ccounter + 1) + ".png")
    plt.clf()

############################################################################ States

    color_map = []
    for node in G:
        match = False

        if node == "REFERE":
            color_map.append('Yellow')
            match = True
        elif node == "ENGLAN":
            color_map.append('Gold')
            match = True
        elif node == "SOUTHA":
            color_map.append('khaki')
            match = True
        elif node == "BRAZIL":
            color_map.append('darkkhaki')
            match = True

        for x in range(len(names)):
            if node in names[x]:
                if "California" in states[x]:
                        color_map.append("Red")
                elif "Montana" in states[x]:
                        color_map.append("lawngreen")
                elif "New York" in states[x]:
                        color_map.append("Purple")
                elif "Puerto Rico" in states[x]:
                        color_map.append("Pink")
                elif "South Carolina" in states[x]:
                        color_map.append("Blue")
                elif "Louisiana" in states[x]:
                        color_map.append("Cyan")
                elif "Texas" in states[x]:
                        color_map.append("Brown")
                else:
                        color_map.append("darkGreen")
                match = True

        if match == False:
            color_map.append('black')

    plt.figure(figsize=(16,9))
    nx.draw(G, pos = posa, node_color=color_map, node_size=size_map, font_weight='bold')
    plt.savefig("clustering/states/cluster" + str(ccounter + 1) + ".png")
    plt.clf()

############################################################################ Dates

    color_map = []
    for node in G:
        match = False

        if node == "REFERE":
            color_map.append('Yellow')
            match = True
        elif node == "ENGLAN":
            color_map.append('Gold')
            match = True
        elif node == "SOUTHA":
            color_map.append('khaki')
            match = True
        elif node == "BRAZIL":
            color_map.append('darkkhaki')
            match = True

        for x in range(len(names)):
            if node in names[x]:
                if 2 == dates[x] or 3 == dates[x]:
                        color_map.append("purple")
                elif 4 == dates[x]:
                        color_map.append("lawngreen")
                elif 5 == dates[x] or 6 == dates[x]:
                        color_map.append("lightcoral")
                elif 7 == dates[x]:
                        color_map.append("Pink")
                elif 8 == dates[x] or 9 == dates[x]:
                        color_map.append("blue")
                elif 10 == dates[x] or 11 == dates[x]:
                        color_map.append("cyan")
                elif 12 == dates[x]:
                        color_map.append("fuchsia")
                elif 13 == dates[x] or 14 == dates[x]:
                        color_map.append("darkred")
                else:
                        color_map.append("darkgreen")
                match = True

        if match == False:
            color_map.append('black')  
    
    plt.figure(figsize=(16,9))
    nx.draw(G, pos = posa, node_color=color_map, node_size=size_map, font_weight='bold')
    plt.savefig("clustering/dates/cluster" + str(ccounter + 1) + ".png")
    plt.clf()