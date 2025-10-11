import matplotlib.pyplot as plt
import networkx as nx

with open("Act11Stats/normalizedLinkageTable.txt", "r") as normalizedLinkageTable, open("Hist1Region.txt") as Hist1Region, open("Hist1_region_features.csv", "r") as features:
    # read the first line, everything after "chrom, start, stop"
    nuclearProfiles = Hist1Region.readline().split("\t")[3:]
    hist1Region = []
    rows = []

    # create the chromosome array
    for row in Hist1Region:
        line = row.split("\t")
        # " chr13:#-# "
        tempChr13 = line[:3]

        # includes the binary values after " chr13:#-# " and deletes the newline
        data = line[3:-1]

        # flattens the array of " chr13:#-# " to a single string
        constructString = ""
        for s in tempChr13:
            constructString += s

        # appends the chromosome
        hist1Region.append(constructString)
        # appends the data for that row
        rows.append(data)

    # initialize empty arrays to import the normalized linkage table from the previous activity
    normalizedLinkageTableArr = []
    linkageArr = []

    for row in normalizedLinkageTable:
        # splits the data after tab delimiter, returns an array
        rowSplit = row.split("\t")

        # add an array of binary values to pull the 2D array back from the import
        normalizedLinkageTableArr.append(rowSplit)
        for col in rowSplit:
            # flatten the array, adding each individual value to a single linkage array
            linkageArr.append(col)

    # sort the values as per the activity
    linkageArr.sort()

    # find the 75th percentile
    _75Percentile = float (linkageArr[int(0.75*len(linkageArr))])
    # print(_75Percentile)

    # initialize the graph edges
    graphEdges = []

    # loop through the linkage table's rows and columns (unique pairs)
    for i in range(len(normalizedLinkageTableArr)):
        for j in range(i + 1, len(normalizedLinkageTableArr)):
            # if the normalized linkage is greater than the 75th percentile
            if float (normalizedLinkageTableArr[i][j]) > _75Percentile:
                # append the i and j unique pair
                graphEdges.append((i, j))

    # initialize a centrality array, containing all zeros over the length of the flattened array
    centrality = [0 for _ in range(len(normalizedLinkageTableArr))]


    for edge in graphEdges:
        # each edge contains two nodes, add to the count for that window
        centrality[edge[0]] += 1
        centrality[edge[1]] += 1

    # create a tuple array containing the window chr13:#-#, the index of that window, and the degree centrality (the count in the centrality array divided by the total number of genomic windows)
    degreeCentrality = [(hist1Region[i], i, (centralityDegree / len(hist1Region) -1)) for i, centralityDegree in enumerate(centrality)]

    # sort on the second parameter (the degree centrality) of the array
    degreeCentrality.sort(key=lambda x: x[2])

    # top five largest degrees
    fiveLargestDegree = degreeCentrality[-5:]

    # initialize an empty map, (window -> feature data)
    chromFeaturesMap = {}

    # to house the chr13:#-#, already done but more explicit
    chromosomes = []

    # extract the headers, (hist1, lad, ctcf, etc...)
    header = features.readline().split(',')[1:]

    for line in features:
        # grab the data on the delimieter of comma
        lineSplit = line.split(",")

        # the chromosome is the first string
        chromForm = lineSplit[0]
        # some formatting things to properly index
        chromForm = chromForm.replace(':', "")
        chromForm = chromForm.replace('-', "")

        # remove the newline, if it exists
        lineSplit[len(lineSplit)-1] = lineSplit[len(lineSplit)-1].replace("\n", "")

        # map the data
        chromFeaturesMap[chromForm] = lineSplit[1:]

        # append the chromosome (chr13:#-#)
        chromosomes.append(chromForm)

    # for chrom in chromFeaturesMap:
    #     print(chrom, chromFeaturesMap[chrom])

    # create a new hubs map, (node -> [node])
    hubs = {fiveLargestDegree[i]: [] for i in range(len(fiveLargestDegree))}

    # loop over the edges and hubs
    for edge in graphEdges:
        for hub in hubs:
            # if the index is not part of this edge, no need to do anything
            if hub[1] not in edge:
                continue

            # if it is, append the matching value
            hubs[hub].append(edge[0] if edge[1] == hub[1] else edge[1])

    print('\n')

    # index out the values of hist1 and LAD from the header array
    hist1Index = header.index('Hist1')
    LADIndex = header.index('LAD')

    for hub in hubs:
        # print each hub and its community of nodes
        print(hub, hubs[hub])

        # hubIndexTranslate = [chromosomes[i] for i in hubs[hub]]
        # print(hub, hubIndexTranslate)

        print(f"hub size: {len(hubs[hub])}")
        hist1Count = 0
        LADCount = 0

        for chromIndex in hubs[hub]:

            # print(hist1Index)
            if int(chromFeaturesMap[chromosomes[chromIndex]][hist1Index]) == 1:
                hist1Count += 1

            if int(chromFeaturesMap[chromosomes[chromIndex]][LADIndex]) == 1:
                LADCount += 1

        print(f"percent hist1region: {round(hist1Count / len(hubs[hub])*100, 2)}%")
        print(f"percent LAD: {round(LADCount / len(hubs[hub])*100, 2)}%\n")

        # nodes are the hub index plus all other community indices
        communityNodes = [hub[1]] + hubs[hub]

        # filter edges: include only those edges where both nodes are in the community nodes list
        communityEdges = [edge for edge in graphEdges if edge[0] in communityNodes and edge[1] in communityNodes]

        # create a new network graph
        community = nx.Graph()

        # add the edges from the array of edges
        community.add_edges_from(communityEdges)

        # create node sizes based on the degree centrality, just multiply by an arbitrary 1000 for size
        nodeSizes = {index: centrality * 1000 for _, index, centrality in degreeCentrality}

        # create a new figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))

        # create a new empty (zeros) matrix
        communityMatrix = [[0 for _ in range(81)] for _ in range(81)]

        # loop through all edges, and directly put the values of the edge into the matrix by activating the value at that intersection
        for edge in communityEdges:
            i, j = edge
            communityMatrix[i][j] = 1
            # symmetry
            communityMatrix[j][i] = 1

        fig.suptitle(f"Community Analysis for Window {hub[1]}")

        # plot network graph on ax1
        position = nx.spring_layout(community)
        nx.draw(
            community,
            pos=position,
            ax=ax1,
            font_size=8,
            with_labels=True,
            font_color="black",
            linewidths=1,
            edgecolors="black",
            edge_color="green",
            node_color="lightblue",
            node_size=[nodeSizes[GW] for GW in communityNodes],
        )
        ax1.set_title("Network Graph")

        # plot heatmap on ax2
        heatmap = ax2.imshow(communityMatrix, cmap='viridis')
        ax2.set_title("Community Heatmap")
        ax2.set_xlabel("Genomic Windows")
        ax2.set_ylabel("Genomic Windows")
        fig.colorbar(heatmap, ax=ax2, label="Edge Presence")

        # show both plots together
        plt.tight_layout()
        plt.show()

