import networkx as nx
import matplotlib.pyplot as plt

with open("Act11Stats/normalizedLinkageTable.txt", "r") as normalizedLinkageTable, open("Hist1Region.txt") as Hist1Region:
    nuclearProfiles = Hist1Region.readline().split("\t")[3:]
    hist1Region = []
    rows = []

    # create the chromosome array
    for row in Hist1Region:
        line = row.split("\t")
        tempChr13 = line[:3]
        data = line[3:-1]

        constructString = ""
        for s in tempChr13:
            constructString += s

        hist1Region.append(constructString)
        rows.append(data)

    # initialize empty arrays to import the normalized linkage table from the previous activity
    normalizedLinkageTableArr = []
    linkageArr = []

    for row in normalizedLinkageTable:
        rowSplit = row.split("\t")
        # add an array of binary values to pull the 2D array back from the import
        normalizedLinkageTableArr.append(rowSplit)
        for col in rowSplit:
            # flatten the array, adding each individual value to a single array
            linkageArr.append(col)

    # sort the values
    linkageArr.sort()

    # find the 75th percentile
    _75Percentile = float (linkageArr[int(0.75*len(linkageArr))])
    print(_75Percentile)

    #
    networkGraph = nx.Graph()
    graphEdges = []

    for i in range(len(normalizedLinkageTableArr)):
        for j in range(i + 1, len(normalizedLinkageTableArr)):
            if float (normalizedLinkageTableArr[i][j]) > _75Percentile:
                graphEdges.append((i, j))

    centrality = [0 for _ in range(len(normalizedLinkageTableArr))]

    for edge in graphEdges:
        centrality[edge[0]] += 1
        centrality[edge[1]] += 1

    degreeCentrality = [(hist1Region[i], i, (centralityDegree / len(hist1Region))) for i, centralityDegree in enumerate(centrality)]

    degreeCentrality.sort(key=lambda x: x[2])

    minDegree = min(degreeCentrality, key=lambda x: x[2])
    maxDegree = max(degreeCentrality, key=lambda x: x[2])
    avgDegree = sum(centralityDegree[2] for centralityDegree in degreeCentrality) / len(degreeCentrality)

    print(minDegree, maxDegree, avgDegree)

    print("\n")

    for i in range(len(degreeCentrality)):
        print(f"{i}: {degreeCentrality[i]}")

    networkGraph.add_edges_from(graphEdges)

    positioning = nx.spring_layout(networkGraph)

    plt.figure(figsize=(8, 8))
    nx.draw(
        networkGraph,
        pos=positioning,
        with_labels=True,
        node_color="green",
        edge_color="gray",
        node_size=200,
        font_size=10,
        font_color="white"
    )
    plt.show()