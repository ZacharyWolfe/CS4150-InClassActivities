import random
import math
import matplotlib.pyplot as plt

def main():
    with open("copy2.txt", "r") as fileRead:

        # 1, 1
        w = {}
        # 0, 1
        x = {}
        # 1, 0
        y = {}

        pairs = {}
        header = fileRead.readline()
        headerSplit = header.split('\t')[3:]
        headerSplit.pop(len(headerSplit) - 1)
        nuclearProfileNames = headerSplit

        extractedNPColumns = [[] for _ in range (len(nuclearProfileNames))]

        for line in fileRead:
            lineSplit = line.split('\t')[3:]
            # remove the newline character at the end.
            lineSplit.pop(len(lineSplit) - 1)

            for i in range(len(lineSplit)):
                extractedNPColumns[i].append(lineSplit[i])
                for j in range(i, len(lineSplit)):
                    if int(lineSplit[i]) and int(lineSplit[j]):
                        w.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): w.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif not int(lineSplit[i]) and int(lineSplit[j]):
                        x.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): x.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif int(lineSplit[i]) and not int(lineSplit[j]):
                        y.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): y.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})

    with (
            open("Act7Stats/Pairs.txt", "w") as pairsFileWrite,
            open("Act7Stats/NonDiagonalPairs.txt", "w") as nonDiagonalPairsFileWrite,
            open("Act7Stats/matrixExport.txt", "w") as matrixFileWrite
        ):
        countDiagonal = 0

        # loop over the matrix to calculate the new equation provided in part1 of clustering
        for i in range(len(nuclearProfileNames)):
            for j in range(len(nuclearProfileNames)):
                key = nuclearProfileNames[i] + ", " + nuclearProfileNames[j]

                # for mirroring over the diagonal
                reverse_key = nuclearProfileNames[j] + ", " + nuclearProfileNames[i]

                # represents the total number of attributes where the attribute of A is 1
                LenA = w.get(key, 0) + y.get(key, 0)

                # represents the total number of attributes where the attribute of B is 1
                LenB = w.get(key, 0) + x.get(key, 0)

                # M11
                Numerator = w.get(key, 0) if not None else w.get(reverse_key, 0)

                # M11 / min {|A|, |B|}
                J = Numerator / min(LenA, LenB) if min(LenA, LenB) != 0 else 0

                # write the similarity to the pair
                pairs[key] = J
                pairsFileWrite.write(f"{key}\t{J}\n")

                # split the key on a comma so we can get both nuclear profiles that make up the key
                splitKey = key.split(', ')

                # if F10A3 == F10A3 essentially
                if splitKey[0] == splitKey[1]:
                    countDiagonal += 1
                else:
                    nonDiagonalPairsFileWrite.write(f"{key}\t{J}\n")

        for i in range(len(nuclearProfileNames)):
            matrixWrite = ""
            for j in range(len(nuclearProfileNames)):
                pair = pairs.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j])

                if pair is None and pairs.get(nuclearProfileNames[j] + ", " + nuclearProfileNames[i]) is None:
                    print(f"{nuclearProfileNames[i]} and {nuclearProfileNames[j]}")
                elif pair is None and pairs.get(nuclearProfileNames[j] + ", " + nuclearProfileNames[i]) is not None:
                    pair = pairs.get(nuclearProfileNames[j] + ", " + nuclearProfileNames[i])

                matrixWrite += f"{pair}\t"
            matrixFileWrite.write(f"{matrixWrite}\n")
            # if pairs.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) is None:
            #     print(pairs.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0))

    matrix = []
    with open("Act7Stats/matrixExport.txt", "r") as matrixFileRead:
        for line in matrixFileRead:
            matrix.append(line.split('\t')[:-1])

        for i in range(len(matrix)):
            for j in range(len(matrix)):
                matrix[i][j] = float(matrix[i][j])

    kSeededRandomNPs = []
    rand = random.randint(0, len(nuclearProfileNames) - 1)

    # grab three different randoms in the range of the nuclear profiles
    while rand not in kSeededRandomNPs and len(kSeededRandomNPs) < 3:
        kSeededRandomNPs.append(rand)
        rand = random.randint(0, len(nuclearProfileNames) - 1)

    print(kSeededRandomNPs)

    # start the simScoreClustering with the matrix we wrote, and the three starting kSeededRandoms
    KLists = simScoreClustering(matrix, nuclearProfileNames, kSeededRandomNPs)
    iterations = [[]]
    iterationNum = 0
    newCentroids = []
    iterationLists = []

    for i in range(100):
        while True:
            while len(newCentroids) < 3:
                for kList in KLists.values():
                    if len(kList) == 0:
                        continue
                    kListVals = [val[1] for val in kList]
                    avgKList = sum(kListVals) / len(kList)

                    diffList = {}
                    for i in range(len(kList)):
                        diff = abs(avgKList - kList[i][1])
                        diffList[nuclearProfileNames[i]] = diff

                    min_diff = min(diffList.values())
                    resNP = [name for name, diff in diffList.items() if diff == min_diff]

                    if len(resNP) > 1:
                        for name in resNP:
                            if name not in newCentroids and len(newCentroids) < 3:
                                newCentroids.append(nuclearProfileNames.index(name))
                    elif len(newCentroids) < 3:
                        newCentroids.append(nuclearProfileNames.index(resNP[0]))
            KLists = simScoreClustering(matrix, nuclearProfileNames, newCentroids)

            if newCentroids in iterations:
                centroidNPs = [nuclearProfileNames[k] for k in newCentroids]
                print(
                    f"\nConverged with {newCentroids}, which correspond to {centroidNPs}, found in iteration: {iterations.index(newCentroids)}")
                break

            iterationNum += 1
            print(f"\niteration: {iterationNum}")
            clusterVariations = []
            for p, klist in KLists.items():
                kListVals = [val[1] for val in klist]
                kListNPs = [val[0] for val in klist]
                iterationLists.append(kListNPs)
                print(f"{p}:\t{kListNPs}")
                sampleSum = 0

                avg = sum(kListVals) / len(klist) if len(klist) > 0 else 0
                for val in klist:
                    sampleSum += ((val[1] - avg) ** 2)

                clusterVariations.append(math.sqrt(sampleSum / (len(klist) - 1)) if len(klist) - 1 > 0 else 0)

            for variation in clusterVariations:
                print(f"cluster variation: {variation}")

            print(f"Total variation: {sum(clusterVariations)}")
            iterations.append(newCentroids)

            newCentroids = []
        iterations = [[]]
    bestCentroidGroupings = iterationLists[-3:]

    # convert to a list of usable indices
    for i in range(len(bestCentroidGroupings)):
        for j in range(len(bestCentroidGroupings[i])):
            bestCentroidGroupings[i][j] = nuclearProfileNames.index(bestCentroidGroupings[i][j])

    for kList in bestCentroidGroupings:
        print(kList)

    clusteringMatrices = []

    for index, kList in enumerate(bestCentroidGroupings):
        plotMatrix = []
        for npIndex in kList:
            npColumn = [int(extractedNPColumns[npIndex][k]) for k in range(len(extractedNPColumns[npIndex]))]
            plotMatrix.append(npColumn)
        clusteringMatrices.append(plotMatrix)

    # open the csv file provided
    with open("Hist1_region_features.csv", "r") as hist1Region:

        # grab the header line
        header = hist1Region.readline()

        # remove "name" column and grab all the others
        header = header.split(",")[1:]

        # remove the endline on the last header
        header[len(header) - 1] = header[len(header) - 1][:-1]

        # create an empty list for each header (to grab the columns)
        headerColumns = {header[i]: [] for i in range(len(header))}

        # grab the genomic windows
        genomicWindows = []

        for line in hist1Region:
            lineSplit = line.split(",")
            genomicWindows.append(lineSplit[0])
            for i in range(1, len(lineSplit)):

                # remove any endlines
                if lineSplit[i].endswith("\n"):
                    lineSplit[i] = lineSplit[i][:-1]

                #
                headerColumns[header[i - 1]].append(lineSplit[i])

    # 5 = stronglyEquatorial
    # 4 = somewhatEquatorial
    # 3 = neitherEquatorialNorApical
    # 2 = somewhatApical
    # 1 = stronglyApical

    radialPosits = {}

    listRadialPosits = ["stronglyApical", "somewhatApical", "neitherEquatorialNorApical", "somewhatEquatorial", "stronglyEquatorial"]

    for radialPosit in listRadialPosits:
        radialPosits[radialPosit] = []

    with open("Act2Stats/radialPositionExport.txt", "r") as radialPositions:
        for line in radialPositions:
            lineSplit = line.split("\t")

            if lineSplit[0] not in nuclearProfileNames:
                continue

            for i in range(5):
                val = int (lineSplit[1]) - 1
                if val == i:
                    radialPosits[listRadialPosits[i]].append(lineSplit[0])

    clusterByFiveRadialPosits = [[[]]]

    for index, cluster in enumerate(bestCentroidGroupings):
        for radialPositVal, radialPosit in radialPosits.items():
            # if radialPosit in bestCentroidGroupings[index]:
            clusterByFiveRadialPosits[index][listRadialPosits.index(radialPositVal)].append(radialPosit)


    for threeByFiveCluster in clusterByFiveRadialPosits:
        print(threeByFiveCluster)

def simScoreClustering(matrix, nuclearProfileNames, centroids):

    # create 3 maps all with empty lists (nuclear profile name mapped to a list of the clustered NPs)
    print(f"length of centroids: {len(centroids)}")
    for centroid in centroids:
        print(centroid)
    KLists = {nuclearProfileNames[centroids[i]]: [] for i in range(3)}
    for i in range(len(nuclearProfileNames)):
        simScore = []
        for k in centroids:
            # grab the score from the matrix using k from the centroids list passed
            score = matrix[i][k]

            if score == 0 and matrix[k][i] != 0:
                # mirrored value over diagonal
                simScore.append(matrix[k][i])
            else:
                simScore.append(score)

        nearestClusterIndex = simScore.index(max(simScore))
        setNPSimScore = [nuclearProfileNames[i], simScore[nearestClusterIndex]]
        KLists[nuclearProfileNames[centroids[nearestClusterIndex]]].append(setNPSimScore)

    return KLists

if __name__ == "__main__":
    main()