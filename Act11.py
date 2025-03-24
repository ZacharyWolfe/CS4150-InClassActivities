import matplotlib.pyplot as plt

with open ("Hist1Region.txt", "r") as fileRead:
    header = fileRead.readline()
    header = header.split("\t")
    header = header[3:-1]

    nuclearProfiles = header
    nuclearProfileData = {nuclearProfiles[i]: [] for i in range(len(nuclearProfiles))}

    hist1Region = []
    cosegList = []
    linkages = []
    rows = []

    for row in fileRead:
        line = row.split("\t")
        tempChr13 = line[:3]
        data = line[3:-1]

        constructString = ""
        for s in tempChr13:
            constructString += s

        hist1Region.append(constructString)

        rows.append(data)

    for i in range(len(rows)):
        linkageLine = []
        for j in range(len(rows)):
            cosegregation = 0
            detectionFrequencyRowA = 0
            detectionFrequencyRowB = 0
            for k in range(len(rows[i])):
                if int (rows[i][k]) == 1:
                    detectionFrequencyRowA += 1

                if int (rows[j][k]) == 1:
                    detectionFrequencyRowB += 1

                if int (rows[i][k]) == 1 and int (rows[j][k]) == 1:
                    cosegregation += 1

            detectionFrequencyRowAPercentage = float(detectionFrequencyRowA / len(rows[i]))
            detectionFrequencyRowBPercentage = float(detectionFrequencyRowB / len(rows[i]))

            cosegregationPercentage = float(cosegregation / len(rows[i]))
            linkage = cosegregationPercentage - (detectionFrequencyRowAPercentage * detectionFrequencyRowBPercentage)

            max_qlinkage = 1

            if linkage < 0:
                max_linkage = min(detectionFrequencyRowAPercentage * detectionFrequencyRowBPercentage, (1 - detectionFrequencyRowAPercentage) * (1 - detectionFrequencyRowBPercentage))
            elif linkage > 0:
                max_linkage = min(detectionFrequencyRowAPercentage * (1 - detectionFrequencyRowBPercentage), (1 - detectionFrequencyRowAPercentage) * detectionFrequencyRowBPercentage)
            else:
                max_linkage = linkage

            linkageLine.append((linkage / max_linkage) if max_linkage > 0 else 0)
        linkages.append(linkageLine)
        # print(linkageLine)

    for row in range(len(linkages)):
        for col in range(len(linkages)):
            print(linkages[row][col])

    plt.figure(figsize=(10, 8))
    plt.title("Normalized Linkage of Two Genomic Windows")
    plt.imshow(linkages, cmap='viridis')
    plt.colorbar()
    plt.xlabel("Genomic Windows")
    plt.ylabel("Genomic Windows")
    plt.show()