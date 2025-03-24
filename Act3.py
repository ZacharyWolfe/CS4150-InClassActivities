import os

def main():
    names = []
    nuclearProfiles = []
    region = []
    genomicRows = []
    GWStartStop = []

    largestGWCol = 0
    smallestGWCol = -1

    largestGWRow = 0
    smallestGWRow = -1

    smallestGWIdentifierRow = 0
    largestGWIdentifierRow = 0

    largestGWIdentifierCol = 0
    smallestGWIdentifierCol = 0

    genomicWindows = 0

    if os.path.exists("Hist1Region.txt"):
        os.remove("Hist1Region.txt")

    with open("GSE64881_segmentation_at_30000bp.passqc.multibam.txt", "r") as readFile, open("Hist1Region.txt", "w") as fileWrite:
        header = readFile.readline()
        names = header.split('\t')[3:]
        nuclearProfiles = [0] * len(names)

        for line in readFile:
            lineBefore = line.split('\t')[:3]
            if lineBefore[0] == "chr13" and (21700000 <= int(lineBefore[1]) <= 24100000 or 21700000 <= int(lineBefore[2]) <= 24100000):
                GWStartStop.append(lineBefore)
                region.append(line)

        for reg in region:
            genomicWindows += 1
            rowResult = 0
            splitLine = reg.split('\t')[3:]
            for k in range(len(splitLine)):
                toInt = int(splitLine[k])
                nuclearProfiles[k] += toInt
                rowResult += toInt

            if rowResult < smallestGWRow or smallestGWRow == -1:
                smallestGWRow = rowResult
                smallestGWIdentifierRow = genomicWindows

            # calculate a new MAX and set the identifier to index the name associated to this row
            if rowResult > largestGWRow:
                largestGWRow = rowResult
                largestGWIdentifierRow = genomicWindows

            genomicRows.append(rowResult)

        removeCols = []
        newRegions = []
        newProfiles = []
        newNames = []

        firstThree = header.split('\t')[:3]
        header = ""

        for chars in firstThree:
            header += chars + '\t'

        for i in range(len(nuclearProfiles)):
            if nuclearProfiles[i] == 0:
                removeCols.append(i + 3)
            else:
                newProfiles.append(nuclearProfiles[i])
                newNames.append(names[i])

        for name in newNames:
            header += name + '\t'

        for reg in region:
            newString = ""
            splitReg = reg.split('\t')
            for i in range(len(splitReg)):
                if i not in removeCols:
                    newString += splitReg[i] + "\t"

            newRegions.append(newString)

        # formatting/writing the new file
        fileWrite.write(header + "\n")
        for i, newReg in enumerate(newRegions):
            if i != len(newRegions) - 1:
                newReg += "\n"
            fileWrite.write(newReg)

    readFile.close()
    fileWrite.close()

    sumOnes = 0

    for i, val in enumerate(nuclearProfiles):
        sumOnes += val

        # if there's a column greater than the current largest
        if val > largestGWCol:
            # set the new largest column and its identifier (the index)
            largestGWCol = val
            largestGWIdentifierCol = i

        # if there's a column smaller than the current smallest
        if val < smallestGWCol or smallestGWCol == -1:
            # set the new smallest column and its identifier (the index)
            smallestGWCol = val
            smallestGWIdentifierCol = i

    radialPositionEstimate = ["strongly apical", "somewhat apical", "neither apical nor equatorial",
                                      "somewhat equatorial", "strongly equatorial"]

    print(f"Q1: Nuclear Profiles: {len(newNames)}\n")
    print(f"Q2: Genomic Windows: {genomicWindows}\n")
    print(f"Q3: Average GW per NP: {round(sumOnes / len(nuclearProfiles), 3)}\n")
    print(f"Q4A: Smallest #GW in any NP: {smallestGWCol, names[smallestGWIdentifierCol]}")
    print(f"Q4B: Largest #GW in any NP: {largestGWCol, names[largestGWIdentifierCol]}\n")
    print(f"Q5A: Average NPs in which a GW is detected: {round(sumOnes / genomicWindows, 3)}")
    print(f"Q5B: Smallest #NPs in which a GW is detected: {smallestGWRow, GWStartStop[smallestGWIdentifierRow]}")
    print(f"Q5C: Largest #NPs in which a GW is detected: {largestGWRow, GWStartStop[largestGWIdentifierRow]}")

    profilesAndClassifications = radialPositionImport(newNames)
    print(f"\nQ6: Most common radial position: {radialPositionEstimate[computeCommonRadialPos(profilesAndClassifications)]}")
    chromosomesAndOccupied = GWSummaryImport(GWStartStop)
    print(f"Q7: Typical compactions of the windows within the HIST1 region: {computeTypicalCompaction(chromosomesAndOccupied)}%")

    readFile.close()
    fileWrite.close()

def radialPositionImport(nuclearProfiles):
    intersection = {}
    with open("Act2Stats/radialPositionExport.txt", "r") as fileRead:
        for line in fileRead:
            tempLine = line.split('\t')[0]
            if tempLine in nuclearProfiles:
                intersection.update({tempLine: int(line.split('\t')[1])})

    fileRead.close()
    return intersection

def computeCommonRadialPos(profilesAndClassifications):
    commonality = [0] * 5
    for key, val in profilesAndClassifications.items():
        commonality[int (val - 1)] += 1

    return commonality.index(max(commonality))

def GWSummaryImport(GWStartStop):
    intersection = {}

    with open ("Act2Stats/GWSummaryUnformatted.txt", "r") as fileRead:
        for line in fileRead:
            lineTemp = line.split('\t')

            beforeList = lineTemp[:3]
            before = ""
            for k in range(len(beforeList)):
                before += beforeList[k] + " "

            if beforeList in GWStartStop:
                intersection.update({before: lineTemp[3]})

    fileRead.close()
    return intersection

def computeTypicalCompaction(chromosomesAndOccupied):
    commonality = {}
    for key, val in chromosomesAndOccupied.items():
        originalCommonality = commonality.get(val, 0)
        commonality.update({val: originalCommonality + 1})

    res = 0
    currBestKey = 0
    for key, val in commonality.items():
        if max(res, val) == val:
            res = val
            currBestKey = key

    return currBestKey

if __name__ == "__main__":
    main()