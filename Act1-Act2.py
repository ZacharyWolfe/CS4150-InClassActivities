#
# @file     Act1-Act2
# @author   Zachary Wolfe (zw224021@ohio.edu)
# @brief    This program reads in a file, counts all NP,
#           genome windows, and other computations as described by Dr. Welch.
# @date     January 16, 2025
# @version  1.0
#

import numpy as np
import os

# Creates a copy of the original file as to prevent accidental overwrites
def copyFile(readFile):
    # open the original file
    with open(readFile, "r") as originalFile, open("copy.txt", "w") as fileWrite:
        # write all the lines in the original to the copy
        for line in originalFile:
            fileWrite.write(line)

def printRadialDict(nameDict, distFromCenter, filePath):
    newPath = "Act2Stats/" + filePath
    with open(newPath, "w") as fileWrite:
        for key in nameDict:
            fileWrite.write(f"{key} at a distance of {round(distFromCenter[key], 3)}%\tfrom the center")

def main():
    file = "GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

    if os.path.exists("copy.txt"):
        os.remove("copy.txt")

    print("\n")

    # copy the file to ensure
    copyFile(file)

    NPGWCount = {}
    NPGWNames = []
    GWStartStop = []
    row1 = True
    sumOnes = 0
    NPWith1GW = 0

    largestGWCol = 0
    smallestGWCol = -1

    largestGWRow = 0
    smallestGWRow = -1

    smallestGWIdentifierRow = 0
    largestGWIdentifierRow = 0

    largestGWIdentifierCol = 0
    smallestGWIdentifierCol = 0

    genomicWindows = 0
    genomicRows = []

    # open the copy we created
    with open("copy.txt", "r") as copyFileRead:
        # read each line
        for line in copyFileRead:
            rowResult = 0

            # if reading the first line, grab all the names, minus the first three
            if row1:
                NPGWNames = line.split('\t')[3:]

                # remove the newline character for F9G5 everything but -1 characters
                NPGWNames[len(NPGWNames) - 1] = NPGWNames[len(NPGWNames) - 1][:-1]

                # creates a new dictionary where the key is i, and the value is zero
                NPGWCount = {i: 0 for i in range(len(NPGWNames))}

            # if not reading the first line
            elif not row1:
                # split the line at each delimiter
                tempLine = line.split('\t')

                # the first 3 lines of the file (chrom, start, stop)
                GWStartStop.append(tempLine[:3])

                # the rest of the line after chrom, start, and stop
                gwLine = tempLine[3:]

                # iterate over the rest of the line
                for k in range(len(gwLine)):
                    # grab the integer at that position in the line
                    toInt = int (gwLine[k])

                    # add to necessary calculations
                    NPGWCount[k] += toInt
                    rowResult += toInt

                # calculate a new MIN and set the identifier to index the name associated to this row
                if rowResult < smallestGWRow or smallestGWRow == -1:
                    smallestGWRow = rowResult
                    smallestGWIdentifierRow = genomicWindows

                # calculate a new MAX and set the identifier to index the name associated to this row
                if rowResult > largestGWRow:
                    largestGWRow = rowResult
                    largestGWIdentifierRow = genomicWindows

                genomicRows.append(rowResult)

            if not row1:
                genomicWindows += 1
            else:
                row1 = False


        # loop over the map
        for key, val in NPGWCount.items():
            # sum all the values
            sumOnes += val

            # if there's a minimum of one, add it to the number of NP's with at least 1 GW
            if val > 0:
                NPWith1GW += 1

            # if there's a column greater than the current largest
            if val > largestGWCol:
                # set the new largest column and its identifier (the index)
                largestGWCol = val
                largestGWIdentifierCol = key

            # if there's a column smaller than the current smallest
            if val < smallestGWCol or smallestGWCol == -1:
                # set the new smallest column and its identifier (the index)
                smallestGWCol = val
                smallestGWIdentifierCol = key

        # close the file
        copyFileRead.close()

    detectionFrequencies = []

    for k, v in NPGWCount.items():
        detectionFrequencies.append((v / sumOnes) * 100)

    meanFreq = np.mean(detectionFrequencies)
    stdFreq = np.std(detectionFrequencies)

    threshold = (stdFreq * 2)

    possibleOutliers = []

    radialPositionEstimate = ["strongly apical", "somewhat apical", "neither apical nor equatorial", "somewhat equatorial", "strongly equatorial"]
    radialPositions = {}
    radialDistanceFromCenter = {}

    # task 1
    with open("Act2Stats/OutlierPrediction.txt", "w") as fileWrite:
        for i in range(len(detectionFrequencies)):
            if abs(meanFreq - detectionFrequencies[i]) / stdFreq > threshold:
                possibleOutliers.append(i)
                fileWrite.write(f"Possible Outlier: {NPGWNames[i]}\n\t\tdetection frequency: {round(detectionFrequencies[i] * 100, 3)}%\n\t\tdistance mean: {round(abs(meanFreq - detectionFrequencies[i]), 3)}\n\t\tz-score: {round(abs(meanFreq - detectionFrequencies[i]) / stdFreq, 3)}\n")

    # task 2
    minFreq = min(detectionFrequencies)
    maxFreq = max(detectionFrequencies)

    with open("Act2Stats/radialPositionExport.txt", "w") as fileWrite:
        for i, freq in enumerate(detectionFrequencies):
            normalizedFreq = (detectionFrequencies[i] - minFreq) / (maxFreq - minFreq)
            radialDistanceFromCenter[NPGWNames[i]] = round(normalizedFreq * 100, 3)
            radialPositions[NPGWNames[i]] = int (round(normalizedFreq * 4)) + 1

            fileWrite.write(f"{NPGWNames[i]}\t{round(normalizedFreq * 4) + 1}\n")
            # 0, 1, 2, 3, 4, 5 ( min of 0, max of 1 because of normalization )

    fileWrite.close()

    # radialPositionsOuterEdge = [ radialPositions[key] for key in radialPositions if radialPositions[key] >= 3]
    stronglyEquatorial = {}
    somewhatEquatorial = {}
    neitherEquatorialNorApical = {}
    somewhatApical = {}
    stronglyApical = {}

    # 5 = stronglyEquatorial
    # 4 = somewhatEquatorial
    # 3 = neitherEquatorialNorApical
    # 2 = somewhatApical
    # 1 = stronglyApical

    for key in radialPositions:
        if radialPositions[key] == 5:
            stronglyEquatorial.update({key: radialPositions[key]})

        if radialPositions[key] == 4:
            somewhatEquatorial.update({key: radialPositions[key]})

        if radialPositions[key] == 3:
            neitherEquatorialNorApical.update({key: radialPositions[key]})

        if radialPositions[key] == 2:
            somewhatApical.update({key: radialPositions[key]})

        if radialPositions[key] == 1:
            stronglyApical.update({key: radialPositions[key]})

    printRadialDict(stronglyEquatorial, radialDistanceFromCenter, "StronglyEquatorial.txt")
    printRadialDict(somewhatEquatorial, radialDistanceFromCenter, "SomewhatEquatorial.txt")
    printRadialDict(neitherEquatorialNorApical, radialDistanceFromCenter, "NeitherEquatorialNorApical.txt")
    printRadialDict(somewhatApical, radialDistanceFromCenter, "SomewhatApical.txt")
    printRadialDict(stronglyApical, radialDistanceFromCenter, "StronglyApical.txt")

    with open("Act2Stats/OutlierStronglySomewhatApical.txt", "w") as fileWrite:
        fileWrite.write(">>> Intersection of possible outliers and strongly/somewhat apical: \n")
        fileWrite.write(f"Num outliers: {len(possibleOutliers)}\n")
        possibleOutliers = [NPGWNames[i] for i in possibleOutliers]
        intersection = set(possibleOutliers).intersection(set.union(set(stronglyApical.keys()), set(somewhatApical.keys())))
        fileWrite.write(f"Num intersection: {len(intersection)}\n")
        fileWrite.write(f"Percent outliers in intersection: {round(len(intersection) / len(possibleOutliers) * 100, 3)}%\n")

        fileWrite.write("\nIntersections:\n")
        for i in intersection:
            fileWrite.write(f"{i}\n")

    # task 3
    sortedGenomicRows = sorted(genomicRows)

    with open ("Act2Stats/GWSummary.txt", "w") as fileWrite, open("Act2Stats/GWSummaryUnformatted.txt", "w") as fileWriteUnformatted:
        for i in range(len(genomicRows)):
            occupied = sortedGenomicRows[i] / len(NPGWCount)
            classification = round(occupied * 10)
            if not classification:
                classification = 1
            fileWrite.write(f"{GWStartStop[i][0]}\tstart:\t{GWStartStop[i][1]}\tstop:\t{GWStartStop[i][2]}\n\t\toccupied:\t{round(occupied * 100, 3)}%\tclass: {classification}\n\n")
            fileWriteUnformatted.write(f"{GWStartStop[i][0]}\t{GWStartStop[i][1]}\t{GWStartStop[i][2]}\t{round(occupied * 100, 3)}\t{classification}\n")
    fileWrite.close()
    fileWriteUnformatted.close()

    # cleanups
    if os.path.exists("copy.txt"):
        os.remove("copy.txt")

    print(f"Q1: Nuclear Profiles: {len(NPGWCount)}\n")
    print(f"Q2: Genomic Windows: {genomicWindows}\n")
    print(f"Q3: Average GW per NP: {round(sumOnes / len(NPGWCount), 3)}\n")
    print(f"Q4A: Smallest #GW in any NP: {smallestGWCol, NPGWNames[smallestGWIdentifierCol]}")
    print(f"Q4B: Largest #GW in any NP: {largestGWCol, NPGWNames[largestGWIdentifierCol]}\n")
    print(f"Q5A: Average NPs in which a GW is detected: {round(sumOnes / genomicWindows, 3)}")
    print(f"Q5B: Smallest #NPs in which a GW is detected: {smallestGWRow, GWStartStop[smallestGWIdentifierRow]}")
    print(f"Q5C: Largest #NPs in which a GW is detected: {largestGWRow, GWStartStop[largestGWIdentifierRow]}")

if __name__ == "__main__":
    main()