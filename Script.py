#
# @file     Script.py
# @author   Zachary Wolfe (zw224021@ohio.edu)
# @brief    This program reads in a file, counts all NP,
#           genome windows, and other computations as described by Dr. Welch.
# @date     January 16, 2025
# @version  1.0
#

import os

# Creates a copy of the original file as to prevent accidental overwrites
def copyFile(readFile):
    # open the original file
    with open(readFile, "r") as originalFile, open("copy.txt", "w") as fileWrite:
        # write all the lines in the original to the copy
        for line in originalFile:
            fileWrite.write(line)

def main():
    file = "GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

    if os.path.exists("copy.txt"):
        os.remove("copy.txt")

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
    genomicRowTotal = 0

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

                genomicRowTotal += rowResult

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

        # cleanup
        if os.path.exists("copy.txt"):
            os.remove("copy.txt")

        print(f"\nQ1: Nuclear Profiles: {len(NPGWCount)}\n")
        print(f"Q2: Genomic Windows: {genomicWindows}\n")
        print(f"Q3: Average GW per NP: {round((sumOnes / len(NPGWCount)), 3)}\n")
        print(f"Q4A: Smallest #GW in any NP: {smallestGWCol, NPGWNames[smallestGWIdentifierCol]}")
        print(f"Q4B: Largest #GW in any NP: {largestGWCol, NPGWNames[largestGWIdentifierCol]}\n")
        print(f"Q5A: Average NPs in which a GW is detected: {round(genomicRowTotal / genomicWindows, 3)}")
        print(f"Q5B: Smallest #NPs in which a GW is detected: {smallestGWRow, GWStartStop[smallestGWIdentifierRow]}")
        print(f"Q5C: Largest #NPs in which a GW is detected: {largestGWRow, GWStartStop[largestGWIdentifierRow]}")

if __name__ == "__main__":
    main()