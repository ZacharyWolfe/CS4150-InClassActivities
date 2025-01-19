import os

# Creates a copy of the original file as to prevent accidental overwrites
def copyFile(readFile):
    with open(readFile, "r") as originalFile, open("copy.txt", "w") as fileWrite:
        for line in originalFile:
            fileWrite.write(line)

def main():
    file = "GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

    if os.path.exists("copy.txt"):
        os.remove("copy.txt")

    copyFile(file)

    NPGWCount = {}
    NPGWNames = []
    GWStartStop = []
    row1 = True
    sumOnes = 0
    NPWith1GW = 0

    largestGW = 0
    smallestGW = -1

    largestGWRow = 0
    smallestGWRow = -1

    smallestGWIdentifierRow = 0
    largestGWIdentifierRow = 0

    genomicWindows = 0
    genomicRowTotal = 0

    largestGWIdentifierCol = 0
    smallestGWIdentifierCol = 0

    with open("GSE64881_segmentation_at_30000bp.passqc.multibam.txt", "r") as copyFileRead:
        for line in copyFileRead:
            rowResult = 0

            if row1:
                NPGWNames = line.split('\t')[3:]

                # Remove the newline character for F9G5 everything but -1 characters
                NPGWNames[len(NPGWNames) - 1] = NPGWNames[len(NPGWNames) - 1][:-1]
                NPGWCount = {i: 0 for i in range(len(NPGWNames))}

            elif not row1:
                tempLine = line.split('\t')
                GWStartStop.append(tempLine[:3])

                gwLine = tempLine[3:]
                for k in range(len(gwLine)):
                    toInt = int (gwLine[k])

                    NPGWCount[k] += toInt
                    sumOnes += toInt
                    rowResult += toInt

                if rowResult < smallestGWRow or smallestGWRow == -1:
                    smallestGWRow = rowResult
                    smallestGWIdentifierRow = genomicWindows

                if rowResult > largestGWRow:
                    largestGWRow = rowResult
                    largestGWIdentifierRow = genomicWindows

                genomicRowTotal += rowResult

            if not row1:
                genomicWindows += 1
            else:
                row1 = False


        for key, val in NPGWCount.items():
            if val > 0:
                NPWith1GW += 1

            if val > largestGW:
                largestGW = val
                largestGWIdentifierCol = key

            if val < smallestGW or smallestGW == -1:
                smallestGW = val
                smallestGWIdentifierCol = key

        copyFileRead.close()

        print(f"\nQ1: Nuclear Profiles: {len(NPGWCount)}\n")
        print(f"Q2: Genomic Windows: {genomicWindows}\n")
        print(f"Q3: Average GW per NP: {round((sumOnes / len(NPGWCount)), 3)}\n")
        print(f"Q4A: Smallest #GW in any NP: {smallestGW, NPGWNames[smallestGWIdentifierCol]}")
        print(f"Q4B: Largest #GW in any NP: {largestGW, NPGWNames[largestGWIdentifierCol]}\n")
        print(f"Q5A: Average NPs in which a GW is detected: {round(genomicRowTotal / genomicWindows, 3)}")
        print(f"Q5B: Smallest #NPs in which a GW is detected: {smallestGWRow, GWStartStop[smallestGWIdentifierRow]}")
        print(f"Q5C: Largest #NPs in which a GW is detected: {largestGWRow, GWStartStop[largestGWIdentifierRow]}")

if __name__ == "__main__":
    main()