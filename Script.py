def modifyCopyFile(readFile):
    with open(readFile, "r") as originalFile, open("copy.txt", "w") as copyFile:
        for line in originalFile:
            delimiters = 0

            for i in range(len(line)):
                if line[i] == '\t':
                    delimiters += 1

                if delimiters == 3 and i + 1 < len(line):
                    copyFile.write(line[i + 1:])
                    break

    copyFile.close()
    originalFile.close()


def main():
    file = "GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

    modifyCopyFile(file)

    NPGWCount = {}
    NPGWNames = []
    line = ""
    row1 = True
    i = 0
    j = 0
    sumOnes = 0
    NPWith1GW = 0
    largestGW = 0
    smallestGW = 0
    columnCount = 0
    largestGWRow = 0
    smallestGWRow = 0
    genomicWindows = 0

    with open("copy.txt", "r") as copyFileRead:
        for line in copyFileRead:
            while i < len(line):

                if line[i] == 'F':
                    NPGWCount[j] = 0
                    k = 0

                    while i + k < len(line) and line[i + k] != '\t':
                        k += 1

                    NPGWNames.append(line[i: i + k])
                    i += k
                    j += 1

                if i < len(line) and line[i] == '\t':
                    columnCount += 1

                elif not row1 and line[i] == '1':
                    sumOnes += 1
                    NPGWCount[columnCount] += 1

                i += 1

            columnCount = 0
            i = 0

            if not row1:
                genomicWindows += 1
            else:
                row1 = False

        for count in NPGWCount.values():
            if count > 0:
                NPWith1GW += 1

            largestGW = max(largestGW, count)
            smallestGW = min(smallestGW, count) if smallestGW != 0 else count

        copyFileRead.close()

        # Output the results
        print(f"Q1: Nuclear Profiles: {len(NPGWCount)}\n")
        print(f"Q2: Genomic Windows: {genomicWindows}\n")
        print(f"Q3: Average GW per NP: {sumOnes / len(NPGWCount)}\n")
        print(f"Q4A: Smallest #GW in any NP: {smallestGW}")
        print(f"Q4B: Largest #GW in any NP: {largestGW}\n")
        print(f"Q5A: Average NPs in which a GW is detected: {(NPWith1GW / len(NPGWCount)) * 100}\n")
        print(f"Q5B: Smallest #NPs in which a GW is detected: {smallestGW}")
        print(f"Q5C: Largest #NPs in which a GW is detected: {largestGW}\n")


if __name__ == "__main__":
    main()

# with open ("GSE64881_segmentation_at_30000bp.passqc.multibam.txt", "r") as file:
#     contents = file.read()



