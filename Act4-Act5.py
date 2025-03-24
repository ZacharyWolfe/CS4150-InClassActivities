import matplotlib.pyplot as plt

def main():
    with open("Hist1Region.txt", "r") as fileRead:
        nuclearProfileNames = []

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

        for line in fileRead:
            lineSplit = line.split('\t')[3:]
            # remove the newline character at the end.
            lineSplit.pop(len(lineSplit) - 1)

            for i in range(len(lineSplit)):
                for j in range(i, len(lineSplit)):
                    # print(f"{nuclearProfileNames[i]} and {nuclearProfileNames[j]}")

                    if int(lineSplit[i]) and int(lineSplit[j]):
                        w.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): w.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif not int(lineSplit[i]) and int(lineSplit[j]):
                        x.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): x.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif int(lineSplit[i]) and not int(lineSplit[j]):
                        y.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): y.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
    fileRead.close()

    with open("Act4Stats/Pairs.txt", "w") as pairsFileWrite, open("Act4Stats/DiagonalPairs.txt", "w") as diagonalPairsFileWrite, open("Act4Stats/NonDiagonalPairs.txt", "w") as nonDiagonalPairsFileWrite, open("Act4Stats/CSVExport.txt", "w") as csvFileWrite, open("Act4Stats/matrixExport.txt", "w") as matrixFileWrite:
        countDiagonal = 0
        csvFileWrites = []

        for i in range(len(nuclearProfileNames)):
            for j in range(len(nuclearProfileNames)):
                key = nuclearProfileNames[i] + ", " + nuclearProfileNames[j]

                wVal = w.get(key, 0)
                xVal = x.get(key, 0)
                yVal = y.get(key, 0)

                denom = wVal + xVal + yVal
                newDenom = 0
                writeVal = 0

                if denom != 0:
                    if wVal == 0:
                        print(f"wVal is zero for {i}: {nuclearProfileNames[i]} and {j}: {nuclearProfileNames[j]}")
                    writeVal = wVal / denom
                else:
                    newKey = nuclearProfileNames[j] + ", " + nuclearProfileNames[i]
                    wVal = w.get(newKey, 0)
                    xVal = x.get(newKey, 0)
                    yVal = y.get(newKey, 0)

                    if wVal == 0:
                        print(f"wVal is zero for {i}: {nuclearProfileNames[i]} and {j}: {nuclearProfileNames[j]}")
                    newDenom = wVal + xVal + yVal

                if newDenom != 0:
                    writeVal = wVal / newDenom

                # print(f"w:{val},\tx:{xVal},\ty:{yVal},\tpair:{key}")

                pairs.update({key: writeVal})
                pairsFileWrite.write(f"{key}\t{writeVal}\n")

                splitKey = key.split(', ')

                if splitKey[0] == splitKey[1]:
                    countDiagonal += 1
                    diagonalPairsFileWrite.write(f"{key}\t{writeVal}\n")
                else:
                    nonDiagonalPairsFileWrite.write(f"{key}\t{writeVal}\n")
                    csvFileWrites.append(f"{writeVal},{splitKey[0]}-{splitKey[1]}\n")

                # print(f"key indices: {nuclearProfileNames.index(splitKey[0])}, {nuclearProfileNames.index(splitKey[1])}")

        sortedCSVWrites = sorted(csvFileWrites, key=lambda x: float(x.split(',')[0]), reverse=True)
        csvFileWrite.writelines(sortedCSVWrites)

        # for pair in pairs:
        #     print(f"{pair}:\t{pairs[pair]}")

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

    pairsFileWrite.close()
    diagonalPairsFileWrite.close()
    nonDiagonalPairsFileWrite.close()
    csvFileWrite.close()

    with open("Act4Stats/matrixExport.txt", "r") as matrixFileRead:
        matrix = []
        for line in matrixFileRead:
            matrix.append(line.split('\t')[:-1])
            print(line.split('\t')[:-1])

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = float(matrix[i][j])
    matrixFileRead.close()


    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap='viridis')
    plt.colorbar()
    plt.title("Jaccard Heatmap")
    plt.show()

    newMatrix = [[0] * len(matrix) for _ in range(len(matrix))]

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            newMatrix[i][j] = 1 - float (matrix[i][j])

    plt.figure(figsize=(10, 8))
    plt.imshow(newMatrix, cmap='viridis')
    plt.colorbar()
    plt.title("Jaccard Distance")
    plt.show()
    # for key, val in pairs.items():
        # print(f"{key}:\t{val}")

if __name__ == "__main__":
    main()