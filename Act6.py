import matplotlib.pyplot as plt
import random

def main():
    with open("copy2.txt", "r") as fileRead:
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
                    if int(lineSplit[i]) and int(lineSplit[j]):
                        w.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): w.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif not int(lineSplit[i]) and int(lineSplit[j]):
                        x.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): x.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif int(lineSplit[i]) and not int(lineSplit[j]):
                        y.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): y.get(
                            nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})

    fileRead.close()

    with open("Act5Stats/Pairs.txt", "w") as pairsFileWrite, open("Act5Stats/NonDiagonalPairs.txt", "w") as nonDiagonalPairsFileWrite, open("Act5Stats/matrixExport.txt", "w") as matrixFileWrite:
        countDiagonal = 0

        for i in range(len(nuclearProfileNames)):
            for j in range(len(nuclearProfileNames)):
                key = nuclearProfileNames[i] + ", " + nuclearProfileNames[j]
                reverse_key = nuclearProfileNames[j] + ", " + nuclearProfileNames[i]

                LenA = w.get(key, 0) + y.get(key, 0)
                LenB = w.get(reverse_key, 0) + x.get(reverse_key, 0)
                Numerator = w.get(key, 0) if not None else w.get(reverse_key, 0)

                J = Numerator / min(LenA, LenB) if min(LenA, LenB) != 0 else 0

                wVal = w.get(key, 0) + w.get(reverse_key, 0)
                xVal = x.get(key, 0) + x.get(reverse_key, 0)
                yVal = y.get(key, 0) + y.get(reverse_key, 0)

                denom = wVal + xVal + yVal
                writeVal = wVal / denom if denom != 0 else 0

                pairs[key] = writeVal
                pairsFileWrite.write(f"{key}\t{writeVal}\n")
                splitKey = key.split(', ')

                if splitKey[0] == splitKey[1]:
                    countDiagonal += 1
                else:
                    nonDiagonalPairsFileWrite.write(f"{key}\t{writeVal}\n")

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
    nonDiagonalPairsFileWrite.close()
    matrixFileWrite.close()

    matrix = []
    with open("Act5Stats/matrixExport.txt", "r") as matrixFileRead:
        for line in matrixFileRead:
            matrix.append(line.split('\t')[:-1])

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = float(matrix[i][j])
    matrixFileRead.close()

    # random.seed(42)
    # i, j, i, j, i, j

    kSeededRandomIJPairs = []
    for k in range(3):
        i = random.randint(0, len(nuclearProfileNames) - 1)
        j = random.randint(0, len(nuclearProfileNames) - 1)
        while pairs[nuclearProfileNames[i] + ", " + nuclearProfileNames[j]] == 0.0:
            i = random.randint(0, len(nuclearProfileNames) - 1)
            j = random.randint(0, len(nuclearProfileNames) - 1)
        kSeededRandomIJPairs.append(i)
        kSeededRandomIJPairs.append(j)

    KLists = {(nuclearProfileNames[kSeededRandomIJPairs[i]], nuclearProfileNames[kSeededRandomIJPairs[i + 1]]): [] for i in range(0, len(kSeededRandomIJPairs), 2)}

    for pair in KLists:
        print(f"{pair}: val: {pairs[pair[0] + ', ' + pair[1]]}")
    for i in range(len(nuclearProfileNames)):
        for j in range(len(nuclearProfileNames)):
            simScore = [abs(float (matrix[kSeededRandomIJPairs[k]][kSeededRandomIJPairs[k + 1]]) - float (matrix[i][j])) for k in range (0, len(kSeededRandomIJPairs), 2)]

            # 0, 1, 2 vs (0, 1), (2, 3), (4, 5)
            offsetIndex = 2 * (simScore.index(min(simScore)))
            kSeededI = kSeededRandomIJPairs[offsetIndex]
            kSeededJ = kSeededRandomIJPairs[offsetIndex + 1]
            nearestCluster = nuclearProfileNames[kSeededI], nuclearProfileNames[kSeededJ]
            KLists[nearestCluster].append(nuclearProfileNames[i] + ", " + nuclearProfileNames[j])

    for key, val in KLists.items():
        print(len(val))

if __name__ == "__main__":
    main()