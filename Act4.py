def main():
    with open("copy2.txt", "r") as fileRead:
        header = fileRead.readline()
        headerSplit = header.split('\t')[3:]

        nuclearProfileNames = headerSplit

        # 1, 1
        w = {}
        # 0, 1
        x = {}
        # 1, 0
        y = {}

        pairs = {}

        for line in fileRead:
            lineSplit = line.split('\t')[3:]
            # remove the newline character at the end.
            lineSplit.pop(len(lineSplit) - 1)

            for i in range(len(lineSplit)):
                for j in range(i, len(lineSplit)):

                    if int(lineSplit[i]) and int(lineSplit[j]):
                        w.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): w.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif not int(lineSplit[i]) and int(lineSplit[j]):
                        x.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): x.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
                    elif int(lineSplit[i]) and not int(lineSplit[j]):
                        y.update({(nuclearProfileNames[i] + ", " + nuclearProfileNames[j]): y.get(nuclearProfileNames[i] + ", " + nuclearProfileNames[j], 0) + 1})
    fileRead.close()

    with open("Act4Stats/DiagonalPairs.txt", "w") as diagonalPairsFileWrite, open("Act4Stats/NonDiagonalPairs.txt", "w") as nonDiagonalPairsFileWrite, open("Act4Stats/CSVExport.txt", "w") as csvFileWrite:
        countDiagonal = 0
        csvFileWrites = []

        for i, (key, val) in enumerate(w.items()):
            xVal = x.get(key, 0)
            yVal = y.get(key, 0)

            print(f"w:{val},\tx:{xVal},\ty:{yVal},\tpair:{key}")

            pairs.update({key: (val / (val + xVal + yVal))})

            splitKey = key.split(', ')

            if splitKey[0] == splitKey[1]:
                countDiagonal += 1
                diagonalPairsFileWrite.write(f"{key}\t{val / (val + xVal + yVal)}\n")
            else:
                nonDiagonalPairsFileWrite.write(f"{key}\t{val / (val + xVal + yVal)}\n")
                csvFileWrites.append(f"{val / (val + xVal + yVal)},{splitKey[0]}-{splitKey[1]}\n")


        print(f"number of diagonals {countDiagonal}")

        sortedCSVWrites = sorted(csvFileWrites, key=lambda x: float(x.split(',')[0]), reverse=True)
        csvFileWrite.writelines(sortedCSVWrites)

    diagonalPairsFileWrite.close()
    nonDiagonalPairsFileWrite.close()
    csvFileWrite.close()
    # for key, val in pairs.items():
        # print(f"{key}:\t{val}")

if __name__ == "__main__":
    main()