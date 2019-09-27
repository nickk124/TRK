result = [1.0, 2.0, 0.2, 0.3, -1, -2, -3, 1, 2, 3, -2, -4, -6, 2, 4, 6, 0.4, 0.2, 0.8, 4.5, 5, 6, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 9, 10, 11, 12, 13, 14, 15]

M = 2



hists = []
edges = []
bincounts = []

for i in range(M):
    bincounts.append(int(result[M + 2 + 6 * M + 4 + i]))

totalBins = sum(bincounts)
firstHistIndex = M + 2 + 6 * M + 4 + M

paramHistogramData = []

for i in range(M):
    hist = []
    edge = []
    bincount = bincounts[i]

    for k in range(bincount):
        hist.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + k])
        edge.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + k + bincount])
    
    edge.append(result[firstHistIndex + 2*sum(bincounts[0:i]) + i + bincount + bincount])

    hist = [str(i) for i in hist]
    edge = [str(i) for i in edge]

    hist = " ".join(hist)
    edge = " ".join(edge)

    paramHistogramData.append(hist + "\t" + edge)


paramHistogramData = "\n".join(paramHistogramData)


print paramHistogramData