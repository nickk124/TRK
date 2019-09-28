import TRKwebpageutils as TRK

print '----------------------------------'
mainData = "1	0.3	2	1\n2	0.3	3	1\n3	0.3	4	1\n4	0.3	5	1\n5	0.3	6	1\n6	0.3	7	1\n7	0.3	8	1\n8	0.3	9	1\n9	0.3	10	1"
guess = "1\n1\n0.3\n0.3"
weighted = False
functionType = "Linear"
priors = ""
hasPriors = False
findPivot = False
optimizeScale = False
findUncertainties = False
scale = 0.25
print 'request received'

scale = float(scale)

#data and guess request strings to python 

mainData = mainData.split("\n")
guess = guess.split("\n")
priorsData = priors.split("\n")

M = len(guess) - 2
N = len(mainData)

x = []
sx = []
y = []
sy = []
w = []

allparamsguess = []

for i in mainData:
    print i
    if len(i) == 0:
        continue
    i.replace("\t"," ") #replaces any potential tab delimiter with a space
    dataPoint = i.split() #list made from the row of inputData

    x.append(float(dataPoint[0]))
    sx.append(float(dataPoint[1]))
    y.append(float(dataPoint[2]))
    sy.append(float(dataPoint[3]))

    if weighted:
        w.append(float(dataPoint[4]))
    else:
        w.append(1.0)

for i in guess:
    if len(i) == 0:
        continue
    allparamsguess.append(float(i))

#priors request strings to python 

lb_values = []
ub_values = []
mu_values = []
dev_values = []

hasLb = []
hasUb = []
hasMu = []
hasDev = []

if hasPriors: #creates list of each of the four prior params
    for i in priorsData:
        if len(i) == 0:
            continue
        i.replace("\t"," ") #replaces any potential tab delimiter with a space
        dataPoint = i.split() #list made from the row of inputData

        iList = []

        lb_s = dataPoint[0] #individual strings
        ub_s = dataPoint[1]
        mu_s = dataPoint[2]
        dev_s = dataPoint[3]

        if lb_s == "x" or lb_s == "X" :
            hasLb.append(0)
            lb_values.append(0)
        else:
            hasLb.append(1)
            lb_values.append(float(lb_s))

        if ub_s == "x" or ub_s == "X":
            hasUb.append(0)
            ub_values.append(0)
        else:
            hasUb.append(1)
            ub_values.append(float(ub_s))

        if mu_s == "x" or mu_s == "X":
            hasMu.append(0)
            mu_values.append(0)
        else:
            hasMu.append(1)
            mu_values.append(float(mu_s))

        if dev_s == "x" or dev_s == "X":
            hasDev.append(0)
            dev_values.append(0)
        else:
            hasDev.append(1)
            dev_values.append(float(dev_s))

onlyBounded = False
onlyGaussian = False
both = False

whichPrior = 0

CountArr = [0, 0, 0, 0]
hasArray = [hasLb, hasUb, hasMu, hasDev]

ct = 0
for list in hasArray:
    for i in list:
        if i==1:
            CountArr[ct] += 1
    ct += 1

if (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] == 0 and CountArr[3] == 0):
    onlyBounded = True
    whichPrior = 2
elif (CountArr[0] == 0 and CountArr[1] == 0) and (CountArr[2] > 0 and CountArr[3] > 0):
    onlyGaussian = True
    whichPrior = 1
elif (CountArr[0] > 0 or CountArr[1] > 0) and (CountArr[2] > 0 and CountArr[3] > 0):
    both = True
    whichPrior = 3

#create int that defines the func type
funcTypeInt = 0

if functionType == 'Linear':
    funcTypeInt = 1
elif functionType == 'Quadratic':
    funcTypeInt = 2
elif functionType == 'Cubic':
    funcTypeInt = 3
elif functionType == 'PowerLaw':
    funcTypeInt = 4
elif functionType == 'Exponential':
    funcTypeInt = 5
elif functionType == 'Logarithmic':
    funcTypeInt = 6

xVec = TRK.DoubleVector(N)
sxVec = TRK.DoubleVector(N)
yVec = TRK.DoubleVector(N)
syVec = TRK.DoubleVector(N)
wVec = TRK.DoubleVector(N)

priorsParamsVec = TRK.DoubleVector(M * 4) #only gaussian or only bounded will just ust the first half of this array
hasPriorsVec = TRK.IntVector(M * 4)

allparamsguessVec = TRK.DoubleVector(M+2)

for i in range(N):
    xVec[i] = x[i]
    sxVec[i] = sx[i]
    yVec[i] = y[i]
    syVec[i] = sy[i]
    wVec[i] = w[i]

for j in range(M+2):
    allparamsguessVec[j] = allparamsguess[j]

if onlyBounded:
    for j in range(0,paramCount):
        priorsParamsVec[2*j] = lb_values[j]
        priorsParamsVec[2*j+1] = ub_values[j]

        hasPriorsVec[2*j] = hasLb[j]
        hasPriorsVec[2*j+1] = hasUb[j]

elif onlyGaussian:
    for j in range(0,paramCount):
        priorsParamsVec[2*j] = mu_values[j]
        priorsParamsVec[2*j+1] = dev_values[j]

        hasPriorsVec[2*j] = hasMu[j]
        hasPriorsVec[2*j+1] = hasDev[j]

elif both:
    for j in range(0,paramCount):
        priorsParamsVec[2*j] = mu_values[j]
        priorsParamsVec[2*j+1] = dev_values[j]

        hasPriorsVec[2*j] = hasMu[j]
        hasPriorsVec[2*j+1] = hasDev[j]

    for j in range(0,paramCount):
        priorsParamsVec[paramCount*2 + 2*j] = lb_values[j]
        priorsParamsVec[paramCount*2 + 2*j+1] = ub_values[j]

        hasPriorsVec[paramCount*2 + 2*j] = hasLb[j]
        hasPriorsVec[paramCount*2 + 2*j+1] = hasUb[j]

result = TRK.DoubleVector()        
print 'performing TRK'

priorsCheckInt = 0
pivotCheckInt = 0
opScaleInt = 0
doMCMCInt = 0

if hasPriors:
    priorsCheckInt = 1
if findPivot:
    pivotCheckInt = 1
if optimizeScale:
    opScaleInt = 1
if findUncertainties:
    doMCMCInt = 1

result = TRK.requestHandler(funcTypeInt, xVec, yVec, wVec, sxVec, syVec, allparamsguessVec, N, pivotCheckInt, priorsCheckInt, priorsParamsVec, hasPriorsVec, opScaleInt, doMCMCInt, scale)
print result

#result = {best fit params, slop, - 1 2 3, + 1 2 3 sigmas, s0, a, b, pivot, bincount1, bincount2 ... , hist1, edges1, hist2, edges2 ...

"""
needed outputs:
scalesData = "s0\na\nb"
pivotData = str(float(pivot))
paramHistogramData = "hist11 hist12 ... hist1n\tedges11 edges12 ... edges1n+1\nhist21 hist 22 ..." = "[hist1]\t[edges1]\n[hist2]..."
bestFit_123SigmasData = "-1 -2 -3\t1 2 3\n-1 -2 -3\t1 2 3\n..." where +/- 1, 2, 3 are +/- 1, 2, 3 standard deviations of each successive param
finalParametersData = "a0\na1\na2 ... \nslopx\nslopy"

"""

#best fit params

finalParametersData = []

for j in range(M + 2):
    finalParametersData.append(result[j])

finalParametersData = [str(i) for i in finalParametersData]
finalParametersData = "\n".join(finalParametersData)

#best fit uncertainties

bestFit_123SigmasData = []

for j in range(M):
    minusSigmas = []
    plusSigmas = []
    for i in range(3):
        minusSigmas.append(str(result[M + 2 + i + j*6]))
    for i in range(3,6):
        plusSigmas.append(str(result[M + 2 + i + j*6]))
    minusSigmas = " ".join(minusSigmas)
    plusSigmas = " ".join(plusSigmas)


    bestFit_123SigmasData.append(minusSigmas + "\t" + plusSigmas)
    
bestFit_123SigmasData = "\n".join(bestFit_123SigmasData)


#scale data
scales = [result[M + 2 + 6 * M], result[M + 2 + 6 * M + 1], result[M + 2 + 6 * M + 2]]
scalesData = "\n".join([str(s) for s in scales])

#pivot data

pivotData = str(result[M + 2 + 6 * M + 3])

#histogram data

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

print '=-=-=-=-=-='
jsonData = dict(finalParametersData = finalParametersData, bestFit_123SigmasData = bestFit_123SigmasData, paramHistogramData = paramHistogramData, pivotData = pivotData, scales = scalesData, )
#from gluon.serializers import json
print jsonData 