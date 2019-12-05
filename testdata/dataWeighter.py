import numpy as np
import pandas as pd
from scipy.integrate import quad

def p(c2, c2n, ec2n):
    return np.exp(-((c2 - c2n)**2.0)/(2.0*(ec2n**2.0)))/np.sqrt(2.0*np.pi*(ec2n**2.0))

def integrand(c2, data, c2n, ec2n):
    N = data.shape[0]
    sum = 0

    for i in range(N):
        c2i = data["c2"][i]
        ec2i = data["ec2"][i]
        sum += p(c2, c2i, ec2i)

    return p(c2, c2n, ec2n) * sum


def makeWeight(data, c2n, ec2n):

    intt = quad(integrand, -np.inf, +np.inf, args=(data, c2n, ec2n))

    intt = intt[0]

    return (1.0/intt)



def main(filename):
    data = pd.read_csv(filename)

    for n in range(data.shape[0]): #make weights
        print("Weighting datapoint %i..." % n)
        c2n = data["c2"][n]
        ec2n = data["ec2"][n]

        data["Weight"][n] = makeWeight(data, c2n, ec2n)

    print("Weights: \n")
    
    for j in range(data.shape[0]):
        print(data["Weight"][j])

    data.to_csv(filename)

    minWeight = np.amin(data["Weight"])

    for n in range(data.shape[0]): #normalize weights
        data["Weight"][n] = data["Weight"][n] / minWeight


    data.to_csv(filename)

main('newData.csv')