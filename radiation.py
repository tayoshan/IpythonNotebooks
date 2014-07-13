import pandas as pd
import numpy as np
import math

def calcS(data, origins, destinations, sep, pop, i, j):
    subset = data[data[origins] == i]
    compDest = subset[subset[sep] > subset[sep][((subset[origins] == i) & (subset[destinations] == j))].values[0]][destinations].values
    sVals.append(data[((data[origins].isin(compDest)) & (data[destinations] == i))][pop].sum())

def run(data, origins, destinations, sep, pop, trips):
    sVals = []
    data.apply(lambda x: s(data, origins, destinations, sep, pop, x.Origin, x.Destination), axis=1)
    data['s'] = sVals
