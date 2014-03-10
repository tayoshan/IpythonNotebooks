import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats.stats import pearsonr
import patsy
import numpy as np
import pandas as pd
from types import MethodType
import time
import numdifftools as ndt

###Regression code###

def regression(data, trips, cost, sep, factors, constraints):
    regstr = str(trips.name) + '~'
    first = True

    if factors != None:
        for factor in factors:
            if first == True:
                regstr = regstr + str(factor.name)
                first = False
            else:
                regstr = regstr + '+' + str(factor.name)

    if constraints != 'unconstrained':
       for constraint in constraints:
        if first != True:
            regstr = regstr + '+' + str(constraints[constraint].name)
        else:
            regstr = regstr + str(constraints[constraint].name)
            first = False

    if first != True:
        regstr = regstr + '+'

    if cost == 'invpow':
        regstr = regstr + 'np.log(' + str(sep.name) + ')'
    else:
        regstr = regstr + str(sep.name)

    results = smf.glm(regstr, data=data, family=sm.families.Poisson(link=sm.families.links.log)).fit()
    cor = pearsonr(results.fittedvalues, trips)[0]
    results.rsquared = cor**2
    results.regstr = regstr
    #Why does results.rsquared work but it isnt tabbale in Ipython?
    return results

###MLE Code###

def setup(data):
    observed = np.sum(data["Data"]*np.log(data["Dij"]*data["Oi"]*data["Dj"]))
    data["beta"] = .001
    data["mu"] = 1
    data["sigma"] = 1
    data["Bj"] = 1.0
    data["Ai"] = 1.0
    data["OldAi"] = 10.0000000000
    data["OldBj"] = 10.000000000
    data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])
    return observed, data

def calcAi(data):
    data["Ai"] = data["Bj"]*(data["Dj"]**data["sigma"])*np.exp(data["Dij"]*data["beta"])

#Step 3: Function to Calculate Bj values
def calcBj(data):
    data["Bj"] = data["Ai"]*(data["Oi"]**data["mu"])*np.exp(data["Dij"]*data["beta"])


#Step 4: Function to check if Ai and Bj have stabilised, if not return to step 2
def balanceFactors(data):
    its = 0
    cnvg = 1

    while cnvg > .0001:
        its = its + 1
        calcAi(data)
        AiBF = (data.groupby("Origin").aggregate({"Ai": np.sum}))
        AiBF["Ai"] = 1/AiBF["Ai"]
        updates = AiBF.ix[pd.match(data["Origin"], AiBF.index), "Ai"]
        data["Ai"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Ai"]
        if its == 1:
            data["OldAi"] = data["Ai"]
        else:
            data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])
            data["OldAi"] = data["Ai"]

        calcBj(data)
        BjBF = data.groupby("Destination").aggregate({"Bj": np.sum})
        BjBF["Bj"] = 1/BjBF["Bj"]
        updates = BjBF.ix[pd.match(data["Destination"], BjBF.index), "Bj"]
        data["Bj"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Bj"]
        if its == 1:
            data["OldBj"] = data["Bj"]
        else:
            data["diff"] = abs((data["OldBj"] - data["Bj"])/data["OldBj"])
            data["OldBj"] = data["Bj"]
        cnvg = np.sum(data["diff"])
        #print its, cnvg, data.beta.ix[0]
    return data

#Step 5: Function to Calculate Tij' (flow estimates)
def estimateFlows(data):
    data["SIM_Estimates"] = (data["Oi"]**data["mu"])*data["Ai"]*(data["Dj"]**data["sigma"])*data["Bj"]*np.exp(data["Dij"]*data["beta"])
    return data

#Step 6: Function to Calculate Sum of all products of Tij' and log distances
def estimateCum(data):
    return np.sum(data["SIM_Estimates"]*np.log(data["Dij"]*data["Oi"]*data["Dj"]))

#Function to compute the Jacobian matrix of parameter values
def Jac(initParams, data):
    #N by 1 array of functions representing the constraints which need to be met to obtain parameters
    functionBk = np.array([computef1(initParams, data), computef2(initParams, data), computef3(initParams, data)])

    #Three N by 1 arrays of Jacobian terms which are stacked into an N by N Matrix (jac)
    f1Jac = ndt.Jacobian(lambda x: np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Dij))- (np.sum(data.Data*np.log(data.Dij))))
    f2Jac = ndt.Jacobian(lambda x: np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Oi))- (np.sum(data.Data*np.log(data.Oi))))
    f3Jac = ndt.Jacobian(lambda x: np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Dj))- (np.sum(data.Data*np.log(data.Dj))))
    f1 = f1Jac(initParams).flatten()
    f2 = f2Jac(initParams).flatten()
    f3 = f3Jac(initParams).flatten()
    jac = np.array([f1,f2,f3])

    #N by 1 array of current parameter values
    Bk = np.array([initParams[0], initParams[1], initParams[2]])

    #Get new parameter estimates by subtracting from the original estimates the inverse of the product of the jacobian matrix and array of constraints
    return Bk - np.dot(np.linalg.inv(jac),functionBk)


#Commputes the value of one
def computef1(x, data):
    return (np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Dij))- (np.sum(data.Data*np.log(data.Dij))))

def computef2(x, data):
    return (np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Oi))- (np.sum(data.Data*np.log(data.Oi))))

def computef3(x, data):
    return (np.sum((data.Oi**x[1])*data.Ai*data.Bj*(data.Dj**x[2])*np.exp(data.Dij*x[0])*np.log(data.Dj))- (np.sum(data.Data*np.log(data.Dj))))

def dConstrain(data):
    its = 0
    observed, data = setup(data)
    data = balanceFactors(data)
    data = estimateFlows(data)
    estimates = estimateCum(data)
    print abs(estimates-observed)
    while abs(estimates - observed) > .0001:
        jac =  Jac([data["beta"].ix[0],data["mu"].ix[0],data["sigma"].ix[0]], data)
        data["beta"], data["mu"], data["sigma"] = jac[0], jac[1], jac[2]
        #data["beta"] = newton(function, data["beta"].ix[0], args=(data.Data, data.Dij, data.Oi, data.Ai, data.Dj, data.Bj))
        data = balanceFactors(data)
        data = estimateFlows(data)
        estimates = estimateCum(data)
        its += 1

    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])
    return data