import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats.stats import pearsonr
import patsy
import numpy as np
import pandas as pd
import time
import numdifftools as ndt
from scipy.optimize import newton
import sys
import inspect
from scipy.optimize import fsolve
from scipy.optimize import fmin
import math


###Regression code###

def regression(data, trips, cost, sep, factors, constraints):

    regstr = str(data[trips].name) + '~'
    first = True

    if factors != None:
        for key in factors.keys():
            for factor in factors[key]:
                if first == True:
                    regstr = regstr + str(factor)
                    first = False
                else:
                    regstr = regstr + '+' + str(factor)

    if constraints != 'unconstrained':
       for constraint in constraints:
        if first != True:
            regstr = regstr + '+' + str(data[constraints[constraint]].name)
        else:
            regstr = regstr + str(data[constraints[constraint]].name)
            first = False

    if first != True:
        regstr = regstr + '+'

    if cost == 'invpow':
        regstr = regstr + 'np.log(' + str(data[sep].name) + ')'
    else:
        regstr = regstr + str(data[sep].name)

    results = smf.glm(regstr, data=data, family=sm.families.Poisson(link=sm.families.links.log)).fit()
    cor = pearsonr(results.fittedvalues, data[trips])[0]
    results.rsquared = cor**2
    results.regstr = regstr
    #Why does results.rsquared work but it isnt tabbale in Ipython?
    return results

###MLE Code###

#Calculate descriptive statistics
def sysDesc(data, trips, sep):
    numObserved = len(data)
    avgDist = np.sum(data[trips])/numObserved
    avgDistTrav = np.sum(data[trips]*data[sep])/np.sum(data[trips])
    #Skipped asymetry index

#Calculate parameter estimate statistics
def peStats(PV, data, params, factors, trips, sep, cost, model, constraints, knowns, estimates):
    firstD = buildLLFunctions([PV], data, params, factors, trips, sep, cost, model, constraints, knowns)
    recalc = buildLLFunctions([PV+.001], data, params, factors, trips, sep, cost, model, constraints, knowns)
    #print firstD, recalc
    diff = firstD[0]-recalc[0]
    #print diff
    secondD = -(1/(diff/.001))
    #print secondD
    print math.sqrt(secondD)
    data[params[0]] = PV



#Check to make sure factors all have been assigned initial parameters
def checkParams(factors, initParams):
    variables = []
    for key in factors.keys():
        if key not in ['origins', 'destinations']:
            sys.exit('Only acceptable keys for factors are "origns" and/or "destinations"')
        for factor in factors[key]:
            variables.append(factor)
    factors = set(variables)
    params = set(initParams.keys())
    params.discard('beta')
    if len(factors.symmetric_difference(params)) > 0:
        sys.exit('The initial paramter keys and the factor names must be symmetrical (excluding beta)')

#Set up balancing factors, total in/out flows, and parameters
def setup(data, trips, sep, cost, factors, constraints, prodCon, attCon, initialParams):

    #For doubly constrained model
    if (prodCon == True) & (attCon == True):

        #Variables for constants and deriving them
        data["Bj"] = 1.0
        data["Ai"] = 1.0
        data["OldAi"] = 10.000000000
        data["OldBj"] = 10.000000000
        data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])

        #Calc total outflows and inflows
        Oi = data.groupby(data[constraints['production']]).aggregate({trips: np.sum})
        data["Oi"] = Oi.ix[pd.match(data[constraints['production']], Oi.index)].reset_index()[trips]
        Dj = data.groupby(data[constraints['attraction']]).aggregate({trips: np.sum})
        data["Dj"] = Dj.ix[pd.match(data[constraints['attraction']], Dj.index)].reset_index()[trips]



    #For Production Constrained model
    if (prodCon == True) & (attCon == False):

        #Calc total outflows
        Oi = data.groupby(data[constraints['production']]).aggregate({trips: np.sum})
        data["Oi"] = Oi.ix[pd.match(data[constraints['production']], Oi.index)].reset_index()[trips]

    #For Attraction Constrained model
    if (prodCon == False) & (attCon == True):

        #Calc total inflows
        Dj = data.groupby(data[constraints['attraction']]).aggregate({trips: np.sum})
        data["Dj"] = Dj.ix[pd.match(data[constraints['attraction']], Dj.index)].reset_index()[trips]

    #For Unconstrained Model
    if (prodCon == False) & (attCon == False):
        pass

    #The following setup is for within all models

    #There is always a beta parameter so set it to user's initial value and add to param list
    data['beta'] = initialParams['beta']
    params = ['beta']

    #This is the observed data for which we want to derive parameters
    if cost == 'negexp':
        knowns = np.exp(data[sep])
    elif cost == 'invpow':
        knowns = data[sep]
    else:
        sys.exit(sys.exit("The distance/cost function must be either 'invpow' or 'negexp'."))

    #If there are additional factors we will include that observed data, add it to param list, and add a data vector for the param
    if factors != None:
        if attCon != False:
            for factor in factors['origins']:
                #Include that information in the model
                knowns = knowns*data[factor]
                #Add to params list
                params.append(str(factor))
                #variable param vector
                data[str(factor) + 'Param'] = initialParams[factor]
        if prodCon != False:
            for factor in factors['destinations']:
                #Include that informatio in the model
                knowns = knowns*data[factor]
                #Add to params list
                params.append(str(factor))
                #variable param vector
                data[str(factor) + 'Param'] = initialParams[factor]

    #Observed information is sum of trips multiplied by the log of known information
    observed = np.sum(data[trips]*np.log(knowns))

    #return observed info, data, knownn info, and params list
    return observed, data, knowns, params



#Function to calculate Ai values
def calcAi(data, sep, cost, factors, model):
    if cost == 'negexp':
        Ai = np.exp(data[sep]*data["beta"])
    elif cost == 'invpow':
        Ai = (data[sep]**data["beta"])
    else:
        sys.exit("The distance/cost function must be either 'invpow' or 'negexp'.")

    if model == 'dConstrained':
        Ai = Ai*data["Bj"]*data["Dj"]

    if factors != None:
        for factor in factors['destinations']:
            Ai = Ai*(data[factor]**data[factor + 'Param'])

    data["Ai"] = Ai



#Function to Calculate Bj values
def calcBj(data, sep, cost, factors, model):
    if cost == 'negexp':
        Bj = np.exp(data[sep]*data["beta"])
    elif cost == 'invpow':
        Bj = (data[sep]**data["beta"])
    else:
        sys.exit("The distance/cost function must be either 'invpow' or 'negexp'.")


    if model == 'dConstrained':
        Bj = Bj*data["Ai"]*data["Oi"]

    if factors != None:
        for factor in factors['origins']:
            Bj = Bj*(data[factor]**data[factor + 'Param'])

    data["Bj"] = Bj



#Function to check if Ai and Bj have stabilised, if not return to step 2
def balanceFactors(data, sep, cost, factors, constraints, model):
    its = 0
    cnvg = 1
    while cnvg > .01:
        its = its + 1
        if model != 'attConstrained':
            calcAi(data, sep, cost, factors, model)
            AiBF = (data.groupby(data[constraints['production']].name).aggregate({"Ai": np.sum}))
            AiBF["Ai"] = 1/AiBF["Ai"]
            updates = AiBF.ix[pd.match(data[constraints['production']], AiBF.index), "Ai"]
            data["Ai"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Ai"]
            if model == 'prodConstrained':
                break
            if its == 1:
                data["OldAi"] = data["Ai"]
            else:
                data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])
                data["OldAi"] = data["Ai"]

        if model != 'prodConstrained':
            calcBj(data, sep, cost, factors, model)
            BjBF = data.groupby(data[constraints['attraction']].name).aggregate({"Bj": np.sum})
            BjBF["Bj"] = 1/BjBF["Bj"]
            updates = BjBF.ix[pd.match(data[constraints['attraction']], BjBF.index), "Bj"]
            data["Bj"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Bj"]
            if its == 1:
                if model == 'attConstrained':
                    break
                data["OldBj"] = data["Bj"]
            else:
                data["diff"] = abs((data["OldBj"] - data["Bj"])/data["OldBj"])
                data["OldBj"] = data["Bj"]
        cnvg = np.sum(data["diff"])
        #print cnvg, its
    return data

#Step 5: Function to Calculate Tij' (flow estimates)
def estimateFlows(data, sep, cost, model, factors):

    if cost == 'negexp':
        decay = np.exp(data[sep]*data['beta'])
    elif cost == 'invpow':
        decay = (data[sep]**data['beta'])
    else:
        sys.exit("The distance/cost function must be either 'invpow' or 'negexp'.")


    if model == 'dConstrained':
        data["SIM_Estimates"] = data["Oi"]*data["Ai"]*data["Dj"]*data["Bj"]*decay

        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    elif model == 'prodConstrained':
        data["SIM_Estimates"] = data["Oi"]*data["Ai"]*decay
        if factors != None:
            for factor in factors['destinations']:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    elif model == 'attConstrained':
        data["SIM_Estimates"] = data["Dj"]*data["Bj"]*decay
        if factors != None:
            for factor in factors['origins']:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    elif model == 'unConstrained':
        data["SIM_Estimates"] = decay
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])

    return data

#Step 6: Function to Calculate Sum of all products of Tij' and log distances
def estimateCum(data, knowns):
    return np.sum(data["SIM_Estimates"]*np.log(knowns))


def buildLLFunctions(PV, data, params, factors, trips, sep, cost, model, constraints, knowns):
    for x, param in enumerate(params):
        if param != 'beta':
            data[str(param) + 'Param'] = PV[x]
        else:
            data[param] = PV[x]
    data = balanceFactors(data, sep, cost, factors, constraints, model)
    data = estimateFlows(data, sep, cost, model, factors)
    estimates = estimateCum(data, knowns)

    def buildFunction(common, data, trips, param, factors, beta=False):

        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    last = len(factors)
                    for count, factor in enumerate(factors[key]):
                        if count+1 != last:
                            f += 'data["'+ str(factor) + '"]**x[' + str(count+1) + ']*'
                        else:
                            f = data[str(factor) ]**PV[1]

        if factors != None:


            if cost == 'negexp':
                decay = np.exp(data[sep]*x[0])
            elif cost == 'invpow':
                decay = (data[sep]**PV[0])
            else:
                sys.exit("The distance/cost function must be either 'invpow' or 'negexp'.")


            if beta == True:
                if cost == 'negexp':
                    return np.sum(common*f*decay*np.log(np.exp(data[param]))) - np.sum(data[trips]*np.log(np.exp(data[param])))
                else:
                    return np.sum(common*f*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))
            else:
                return np.sum(common*f*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))

        else:

            if cost == 'negexp':
                decay = np.exp(data[sep]*x)
            elif cost == 'invpow':
                decay = (data[sep]**PV[0])
            else:
                sys.exit("The distance/cost function must be either 'invpow' or 'negexp'.")


            if beta == True:
                if cost == 'negexp':
                    return np.sum(common*decay*np.log(np.exp(data[param]))) - np.sum(data[trips]*np.log(np.exp(data[param])))
                else:
                    return np.sum(common*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))
            else:
                return np.sum(common*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))


    functions = []

    if model == 'dConstrained':
        common = data['Ai']*data['Oi']*data['Bj']*data['Dj']
        func = buildFunction(common, data, trips, sep, factors, beta=True)
        functions.append(func)
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    func = buildFunction(common, data, trips, factor, factors)
                    functions.append(func)


    if model == 'prodConstrained':
        common = data['Ai']*data['Oi']
        func = buildFunction(common, data, trips, sep, factors, beta=True)
        functions.append(func)
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    func = buildFunction(common, data, trips, factor, factors)
                    functions.append(func)

    if model == 'attConstrained':
        common = data['Bj']*data['Dj']
        func = buildFunction(common, data, trips, sep, factors, beta=True)
        functions.append(func)
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    func = buildFunction(common, data, trips, factor, factors)
                    functions.append(func)

    if model == 'unConstrained':
        common = 1
        func = buildFunction(common, data, trips, sep, factors, beta=True)
        functions.append(func)
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    func = buildFunction(common, data, trips, factor, factors)
                    functions.append(func)


    return functions


def run(observed, data, knowns, params, trips, sep, cost, factors, constraints, model):
    print 'Model selected: ' + model
    data = balanceFactors(data, sep, cost, factors, constraints, model)
    data = estimateFlows(data, sep, cost, model, factors)
    estimates = estimateCum(data, knowns)
    its = 0

    while abs(estimates - observed) > .001:
        paramSingle = []
        for param in params:
            if param != 'beta':
                paramSingle.append(data[str(param) + 'Param'].ix[0])
            else:
                paramSingle.append(data[param].ix[0])

        updates = fsolve(buildLLFunctions, paramSingle, (data, params, factors, trips, sep, cost, model, constraints, knowns))
        print updates, abs(estimates - observed)

        data = balanceFactors(data, sep, cost, factors, constraints, model)
        data = estimateFlows(data, sep, cost, model, factors)
        estimates = estimateCum(data, knowns)
        its += 1
        if its > 25:
            break
    variance = peStats(updates[0], data, params, factors, trips, sep, cost, model, constraints, knowns, estimates)
    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])
    cor = pearsonr(data.SIM_Estimates, data.Data)[0]
    return data, cor



