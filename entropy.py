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

    if cost == 'pow':
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
    origins = len(data['Origin'].unique())
    destinations = len(data['Destination'].unique())
    pairs = len(data)
    obsInt = np.sum(data[trips])
    predInt = np.sum(data['SIM_Estimates'])
    avgDist = round(np.sum(data[sep])*1.00000/pairs)
    avgDistTrav = round((np.sum(data[trips]*data[sep]))*1.00000/np.sum(data[trips])*1.00000)
    obsMeanTripLen = (np.sum(data[trips]*data[sep]))*1.00000/obsInt*1.000000
    predMeanTripLen = (np.sum(data['SIM_Estimates']*data[sep]))/predInt
    aSymSum = 0
    for o in data['Origin'].unique():
        for d in data['Destination'].unique():
            if o != d:
                aSymSum += abs((data[trips][(data['Origin'] == o) & (data['Destination'] == d)].values) - (data[trips][(data['Origin'] == d) & (data['Destination'] == o)].values))
    if origins == destinations:
        aSymInd = 50.0000*aSymSum[0]/(np.sum(data[trips]))
    else: aSymInd = 'N/A'
    #three likelihood statistics
    percentDev = round(((np.sum(abs(data[trips]-data['SIM_Estimates'])))/np.sum(data[trips]))*100, 3)
    intMean = round(np.sum(data[trips]/pairs), 1)
    percentDevMean = round(((np.sum(abs(data[trips]-intMean)))/np.sum(data[trips]))*100, 3)
    percentDevRed = (abs((percentDev-percentDevMean))/percentDevMean)*100
    pij = data[trips]/np.sum(data[trips])
    phatij = data['SIM_Estimates']/np.sum(data[trips])
    infoGain = np.sum(pij*np.log((pij/phatij)))
    sij = (pij+phatij)/2
    psiStat = np.sum(pij*np.log(pij/sij)) + np.sum(phatij*np.log(phatij/sij))
    #why is MDI only calculated once? skipped
    srmse = ((np.sum((data[trips]-data['SIM_Estimates'])**2)/pairs)**.5)/(np.sum(data[trips])/pairs)
    maxEntropy = round(np.log(pairs), 4)
    predEntropy = round(-np.sum(phatij*np.log(phatij)), 4)
    obsEntropy = round(-np.sum(pij*np.log(pij)), 4)
    diffPredEnt = round(maxEntropy - predEntropy, 4)
    diffObsEnt = round(maxEntropy - obsEntropy, 4)
    diffEntropy = round(predEntropy - obsEntropy, 4)
    entropyRS = round(diffPredEnt/diffObsEnt, 4)
    varPredEnt = round(((np.sum(phatij*(np.log(phatij)**2))-predEntropy**2)/obsInt) + ((pairs-1)/(2*obsInt**2)), 11)
    varObsEnt = round(((np.sum(pij*np.log(pij)**2)-obsEntropy**2)/obsInt) + ((pairs-1)/(2*obsInt**2)), 11)
    tStatEnt = round((predEntropy-obsEntropy)/((varPredEnt+varObsEnt)**.5), 4)
    return origins, destinations, pairs, obsInt, predInt, avgDist, avgDistTrav, obsMeanTripLen, predMeanTripLen, aSymInd, percentDev, percentDevMean, percentDevRed, pij, phatij, infoGain, psiStat, srmse, maxEntropy, predEntropy, obsEntropy, diffPredEnt, diffObsEnt, diffEntropy, entropyRS, varPredEnt, varObsEnt, tStatEnt

#Calculate parameter estimate statistics
def peStats(PV, data, params, factors, trips, sep, cost, model, constraints, knowns, estimates):
    if len(PV) == 1:
        firstD = buildLLFunctions(PV, data, params, factors, trips, sep, cost, model, constraints, knowns)
        recalc = buildLLFunctions(PV+.001, data, params, factors, trips, sep, cost, model, constraints, knowns)
        diff = firstD[0]-recalc[0]
        secondD = -(1/(diff/.001))
        data[params[0]] = PV[0]
        return [math.sqrt(secondD)]
    elif len(PV) > 1:
        counter = 0
        varMatrix = np.zeros((len(PV),len(PV)))

        for x, param in enumerate(PV):
            varParams = list(PV)
            varParams[x] = varParams[x] + .1

            varMatrix[x] = buildLLFunctions(varParams, data, params, factors, trips, sep, cost, model, constraints, knowns)







        return np.sqrt(np.linalg.inv(np.transpose(varMatrix)).diagonal())











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
def setup(data, trips, sep, cost, factors, constraints, prodCon, attCon, initialParams, totalFlows):

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
        if factors == None:
            Dj = data.groupby(data[totalFlows]).aggregate({trips: np.sum})
            data["Dj"] = Dj.ix[pd.match(data[totalFlows], Dj.index)].reset_index()[trips]

    #For Attraction Constrained model
    if (prodCon == False) & (attCon == True):

        #Calc total inflows
        Dj = data.groupby(data[constraints['attraction']]).aggregate({trips: np.sum})
        data["Dj"] = Dj.ix[pd.match(data[constraints['attraction']], Dj.index)].reset_index()[trips]
        if factors == None:
            Oi = data.groupby(data[totalFlows]).aggregate({trips: np.sum})
            data["Oi"] = Oi.ix[pd.match(data[totalFlows], Oi.index)].reset_index()[trips]

    #For Unconstrained Model
    if (prodCon == False) & (attCon == False):
        pass

    #The following setup is for within all models

    #There is always a beta parameter so set it to user's initial value and add to param list
    data['beta'] = initialParams['beta']
    params = ['beta']

    #This is the observed data for which we want to derive parameters
    if cost == 'exp':
        knowns = data[sep]
    elif cost == 'pow':
        knowns = np.log(data[sep])
    else:
        sys.exit(sys.exit("The distance/cost function must be either 'pow' or 'exp'."))

    #If there are additional factors we will include that observed data, add it to param list, and add a data vector for the param
    if factors != None:
        if attCon != False:
            for factor in factors['origins']:
                #Include that information in the model
                knowns = knowns+np.log(data[factor])
                #Add to params list
                params.append(str(factor))
                #variable param vector
                data[str(factor) + 'Param'] = initialParams[factor]
        if prodCon != False:
            for factor in factors['destinations']:
                #Include that informatio in the model
                knowns = knowns+np.log(data[factor])
                #Add to params list
                params.append(str(factor))
                #variable param vector
                data[str(factor) + 'Param'] = initialParams[factor]

    #Observed information is sum of trips multiplied by the log of known information
    observed = np.sum(data[trips]*knowns)


    #return observed info, data, knownn info, and params list
    return observed, data, knowns, params



#Function to calculate Ai values
def calcAi(data, sep, cost, factors, model):
    if cost == 'exp':
        Ai = np.exp(data[sep]*data["beta"])
    elif cost == 'pow':
        Ai = (data[sep]**data["beta"])
    else:
        sys.exit("The distance/cost function must be either 'pow' or 'exp'.")

    if factors != None:
        for factor in factors['destinations']:
            Ai = Ai*(data[factor]**data[factor + 'Param'])

    else:
        Ai = Ai*data['Dj']

    if model == 'dConstrained':
        Ai = Ai*data["Bj"]


    data["Ai"] = Ai



#Function to Calculate Bj values
def calcBj(data, sep, cost, factors, model):
    if cost == 'exp':
        Bj = np.exp(data[sep]*data["beta"])
    elif cost == 'pow':
        Bj = (data[sep]**data["beta"])
    else:
        sys.exit("The distance/cost function must be either 'pow' or 'exp'.")

    if factors != None:
        for factor in factors['origins']:
            Bj = Bj*(data[factor]**data[factor + 'Param'])

    else:
        Bj = Bj*data['Oi']

    if model == 'dConstrained':
        Bj = Bj*data["Ai"]



    data["Bj"] = Bj



#Function to check if Ai and Bj have stabilised, if not return to step 2
def balanceFactors(data, sep, cost, factors, constraints, model):
    its = 0
    cnvg = 1
    while cnvg > .0001:
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

    if cost == 'exp':
        decay = np.exp(data[sep]*data['beta'])
    elif cost == 'pow':
        decay = (data[sep]**data['beta'])
    else:
        sys.exit("The distance/cost function must be either 'pow' or 'exp'.")


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
        else:
            data["SIM_Estimates"] = data["SIM_Estimates"]*data['Dj']

    elif model == 'attConstrained':
        data["SIM_Estimates"] = data["Dj"]*data["Bj"]*decay
        if factors != None:
            for factor in factors['origins']:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
        else:
            data["SIM_Estimates"] = data["SIM_Estimates"]*data['Oi']


    elif model == 'unConstrained':
        data["SIM_Estimates"] = decay
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])

    return data

#Step 6: Function to Calculate Sum of all products of Tij' and log distances
def estimateCum(data, knowns):
    return np.sum(data["SIM_Estimates"]*knowns)


def buildLLFunctions(PV, data, params, factors, trips, sep, cost, model, constraints, knowns, peM=False):
    for x, param in enumerate(params):
        if param != 'beta':
            data[str(param) + 'Param'] = PV[x]
        else:
            data[param] = PV[x]

    if peM == False:
        data = balanceFactors(data, sep, cost, factors, constraints, model)



    def buildFunction(common, data, trips, param, factors, beta=False):

        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    last = len(factors)
                    for count, factor in enumerate(factors[key]):
                        if count+1 != last:
                            f += 'data["'+ str(factor) + '"]**PV[' + str(count+1) + ']*'
                        else:
                            f = 'data["'+ str(factor) + '"]**PV[' + str(count+1) + ']'

        if factors != None:


            if cost == 'exp':
                decay = np.exp(data[sep]*PV[0])
            elif cost == 'pow':
                decay = (data[sep]**PV[0])
            else:
                sys.exit("The distance/cost function must be either 'pow' or 'exp'.")


            if beta == True:
                if cost == 'exp':
                    return np.sum(common*eval(f)*decay*data[param]) - np.sum(data[trips]*data[param])
                else:
                    return np.sum(common*eval(f)*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))
            else:
                return np.sum(common*eval(f)*decay*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))

        else:

            if cost == 'exp':
                decay = np.exp(data[sep]*PV[0])
            elif cost == 'pow':
                decay = (data[sep]**PV[0])
            else:
                sys.exit("The distance/cost function must be either 'pow' or 'exp'.")


            if beta == True:
                if cost == 'exp':
                    return np.sum(common*decay*data[param]) - np.sum(data[trips]*data[param])
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
        if factors == None:
            common = common*data['Dj']
        func = buildFunction(common, data, trips, sep, factors, beta=True)
        functions.append(func)
        if factors != None:
            for key in factors.keys():
                for factor in factors[key]:
                    func = buildFunction(common, data, trips, factor, factors)
                    functions.append(func)

    if model == 'attConstrained':
        common = data['Bj']*data['Dj']
        if factors == None:
            common = common*data['Oi']
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
    #print 'Model selected: ' + model
    data = balanceFactors(data, sep, cost, factors, constraints, model)
    data = estimateFlows(data, sep, cost, model, factors)
    estimates = estimateCum(data, knowns)
    its = 0


    while abs(estimates - observed) > 1:
        paramSingle = []
        for param in params:
            if param != 'beta':
                paramSingle.append(data[str(param) + 'Param'].ix[0])
            else:
                paramSingle.append(data[param].ix[0])

        updates = fsolve(buildLLFunctions, paramSingle, (data, params, factors, trips, sep, cost, model, constraints, knowns))
        #print updates, abs(estimates - observed)
        for x, each in enumerate(params):
            updates[x] = round(updates[x], 7)

        data = balanceFactors(data, sep, cost, factors, constraints, model)
        data = estimateFlows(data, sep, cost, model, factors)
        estimates = estimateCum(data, knowns)


        its += 1
        if its > 100:
            break

    if 'Ai' in data.columns.names:
        Ai = data.Ai.values.copy()
    if 'Bj' in data.columns.names:
        Bj = data.Bj.values.copy()
    ests = data.SIM_Estimates.values.copy()
    for x, each in enumerate(params):
        updates[x] = round(updates[x], 7)

    variance = peStats(updates, data, params, factors, trips, sep, cost, model, constraints, knowns, estimates)
    origins, destinations, pairs, obsInt, predInt, avgDist, avgDistTrav, obsMeanTripLen, predMeanTripLen, aSymInd, percentDev, percentDevMean, percentDevRed, pij, phatij, infoGain, psiStat, srmse, maxEntropy, predEntropy, obsEntropy, diffPredEnt, diffObsEnt, diffEntropy, entropyRS, varPredEnt, varObsEnt, tStatEnt = sysDesc(data, trips, sep)

    if 'Ai' in data.columns.names:
        data.Ai = Ai
    if 'Bj' in data.columns.names:
        data.Bj = Bj
    data.SIM_Estimates = ests

    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])
    cor = pearsonr(data.SIM_Estimates, data.Data)[0]
    sumStr = ''
    sumStr += '\n'
    sumStr += '\nModel type: ' + str(model)
    sumStr += '\nWith ' + str(origins) + ' origins and ' + str(destinations) + ' destinations.'
    sumStr += '\n'
    sumStr += '\nThe observed mean trip length: ' +  str(obsMeanTripLen)
    sumStr += '\nThe predicted mean trip length: ' + str(predMeanTripLen)
    sumStr += '\n'
    for x,param in enumerate(params):
        sumStr += '\nMaximum Likelihood ' + param + ' parameter: ' + str(updates[x])
        sumStr += '\n'
    sumStr += '\nAfter ' + str(its) + 'iterations of the calibration routine with a cost/distance function of: ' +  str(cost)
    sumStr += '\n'
    sumStr += '\nThe number of origin-destination pairs considered = ' + str(pairs)
    sumStr += '\n'
    sumStr += '\nThe total interactions observed: ' + str(obsInt)
    sumStr += '\nThe total interactions predicted: ' + str(predInt)
    sumStr += '\n'
    sumStr += '\nThe Asymmetry Index for this interaction data: ' + str(aSymInd)
    sumStr += '\n'
    sumStr += '\nRegressing the observed interactions on the predicted interactions yields and r-squared value of: ' + str(cor*cor)
    sumStr += '\n'
    sumStr += '\nT statistic for regression: still needs to be computed'
    sumStr += '\n'
    sumStr += '\nPercentage deviation of observed interaction from the mean: ' + str(percentDevMean)
    sumStr += '\n'
    sumStr += '\nPercentage deviation of predicted interaction from the observed interaction: ' + str(percentDev)
    sumStr += '\n'
    sumStr += '\nPercentage reduction in deviation: ' + str(percentDevRed)
    sumStr += '\n'
    sumStr += '\nAyeni S Information Statistic (psi) = ' + str(psiStat)
    sumStr += '\n'
    sumStr += '\nMinimum Discriminant Information Statistic: still needs to be computed'
    sumStr += '\n'
    sumStr += '\nThe standardied root mean square error statistic: ' + str(srmse)
    sumStr += '\n'
    sumStr += '\nThe maximum entropy for ' + str(pairs) + ' cases: ' + str(maxEntropy)
    sumStr += '\nThe entropy of the predicted interactions: ' + str(predEntropy)
    sumStr += '\nThe entropy of the observed interactions: ' + str(obsEntropy)
    sumStr += '\n'
    sumStr += '\nMaximum entropy - entropy of predicted interactions: ' + str(diffPredEnt)
    sumStr += '\n'
    sumStr += '\nEntropy of predicted interactions - entropy of observed interactions: ' + str(diffEntropy)
    sumStr += '\n'
    sumStr += '\nEntropy ratio statistic: ' + str(entropyRS)
    sumStr += '\n'
    sumStr += '\nVariance of the entropy of predicted interactions: ' + str(varPredEnt)
    sumStr += '\n'
    sumStr += '\nVariance of the entropy of observed interactions: ' + str(varObsEnt)
    sumStr += '\n'
    sumStr += '\nT statistic for the absolute entropy difference: ' + str(tStatEnt)
    sumStr += '\n'
    sumStr += '\nInformation gain statistic: ' + str(infoGain)
    sumStr += '\n'
    sumStr += '\nAverage distance traveled in system: ' + str(avgDistTrav)
    sumStr += '\nAverage origin-destination separation: ' + str(avgDist)
    sumStr += '\n'
    for x, param, in enumerate(params):
        sumStr += '\nStandard error of the ' + str(param) + ' parameter: ' + str(variance[x])
        sumStr += '\n'
    return data, cor, sumStr



