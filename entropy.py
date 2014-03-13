import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats.stats import pearsonr
import patsy
import numpy as np
import pandas as pd
from types import MethodType
import time
import numdifftools as ndt
from scipy.optimize import newton

###Regression code###

def regression(data, trips, cost, sep, factors, constraints):
    regstr = str(trips.name) + '~'
    first = True

    if factors != None:
        for ODs in factors.values():
            for factor in ODs:
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

def setup(data, trips, sep, cost, factors, constraints, destCon, attCon, initialParams):
    #For doubly constrained model
    if destCon == True & attCon == True:
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
        #There is always a beta parameter so set it to user's initial value and add to param list
        data['beta'] = initialParams['beta']
        params = ['beta']
        #Represents all known information in the system to include in the model
        knowns = data[sep]*data['Oi']*data[trips]
        #If there are additional factors
        if factors != None:
            for factor in factors:
                #Include that informatio in the model
                knowns = knowns*data[factor + 'param']
                #Add to params list
                params.append(factor + 'param')
        #Observed information is sum of trips multiplied by the log of known information
        observed = np.sum(data[trips]*np.log(knowns))
    #return observed info, data, knownn info, and params list

    if destCon == True & attCon == False:
        pass

    if destCon == False & attCon == True:
        pass

    if destCon == False& attCon == False:
        pass

    return observed, data, knowns, params


def calcAi(data, sep, factors):
    Ai = data["Bj"]*data["Dj"]
    if factors != None:
        for factor in factors['origins']:
            Ai = Ai*factor**data[factor.name + 'param']
    data["Ai"] = Ai*np.exp(data[sep]*data["beta"])

#Step 3: Function to Calculate Bj values
def calcBj(data, sep, factors):
    Bj = data["Ai"]*data["Oi"]
    if factors != None:
        for factor in factors['destinations']:
            Bj = Bj*factor**data[factor.name + 'param']
    data["Bj"] = Bj*np.exp(data[sep]*data["beta"])

#Step 4: Function to check if Ai and Bj have stabilised, if not return to step 2
def balanceFactors(data, sep, factors, constraints):
    its = 0
    cnvg = 1
    while cnvg > .0001:
        its = its + 1
        calcAi(data, sep, factors)
        AiBF = (data.groupby(data[constraints['production']].name).aggregate({"Ai": np.sum}))
        AiBF["Ai"] = 1/AiBF["Ai"]
        updates = AiBF.ix[pd.match(data[constraints['production']], AiBF.index), "Ai"]
        data["Ai"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Ai"]
        if its == 1:
            data["OldAi"] = data["Ai"]
        else:
            data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])
            data["OldAi"] = data["Ai"]

        calcBj(data, sep, factors)
        BjBF = data.groupby(data[constraints['attraction']].name).aggregate({"Bj": np.sum})
        BjBF["Bj"] = 1/BjBF["Bj"]
        updates = BjBF.ix[pd.match(data[constraints['attraction']], BjBF.index), "Bj"]
        data["Bj"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Bj"]
        if its == 1:
            data["OldBj"] = data["Bj"]
        else:
            data["diff"] = abs((data["OldBj"] - data["Bj"])/data["OldBj"])
            data["OldBj"] = data["Bj"]
        cnvg = np.sum(data["diff"])
    return data

#Step 5: Function to Calculate Tij' (flow estimates)
def estimateFlows(data, sep):
    data["SIM_Estimates"] = data["Oi"]*data["Ai"]*data["Dj"]*data["Bj"]*np.exp(data["Dij"]*data[sep])
    return data

#Step 6: Function to Calculate Sum of all products of Tij' and log distances
def estimateCum(data, knowns):
    return np.sum(data["SIM_Estimates"]*np.log(knowns))


def buildLLFunctions(data, params, factors, trips, sep, model):
    if model == 'dconstrained':
        common = data['Ai']*data['Oi']*data['Bj']*data['Dj']
    if model == 'destConstrained':
        common = data['Ai']*data['Oi']
    if model == 'attConstrained':
        common = data['Bj']*data['Dj']
    if model == 'unConstrained':
        pass

    paramCount = 1
    if params > 1:
        if factors != None:
            for factor in factors['origins'] | factors['destinations']:
                common = common*(data[factor]**x[paramCount])

    def buildFunction(common, data, trips, sep, param):
        return lambda x: np.sum(common*np.exp(data[sep]*x[0])*np.log(data[[param]])) - np.sum(data[trips]*data[param])

    functions = []
    for param in params:
        functions.append(buildFunction(common, data, trips, sep, param))

    return functions



def new(paramSingle, data, params, sep, trips, functions):
    if len(functions) > 1:
        return Jac(paramSingle, data, params, sep, trips, functions)
    else:
        return [newton(functions[0], paramSingle[0])]



#Function to compute the Jacobian matrix of parameter values
def Jac(PV, data, params, sep, trips):
    fBk = []
    common = computefBk(data, sep)
    for param in PV:
        fBk.append(np.sum(common(PV)*np.log(param))-np.sum(data[trips]*np.log(param)))
    #N by 1 array of functions representing the constraints which need to be met to obtain parameters
    fBk = np.array(fBk)

    #Three N by 1 arrays of Jacobian terms which are stacked into an N by N Matrix (jac)
    f1Jac = ndt.Jacobian(lambda x: np.sum((data.Oi**x[1])*data.Ai*data.Bj*data.Dj*np.exp(data.Dij*x[0])*np.log(data.Dij))- (np.sum(data.Data*np.log(data.Dij))))
    f1 = f1Jac(PV).flatten()
    #f1Jac = ndt.Jacobian(lambda x: np.sum(data.Oi*data.Ai*data.Bj*data.Dj*np.exp(data.Dij*x[0])*np.log(data.Dij))- (np.sum(data.Data*np.log(data.Dij))))

    jac = np.array([f1])
    #N by 1 array of current parameter values
    Bk = np.array(PV)

    print Bk, jac, fBk
    #Get new parameter estimates by subtracting from the original estimates the inverse of the product of the jacobian matrix and array of constraints
    return Bk - np.dot(np.linalg.inv(jac),fBk)


#Commputes the value of one
def computefBk(data, sep):
    paramCount = 1
    common = lambda PV: data.Oi*data.Ai*data.Bj*data.Dj*np.exp(data[sep]*PV[0])
    return common





def dConstrain(observed, data, knowns, params, trips, sep, factors, constraints):
    model = 'dconstrained'
    its = 0
    data = balanceFactors(data, sep, factors, constraints)
    data = estimateFlows(data, sep)
    estimates = estimateCum(data, knowns)
    functions = buildLLFunctions(data, params, factors, trips, sep, model)
    while abs(estimates - observed) > .0001:
        paramSingle = []
        for param in params:
            paramSingle.append(data[param].ix[0])
        jac = new(paramSingle, data, params, sep, trips, functions)
        print params
        for x, param in enumerate(params):
            print jac[x]
            data[param] = jac[x]
        data = balanceFactors(data, sep, factors, constraints)
        data = estimateFlows(data, sep)
        estimates = estimateCum(data, knowns)
        its += 1

    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])
    return data

def attConstrain(data, observed, knowns):
    return data

def destConstrain(data, observed, knowns):
    return data

def unConstrain(data, observed, knowns):
    return data
