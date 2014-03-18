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

    regstr = str(data[trips].name) + '~'
    first = True

    if len(factors) != 0:
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

def setup(data, trips, sep, cost, factors, constraints, prodCon, attCon, initialParams):

    #For doubly constrained model
    if prodCon == True & attCon == True:

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
    if prodCon == True & attCon == False:

        #Calc total outflows
        Oi = data.groupby(data[constraints['production']]).aggregate({trips: np.sum})
        data["Oi"] = Oi.ix[pd.match(data[constraints['production']], Oi.index)].reset_index()[trips]

    #For Attraction Constrained model
    if prodCon == False & attCon == True:

        #Calc total inflows
        Dj = data.groupby(data[constraints['attraction']]).aggregate({trips: np.sum})
        data["Dj"] = Dj.ix[pd.match(data[constraints['attraction']], Dj.index)].reset_index()[trips]

    #For Unconstrained Model
    if prodCon == False& attCon == False:
        pass

    #The following setup is for within all models

    #There is always a beta parameter so set it to user's initial value and add to param list
    data['beta'] = initialParams['beta']
    params = ['beta']

    #This is the observed data for which we want to derive parameters
    knowns = data[sep]

    #If there are additional factors we will include that observed data, add it to param list, and add a data vector for the param
    if len(factors) != 0:
        if attCon != False:
            for factor in factors['origins']:
                #Include that informatio in the model
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
def calcAi(data, sep, factors, model):
    if model == 'dConstrained':
        Ai = data["Bj"]*data["Dj"]
    else:
        Ai = data['Dij']
    if len(factors) != 0:
        for factor in factors['destinations']:
            Ai = Ai*(data[factor]**data[factor + 'Param'])
    data["Ai"] = Ai*np.exp(data[sep]*data["beta"])



#Function to Calculate Bj values
def calcBj(data, sep, factors, model):
    if model == 'dConstrained':
        Bj = data["Ai"]*data["Oi"]
    else:
        Bj = data['Dij']
    if len(factors) != 0:
        for factor in factors['origins']:
            Bj = Bj*(data[factor]**data[factor + 'Param'])
    data["Bj"] = Bj*np.exp(data[sep]*data["beta"])



#Function to check if Ai and Bj have stabilised, if not return to step 2
def balanceFactors(data, sep, factors, constraints, model):
    its = 0
    cnvg = 1
    while cnvg > .0001:
        its = its + 1
        if model != 'attConstrained':
            calcAi(data, sep, factors, model)
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
            calcBj(data, sep, factors, model)
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
    return data

#Step 5: Function to Calculate Tij' (flow estimates)
def estimateFlows(data, sep, model, factors):
    if model == 'dConstrained':
        data["SIM_Estimates"] = data["Oi"]*data["Ai"]*data["Dj"]*data["Bj"]*np.exp(data[sep]*data['beta'])
        for key in factors.keys():
            for factor in factors[key]:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    elif model == 'prodConstrained':
        data["SIM_Estimates"] = data["Oi"]*data["Ai"]*np.exp(data[sep]*data['beta'])
        if len(factors) != 0:
            for factor in factors['destinations']:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    elif model == 'attConstrained':
        data["SIM_Estimates"] = data["Dj"]*data["Bj"]*np.exp(data[sep]*data['beta'])
        if len(factors) != 0:
            for factor in factors['origins']:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])
    else:
        data["SIM_Estimates"] = np.exp(data[sep]*data['beta'])
        for key in factors.keys():
            for factor in factors[key]:
                data["SIM_Estimates"] = data["SIM_Estimates"]*(data[factor]**data[str(factor) + 'Param'])


    return data

#Step 6: Function to Calculate Sum of all products of Tij' and log distances
def estimateCum(data, knowns):
    return np.sum(data["SIM_Estimates"]*np.log(knowns))


def buildLLFunctions(data, params, factors, trips, sep, model, PV):


    def buildFunction(PV, common, data, trips, param, factors):

        for key in factors.keys():
            for factor in factors[key]:
                f = ''
                last = len(factors)
                for count, factor in enumerate(factors['destinations']):#
                    print count, factor
                    if count != last:
                        f += 'data["'+ str(factor) + '"]**x[' + str(count) + ']*'
                    else:
                        f += 'data["'+ str(factor) + '"]**x[' + str(count) + ']'


        def function(x):

            return np.sum(common*eval(f)*np.exp(data[sep]*x)*np.log(data[param])) - np.sum(data[trips]*np.log(data[param]))
        return function

    functions = []

    if model == 'dConstrained':
        common = data['Ai']*data['Oi']*data['Bj']*data['Dj']
        func = buildFunction(PV, common, data, trips, sep, factors)
        functions.append(func)
        for key in factors.keys():
            for factor in factors[key]:
                func = buildFunction(PV, common, data, trips, factor, factors)
                functions.append(func)


    if model == 'prodConstrained':
        common = data['Ai']*data['Oi']
        func = buildFunction(PV, common, data, trips, sep, factors)
        functions.append(func)
        for key in factors.keys():
            for factor in factors[key]:
                func = buildFunction(PV, common, data, trips, factor, factors)
                functions.append(func)

    if model == 'attConstrained':
        common = data['Bj']*data['Dj']
        func = buildFunction(PV, common, data, trips, sep, factors)
        functions.append(func)
        for key in factors.keys():
            for factor in factors[key]:
                func = buildFunction(PV, common, data, trips, factor, factors)
                functions.append(func)

    if model == 'unConstrained':
        func = buildFunction(PV, common, data, trips, sep, factors)
        functions.append(func)
        for key in factors.keys():
            for factor in factors[key]:
                func = buildFunction(PV, common, data, trips, factor, factors)
                functions.append(func)






    return functions



def new(paramSingle, data, params, sep, trips, functions):
    if len(functions) > 1:
        return Jac(paramSingle, data, params, sep, trips, functions)
    else:

        return [newton(functions[0], paramSingle[0])]



#Function to compute the Jacobian matrix of parameter values
def Jac(PV, data, params, sep, trips, functions):
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

    #Get new parameter estimates by subtracting from the original estimates the inverse of the product of the jacobian matrix and array of constraints
    return Bk - np.dot(np.linalg.inv(jac),fBk)


#Commputes the value of one
def computefBk(data, sep):
    paramCount = 1
    common = lambda PV: data.Oi*data.Ai*data.Bj*data.Dj*np.exp(data[sep]*PV[0])
    return common





def dConstrain(observed, data, knowns, params, trips, sep, factors, constraints):
    model = 'dConstrained'
    its = 0
    data = balanceFactors(data, sep, factors, constraints, model)
    data = estimateFlows(data, sep, model, factors)
    estimates = estimateCum(data, knowns)
    while abs(estimates - observed) > .0001:
        paramSingle = []
        for param in params:
            paramSingle.append(data[param].ix[0])
        functions = buildLLFunctions(data, params, factors, trips, sep, model, paramSingle)
        jac = new(paramSingle, data, params, sep, trips, functions)
        print jac, estimates, observed, its
        for x, param in enumerate(params):
            data[param] = jac[x]
        data = balanceFactors(data, sep, factors, constraints, model)
        data = estimateFlows(data, sep, model, factors)
        estimates = estimateCum(data, knowns)
        its += 1

    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])
    cor = pearsonr(data.SIM_Estimates, data.Data)[0]
    return data, cor

def prodConstrain(observed, data, knowns, params, trips, sep, factors, constraints):
    model = 'prodConstrained'
    print 'production constrained model chosen'

    its = 0
    data = balanceFactors(data, sep, factors, constraints, model)
    data = estimateFlows(data, sep, model, factors)
    estimates = estimateCum(data, knowns)
    while abs(estimates - observed) > .0001:
        paramSingle = []
        for param in params:
            paramSingle.append(data[param].ix[0])
        functions = buildLLFunctions(data, params, factors, trips, sep, model, paramSingle)
        jac = new(paramSingle, data, params, sep, trips, functions)
        print jac, estimates, observed, its
        for x, param in enumerate(params):
            data[param] = jac[x]
        data = balanceFactors(data, sep, factors, constraints, model)
        data = estimateFlows(data, sep, model, factors)
        estimates = estimateCum(data, knowns)
        its += 1

    print "After " + str(its) + " runs, beta is : " + str(data["beta"].ix[0])

    cor = pearsonr(data.SIM_Estimates, data.Data)[0]
    return data, cor


def attConstrain(observed, data, knowns, params, trips, sep, factors, constraints):
    model = 'attConstrained'
    print 'attraction constrained model chosen'
    cor = pearsonr(data.results.SIM_Estimates, data.Data)[0]
    return data, cor


def unConstrain(observed, data, knowns, params, trips, sep, factors, constraints):
    model = 'unConstrained'
    print 'Unconstrained model chosen'
    cor = pearsonr(data.results.SIM_Estimates, data.Data)[0]
    return data, cor
