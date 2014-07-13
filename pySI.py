# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import pandas as pd
import numpy as np
import entropy
import sys




class calibrate:
    """
    Calibration class is a set calibration routines for estimating parameters of SI models

    """

    def __init__(self, data, origins, destinations, trips, sep, dataForm='adj', diagFilter = True, cost='pow', factors=None, constraints={}, Oi=None, Dj=None, totalFlows=None):
        """
        initialize a model to be calibrated

        Parameters
        ----------
        data :          pandas data frame or csv with header composed of columns of necessary data
        origins:        str, value of column name for origins
        destinations:   str, value of column name for destinations
        trips:          str, value of column name for observed flows in data object
        sep:            str, deafults to 'pow', value of column name for distance values in data object
        dataForm:       not yet being used but will be to allow matrix data format
        diagFilter:     boolean, defaults to True, to filter intra-zonal flows
        cost:           str, either 'exp' for exponential or 'pow' for power distance function
        factors:        dict, defaults to None,  whose keys can be only 'origins' and/or 'destinations' and whose values are each a list of strings for names of columns in data object
        constraints:    dict, deaults to empty dict, whose keys can only 'production' and/or 'attraction' and whose values strings for the name of columns representing origins and/or destinations respectively
        Oi:             str, defaults to None, to provide total outflow column rather than calculate from data
        Dj:             str, defaults to None, to provide totoal inflow column rather than calculate from the data
        totalFlows:     str, defaults to None, to provide a field to tally total in/out flows if none is provided (only for GUI extension)
        """
        self.data = data
        self.origins = origins
        self.destinations = destinations
        self.cost = cost
        self.dataForm = dataForm
        self.constraints = constraints
        self.factors = factors
        self.trips = trips
        self.sep = sep
        self.Oi = Oi
        self.Dj = Dj
        self.totalFlows = totalFlows
        self.status = False
        self.data[[origins, destinations]] = self.data[[origins, destinations]].astype(str)

        if diagFilter == True:
            if self.dataForm == 'adj':
                self.data = self.data[self.data[self.origins] != self.data[self.destinations]].reset_index(level = 0, drop = True)
            else:
                print 'Need to implement method to filter matrix for inter-zone data'

        self.prodCon = False
        self.attCon = False
        if 'production' in self.constraints.keys():
            self.prodCon = True
        if 'attraction' in self.constraints.keys():
            self.attCon = True

    def summary(self):
        if self.method == 'regression':
            print self.results.summary()
        elif self.method == 'mle':
            print self.results.summaryString

    def rsquared(self):
        if self.method == 'regression':
            print self.results.rsquared
        elif self.method == 'mle':
            print self.results.rsquared

    def regression(self):
        """
        calibrate max-entropy model using regression parameter estimation

        Parameters
        ---------
        None
        """
        self.method = 'regression'
        self.results = entropy.regression(self.data, self.trips, self.cost, self.sep, self.factors, self.constraints)
        self.status = True
        return self

    def mle(self, initialParams):
        """
        calibrate max-entropy gravity model using maximum lilelihood parameter estimation

        Parameters
        ----------
        initialParams: dict, whose keys must be either 'beta' or equal to factors which have been specificed for parameter estimation and whose values are an integer
        """
        self.method = 'mle'
        self.initialParams = initialParams

        observed, data, knowns, params = entropy.setup(self.data, self.trips, self.sep, self.cost, self.factors,self.constraints, self.prodCon, self.attCon, self.initialParams, self.Oi, self.Dj, self.totalFlows)

        if (self.prodCon == True) & (self.attCon == True):
            self.model = 'dConstrained'
        elif (self.prodCon == True) & (self.attCon == False):
            self.model = 'prodConstrained'
        elif (self.prodCon == False) & (self.attCon == True):
            self.model = 'attConstrained'
        elif (self.prodCon == False) & (self.attCon == False):
            self.model = 'unConstrained'

        self.results, cor, summary, finalParams = entropy.run(observed, data, self.origins, self.destinations, knowns, params, self.trips, self.sep, self.cost, self.factors, self.constraints, self.model, self.initialParams)
        self.results.rsquared = cor**2
        self.results.summaryString = summary
        self.status = True
        self.finalParams = finalParams
        return self


'''
class simulate():
    """
    Simulation class is a set of functions for simulating flows given a calibrated model or a model without parameters
    """
    #Function to calculate Ai balancing factor values
    def calcAi(data, sep, cost, factors, model):
        """
        calculate Ai balancing factor
        """

        #add distance data with appropriate functional form
        if cost == 'exp':
            Ai = np.exp(data[sep]*data["beta"])
        elif cost == 'pow':
            Ai = (data[sep]**data["beta"])
        else:
            sys.exit("The distance/cost function must be either 'pow' or 'exp'.")

        #Add factors
        if factors != None:
            for factor in factors['destinations']:
                Ai = Ai*(data[factor]**data[factor + 'Param'])

        else:
            Ai = Ai*data['Dj']

        #If model is doubly constrained add destination balancing factor
        if model == 'dConstrained':
            Ai = Ai*data["Bj"]


        data["Ai"] = Ai



    #Function to Calculate Bj values
    def calcBj(data, sep, cost, factors, model):
        """
        calculate Bj balancing factor
        """

        #add distance data with appropriate functional form
        if cost == 'exp':
            Bj = np.exp(data[sep]*data["beta"])
        elif cost == 'pow':
            Bj = (data[sep]**data["beta"])
        else:
            sys.exit("The distance/cost function must be either 'pow' or 'exp'.")

        #Add factors
        if factors != None:
            for factor in factors['origins']:
                Bj = Bj*(data[factor]**data[factor + 'Param'])

        else:
            Bj = Bj*data['Oi']

        #If model is doubly constrained add origin balancing factor
        if model == 'dConstrained':
            Bj = Bj*data["Ai"]



        data["Bj"] = Bj


    #Function to check if Ai and Bj have stabilised, if not return to step 2
    #Only get called for att, prod, and doubly constrained - not unconstrained
    def balanceFactors(data, sep, cost, factors, constraints, model):
        """
        calculate balancing factors and balance the balancing factors if doubly constrained model
        """
        its = 0
        cnvg = 1
        while cnvg > .001:
            its = its + 1
            #If model is prod or doubly constrained
            if model != 'attConstrained':
                calcAi(data, sep, cost, factors, model)
                AiBF = (data.groupby(data[constraints['production']].name).aggregate({"Ai": np.sum}))
                AiBF["Ai"] = 1/AiBF["Ai"]
                updates = AiBF.ix[pd.match(data[constraints['production']], AiBF.index), "Ai"]
                data["Ai"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Ai"]
                #If model is prod constrained stop here - dont need to balance
                if model == 'prodConstrained':
                    break
                if its == 1:
                    data["OldAi"] = data["Ai"]
                else:
                    data["diff"] = abs((data["OldAi"] - data["Ai"])/data["OldAi"])
                    data["OldAi"] = data["Ai"]
            #If model is att or doubly constrained
            if model != 'prodConstrained':
                calcBj(data, sep, cost, factors, model)
                BjBF = data.groupby(data[constraints['attraction']].name).aggregate({"Bj": np.sum})
                BjBF["Bj"] = 1/BjBF["Bj"]
                updates = BjBF.ix[pd.match(data[constraints['attraction']], BjBF.index), "Bj"]
                data["Bj"] = updates.reset_index(level=0, drop=True) if(updates.notnull().any()) else data["Bj"]
                if its == 1:
                    #If model is att constrained stop here - dont need to balance
                    if model == 'attConstrained':
                        break
                    data["OldBj"] = data["Bj"]
                else:
                    data["diff"] = abs((data["OldBj"] - data["Bj"])/data["OldBj"])
                    data["OldBj"] = data["Bj"]
            cnvg = np.sum(data["diff"])
            #print cnvg, its

    return data


    def setDistance():
        pass

    def setFactor():
        pass

    def setOi():
        pass

    def setDj():
        pass

    def setParam():
        pass


    def isCalibrated(self):
        if self.calibratedModel != None:
            if self.calibratedModel.status == True:
                self.data = self.calibratedModel.data
                self.origins= self.calibratedModel.origins
                self.destinations = self.calibratedModel.destinations
                self.cost = self.calibratedModel.cost
                self.dataForm = self.calibratedModel.dataForm
                self.constraints = self.calibratedModel.constraints
                self.factors = self.calibratedModel.factors
                self.trips = self.calibratedModel.trips
                self.sep = self.calibratedModel.sep
                self.Oi = self.calibratedModel.Oi
                self.Dj = self.calibratedModel.Dj
                self.totalFlows = self.calibratedModel.totalFlows
                self.model = self.calibratedModel.model

            else:
                print 'calibrate object has not been calibrated yet'

        else:
                self.factors = None
                for input in self.inputs:
                    setattr(self, input, self.inputs[input])
                self.data["Bj"] = 1.0
                self.data["Ai"] = 1.0
                self.data["OldAi"] = 10.000000000
                self.data["OldBj"] = 10.000000000
                self.data["diff"] = abs((self.data["OldAi"] - self.data["Ai"])/self.data["OldAi"])


    def __init__(self, calibratedModel=None, **inputs):#data, origins, destinations, trips, sep, dataForm='adj', diagFilter = True, cost='pow', factors=None, constraints={}, Oi=None, Dj=None, totalFlows=None):
        self.calibratedModel = calibratedModel
        self.inputs = inputs
        self.isCalibrated()

    def entropy(self, removeNode=[], setParams={}, setFactors={}):
        """
        using max entropy model with parameters derived from pySI.calibrate (regression, mle)
        """
        #Remove a node from the interaction network (removes inflows and outflows from the node)
        if len(removeNode) > 0:
            self.data = self.data[~((self.data[self.origins].isin(removeNode)) | (self.data[self.destinations].isin(removeNode)))]
            self.data = self.data.reset_index()
            self.data = self.data.drop('index', axis=1)

        #This is to add a new node and is unfinished as it will require distances to be computed on the fly
        #Will need to add in the addNode=[] parameter into the function
        #if len(addNode) > 0:
            #for node in addNode:
                #origins = set(self.data[self.origins].unique())
                #destinations = set(self.data[self.destinations].unique())
                #inflows = it.product(origins, [node])
                #outflows = it.product([node], destinations)

        #Check for model type based on constraints
        self.prodCon = False
        self.attCon = False
        if 'production' in self.constraints.keys():
            self.prodCon = True
        if 'attraction' in self.constraints.keys():
            self.attCon = True

        if (self.prodCon == True) & (self.attCon == True):
            self.model = 'dConstrained'
        elif (self.prodCon == True) & (self.attCon == False):
            self.model = 'prodConstrained'
        elif (self.prodCon == False) & (self.attCon == True):
            self.model = 'attConstrained'
        elif (self.prodCon == False) & (self.attCon == False):
            self.model = 'unConstrained'

        if self.model == 'prodConstrained':
            Oi = self.data.groupby(self.data[self.constraints['production']]).aggregate({self.trips: np.sum})
            self.data["Oi"] = Oi.ix[pd.match(self.data[self.constraints['production']], Oi.index)].reset_index()[self.trips]

        if self.model == 'attConstrained':
            Dj = self.data.groupby(self.data[self.constraints['attraction']]).aggregate({self.trips: np.sum})
            self.data["Dj"] = Dj.ix[pd.match(self.data[self.constraints['attraction']], Dj.index)].reset_index()[self.trips]

        if self.model == 'dConstrained':
            Oi = self.data.groupby(self.data[self.constraints['production']]).aggregate({self.trips: np.sum})
            self.data["Oi"] = Oi.ix[pd.match(self.data[self.constraints['production']], Oi.index)].reset_index()[self.trips]

            Dj = self.data.groupby(self.data[self.constraints['attraction']]).aggregate({self.trips: np.sum})
            self.data["Dj"] = Dj.ix[pd.match(self.data[self.constraints['attraction']], Dj.index)].reset_index()[self.trips]

        for factor in setParams:
            if factor == 'beta':
                self.data['beta'] = setParams[factor]
            self.data[factor + 'Param'] = setParams[factor]

        for factor in setFactors:
            if type(setFactors[factor]) == list:
                for each in setFactors[factor]:
                    association, location, val = each
                    if association.lower() == 'origins':
                        association = self.origins
                    elif association.lower() == 'destinations':
                        association = self.destinations
                    else:
                        print 'Please specify whether you want to change a variable associated with an origin(s) or a destination(s) using "origins" or "destinations"'
                    self.data[factor][self.data[association] == location] = val
            else:
                association, location, val = setFactors[factor]
                if association.lower() == 'origins':
                    association = self.origins
                elif association.lower() == 'destinations':
                    association = self.destinations
                else:
                    print 'Please specify whether you want to change a variable associated with an origin(s) or a destination(s) using "origins" or "destinations"'
                self.data[factor][self.data[association] == location] = val


        self.data = entropy.balanceFactors(self.data, self.sep, self.cost, self.factors, self.constraints, self.model)

        self.estimates = entropy.estimateFlows(self.data, self.sep, self.cost, self.model, self.factors)

    #To be implemented
    '''
    def radiation(self):
        """
        Radiation model implementation
        """
        self.method = 'radiation'
        return self
<<<<<<< HEAD
    '''
=======
'''
>>>>>>> origin/master
