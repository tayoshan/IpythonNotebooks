# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import pandas as pd
import numpy as np
import entropy
import sys



class calibrate:

    def checkCols(self, data, trips, sep, factors, constraints, totalFlows):
            userInput = [trips, sep]
            if len(totalFlows) > 0:
                userInput.append(totalFlows)
            if factors != None:
                for key in factors.keys():
                    for factor in factors[key]:
                        userInput.append(factor)
            if len(constraints) > 0:
                for key in constraints.keys():
                    userInput.append(constraints[key])
                userInput = set(userInput)
                cols = set(data.columns)
            if userInput.issubset(cols) != True:
                print userInput
                print cols
                print userInput.issubset(cols)
                sys.exit('Not all input data exists in the dataset - check spelling to ensure columns and input match')

    def checkFactors(self, prodCon, attCon, factors):


        if (prodCon == False) & (attCon == False):

            try:
                if (len(factors['destinations']) < 1) | (len(factors['origins']) < 1):
                    sys.exit('In an un-constrained model there must be at least one variable for the origin propulsiveness ("origins" key in factors dict) and one for the destination attractiveness ("destinations" key in factors dict)')
            except TypeError:
                sys.exit('In an un-constrained model there must be at least one variable for the origin propulsiveness ("origins" key in factors dict) and one for the destination attractiveness ("destinations" key in factors dict)')
            except TypeError:
                sys.exit('In an un-constrained model there must be at least one variable for the origin propulsiveness ("origins" key in factors dict) and one for the destination attractiveness ("destinations" key in factors dict)')

    #def checkLen(self.trips, self.sep, self.factors, self):

    def __init__(self, data, trips, sep, dataForm='adj', diagFilter = True, cost='exp', factors=None, constraints={}, totalFlows=''):
        self.data = data
        self.cost = cost
        self.dataForm = dataForm
        self.constraints = constraints
        self.factors = factors
        self.trips = trips
        self.sep = sep
        self.totalFlows = totalFlows




        self.checkCols(self.data, self.trips, self.sep, self.factors, self.constraints, self.totalFlows)

        if diagFilter == True:
            if self.dataForm == 'adj':
                self.data = self.data[self.data['Origin'] != self.data['Destination']].reset_index(level = 0, drop = True)
            else:
                print 'Need to implement method to filter matrix for inter-zone data'

        self.prodCon = False
        self.attCon = False
        if 'production' in self.constraints.keys():
            self.prodCon = True
        if 'attraction' in self.constraints.keys():
            self.attCon = True

        self.checkFactors(self.prodCon, self.attCon, self.factors)




    def gravity(self):
        self.method = 'gravity'
        return self

    def regression(self):
        self.method = 'regression'
        self.results = entropy.regression(self.data, self.trips, self.cost, self.sep, self.factors, self.constraints)
        #Why does results.rsquared work but it isnt tabbale in Ipython?
        return self


    def mle(self, initialParams):
        self.method = 'mle'
        self.initialParams = initialParams



        observed, data, knowns, params = entropy.setup(self.data, self.trips, self.sep, self.cost, self.factors,self.constraints, self.prodCon, self.attCon, self.initialParams, self.totalFlows)

        if self.factors != None:
            entropy.checkParams(self.factors, self.initialParams)

        if (self.prodCon == True) & (self.attCon == True):
            self.model = 'dConstrained'
        elif (self.prodCon == True) & (self.attCon == False):
            self.model = 'prodConstrained'
        elif (self.prodCon == False) & (self.attCon == True):
            self.model = 'attConstrained'
        elif (self.prodCon == False) & (self.attCon == False):
            self.model = 'unConstrained'

        self.results, cor, sumStr = entropy.run(observed, data, knowns, params, self.trips, self.sep, self.cost, self.factors, self.constraints, self.model)
        self.results.rsquared = cor**2
        self.results.sumStr = sumStr


        return self


    def choice(self):
        self.method = 'choice'
        return self

    def compDest(self):
        self.method = 'compDest'
        return self


    def radiation(self):
        self.method = 'radiation'
        return self



