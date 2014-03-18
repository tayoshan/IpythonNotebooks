# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import pandas as pd
import numpy as np
import entropy



class calibrate:
    def __init__(self, data, trips, sep, dataForm='adj', diagFilter = True, cost='negexp', factors=None, constraints={}):
        self.data = data
        self.cost = cost
        self.dataForm = dataForm
        self.constraints = constraints
        self.factors = factors
        if diagFilter == True:
            if self.dataForm == 'adj':
                self.data = self.data[self.data['Origin'] != self.data['Destination']].reset_index(level = 0, drop = True)
            else:
                print 'Need to implement method to filter matrix for inter-zone data'
        self.trips = trips
        self.sep = sep


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
        self.prodCon = False
        self.attCon = False
        self.initialParams = initialParams
        if 'production' in self.constraints.keys():
            self.prodCon = True
        if 'attraction' in self.constraints.keys():
            self.attCon = True

        observed, data, knowns, params = entropy.setup(self.data, self.trips, self.sep, self.cost, self.factors,self.constraints, self.prodCon, self.attCon, self.initialParams)

        if (self.prodCon == True) & (self.attCon == True):
            self.results, cor = entropy.dConstrain(observed, data, knowns, params, self.trips, self.sep, self.factors, self.constraints)

        elif (self.prodCon == True) & (self.attCon == False):
            self.results, cor = entropy.prodConstrain(observed, data, knowns, params, self.trips, self.sep, self.factors, self.constraints)

        elif (self.prodCon == False) & (self.attCon == True):
            self.results, cor = entropy.attConstrain(observed, data, knowns, params, self.trips, self.sep, self.factors, self.constraints)

        elif (self.prodCon == False) & (self.attCon == False):
            self.results, cor = entropy.unConstrain(observed, data, knowns, params, self.trips, self.sep, self.factors, self.constraints)


        self.results.rsquared = cor**2

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



