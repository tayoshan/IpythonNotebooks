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
        self.trips = self.data[trips]
        self.sep = self.data[sep]


    def gravity(self):
        self.method = 'gravity'
        return self

    def regression(self):
        self.method = 'regression'
        self.results = entropy.regression(self.data, self.trips, self.cost, self.sep, self.factors, self.constraints)
        #Why does results.rsquared work but it isnt tabbale in Ipython?
        return self


    def mle(self):
        self.method = 'mle'
        self.destCon = False
        self.attCon = False
        if 'attraction' in self.constraints.keys():
            self.destCon = True
        if 'production' in self.constraints.keys():
            self.attCon = True
        observed, data, knowns, params = entropy.setup(self.data, self.trips, self.sep, self.cost, self.factors,self.constraints, self.destCon, self.attCon)
        if self.destCon == True & self.attCon == True:
            data = entropy.dConstrain(observed, data, knowns, params, self.sep, self.factors, self.constraints)


        if self.destCon == True & self.attCon == False:
            data = entropy.attConstrain(observed, data, knowns, self.factors, self.constraints)
        if self.destCon == False & self.attCon == True:
            data = entropy.destConstrain(observed, data, knowns, self.factors, self.constraint)
        if self.destCon == False & self.destCon == False:
            data = entropy.unConstrain(observed, data, knowns, self.factors, self.constraints)

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



