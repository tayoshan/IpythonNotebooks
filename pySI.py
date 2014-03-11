# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import pandas as pd
import numpy as np
import entropy



class calibrate:
    def __init__(self, data, trips, sep, dataForm='adj', diagFilter = True, cost='negexp', factors=None, constraints='unconstrained'):
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

    def mle(self):
        self.method = 'mle'
        return self

    def regression(self):
        self.method = 'regression'
        self.results = entropy.regression(self.data, self.trips, self.cost, self.sep, self.factors, self.constraints)
        #Why does results.rsquared work but it isnt tabbale in Ipython?
        return self

    def radiation(self):
        self.method = 'radiation'
        return self

    def choice(self):
        self.method = 'choice'
        return self

    def compDest(self):
        self.method = 'compDest'
        return self

    def check(self):
        return self.method, self.data, self.cost, self.model




