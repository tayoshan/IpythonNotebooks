# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

class calibrate:
    def __init__(self, data, cost, model):
        self.data = data
        self.cost = cost
        self.model = model

    def mle(self):
        self.method = 'mle'
        return self




    def regression(self):
        self.method = 'regression'
        return self


    def check(self):
        return self.method, self.data, self.cost, self.model




