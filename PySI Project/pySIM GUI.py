#GUI for pySIM - Framework for spatial interaction modelling

import ttk
from Tkinter import *

import Tix
import tkFileDialog
import tkMessageBox
#import os
import sys
import FileIO


import numpy as np
import pandas as pd
import pySI as SI
import csv
import shapely
from math import sqrt, exp
from scipy import stats
from scipy.stats.mstats import mquantiles
from copy import copy
from datetime import datetime

FILE_TYPE = {"csv":0, "dbf":1, "txt":2, "shp":3, "ctl":4}
OPT_CRITERIA = {0: "AICc", 1: "AIC", 2: "BIC", 3: "CV"}
OPTION = {0:"OFF",1:"ON"}
DIST_TYPE = {0: "Euclidean distance", 1: "Spherical distance"}
models = {0:'Unconstrained', 1:'Production-Constrained', 2:'Attraction-Constrained', 3:'Doubly-Constrained'}


class Info(Toplevel):
    """
    output window
    """
    def __init__(self, master=None,summary=''):
        Toplevel.__init__(self, master)
        self.transient(master)
        #Frame.__init__(self, master)
        #self.pack()

        # make window stretchable
        top=self.winfo_toplevel()
        top.rowconfigure(0, weight=1)
        top.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        self.title("Summary")
        self.txt_summary = ttk.Tkinter.Text(self)
        #self.txt_summary["width"] = 100
        #self.txt_summary["scrollregion"]=(0, 0, 1200, 800)
        self.txt_summary.grid(row=0,column=0,sticky=N+S+E+W)
        self.txt_summary.insert(END,summary)
        #self.txt_summary = ttk.Label(self)
        #self.txt_summary["width"] = 200
        #self.txt_summary.grid(row=0,column=0)
        #self.txt_summary["text"] = summary

        self.scrollY = ttk.Scrollbar(self, orient=VERTICAL, command=self.txt_summary.yview )
        self.scrollY.grid(row=0, column=1, sticky=N+S )
        self.scrollX = ttk.Scrollbar ( self, orient=HORIZONTAL, command=self.txt_summary.xview )
        self.scrollX.grid(row=1, column=0, sticky=E+W )

        self.txt_summary["xscrollcommand"] = self.scrollX.set
        self.txt_summary["yscrollcommand"] = self.scrollY.set

class mainGUI(Frame):
    """
    GUI to implement GWR functions

    """
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()

        self.createWidget()


    def set_modtype(self):
        """
        set model type
        """
        val = self.var_mod.get()
        if val == 1:
            self.lstb_DA["state"] = DISABLED
            self.lstb_DA["bg"] = 'gray'
            self.lstb_OP["state"] = DISABLED
            self.lstb_OP["bg"] = 'gray'
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = NORMAL
            self.var_inFlow.set(1)
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = NORMAL
            self.var_outFlow.set(1)
        elif val == 2:
            self.lstb_DA["state"] = NORMAL
            self.lstb_DA["bg"] = 'white'
            self.lstb_OP["state"] = DISABLED
            self.lstb_OP["bg"] = 'gray'
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = NORMAL
            self.var_outFlow.set(1)
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = DISABLED
            self.var_inFlow.set(0)
        elif val == 3:
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = NORMAL
            self.var_inFlow.set(1)
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = NORMAL
            self.var_outFlow.set(1)
            self.lstb_DA["state"] = DISABLED
            self.lstb_DA["bg"] = 'gray'
            self.lstb_OP["state"] = DISABLED
            self.lstb_OP["bg"] = 'gray'

        elif val == 4:
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = NORMAL
            self.var_inFlow.set(1)
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = DISABLED
            self.var_outFlow.set(0)
            self.lstb_DA["state"] = DISABLED
            self.lstb_DA["bg"] = 'gray'
            self.lstb_OP["state"] = NORMAL
            self.lstb_OP["bg"] = 'white'
        elif val == 5:
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = NORMAL
            self.var_inFlow.set(1)
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = NORMAL
            self.var_outFlow.set(1)
            self.lstb_DA["state"] = DISABLED
            self.lstb_DA["bg"] = 'gray'
            self.lstb_OP["state"] = DISABLED
            self.lstb_OP["bg"] = 'gray'


        else:
            self.txt_Oi["state"] = DISABLED
            self.cbtn_outFlow["state"] = DISABLED
            self.var_outFlow.set(0)
            self.txt_Dj["state"] = DISABLED
            self.cbtn_inFlow["state"] = DISABLED
            self.var_inFlow.set(0)
            self.lstb_DA["state"] = NORMAL
            self.lstb_DA["bg"] = 'white'
            self.lstb_OP["state"] = NORMAL
            self.lstb_OP["bg"] = 'white'


    def checkTotalFlows(self, attribute):
        """
        turn Oi and Dj fields on/off depending on auto-compute check boxes
        """

        if str(attribute["state"]) == 'normal':
            attribute["state"] = DISABLED
        else:
            attribute["state"] = NORMAL




    def changeParams(self):
        self.params = []

        def updateParams(self):
            for param in self.params:
                print param
                if param == 'beta':
                    self.initParams['beta'] = self.txt_distance.get()
                    print self.initParams
                self.initParams[param] = eval('self.'+'txt'+param[4:]).get()
            paramWindow.destroy()


        self.params.append('beta')

        if str(self.lstb_OP["state"]) == 'normal':
            for param in [self.lstb_OP.get(i) for i in range(self.lstb_OP.size())]:
                self.params.append(param)

        if str(self.lstb_DA["state"]) == 'normal':
            for param in [self.lstb_DA.get(i) for i in range(self.lstb_DA.size())]:
                self.params.append(param)

        if any(self.params) < 1:
            tkMessageBox.showwarning('Warning', 'No variables have been selected which require estimation')
            return

        paramWindow = Toplevel(self)
        paramWindow.resizable(FALSE,FALSE)

        self.winFrame= ttk.Frame(paramWindow)
        self.winFrame["width"] = 200
        self.winFrame.grid(row=0,column=1, padx=0,pady=2)



        self.param_Form = ttk.LabelFrame(self.winFrame)
        self.param_Form["text"] = "Initial Paramter Value(s)"
        self.param_Form.grid(row=0,column=0,padx=5, pady=5)

        for x, param in enumerate(self.params):

            if param == 'beta':
                print '1'
                label = 'lbl'+ str('beta')
                txt = 'txt' + str('beta')
                var = 'var' + str('beta')

            else:
                print '2'
                label = 'lbl'+ str(param)
                txt = 'txt' + str(param)
                var = 'var' + str(param)

            setattr(self, label, ttk.Label(self.param_Form))
            if param == 'beta':
                print '1'
                eval('self.'+label)["text"] = 'beta'
            else:
                print '2'
                eval('self.'+label)["text"] = param
            eval('self.'+label).grid(row=x,column=0,padx=2,sticky=W)

            setattr(self, var, StringVar())
            setattr(self, txt, ttk.Entry(self.param_Form))
            eval('self.'+txt)["width"] = 10
            eval('self.'+txt)["textvariable"] = eval('self.'+var)
            if param == 'beta':
                eval('self.'+var).set(str(self.initParams['beta']))
            else:
                print param
                eval('self.'+var).set(self.initParams[param])
            eval('self.'+txt).grid(row=x,column=2,padx=2)


        self.btn_paramOK = ttk.Button(self.winFrame)
        self.btn_paramOK["text"] = "Ok"
        self.btn_paramOK["width"] = 10
        self.btn_paramOK["command"] = lambda arg1=None: updateParams(self) #
        self.btn_paramOK.grid(row=1, column=0,padx=8,pady=3,sticky=W)

        self.btn_paramCan = ttk.Button(self.winFrame)
        self.btn_paramCan["text"] = "Cancel"
        self.btn_paramCan["width"] = 10
        self.btn_paramCan["command"] =  lambda arg1=None: paramWindow.destroy()
        self.btn_paramCan.grid(row=1, column=2,padx=8,pady=3,sticky=W)



    def reset(self):
        """
        reset GUI to default settings
        """

        self.txt_open.delete(0, END)
        #self.lbx_mFormula.delete(0, END)
        self.txt_Observed.delete(0, END)
        self.txt_Oi.delete(0, END)
        self.txt_Dj.delete(0, END)
        self.lstb_OP.delete(0, END)
        self.lstb_DA.delete(0, END)
        self.lstb_allvars.delete(0,END)
        self.txt_flebeta.delete(0,END)
        self.txt_Origins.delete(0,END)
        self.txt_Destinations.delete(0,END)
        self.txt_Distance.delete(0,END)
        self.var_mod.set(0)
        self.initParams = {}
        self.params = []


    def openFile(self,target):
        """
        open file dialog for input data or prediction data
        target: 0--input data for regression; 1--input data for prediction
        """
        options = {}
        if target == 0: # regression data file
            options["filetypes"] = [("CSV", "*.csv"), ("DBF", "*.dbf"), ("TEXT", "*.txt"), ("SHP", "*.shp"), ("CONTROL", "*.ctl")]
            options["title"] = "Open File"
            options['defaultextension'] = '.txt'
            name_openfle = tkFileDialog.askopenfilename(**options)
        if target == 1: # prediction data file
            options["filetypes"] = [("CSV", "*.csv"), ("DBF", "*.dbf"), ("TEXT", "*.txt"), ("SHP", "*.shp")]
            options['defaultextension'] = '.txt'
            options["title"] = "Open File"
            name_openfle = tkFileDialog.askopenfilename(**options)
        if target == 2: # control file
            options["filetypes"] = [("CONTROL", "*.ctl")]
            options['defaultextension'] = '.ctl'
            options["title"] = "Save File"
            name_openfle = tkFileDialog.asksaveasfilename(**options)
        if target == 3:  # summary file
            options["filetypes"] = [("TEXT", "*.txt")]
            options['defaultextension'] = '.txt'
            options["title"] = "Save File"
            name_openfle = tkFileDialog.asksaveasfilename(**options)
        if target == 4 or target == 5: # local estimates and prediction file
            options["filetypes"] = [("CSV", "*.csv"), ("DBF", "*.dbf"), ("TEXT", "*.txt")]
            options['defaultextension'] = '.csv'
            options["title"] = "Save File"
            name_openfle = tkFileDialog.asksaveasfilename(**options)


        # openfile on own code: input data for regression model
        if name_openfle and target == 0:  #read_FILE = {0: read_CSV, 1: read_DBF, 2: read_TXT, 3: read_SHP, }
            fType = FILE_TYPE[name_openfle[-3:]]
            if fType > 3: # open control file
                self.fillGUI(name_openfle)

            else:
                self.reset()
                self.txt_open.insert(0,name_openfle) # show file name
                self.flepath_open = name_openfle
                self.fleType = FILE_TYPE[self.flepath_open[-3:]]

                allData = FileIO.read_FILE[self.fleType](self.flepath_open)
                if self.fleType == 3: # shapefile
                    self.coords = allData[0]
                    self.lstFlds = allData[1]
                    self.dicData = allData[2]
                else:
                    self.lstFlds = allData[0]
                    self.dicData = allData[1]

            # fill the fields list
            if fType <= 3:
                #self.lstb_allvars.delete(0,END)
                nflds = len(self.lstFlds)
                for i in range(nflds):
                    if i < 9:
                        id_fld = '00' + str(i+1)
                    else:
                        if i < 99:
                            id_fld = '0' + str(i+1)
                        else:
                            id_fld = str(i+1)
		    self.lstb_allvars.insert(i,id_fld+' '+self.lstFlds[i])

        # openfile on own code: input data for prediction
        if name_openfle and target == 1:  #read_FILE = {0: read_CSV, 1: read_DBF, 2: read_SHP, 3: read_TXT}
            self.txt_flepred.delete(0,END)
            self.txt_flepred.insert(0,name_openfle) # show file name
            #fleType = FILE_TYPE[name_openfle[-3:]]
            #allData = FileIO.read_FILE[fleType](name_openfle)
            #if fleType == 3: # shapefile
                #self.coords_pred = allData[0]
            #else:
                #self.coords_pred = allData[1]
	    # insert default names for local estimate file
            de_flepredout = name_openfle[:-3] + 'csv'
            self.txt_flepredout.delete(0,END)
            self.txt_flepredout.insert(0,de_flepredout)

        # get name for control file
        if name_openfle and target == 2:
            self.txt_flectrl.delete(0,END)
            self.txt_flectrl.insert(0,name_openfle)
	    # insert default names for summary file and local estimate file
            de_flesum = name_openfle[:-4] + '_summary.txt'
            self.txt_flesum.delete(0,END)
            self.txt_flesum.insert(0,de_flesum)
            de_flebeta = name_openfle[:-4] + '_listwise.csv'
            self.txt_flebeta.delete(0,END)
            self.txt_flebeta.insert(0,de_flebeta)

        # get name for summary file
        if name_openfle and target == 3:
            self.txt_flesum.delete(0,END)
            self.txt_flesum.insert(0,name_openfle)

        # get name for local estimates file
        if name_openfle and target == 4:
            self.txt_flebeta.delete(0,END)
            self.txt_flebeta.insert(0,name_openfle)

        # get name for local estimates file for prediction
        if name_openfle and target == 5:
            self.txt_flepredout.delete(0,END)
            self.txt_flepredout.insert(0,name_openfle)


    def addVars(self,txt_obj=None,lstb_obj=None):
        """
        add variables
        """
        #self.txt_id.insert(0,self.lstb_allvars.get[self.lstb_allvars.curselection()])
        if not txt_obj is None and self.lstb_allvars.curselection():
            if str(txt_obj["state"]) == 'normal':
                curr_id = self.lstb_allvars.curselection()
                curr_var = self.lstb_allvars.get(curr_id)
                if txt_obj.get(): # return the variable to the all variable list
                    pre_var = txt_obj.get()
                    self.lstb_allvars.insert(END,pre_var)
                    txt_obj.delete(0,END)
                txt_obj.insert(0,curr_var)#
                self.lstb_allvars.delete(curr_id)
                if txt_obj == self.txt_Distance:
                    self.initParams['beta'] = 0


        if not lstb_obj is None and self.lstb_allvars.curselection():
            if str(lstb_obj["state"]) == 'normal':
                curr_id = self.lstb_allvars.curselection()
                curr_var = self.lstb_allvars.get(curr_id)
                lstb_obj.insert(END,curr_var)
                self.lstb_allvars.delete(curr_id)
                self.initParams[curr_var[4:]] = 1



    def outVars(self,txt_obj=None,lstb_obj=None):
        """
        remove variables
        """
        if not lstb_obj is None and lstb_obj.curselection():
            curr_id = lstb_obj.curselection()
            curr_fld = lstb_obj.get(curr_id)
            lstb_obj.delete(curr_id)
            self.initParams.pop(curr_fld[4:])
            self.lstb_allvars.insert(END,curr_fld)#int(curr_fld[:3])

        if not txt_obj is None and txt_obj.get():
            if str(txt_obj["state"]) == 'normal':
                curr_fld = txt_obj.get()
                txt_obj.delete(0,len(curr_fld))
                if txt_obj == self.txt_Distance:
                    self.initParams.pop('beta')
                self.lstb_allvars.insert(END,curr_fld)# int(curr_fld[:3])


    def checkVars(self):
        """
        check validity of model setting before run
        """
        # 0 check open file
        if not self.txt_open.get():
            tkMessageBox.showwarning("Warning", "Please open a data file!")

        #Standard model terms
        if not self.txt_Origins.get():
            tkMessageBox.showwarning("Warning", "Please select a variable for the origin attributes!")
        if not self.txt_Destinations.get():
            tkMessageBox.showwarning("Warning", "Please select a variable for the destination attributes")
        if not self.txt_Observed.get():
            tkMessageBox.showwarning("Warning", "Please select a variable for observed flows!")
        if not self.txt_Distance.get():
            tkMessageBox.showwarning("Warning", "Please select a variable for distance/cost values!")

        #Model-specific terms

        if str(self.txt_Oi["state"]) == "normal":
            if not self.txt_Oi.get():
                tkMessageBox.showwarning("Warning", "Please select a variable for total outflow for model type: " + models[self.var_mod.get()] + "!")
        if str(self.txt_Dj["state"]) == "normal":
            if not self.txt_Dj.get():
                tkMessageBox.showwarning("Warning", "Please select a variable for total inflow for model type: " + models[self.var_mod.get()] + "!")

        if str(self.lstb_OP["state"]) == "normal":
            if not self.lstb_OP.get(0):
                tkMessageBox.showwarning("Warning", "Please select at least one variable for origin propulsiveness for model type: " + models[self.var_mod.get()] + "!")
        if str(self.lstb_DA["state"]) == "normal":
            if not self.lstb_DA.get(0):
                tkMessageBox.showwarning("Warning", "Please select at least one variable for destination attractivness for model type: " + models[self.var_mod.get()] + "!")


        #output files
        if not self.txt_flectrl.get():
            tkMessageBox.showwarning("Warning", "Please input the control file name!")

        if not self.txt_flesum.get():
            tkMessageBox.showwarning("Warning", "Please input the summary file name!")

        if not self.txt_flebeta.get():
            tkMessageBox.showwarning("Warning", "Please input the local estimates file name!")


    def saveFle_ctrl(self,fleName, dic):
	"""
	save control file
	Arguments:
	    fleName  : string
	               path of control file (.txt)
	    dic      : dictionary,
	               key: item, value: model settings
	"""
        with open(fleName, 'w') as txtfile:
        # 1 input file
            txtfile.write('Data file:\n') #line 1
            txtfile.write(''.join([dic["flepath_open"],'\n'])) # line 2

            # 2 variable names
            txtfile.write(': '.join(['Model type',str(dic["modtype"])])) # line 13
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of origins column ',self.txt_Origins.get()])) # line 3
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of destinations column ',self.txt_Destinations.get()])) # line 4
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of observed data/trips column ',self.txt_Observed.get()])) # line 5
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of distance/cost/separation column ',str(dic["distance"])])) # line 6
            txtfile.write('\n')
            txtfile.write(': '.join(['Distance/cost/separation function ',str(dic["ctype"])])) # line 9
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of total outflow column',str(dic["Oi"])])) # line 7
            txtfile.write('\n')
            txtfile.write(': '.join(['Name of total inflow column',str(dic["Dj"])])) # line 8
            txtfile.write('\n')

            # Vi
            nflds_unused = self.lstb_OP.size()
            txtfile.write(':'.join(['Origin propulsiveness variable columns',str(nflds_unused)])) # line 23
            txtfile.write('\n')
            for i in range(nflds_unused):
                txtfile.write(''.join([self.lstb_OP.get(i),'\n']))

            # Wj
            nflds_unused = self.lstb_DA.size()
            txtfile.write(':'.join(['Destination attractiveness variable columns',str(nflds_unused)])) # line 23
            txtfile.write('\n')
            for i in range(nflds_unused):
                txtfile.write(''.join([self.lstb_DA.get(i),'\n']))


            #Parameter initial values
            txtfile.write(''.join(['Number of initial parameter estimates:', str(len(self.initParams)), '\n']))
            for param in self.initParams:
                txtfile.write(''.join([param + ':' + str(self.initParams[param]),'\n']))


            # Unused variables
            nflds_unused = self.lstb_allvars.size()
            txtfile.write(':'.join(['Unused variables',str(nflds_unused)])) # line 23
            txtfile.write('\n')
            for i in range(nflds_unused):
                txtfile.write(''.join([self.lstb_allvars.get(i),'\n']))

            # control file
            txtfile.write('Path of control file:\n') # line 24
            txtfile.write(''.join([dic["flepath_ctrl"],'\n']))

            # summary file
            txtfile.write('Path of summary file:\n') # line 25
            txtfile.write(''.join([dic["flepath_sum"],'\n']))

            # listwise file
            txtfile.write('Path of local estimates file:\n') # line 26
            txtfile.write(''.join([dic["flepath_beta"],'\n']))


            txtfile.close()

    def saveFle_beta(self, results, flepath):
	"""
	save local estimates into a .txt file

	Arguments:
	    GWRMod      : GWR model

	    flepath     : string
	                  file path
	"""

        results.data.to_csv(flepath)

    '''
	# 3 write file
	fleType = FILE_TYPE[flepath[-3:]]
	FileIO.write_FILE[fleType](flepath,headers,lstBeta)
    '''


    def fillGUI(self, flepath):
	"""
	fill GUI using control file
	"""
        with open(flepath, 'rb') as txtfile:

            # 1 input file
            txtfile.readline()
            self.txt_open.delete(0,END) # clear the text first
            self.flepath_open = txtfile.readline().strip()
            self.txt_open.insert(0,self.flepath_open)

            # 2 variable names
            modelType = txtfile.readline().strip().split(':')[1].strip()
            self.var_mod.set(modelType)
            self.set_modtype()


            origins = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Origins.delete(0,END)
            self.txt_Origins.insert(0,origins)

            destinations = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Destinations.delete(0,END)
            self.txt_Destinations.insert(0,destinations)

            observed = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Observed.delete(0,END)
            self.txt_Observed.insert(0,observed)

            distances = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Distance.delete(0,END)
            self.txt_Distance.insert(0,distances)

            ctype = int(txtfile.readline().strip().split(':')[1])
            self.cmb_ctype.set(self.cmb_ctype['values'][ctype])

            Oi = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Oi.delete(0,END)
            self.txt_Oi.insert(0,Oi)

            Dj = txtfile.readline().strip().split(':')[1].strip()
            self.txt_Dj.delete(0,END)
            self.txt_Dj.insert(0,Dj)

            # x local
            n_Vi = int(txtfile.readline().strip().split(':')[1])
            self.lstb_OP.delete(0,END)
            for i in range(n_Vi):
                self.lstb_OP.insert(i,txtfile.readline().strip())

            # x local
            n_Wj = int(txtfile.readline().strip().split(':')[1])
            self.lstb_DA.delete(0,END)
            for i in range(n_Wj):
                self.lstb_DA.insert(i,txtfile.readline().strip())

            numparams = int(txtfile.readline().strip().split(':')[1])
            for i in range(numparams):
                paramLine = txtfile.readline().strip().split(':')
                print paramLine
                self.initParams[paramLine[0]] = int(paramLine[1])
                self.params.append(paramLine[0])


            # x unused
            nflds_unused = int(txtfile.readline().strip().split(':')[1])
            self.lstb_allvars.delete(0,END)
            for i in range(nflds_unused):
                self.lstb_allvars.insert(i,txtfile.readline().strip())


            # control file
            txtfile.readline()
            self.txt_flectrl.delete(0,END)
            self.txt_flectrl.insert(0,txtfile.readline().strip())

            # summary file
            txtfile.readline()
            self.txt_flesum.delete(0,END)
            self.txt_flesum.insert(0,txtfile.readline().strip())

            # listwise file
            txtfile.readline()
            self.txt_flebeta.delete(0,END)
            self.txt_flebeta.insert(0,txtfile.readline().strip())




	txtfile.close()

    def run_mod(self,ref):
        """
        run model
        """
        # 1 check model setting
        self.checkVars()

        #---------------------------------------setting information, stored in control file--------------------------------

        dic_ctrl = {}
        # open data file
        if self.flepath_open <> self.txt_open.get():
    	# read data file
            self.flepath_open = self.txt_open.get().strip()
            print self.flepath_open
            self.fleType = FILE_TYPE[self.flepath_open[-3:]]
            allData = FileIO.read_FILE[self.fleType](self.flepath_open)
            if self.fleType == 3: # shapefile
                self.coords = allData[0]
                self.lstFlds = allData[1]
                self.dicData = allData[2]
            else:
                self.lstFlds = allData[0]
                self.dicData = allData[1]

        dic_ctrl["flepath_open"] = self.flepath_open
        dic_ctrl["origins"] = self.txt_Origins.get()[4:].strip()
        dic_ctrl["destinations"] = self.txt_Destinations.get()[4:].strip()
        dic_ctrl["observed"] = self.txt_Observed.get()[4:].strip()
        dic_ctrl["distance"] = self.txt_Distance.get()
        dic_ctrl["Oi"] = self.txt_Oi.get()
        dic_ctrl["Dj"] = self.txt_Dj.get()
        dic_ctrl['ctype'] = self.cmb_ctype.current()


                # 9 get output file
        dic_ctrl["flepath_ctrl"] = self.txt_flectrl.get().strip()
        dic_ctrl["flepath_beta"] = self.txt_flebeta.get().strip()
        dic_ctrl["flepath_sum"] = self.txt_flesum.get().strip()


        # 5 get model type
        dic_ctrl["modtype"] = self.var_mod.get()



        dic_ctrl["Vi"] = []
        dic_ctrl["Wj"] = []
        self.Vi= self.lstb_OP.size()
        self.Wj = self.lstb_DA.size()
        for i in range(self.Vi):
            dic_ctrl["Vi"].append(self.lstb_OP.get(i)[4:].strip())
        for i in range(self.Wj):
            dic_ctrl["Wj"].append(self.lstb_DA.get(i)[4:].strip())



        #------------------------------------------run model--------------------------------------------------------
        begin_t = datetime.now()

        #Get data file and convert to pandas data frame
        data = pd.DataFrame(pd.read_csv(self.flepath_open))
        data[['Origin', 'Destination']] = data[['Origin', 'Destination']].astype(str)

        #Format GUI parameters for use in pySI API
        trips=str(self.txt_Observed.get()[4:])
        sep = str(self.txt_Distance.get()[4:])
        cost = str(self.cmb_ctype.get())
        if cost == 'Exponential':
            cost = 'exp'
        elif cost == 'Power':
            cost = 'pow'
        origin = str(self.txt_Origins.get()[4:])
        destination = str(self.txt_Destinations.get()[4:])
        model_type = self.var_mod.get()
        Oi = self.txt_Oi.get()[4:]
        Dj = self.txt_Dj.get()[4:]
        originFactors = []
        for factor in dic_ctrl["Vi"]:
            originFactors.append(factor)
        destFactors = []
        for factor in dic_ctrl["Wj"]:
            destFactors.append(factor)

        print self.initParams



        #Establish model and then run
        if model_type == 0:
            pass

        elif model_type == 1:
            constraints = {'production':origin}
            if len(Oi) > 0 and len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Oi=Oi, Dj=Dj)
            elif len(Oi) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Oi=Oi,  totalFlows=destination)
            else:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, totalFlows=destination)
            results = model.mle(self.initParams)

        elif model_type == 2:
            constraints = {'production':origin}
            if len(Oi) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, factors={'destinations':destFactors}, constraints=constraints, Oi=Oi)
            else:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, factors={'destinations':destFactors}, constraints=constraints, totalFlows=destination)
            results = model.mle(self.initParams)

        elif model_type == 3:
            constraints = {'attraction':destination}
            if len(Oi) > 0 and len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Oi=Oi, Dj=Dj)
            elif len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Dj=Dj, totalFlows=origin)
            else:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, totalFlows=origin)
            results = model.mle(self.initParams)

        elif model_type == 4:
            constraints = {'attraction':destination}
            if len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, factors={'origins':originFactors}, constraints=constraints, Dj=Dj)
            else:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, factors={'origins':originFactors}, constraints=constraints, totalFlows=origin)
            results = model.mle(self.initParams)
        elif model_type == 5:
            constraints = {'production':origin, 'attraction':destination}
            if len(Oi) > 0 and len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Oi=Oi, Dj=Dj)
            elif len(Oi) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Oi=Oi)
            elif len(Dj) > 0:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints, Dj=Dj)
            else:
                model = SI.calibrate(data=data, trips=trips, sep=sep, cost=cost, constraints=constraints)
            results = model.mle(self.initParams)


        #print trips, sep, cost, origin, destination, originFactors, destFactors, constraints



        # 10 output
        end_t = datetime.now()

        # 10.1 running time
        str_out = "%-21s: %s %s\n\n" % ('Program started at', datetime.date(begin_t), datetime.strftime(begin_t,"%H:%M:%S"))
        str_out += "%-21s: %s %s\n\n" % ('Program terminated at', datetime.date(end_t), datetime.strftime(end_t,"%H:%M:%S"))


        str_out += results.results.sumStr

        # 11 save to file
        # 11.1 save control file


        self.saveFle_ctrl(dic_ctrl["flepath_ctrl"],dic_ctrl)

        # 11.2 save summary file
        with open(dic_ctrl["flepath_sum"], 'w') as sumfile:
            sumfile.write(str_out)
            sumfile.close()

        # 11.3 save listwise file
	    self.saveFle_beta(results, dic_ctrl["flepath_beta"])


        # 12 output window
        self.infoWin = Info(self, str_out)
        self.infoWin.mainloop()

        #if not self.openInfo:
            #self.openInfo = True
            #self.infoWin = Info(self, str_out)
            #self.infoWin.mainloop()

        #if not self.infoWin.winfo_exists():
            #self.infoWin = Info(self, str_out)
            #self.infoWin.mainloop()

        #if self.infoWin.state == "withdrawn": # if the toplevel window is invisible, make it visible
            #self.infoWin.deiconify()



    def createWidget(self):
        """
        create controls

        """
        self.initParams = {}
        self.params = []

        # canvas containing all the widgets
        self.Container = Canvas(self)
        self.Container["height"] = 600
        self.Container.grid(row=0,column=0, columnspan=70, rowspan=80)

        #----------------------------Frame 1----------------------------------------

        self.left= ttk.Frame(self)
        self.left["width"] = 200
        self.left.grid(row=0,column=1, padx=0,pady=2)


        # frame 1: open data file

        self.frm_open = ttk.LabelFrame(self.left)
        self.frm_open["text"] = "Data File"
        self.frm_open["width"] = 200
        self.frm_open.grid(row=0,column=1, padx=5,pady=5)


        # text: open data file
        self.txt_open = ttk.Entry(self.frm_open)#
        self.txt_open["width"] = 24
        self.txt_open.grid(row=0,column=0, sticky=W, padx=5,pady=2)#

        # button: open file
        self.btn_open = ttk.Button(self.frm_open)#
        self.btn_open["width"] = 3
        self.img_open = PhotoImage(file=sys.path[0] + "\\img\\openfolder.gif")
        self.btn_open["image"] = self.img_open
        self.btn_open["command"] = lambda arg1=0: self.openFile(arg1)
        self.btn_open.grid(row=0,column=1,padx=5)#

        #-------------------------------Frame 2--------------------------------------------

        # frame 2: choose the model type

        self.frm_mod = ttk.LabelFrame(self.left)
        self.frm_mod["text"] = "Model Type"
        self.frm_mod.grid(row=1,column=1,padx=5, pady=5)



        # control variable for type of models
        self.var_mod = IntVar()

        # radiobutton: Unconstrained
        self.rbtn_gau = ttk.Radiobutton(self.frm_mod)
        self.rbtn_gau["text"] = "Poisson Unconstrained (Vi*Wj*f(Dij)"
        self.rbtn_gau["variable"] = self.var_mod
        self.rbtn_gau["value"] = 0
        self.rbtn_gau["width"] = 40
        self.rbtn_gau["command"] = self.set_modtype
        self.rbtn_gau.grid(row=0,column=0,padx=5,pady=5,sticky=W)

        # radiobutton: Origin Constrained
        self.rbtn_poson = ttk.Radiobutton(self.frm_mod)
        self.rbtn_poson["text"] = "Production-Constrained (Ai*Oi*Dj*f(Dij))"
        self.rbtn_poson["variable"] = self.var_mod
        self.rbtn_poson["value"] = 1
        self.rbtn_poson["width"] = 40
        self.rbtn_poson["command"] = self.set_modtype
        self.rbtn_poson.grid(row=1,column=0,padx=5,pady=5,sticky=W)

        # radiobutton: Origin Constrained (extra variable)
        self.rbtn_poson = ttk.Radiobutton(self.frm_mod)
        self.rbtn_poson["text"] = "Production-Constrained (Ai*Oi*Wj^u*f(Dij))"
        self.rbtn_poson["variable"] = self.var_mod
        self.rbtn_poson["value"] = 2
        self.rbtn_poson["width"] = 40
        self.rbtn_poson["command"] = self.set_modtype
        self.rbtn_poson.grid(row=2,column=0,padx=5,pady=5,sticky=W)


        # radiobutton: Destination Constrained
        self.rbtn_log = ttk.Radiobutton(self.frm_mod)
        self.rbtn_log["text"] = "Destination-Constrained (Bj*Dj*Oi*f(Dij))"
        self.rbtn_log["variable"] = self.var_mod
        self.rbtn_log["value"] = 3
        self.rbtn_log["width"] = 40
        self.rbtn_log["command"] = self.set_modtype
        self.rbtn_log.grid(row=3,column=0,padx=5,pady=5,sticky=W)


        # radiobutton: Destination Constrained (extra variable)
        self.rbtn_log = ttk.Radiobutton(self.frm_mod)
        self.rbtn_log["text"] = "Destination-Constrained (Bj*Dj*Vi^u*f(Dij))"
        self.rbtn_log["variable"] = self.var_mod
        self.rbtn_log["value"] = 4
        self.rbtn_log["width"] = 40
        self.rbtn_log["command"] = self.set_modtype
        self.rbtn_log.grid(row=4,column=0,padx=5,pady=5,sticky=W)

        # radiobutton: Doubly Constrained
        self.rbtn_log = ttk.Radiobutton(self.frm_mod)
        self.rbtn_log["text"] = "Doubly Constrained (Ai*Oi*Bj*Dj*f(Dij)"
        self.rbtn_log["variable"] = self.var_mod
        self.rbtn_log["value"] = 5
        self.rbtn_log["width"] = 40
        self.rbtn_log["command"] = self.set_modtype
        self.rbtn_log.grid(row=5,column=0,padx=5,pady=5,sticky=W)



        self.var_mod.set(0)

        self.frm_mFormula = ttk.LabelFrame(self.left)
        self.frm_mFormula["text"] = "Standard Model Terms"
        self.frm_mFormula.grid(row=2,column=1,padx=5, pady=5)



        # button: remove Origins attribute field
        self.btn_outOrigins = ttk.Button(self.frm_mFormula)
        self.btn_outOrigins["text"] = '>'
        self.btn_outOrigins["width"] = 2
        self.btn_outOrigins.grid(row=0,column=3)

        # button: add observed interaction
        self.btn_addOrigins = ttk.Button(self.frm_mFormula)
        self.btn_addOrigins["text"] = '<'
        self.btn_addOrigins["width"] = 2
        self.btn_addOrigins.grid(row=0,column=2)

        # label: Observed Interaction
        self.lbl_Origins = ttk.Label(self.frm_mFormula)
        self.lbl_Origins["text"] = 'Origins'
        self.lbl_Origins.grid(row=0,column=1,padx=2,sticky=W)

        # text: Observed Interaction
        self.txt_Origins = ttk.Entry(self.frm_mFormula)
        self.txt_Origins["width"] = 16
        self.txt_Origins.grid(row=0,column=0,padx=2)

        self.btn_outOrigins["command"] = lambda arg1=self.txt_Origins, arg2=None: self.outVars(arg1, arg2)
        self.btn_addOrigins["command"] = lambda arg1=self.txt_Origins, arg2=None: self.addVars(arg1, arg2)

        # button: remove Destinations attribute field
        self.btn_outDestinations = ttk.Button(self.frm_mFormula)
        self.btn_outDestinations["text"] = '>'
        self.btn_outDestinations["width"] = 2
        self.btn_outDestinations.grid(row=1,column=3)

        # button: add Oi
        self.btn_addDestinations = ttk.Button(self.frm_mFormula)
        self.btn_addDestinations["text"] = '<'
        self.btn_addDestinations["width"] = 2
        self.btn_addDestinations.grid(row=1,column=2)

        # label: Oi total Outflow
        self.lbl_Destinations = ttk.Label(self.frm_mFormula)
        self.lbl_Destinations["text"] = 'Destinations'
        self.lbl_Destinations.grid(row=1,column=1,padx=2,sticky=W)#

        # text: Oi total Outflow
        self.txt_Destinations = ttk.Entry(self.frm_mFormula)
        self.txt_Destinations["width"] = 16
        self.txt_Destinations.grid(row=1,column=0,padx=2)

        self.btn_addDestinations["command"] = lambda arg1=self.txt_Destinations, arg2=None: self.addVars(arg1, arg2)
        self.btn_outDestinations["command"] = lambda arg1=self.txt_Destinations, arg2=None: self.outVars(arg1, arg2)


        # button: remove observed interaction
        self.btn_outObserved = ttk.Button(self.frm_mFormula)
        self.btn_outObserved["text"] = '>'
        self.btn_outObserved["width"] = 2
        self.btn_outObserved.grid(row=2,column=3)

        # button: add observed interaction
        self.btn_addObserved = ttk.Button(self.frm_mFormula)
        self.btn_addObserved["text"] = '<'
        self.btn_addObserved["width"] = 2
        self.btn_addObserved.grid(row=2,column=2)

        # label: Observed Interaction
        self.lbl_Observed = ttk.Label(self.frm_mFormula)
        self.lbl_Observed["text"] = 'Observed Interaction'
        self.lbl_Observed.grid(row=2,column=1,padx=2,sticky=W)

        # text: Observed Interaction
        self.txt_Observed = ttk.Entry(self.frm_mFormula)
        self.txt_Observed["width"] = 16
        self.txt_Observed.grid(row=2,column=0,padx=2)

        self.btn_outObserved["command"] = lambda arg1=self.txt_Observed, arg2=None: self.outVars(arg1, arg2)
        self.btn_addObserved["command"] = lambda arg1=self.txt_Observed, arg2=None: self.addVars(arg1, arg2)


        # button: remove distance/cost variable
        self.btn_outDistance = ttk.Button(self.frm_mFormula)
        self.btn_outDistance["text"] = '>'
        self.btn_outDistance["width"] = 2
        self.btn_outDistance.grid(row=3,column=3)

        # button: add distance/cost variable
        self.btn_addDistance = ttk.Button(self.frm_mFormula)
        self.btn_addDistance["text"] = '<'
        self.btn_addDistance["width"] = 2
        self.btn_addDistance.grid(row=3,column=2)

        # label: distance/cost variable
        self.lbl_Distance = ttk.Label(self.frm_mFormula)
        self.lbl_Distance["text"] = 'Distances/Costs'
        self.lbl_Distance.grid(row=3,column=1,padx=2,sticky=W)

        # text: distance/cost variable
        self.txt_Distance = ttk.Entry(self.frm_mFormula)
        self.txt_Distance["width"] = 16
        self.txt_Distance.grid(row=3,column=0,padx=2)

        self.btn_outDistance["command"] = lambda arg1=self.txt_Distance, arg2=None: self.outVars(arg1, arg2)
        self.btn_addDistance["command"] = lambda arg1=self.txt_Distance, arg2=None: self.addVars(arg1, arg2)



        #Cost function
        self.lbl_cost = ttk.Label(self.frm_mFormula)
        self.lbl_cost["text"] = "Cost/Distance Function"
        self.lbl_cost["width"] = 24
        self.lbl_cost.grid(row=4,column=1,padx=2)#,pady=2


        # combox: Cost Function type selection
        self.cmb_ctype = ttk.Combobox(self.frm_mFormula)
        self.cmb_ctype["values"] = ["Power","Exponential"]
        self.cmb_ctype.set("Power")
        self.cmb_ctype["width"] = 13
        self.cmb_ctype.grid(row=4,column=0,columnspan=6,padx=2,sticky=W)


        #--------------------------------Frame 4----------------------------------------------

        # frame 4: list all the variables in the data file
        self.frm_varlist = ttk.LabelFrame(self)
        self.frm_varlist["text"] = "Variables"
        self.frm_varlist["width"] = 200
        self.frm_varlist.grid(row=0,column=2,rowspan=2,pady=2)

        # list box: list all the variables
        self.lstb_allvars = ttk.Tkinter.Listbox(self.frm_varlist)
        self.lstb_allvars["height"] = 25
        self.lstb_allvars["width"] = 23
        self.lstb_allvars.grid(row=0,column=0,padx=2,pady=3)

        # scroll bar
        self.scrollY_allvars = ttk.Scrollbar(self.frm_varlist, orient=VERTICAL, command=self.lstb_allvars.yview )
        self.scrollY_allvars.grid(row=0, column=1, sticky=N+S )
        self.scrollX_allvars = ttk.Scrollbar(self.frm_varlist, orient=HORIZONTAL, command=self.lstb_allvars.xview )
        self.scrollX_allvars.grid(row=1, column=0, sticky=E+W )

        self.lstb_allvars["xscrollcommand"] = self.scrollX_allvars.set
        self.lstb_allvars["yscrollcommand"] = self.scrollY_allvars.set


        #--------------------------------Frame 5----------------------------------------------

        # frame5: Model terms
        self.frm_mTerms = ttk.LabelFrame(self)
        self.frm_mTerms["text"] = "Model-Specific Terms"
        self.frm_mTerms["width"] = 200
        self.frm_mTerms.grid(row=0,column=3,rowspan=3,padx=10,pady=5, columnspan=5)




        # button: remove Oi
        self.btn_outOi = ttk.Button(self.frm_mTerms)
        self.btn_outOi["text"] = '<'
        self.btn_outOi["width"] = 2
        self.btn_outOi.grid(row=1,column=0)

        # button: add Oi
        self.btn_addOi = ttk.Button(self.frm_mTerms)
        self.btn_addOi["text"] = '>'
        self.btn_addOi["width"] = 2
        self.btn_addOi.grid(row=1,column=1)

        # label: Oi total Outflow
        self.lbl_Oi = ttk.Label(self.frm_mTerms)
        self.lbl_Oi["text"] = 'Oi (Total Outflow)'
        self.lbl_Oi.grid(row=1,column=2,padx=2,sticky=W)#

        # text: Oi total Outflow
        self.txt_Oi = ttk.Entry(self.frm_mTerms)
        self.txt_Oi["width"] = 16
        self.txt_Oi["state"] = DISABLED
        self.txt_Oi.grid(row=1,column=3,padx=2)

        self.btn_addOi["command"] = lambda arg1=self.txt_Oi, arg2=None: self.addVars(arg1, arg2)
        self.btn_outOi["command"] = lambda arg1=self.txt_Oi, arg2=None: self.outVars(arg1, arg2)

        # checkButton: compute total outflows from data
        self.var_outFlow = IntVar()
        self.cbtn_outFlow = ttk.Checkbutton(self.frm_mTerms)
        self.cbtn_outFlow["text"] = "Auto-compute total outflows"
        self.cbtn_outFlow["state"] = DISABLED
        self.cbtn_outFlow["variable"] = self.var_outFlow
        self.var_outFlow.set(0)
        self.cbtn_outFlow["command"] = lambda arg1=self.txt_Oi: self.checkTotalFlows(arg1)
        self.cbtn_outFlow.grid(row=2,column=1, columnspan=3, sticky=W, padx=5, pady=1)

        # button: remove Dj
        self.btn_outDj = ttk.Button(self.frm_mTerms)
        self.btn_outDj["text"] = '<'
        self.btn_outDj["width"] = 2
        self.btn_outDj.grid(row=3,column=0)

        # button: add Dj
        self.btn_addDj = ttk.Button(self.frm_mTerms)
        self.btn_addDj["text"] = '>'
        self.btn_addDj["width"] = 2
        self.btn_addDj.grid(row=3,column=1)

        # label: Dj total Inflow
        self.lbl_Dj = ttk.Label(self.frm_mTerms)
        self.lbl_Dj["text"] = 'Dj (Total Inflow)'
        self.lbl_Dj.grid(row=3,column=2,padx=2,sticky=W)

        # text: Dj total Inflow
        self.txt_Dj = ttk.Entry(self.frm_mTerms)
        self.txt_Dj["width"] = 16
        self.txt_Dj["state"] = DISABLED
        self.txt_Dj.grid(row=3,column=3,padx=2)

        self.btn_addDj["command"] = lambda arg1=self.txt_Dj, arg2=None: self.addVars(arg1, arg2)
        self.btn_outDj["command"] = lambda arg1=self.txt_Dj, arg2=None: self.outVars(arg1, arg2)

        # checkButton: compute total inflows from data
        self.var_inFlow = IntVar()
        self.cbtn_inFlow = ttk.Checkbutton(self.frm_mTerms)
        self.cbtn_inFlow["text"] = "Auto-compute total inflows"
        self.cbtn_inFlow["state"] = DISABLED
        self.cbtn_inFlow["variable"] = self.var_inFlow
        self.var_inFlow.set(0)
        self.cbtn_inFlow["command"] = lambda arg1=self.txt_Dj: self.checkTotalFlows(arg1)
        self.cbtn_inFlow.grid(row=4,column=1, columnspan=3, sticky=W, padx=5, pady=2)

        # label: Origin Propulsiveness
        self.lbl_OP = ttk.Label(self.frm_mTerms)
        self.lbl_OP["text"] = 'Vi (Origin Attributes)'
        self.lbl_OP.grid(row=5,column=2,padx=2,pady=2,sticky=W)

        # listbox: Origin Propulsiveness
        self.lstb_OP = ttk.Tkinter.Listbox(self.frm_mTerms)
        self.lstb_OP["height"] = 8
        self.lstb_OP["width"] = 25
        #self.lstb_OP.insert(0,"000 Intercept")
        self.lstb_OP.grid(row=6,column=2,columnspan=2,padx=2,sticky=W)


        # scrollbar: Origin Propulsiveness
        self.scrollY_OP = ttk.Scrollbar(self.frm_mTerms, orient=VERTICAL, command=self.lstb_OP.yview )
        self.scrollY_OP.grid(row=6, column=3, sticky=N+S+E )
        self.scrollX_OP = ttk.Scrollbar(self.frm_mTerms, orient=HORIZONTAL, command=self.lstb_OP.xview )
        self.scrollX_OP.grid(row=6, column=2, columnspan=2,sticky=E+W+S )

        self.lstb_OP["xscrollcommand"] = self.scrollX_OP.set
        self.lstb_OP["yscrollcommand"] = self.scrollY_OP.set

        # button: origin propulsiveness
        self.btn_outOP = ttk.Button(self.frm_mTerms)
        self.btn_outOP["text"] = '<'
        self.btn_outOP["width"] = 2
        self.btn_outOP["command"] = lambda arg1=None, arg2=self.lstb_OP: self.outVars(arg1, arg2)
        self.btn_outOP.grid(row=6,column=0)

        # button: add origin propulsiveness
        self.btn_addOP = ttk.Button(self.frm_mTerms)
        self.btn_addOP["text"] = '>'
        self.btn_addOP["width"] = 2
        self.btn_addOP["command"] = lambda arg1=None, arg2=self.lstb_OP: self.addVars(arg1, arg2)
        self.btn_addOP.grid(row=6,column=1)

        # label: destination attraction
        self.lbl_DA = ttk.Label(self.frm_mTerms)
        self.lbl_DA["text"] = 'Wj (Destination Attributes)'
        self.lbl_DA.grid(row=7,column=2,padx=2,sticky=W)

        # listbox: destination attraction
        self.lstb_DA = ttk.Tkinter.Listbox(self.frm_mTerms)
        self.lstb_DA["height"] = 8
        self.lstb_DA["width"] = 25
        self.lstb_DA["state"] = NORMAL
        self.lstb_DA.grid(row=8,column=2,columnspan=2,padx=2,pady=3,sticky=W)

        # scrollbar: destination attraction
        self.scrollY_DA = ttk.Scrollbar(self.frm_mTerms, orient=VERTICAL, command=self.lstb_DA.yview )
        self.scrollY_DA.grid(row=8, column=3, sticky=N+S+E )
        self.scrollX_DA = ttk.Scrollbar(self.frm_mTerms, orient=HORIZONTAL, command=self.lstb_DA.xview )
        self.scrollX_DA.grid(row=8, column=2, columnspan=2,sticky=E+W+S )

        self.lstb_DA["xscrollcommand"] = self.scrollX_DA.set
        self.lstb_DA["yscrollcommand"] = self.scrollY_DA.set

        # button: remove destination attraction
        self.btn_outDA = ttk.Button(self.frm_mTerms)
        self.btn_outDA["text"] = '<'
        self.btn_outDA["width"] = 2
        self.btn_outDA["command"] = lambda arg1=None, arg2=self.lstb_DA: self.outVars(arg1, arg2)
        self.btn_outDA.grid(row=8,column=0)

        # button: add destination attraction
        self.btn_addDA = ttk.Button(self.frm_mTerms)
        self.btn_addDA["text"] = '>'
        self.btn_addDA["width"] = 2
        self.btn_addDA["command"] = lambda arg1=None, arg2=self.lstb_DA: self.addVars(arg1, arg2)
        self.btn_addDA.grid(row=8,column=1)


        self.btn_initParams = ttk.Button(self.frm_mTerms)
        self.btn_initParams["text"] = "Change initial parameter value(s)"
        self.btn_initParams["width"] = 30
        self.btn_initParams["command"] = lambda arg1=self: self.changeParams()
        self.btn_initParams.grid(row=9, column=0, columnspan=4,padx=8,pady=3,sticky=W+E)

        #--------------------------------Frame 7----------------------------------------------

        # frame: output
        self.frm_output = ttk.LabelFrame(self)
        self.frm_output["text"] = "Outputs"
        #self.frm_mod["width"] = 200
        self.frm_output.grid(row=5,column=1,columnspan=3,padx=5,pady=3,sticky=W)

        # label: control file
        self.lbl_flectrl = ttk.Label(self.frm_output)
        self.lbl_flectrl["text"] = "Control file"
        self.lbl_flectrl.grid(row=0, column=0,padx=2,sticky=W)

        # text: control file
        self.txt_flectrl = ttk.Entry(self.frm_output)
        self.txt_flectrl["width"] = 55
        self.txt_flectrl.grid(row=0, column=1,padx=2,sticky=W)

        # button: control file
        self.btn_flectrl = ttk.Button(self.frm_output)
        self.btn_flectrl["width"] = 3
        self.btn_flectrl["command"] = lambda arg1=2: self.openFile(arg1)
        self.btn_flectrl["image"] = self.img_open
        self.btn_flectrl.grid(row=0, column=2,padx=5)

        # label: summary file
        self.lbl_flesum = ttk.Label(self.frm_output)
        self.lbl_flesum["text"] = "Summary file"
        self.lbl_flesum.grid(row=1, column=0,padx=2,sticky=W)

        # text: summary file
        self.txt_flesum = ttk.Entry(self.frm_output)
        self.txt_flesum["width"] = 55
        self.txt_flesum.grid(row=1, column=1,padx=2,sticky=W)

        # button: summary file
        self.btn_flesum = ttk.Button(self.frm_output)
        self.btn_flesum["width"] = 3
        self.btn_flesum["image"] = self.img_open
        self.btn_flesum["command"] = lambda arg1=3: self.openFile(arg1)
        self.btn_flesum.grid(row=1, column=2,padx=5)



        # label: Output file
        self.lbl_flebeta = ttk.Label(self.frm_output)
        self.lbl_flebeta["text"] = "Output"
        self.lbl_flebeta.grid(row=2, column=0,padx=2,sticky=W)

        # text: Output file
        self.txt_flebeta = ttk.Entry(self.frm_output)
        self.txt_flebeta["width"] = 55
        self.txt_flebeta.grid(row=2, column=1,padx=2,sticky=W)

        # button: output file
        self.btn_flebeta = ttk.Button(self.frm_output)
        self.btn_flebeta["width"] = 3
        self.btn_flebeta["image"] = self.img_open
        self.btn_flebeta["command"] = lambda arg1=4: self.openFile(arg1)
        self.btn_flebeta.grid(row=2, column=2,padx=5)


        #-------------------------------------------run button------------------------------------
        # button: run
        self.btn_run = ttk.Button(self)
        self.btn_run["text"] = "Run"
        self.btn_run["width"] = 10
        self.btn_run["command"] = lambda arg1=1: self.run_mod(1) #
        self.btn_run.grid(row=5, column=4,padx=8,pady=3,sticky=E)


if __name__=='__main__':


    root = Tk()
    root.title("pySI: A Spatial Interaction Modelling Framework")

    app = mainGUI(master=root)
    app.mainloop()
    #root.destroy()