#! /usr/bin/env python

# The tRecX package is free for personal use.
# Any commercial use of the code or parts of it is excluded.
# Restrictions for academic use apply. 
# 
# See terms of use in the LICENSE file included with the source distribution
# 
# Copyright (c) 2015,2016,2017,2018 by Armin Scrinzi
# End of license
 
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.interpolate
from scipy.interpolate import griddata
import copy  as cp
from .tools import *
from .plotTools import *

class Compare:
    def __init__(self,linY):
        self.tolerance=1.e-4
        self.errRange=1.e-5  # plot in [maxErr*errRange,maxErr]
        self.plotRange=1.e-7 # plot in [max*plotRange,max]
        self.linY=linY
        self.args,self.flags=argsAndFlags();
        self.fax,(self.ax1,self.ax2)=plt.subplots(nrows=2,sharex=True)

        gs = gridspec.GridSpec(2,1,height_ratios=[3,3])
        self.ax1=plt.subplot(gs[0])
        self.ax2=plt.subplot(gs[1])

        # axis labels
        self.ax2.set_xlabel("x (a.u.)",fontsize=14)
        self.ax2.set_ylabel("Relative error",fontsize=14)

        self.fmin= 1e10
        self.fmax=-1e10
        self.errmax=-1
        self.errmin=1

    def addPlot(self,xVal,fVal,lab,colr):
        self.fmin=min(self.fmin,np.min(fVal))
        self.fmax=max(self.fmax,np.max(fVal))

        if self.linY: self.ax1.plot(xVal,fVal,label=lab,color=colr)
        else:         self.ax1.semilogy(xVal,fVal,label=lab,color=colr)

    def reference(self,xVal,fVal,lab,colr):
        self.xref=xVal
        self.fref=fVal

        fmax=max(self.fmax,np.max(fVal))
        self.ax1.set_xlim([np.min(xVal),np.max(xVal)])
        self.ax2.set_xlim([np.min(xVal),np.max(xVal)])
        self.xmin,self.xmax=getRange(self.flags,[np.min(xVal),np.max(xVal)],"-xrange",self.errRange)

        self.addPlot(xVal,fVal,lab,'black')

    def compare(self,xVal,fVal,lab,colr):
        f1=self.fref
        fx=griddata(xVal,fVal,(self.xref,),method='linear')
        #plt.suptitle(dir+"*"+which,fontsize=20)
        self.addPlot(self.xref,fx,"",colr)

        self.fmin=min(self.fmin,np.min(fx))
        self.fmax=max(self.fmax,np.max(fx))
        if np.max(abs(fx-self.fref)) !=0:
            ferr=[]
            xerr=[]
            for n in range(len(fx)):
                if (fx[n]!=0 or f1[n]!=0) and fx[n]-self.fref[n]!=0:
                    ferr.append(abs(fx[n]-self.fref[n])/abs(0.5*(fx[n]+f1[n])))
                    xerr.append(self.xref[n])

            self.ax2.semilogy(xerr,ferr,label=lab,color=colr)
            self.errmin=min(self.errmin,np.min(ferr))
            self.errmax=max(self.errmax,np.max(ferr))

            # rms deviations
            relErr=(fx-self.fref)/(self.fref+self.tolerance*np.max(self.fref))
            errRMS=np.sqrt( np.sum(relErr*relErr) * (xVal[1]-xVal[0]) / (np.max(xVal)-np.min(xVal)) )
            print("errRMS:",errRMS,"\t",self.tolerance,"\t",np.max(self.fref))
        else: print("identical")

    def plot(self):
        # remove ticks in pupper panel
        self.fax.subplots_adjust(hspace=0.1)
        plt.setp([a.get_xticklabels() for a in self.fax.axes[:-1]], visible=False)

        errmin=self.errmin
        errmax=self.errmax
        print("Error range: ",errmin,errmax,self.errRange,"set range by -erange=[eMin,eMax]")
        errmin=max(errmax*self.errRange,errmin)
        errmin,errmax=getRange(self.flags,[errmin,errmax],"-erange",self.errRange)
        self.ax2.set_ylim([errmin,errmax])

        fmin=self.fmin
        fmax=self.fmax
        fmin,fmax=getRange(self.flags,[fmin,fmax],"-vrange")
        self.ax1.set_ylim([fmin,fmax])

        #finalize
        self.ax1.legend().draw_frame(False)
        if self.ax2.legend()!=None:
            self.ax2.legend().draw_frame(False)



