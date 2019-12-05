import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.interpolate
from scipy.interpolate import griddata
import copy  as cp
from math import *

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm

from .tools import *
from .plotTools import *
from .compareToFirst import *


def plot2d(fax,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info):
    file=datfil.name
    xVal=datfil.xAxis()
    yVal=datfil.yAxis()
    fVal=datfil.column(col)

    xGrid,yGrid,vMatr=xyzToAxisArray(xVal,yVal,fVal)
    if "-sum" in flags:
        if info==files[0] and col==datfil.datCols[0]:
            xSum=xGrid
            ySum=yGrid
            fSum=vMatr
        else:
            fMatr+=fSum
            if (xSum!=xGrid).any() or (ySum!=yGrid).any():
                exit("axes do no match")

    # transformations:
    # energy-axis in eV
    if "-symX12" in flags:
        tmp=vMatr
        vMatr=tmp+tmp.transpose()
        print("symmetrized 1<-->2")
    if "-eV" in flags: momentumToEnergyEv(xGrid,yGrid,vMatr)
    if "-normalize" in flags: vMatr=vMatr/np.max(abs(vMatr))


    ratio=1.e-9
    if "-linY" in flags: ratio=None
    vMin,vMax=getRange(flags,fVal,"-vrange",ratio)
    xMin,xMax=getRange(flags,xGrid,"-xrange")
    yMin,yMax=getRange(flags,yGrid,"-yrange")

    # extend the plot sizes to present
    plot.extend(np.min(fVal),np.max(fVal),np.min(xGrid),np.max(xGrid),np.min(yGrid),np.max(yGrid))

    if "-linY" in flags or "-linV" in flags: lev,tic=linLevels(vMin,vMax)
    else:                                    lev,tic=logLevels(vMin,vMax)


    if "-compare" in flags:
        if countAll==0:
            refValues=vMatr
            refX=xGrid
            refY=yGrid
        else:
            if len(refX)!=len(xGrid) or len(refY)!=len(yGrid):
                exit("plot dimension differ - cannot compare")

            for i in range(vMatr.shape[0]):
               for j in range(vMatr.shape[1]):
                   if refValues[i,j]!=0.: vMatr[i,j]=abs(vMatr[i,j]-refValues[i,j])/(abs(refValues[i,j])+abs(vMatr[i,j]))*2.
                   else: vMatr[i,j]=max(vMin*10.,vMax*1.e-12)

            vMin,vMax=getRange(flags,vMatr,"-vrange",ratio)
            vMin=max(vMax*1.e-5,vMin)
            print("min/max",np.min(vMatr),np.max(vMatr))
    else:
        vMin,vMax=getRange(flags,plot.dim[2],"-vrange=")

    if "-polar" in flags:
        axn=fax.add_subplot(gs[nPlot],projection='polar')

        # guessing coordinate meaning from file name
        if file.find("/spec")==len(file)-5:
            sys.exit("need \"Eta\" or \"Phi\" in file name for -polar, file-name is: "+file)
        theta=yGrid
        if file.find("Phi")!=-1 or (len(datfil.axisName)==2 and datfil.axisName[1].find("Eta")!=-1):
            theta=np.arccos(yGrid)
            for j in range(len(yGrid)): vMatr[j,:]*=sqrt(1-min(yGrid[j]*yGrid[j],1)) # correct weight for theta

        if "-normalize" in flags: fVal/=np.max(fVal)
        vMin,vMax=getRange(flags,vMatr,"-vrange",ratio)
        if "-logY" in flags or not "-linY" in flags:
            vMin=adjustLogRange(vMin,vMax)
            vMatr+=vMin


        CS=axn.contourf(theta,xGrid,np.transpose(vMatr,[1,0]),lev,cmap=cmap,vmin=vMin,vmax=vMax,locator=ticker.LogLocator())
        axn.set_rmax(float(flagValue(flags,"-rmax",xMax)))
    else:
        axn=fax.add_subplot(gs[nPlot])
        if datfil.name.find("kXkY")==len(file)-4 or "-equalAx" in flags: axn.set_aspect('equal', 'box')
        if doPlot: CS,tic=plotLinOrLog(axn,xGrid,yGrid,vMatr,flags)

    title=info
    if file!=info: title=dir+"*"+which
    if pars.item("I(W/cm2)",0)!=None:
        title+="\n"+pars.item("I(W/cm2)",0)+", "+pars.item("FWHM",0)+", CEO="+pars.item("phiCEO",0)
        axn.set_title(title)

    if doPlot:
       cbar=fax.colorbar(CS,ticks=tic,ax=axn, shrink=0.7)
       cbar.ax.set_ylabel('yield')

    if "-sum" in flags and doPlot:
        print("sum on file",file+"_sum")
        f=open(file+"_sum",'w')
        f.write("#  \n")
        for k in range(len(xVal)):
            if k>0 and xVal[k]!=xVal[k-1]: f.write("\n")
            f.write(str(xVal[k])+", "+str(yVal[k])+", "+str(fVal[k])+"\n")
        f.close()
    if "-peaks" in flags:
        HeIp=2.88738
        if doPlot: twoElectronEnergies(xGrid,pars,axn,flags,HeIp)
        lineoutAllSumEnergy(datfil,xGrid,yGrid,vMatr,HeIp)
