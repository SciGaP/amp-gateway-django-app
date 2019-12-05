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
from math import *

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm

from .tools import *

def saveAndShow(savename=None):
    # get python file name (as run from the command line)
    if savename!=None:
        name=savename
    else: 
        name=sys.argv[0].split(".py")[0]+".png"
    plt.savefig(name,bbox_inches='tight')
    plt.show(block=False)
    answer = input('"u" to update ')
    plt.close()


class Plot:
    def lineoutPosWidth(self,Flag):
        if Flag.find(",")!=-1:
            pos=float(Flag[10:Flag.find(",")])
            wid=float(Flag[Flag.find(",")+1:])
        else:
            pos=float(Flag[10:])
            wid=0.
        return pos,wid

    def setLineout(self,Flags):
        self.atX=None # lineout at X-coordinate
        self.atY=None # lineout at Y-coordinate
        for f in Flags:
            if f.find("-lineoutX=")!=-1: self.atX,self.wid=self.lineoutPosWidth(f)
            if f.find("-lineoutY=")!=-1: self.atY,self.wid=self.lineoutPosWidth(f)
        
    """ collect parameters and modify data for figure """
    def __init__(self,Flags):
        self.setLineout(Flags)
        self.dim=[[1.e10,-1.e10],[1.e10,-1.e10],[1.e10,-1.e10],]

    def layout(self,Nsub):
        """ layout for subplots in figure (rows/cols"""
        if Nsub<3:   nr=1
        elif Nsub<9: nr=2
        else:        nr=3
        nc=(Nsub-1)/nr+1

        return nr,nc

    def extend(self,fmin,fmax,xmin,xmax,ymin=None,ymax=None):
        """ extend the plot size to present """
        def extend1(xmin,xmax,dim):
            dim[0]=min(dim[0],xmin)
            dim[1]=max(dim[1],xmax)
            return dim

        self.dim[0]=extend1(xmin,xmax,self.dim[0])
        self.dim[2]=extend1(fmin,fmax,self.dim[2])
        if ymin!=None: self.dim[1]=extend1(ymin,ymax,self.dim[1])

    def isLineout(self):
        return self.atX!=None or self.atY!=None


    def lineout(Plot,Xcol,Ycol,Zcol,Flags=None):
        """ exctract line-out according to Plot """
        # axes and matrix of values
        x,y,v=xyzToAxisArray(Xcol,Ycol,Zcol)

        if Flags!=None: Plot.setLineout(Flags)

        # extract row/col nearest to atX/atY
        if Plot.atX!=None:
            colMin = (np.abs(x-Plot.atX+Plot.wid/2)).argmin()
            colMax = (np.abs(x-Plot.atX-Plot.wid/2)).argmin()+1
            return x[colMin],y,np.divide(np.sum(v[:,colMin:colMax],axis=1),colMax-colMin)
        elif Plot.atY!=None:
            rowMin = (np.abs(y-Plot.atY+Plot.wid/2)).argmin()
            rowMax = (np.abs(y-Plot.atY-Plot.wid/2)).argmin()+1
            print("sum",rowMin,rowMax,Plot.atY,y[0],y[-1])
            return y[rowMin],x,np.divide(np.sum(v[rowMin:rowMax,:],axis=0),rowMax-rowMin)
        else:
            print("specify lineout point -lineoutX=XVAL or -lineoutY=YVAL")
            exit(1)

def energyPeakPositions(E0,Omega,Emin,Emax):
    ePos =[]
    nPhot=[]
    epos=E0-Omega
    nphot=-1
    while epos<Emax-Omega:
        epos=epos+Omega
        nphot+=1
        if epos>Emin:
            ePos.append(epos)
            nPhot.append(nphot)
    return ePos,nPhot


def lineoutSumEnergy(SumEnergy,Width,xGrid,yGrid,vMatr,Flags=None):
    """ create a lineout si(xGrid) along xGrid+yGrid=SumEnergy"""
    res=[]
    dGrid=[]
    for xE in np.arange(np.min(xGrid),np.max(xGrid)+0.1,0.1):
        yE=SumEnergy-xE
        if Width<xE and xE<SumEnergy-Width and yE>0:
            dGrid.append(2*xE-SumEnergy)
            """ average over square surrounding (xE,yE)"""
            iLow=findFirstLarger(xGrid,xE-Width/2)
            iUpp=findFirstLarger(xGrid,xE+Width/2)
            jLow=findFirstLarger(yGrid,yE-Width/2)
            jUpp=findFirstLarger(yGrid,yE+Width/2)
            cnt=0
            sum=0
            for i in range(iLow,iUpp):
                for j in range(jLow,jUpp):
                    cnt+=1
                    sum+=vMatr[i,j]
            if cnt>0: res.append(sum/cnt)
            else:     res.append(0.)
    return dGrid,res

def lineoutAllSumEnergy(Datfil,xGrid,yGrid,vMatr,Ip):
    pars=RunParameters(Datfil.name[:Datfil.name.rfind('/')])
    omega,up,ip=pulseParameters(pars,Ip)
    esum,nphot=energyPeakPositions(-ip-2*up,omega,0,np.max(xGrid))
    for k in range(len(esum)):
        dgrid,espec=lineoutSumEnergy(esum[k],omega*0.5,xGrid,yGrid,vMatr,Ip)
        if len(dgrid)>0:
            f=open(Datfil.name+"_linoutN"+str(nphot[k]),'w')
            f.write("# lineout at sumEnergy="+str(esum[k])+", witdh="+str(omega*0.5)+"\n")
            for l in range(len(dgrid)):
                f.write(str(dgrid[l])+", "+str(espec[l])+"\n")
            f.close()

            tspec=np.fft.rfft(espec,norm="ortho")
            f=open(Datfil.name+"_linoutT"+str(nphot[k]),'w')
            f.write("# Fourier transform at sumEnergy="+str(esum[k])+", witdh="+str(omega*0.5)+"\n")
            dt=0.5*omega/(np.max(dgrid)-np.min(dgrid))
            for l in range(len(tspec)):
                f.write(str(l*dt)+", "+str(abs(tspec[l]))+"\n")
            f.close()



def normalize(flags,x,y,atX=None):
    for f in flags:
        if f.find("-normalize")!=-1:
           try:
               x0=float(f.split("=")[1])
               print("normalized at x=",x[np.abs(x - x0).argmin()],"is",y[np.abs(x - x0).argmin()])
               y/=y[np.abs(x - x0).argmin()]
           except:
               print("normalized at maximum",np.max(np.abs(y)))
               y/=np.max(np.abs(y))
           return
    if atX!=None:
        print("normalized at x=",x[np.abs(x - atX).argmin()],"is",y[np.abs(x - atX).argmin()])
        y/=y[np.abs(x-atX).argmin()]


# average around values around x-points
def average(flags,x,y,av=None):

    for f in flags:
        if f.find("-average=")!=-1:
            av=float(f.split("=")[1])/2.
    if av==None or av==0.: return x,y

    ya=np.zeros((len(x)))
    for k in range(len(x)):
        na=0
        for l in range(len(x)):
            if x[k]-av<x[l] and x[l]<x[k]+av:
                na+=1
                ya[k]+=y[l]
        ya[k]/=na
    return x,ya


def getRange(flags,data,command,minRatio=None):

    for flag in flags:
        if flag.find(command)==0:
            return floatRange(flag)

    up= np.max(data)
    if minRatio!=None:
        low=up*minRatio
    else:
        low=np.min(data)
    return low,up

def adjustLogRange(vMin,vMax):
    if vMin<=0:
        vMin=vMax*1e-5
        print("lowest value is <= 0, using default vmin="+str(vMin)+", may specify -vrange=[vmin,vmax]")
    return vMin


def plotLinOrLog(axn,xGrid,yGrid,vMatr,flags):
    cmap = plt.cm.get_cmap("gnuplot")
    vMin,vMax=getRange(flags,vMatr,"-vrange=")
    yMin,yMax=getRange(flags,yGrid,"-yrange=")
    xMin,xMax=getRange(flags,xGrid,"-xrange=")

    if "-linV" in flags:
        lev,tic=linLevels(vMin,vMax)
        CS=axn.contourf(xGrid,yGrid,vMatr,lev,cmap=cmap,vmin=vMin,vmax=vMax)
    else:
        vMin=adjustLogRange(vMin,vMax)
        lev,tic=logLevels(vMin,vMax)
        CS=axn.contourf(xGrid,yGrid,vMatr,lev,cmap=cmap,vmin=vMin,vmax=vMax,locator=ticker.LogLocator())
    axn.set_xlim(xMin,xMax)
    axn.set_ylim(yMin,yMax)
    return CS,tic

def xyzToAxisArray(Xcol,Ycol,Zcol):
    """ convert from 3-column to 2dim data """
    fastFirst=Xcol[0]!=Xcol[1];

    if not fastFirst:
        for k in range(len(Xcol)):
            if(Xcol[0]!=Xcol[k]):
                yAxis=Ycol[:k]
                break
        xAxis=np.array([Xcol[k] for  k in range(0,len(Xcol),len(yAxis))])

        vArray=Zcol.reshape(len(xAxis),len(yAxis)).transpose() # numpy is row-wise
    else:
        for k in range(len(Ycol)):
            if(Ycol[0]!=Ycol[k]):
                xAxis=Xcol[:k]
                break
        yAxis=np.array([Ycol[k] for  k in range(0,len(Ycol),len(xAxis))])
        vArray=Zcol.reshape(len(yAxis),len(xAxis)) # numpy is row-wise

    return xAxis,yAxis,vArray

def momentumToEnergyEv(kX,kY,sigmaK):
    au2eV=27.211386
    for i in range(len(kX)):
        for j in range(len(kY)):
            sigmaK[i,j]*=kX[i]*kY[j]/(au2eV*au2eV)
    for k in range(len(kX)): kX[k]*=kX[k]*0.5*au2eV
    for k in range(len(kY)): kY[k]*=kY[k]*0.5*au2eV

def pulseParameters(inputPars,ionizationPot):
    # check wave-length
    lam=inputPars.allItems("lambda(nm)")
    for l in lam:
        if l!=lam[0]: print("multiple wave-length, using first: ",lam)
    wavelength=float(inputPars.item("lambda(nm)",0))

    inte=inputPars.allItems("I(W/cm2)")
    intensity=float(inte[0][0])

    omega=45.5633/wavelength
    up=intensity/(4*omega*omega)/3.50944e16
    ip=ionizationPot

    au2eV=27.211386
    omega*=au2eV
    up*=au2eV
    ip*=au2eV

    return omega,up,ip

def twoElectronEnergies(eGrid,inputPars,axn,flags,ionizationPot,photonNumber=None,colr=None):
    """
    compute peak postitions
    from wavelength, Ip, and intensity
    add lines to plot
    """

    print("Ip",ionizationPot)
    omega,up,ip=pulseParameters(inputPars,ionizationPot)
    title="om="+str(omega)+", Up="+str(up)+", ip="+str(ip)

    # draw lines at n omega -ip-2*up
    emax=np.max(eGrid)
    emin=np.min(eGrid)
    epos,nphot=energyPeakPositions(-ip-2*up,omega,emin,emax)
    actualPlot=-1
    for k in range(len(epos)):
        if photonNumber==None or nphot[k] in photonNumber:
            c="b"
            w=1
            actualPlot+=1
            if colr!=None:
                c=colr[actualPlot%len(colr)]
                w=3
            axn.plot([0,epos[k]],[epos[k],0],'--',color=c,linewidth=w)
            bbox_props = dict(boxstyle="square,pad=0.", fc="white", ec="b", lw=0, alpha=0.8)
            axn.text(epos[k]/2,epos[k]/2-10,str(nphot[k]),color=c, fontsize=14, fontweight='bold',bbox=bbox_props)

    # draw diagonals
    ediag=-omega
    while ediag<emax-omega:
        ediag+=omega
        if photonNumber==None: axn.plot([ediag,emax],[0,emax-ediag],'--',color='y')

    return title


def peakPositions(pltsize,inputPars,ax1,flags):
    """
    compute peak postitions
    from wavelength, Ip, and intensity
    add lines to plot
    """
    title=inputPars.run

    # get wave-length, intensity, ip
    name=""
    # check wave-length
    lam=inputPars.allItems("lambda(nm)")
    for l in lam:
        if l!=lam[0]: print("multiple wave-length, using first: ",lam)
    wavelength=float(inputPars.item("lambda(nm)",0))
    # add up all intensities
    inte=inputPars.allItems("I(W/cm2)")
    intensity=0
#    for i in inte: intensity+=float(i)
    intensity=float(inte[0][0])

    omega=45.5633/wavelength
    up=intensity/(4*omega*omega)/3.50944e16
    ip=0.903

    if "-eV" in flags:
        up*=27.211
        ip*=27.211
        omega*=27.211

    # draw lines at n omega - ip - up
    epos=-ip
    title="om="+str(omega)+", Up="+str(up)+", ip="+str(ip)
    epos=epos-up
    title=title+" Up subtracted"
    plt.title(title)

    nphot=0
    emax=pltsize.dim[0][1]*pltsize.dim[0][1]*0.5
    while epos<min(emax-omega,pltsize.dim[0][1]-omega):
        epos=epos+omega
        nphot=nphot+1
        if epos>0:
            kpos=epos
            if not "-eV" in flags: kpos=sqrt(2*epos)
            ax1.plot([kpos,kpos],pltsize.dim[2],color='b')
            ax1.text(kpos,pltsize.dim[2][1]*0.5,str(nphot),rotation=90)

    return name

def linLevels(vmin,vmax):
    """
    a set of linearly spaced levels in [vmin,vmax]
    ticks are at 16 intervals
    """
    lmin=int(log10(abs(vmin))+1)
    lmax=int(log10(abs(vmax))+1)
    if vmin<0: lmin=-lmin
    if vmax<0: lmax=-lmax

    tics=[]
    fact=(vmax-vmin)/16
    levs=[vmin+fact,]
    for n in range(1,17): levs.append(levs[n-1]+fact)
    for l in range(lmin,lmax): tics.append(pow(10,l))
    return levs,tics

def logLevels(vmin,vmax):
    """
    a set of logarithmically spaced levels in [vmin,vmax]
    ticks are returned at powers of 10
    """
    lmin=int(log10(vmin)+1)
    lmax=int(log10(vmax)+1)

    tics=[]
    fact=pow(vmax/vmin,float(1./16.))
    levs=[vmin*fact,]
    for n in range(1,17): levs.append(levs[n-1]*fact)
    for l in range(lmin,lmax): tics.append(pow(10,l))
    return levs,tics

def minGraph(x,y):
    my=-np.array(y)
    return maxGraph(x,my)


def maxGraph(x,y):
    """
    locate the maximum in a graph, by quadratic interpolation
    """
    kmax=np.argmax(y)
    if kmax==0 or kmax==len(y):
        print("WARNING: no local maximum")
        return x[kmax]
    
    xMat=np.zeros((3,3))
    for i in range(3):
        xij=1
        for j in range(3):
            xMat[i,j]=xij
            xij*=x[kmax-1+i]

    c=np.linalg.solve(xMat,y[kmax-1:kmax+2])
    
    xmax=-c[1]/(2*c[2])
    if x[kmax-1]>xmax or x[kmax+1]<xmax:
        print("parabolic location of maximum failed: ",xm,"points",x[kmax-1:kmax+2])
    
    return xmax



class Legend:
    """
    for creating labels and formating the legend
    - use "-label=PARNAME" on the command line to extract values from input (linp-file) to legend
    """
    def __init__(self,Dir,Postfix,Flags,Pars=None):
        self.previous=""
        self.previousDir=""
        self.dir=Dir
        self.postfix=Postfix
        self.labShort=""
        self.labName=""
        for f in Flags:
            if f.find("-label=")!=-1:
                self.labName=f[7:].strip()
                if self.labName[ 0]=="'" or self.labName[ 0]=='"':self.labName=self.labName[1:]
                if self.labName[-1]=="'" or self.labName[-1]=='"':self.labName=self.labName[:-1]
                if Pars!=None:
                    for n in self.labName.split(","):
                        self.labShort+=","+Pars.short(n,0)
                if(len(self.labShort)!=0): self.labShort=self.labShort[1:]

    def label(self,ax1,col,datfil,info,pars):
        """
        create a new label for a plot, add header info as appropriate
        - if there is a dir or postfix add legend header
        - if there is a new file or info and multiple columns, add legend header
        - put into legend label the info that is not in header
        """

        file=self.dir+"/"+info+"/"+self.postfix

        parVal=""
        # if input parameters are given, add values to label
        if self.labName!="":
            for l in self.labName.split(","):
                parVal+=pars.item(l,0)+","
        parVal=parVal[:-1]

        # append file info or column number to label
        if len(datfil.datCols)==1:
            # if single column, just use file name
            lab=parVal+" ["
            if self.dir+self.postfix=="": lab+=file
            else                        : lab+=info
            lab+="]"
        else:
            # add column names
            lab=str(datfil.cols[col][0])

        # directory or post-fix are given, write into header
        curDir=self.dir+"*"+self.postfix
        if self.dir+self.postfix!="" and curDir!=self.previousDir:
            labHead=""
            if self.labShort!="": labHead=self.labShort+": "+curDir
            else:                 labHead=curDir
            if labHead!="": ax1.plot([1,],[1,],' ',label=labHead)
            self.previousDir=curDir

        # if file name changed and not in label, force into legend
        if self.previous!=info and lab.find(file)==-1:
            self.previous=info
            ax1.plot([1,],[1,],' ',label=parVal+" ["+info.split('[')[0]+"]",alpha=1)

        return lab

