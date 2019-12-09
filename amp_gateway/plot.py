#! /usr/bin/env python

# The tRecX package is free for personal use.
# Any commercial use of the code or parts of it is excluded.
# Restrictions for academic use apply. 
# 
# See terms of use in the LICENSE file included with the source distribution
# 
# Copyright (c) 2015,2016,2017,2018 by Armin Scrinzi
# End of license
 
"""
plot one or several data-files with column-wise storage

columns can be specified in square brackets after the file name
        in the form fileName[2,3,7-9] for plotting columns 2,3,7,8,9
        without square bracket, all columns are plotted
column 0 of each file is interpreted as x-axis by default
       for different x-axis use [3:7,8], which plot columns 7,8 against 3
column headers will be recognized and used as legend labels
"""

import io
import re
import sys

from matplotlib.figure import Figure
import matplotlib.cm as cm
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
from .plot1d import *
from .plot2d import *

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def addPolar(CosX,Y):
    return 1



def findFirstOf(string,chars):
    """ find first position in string of any from as string of characters"""
    pos=100000
    for k in range(len(chars)):
        if string.find(chars[k])!=-1: pos=min(pos,string.find(chars[k]))
    return pos

help=[
"-b=DIR......... prefix input file by DIR"
,"-a=WHICH..... postfix input file by WHICH"
,"-label=PAR1,PAR2... parameter value(s) for PAR1,PAR2,... will be used in plot legend"
,"-linV.... plot on linear value-axis (default is log)"
,"-linY.... plot on linear y-axis (default is log)"
,"-logY.... plot on logarithmic y-axis (default)"
,"-peaks... mark photon peaks at n*omega-Ip-Up or (-2 Up for Helium)"
,"-eV...... x -> x^2/2 *Rydberg"
,"-vrange[vmin,vmax]... function value range for ALL plots "
,"-xrange[xmin,xmax]... plot axis range for All plots "
,"-yrange[ymin,ymax]... plot axis range for All plots "
,"-erange[emin,emax]... range for plotting errors "
,"-polar            ... polar plot of 1d or 2d data"
,"-equalAx          ... 2d with commensurate axes"
,"-rmax=r           ... radius in polar plot"
,"-lineoutX=x[,w]   ... lineout of 2d plot nearest to x or sum [x-w/2,x+w/2], similar for Y"
,"-normalize[=x0]   ... scale to maximal value = 1 (nearest to x0)"
,"-symX12           ... symmetrize 2d plots by 1 <--> 2 (special for He)"
,"-compare          ... compare multiple 2d plots to first"
,"-maxgraph         ... print location of graph's maximum"
,"-mingraph         ... print location of graph's minimum"
,"-sum              ... add files and columns single plot"
]
# (files,flags) = argsAndFlags();

#  check flags
# fail=""
# for f in flags:
#     f0=f
#     if f0.find("=")!=-1:f0=f0[:f0.find("=")]
#     for h in help:
#         h0=h
#         ff=findFirstOf(h," [=.,")
#         if ff!=-1: h0=h[:ff]
#         if f0==h0: break
#     else: fail=f
# if fail!="":
#     print("undocumented flag: "+fail)
#     print("allowed flags:")
#     for h in help: print("  ",h)
#     sys.exit(0)


# if len(files)<1:
#     print("Usage: ")
#     print("   plot.py file1[{xcol:}colrange] file2[{xcol:}colrange] ... {flags}")
#     print("Example:")
#     print("   plot -b=Argon/ -a=/spec[3:22,26,30] 0015 0014 -label='I(W/cm2),lambda(nm)'")
#     print("     will plot column 3 vs columns 22,26,30 of files Argon/0015/spec and Argon/0014/spec")
#     print("     with legend labels including intensity and wavelength")
#     print("Command line flags:")
#     for h in help: print("  ",h)
#     sys.exit(0)

class TRecXPlotViewProvider:
    display_type = 'image'
    name = "TRecX Plot"
    test_output_file = os.path.join(BASE_DIR, "data", "spec_total")
    # plot2d test case
    # test_output_file = os.path.join(BASE_DIR, "data", "plot2d", "spec_total")
    # Modl_RunID test case
    # test_output_file = os.path.join(BASE_DIR, "data", "Modl_RunID", "Modl_RunID")

    def generate_data(self, request, experiment_output, experiment, output_file=None):

        # files = [os.path.basename(output_file.name)]
        # flags = [f"-b={os.path.dirname(output_file.name)}/"]
        file_name = os.path.basename(output_file.name)
        # parse the Modl_RunID file
        if file_name == "Modl_RunID":
            model_runid = output_file.read().decode()
            m = re.match(r"(\S+) (\S+)", model_runid)
            if m is None:
                raise Exception(f"Invalid Modl_RunID file contents: {model_runid}")
            model, run_id = m.group(1, 2)
            files = [os.path.join(os.path.dirname(output_file.name), 
                     "ARCHIVE", model, run_id, "spec_total")]
        else:
            files = [output_file.name]
        flags = []
        dir,which=prePost(flags)

        # get total number of datasets to be plotted
        nSets=0
        for f in files:
            file=dir+f+which
            if file.find("[")==-1:
                datfil=DataFile(file)
                nSets+=len(datfil.datCols)
            else:
                nSets=nSets+len(expandRange(file[max(file.find('['),file.find(':'))+1:file.find(']')]) )



        theta=[] # may need this
        replot=False

        comp=None
        countAll=-1
        fax = None

        # loop through input files
        plot=Plot(flags)
        previous=""
        previousDir=""
        nPlot=-1
        #loop through files
        for info in files:
            file=dir+info+which
            pars=RunParameters(file[:file.rfind('/')])
            # get data into selected columns
            datfil=DataFile(file)

            if info==files[0]:
                # first file --- set up the figure
                if datfil.isMulti and not plot.isLineout():
                    """ two axes -- separate panels """
                    fax=Figure()
                    cmap = cm.get_cmap("gnuplot")
                    nr,nc=plot.layout(nSets)
                    gs = gridspec.GridSpec(nr,nc,height_ratios=[1]*nr)
                else:
                    """ single axis --- all in one panel """
                    gs = gridspec.GridSpec(1,1,height_ratios=[1,])
                    fax=Figure()
                    axn=fax.add_subplot(gs[0])

            if info==files[0]: leg=Legend(dir,which,flags)

            # loop through columns in file
            for col in datfil.datCols:
                nPlot+=1
                countAll+=1
                # sum all data into single plot, plot when arrived at last
                doPlot="-sum" not in flags or (info==files[-1] and col==datfil.datCols[-1])

                if datfil.isMulti and not plot.isLineout(): plot2d(fax,         datfil,flags,pars,plot,col,gs,nPlot,doPlot,info)
                else:                              axn,comp=plot1d(fax,axn,comp,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info,leg,countAll+1==len(files)*len(datfil.datCols))

        if (comp==None and (not datfil.isMulti or plot.isLineout())):
            # single plot, adjust sizes
            xMin,xMax=getRange(flags,plot.dim[0],"-xrange")
            vMin,vMax=getRange(flags,plot.dim[2],"-vrange")

            axn.set_xlim(xMin,xMax)
            axn.set_ylim(vMin,vMax)

            # axis labels
            xlab="x"
            if "-eV" in flags: xlab="eV"
            axn.set_xlabel(xlab,fontsize=14)

            # Print legend
            lg = axn.legend()

            # mark peak positions (if so desired)
            if "-peaks" in flags: peakPositions(plot,pars,axn,flags)

        # Export plot as image buffer
        buffer = io.BytesIO()
        fax.savefig(buffer, format='png')
        image_bytes = buffer.getvalue()
        buffer.close()

        # return dictionary with image data
        return {
            'image': image_bytes,
            'mime-type': 'image/png'
        }

# answer="u"
# while answer=="u":
#     comp=allPlots()
#     plt.show(block=False)
#     answer = input('"u" to update ')

#     # get python file name (as run from the command line)
#     name=sys.argv[0].split(".py")[0]
#     if "-compare" in flags: comp.fax.savefig(name+".png")
#     else: fax.savefig(name+".png")
#     plt.close()
#     replot=True

