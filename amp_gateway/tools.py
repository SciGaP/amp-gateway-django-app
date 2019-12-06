#!/usr/bin/env python

# The tRecX package is free for personal use.
# Any commercial use of the code or parts of it is excluded.
# Restrictions for academic use apply. 
# 
# See terms of use in the LICENSE file included with the source distribution
# 
# Copyright (c) 2015,2016,2017,2018 by Armin Scrinzi
# End of license
 
# -*- coding: utf-8 -*-

import sys
import numpy as np
import os.path

def argsAndFlags():
    """
    separete sys.argv into flags (argumentes starting with "-") and the rest
    """
    args=[]
    flags=[]
    for item in sys.argv[1:]:
        if item[0]=="-": flags.append(item)
        else:            args.append(item)
    return args,flags

def flagValue(flags,command,default=None):

    for flag in flags:
        if flag.find(command)==0 and flag.find("=")!=-1:
            return flag.split("=")[1]
    if default==None:
        print("specify value for flag ",command)
        sys.exit(1)

    return default


def prePost(flags):
    """ resolve flags for pre-and post-fixes by -dir=PreFix and -which=PostFix """
    for flag in flags:
        if flag.find("-dir")==0:
            print("obsolete flag",flag," use \"before\" flag -b=... instead")
        if flag.find("-dir")==0:
            print("obsolete flag",flag," use \"after\" flag -a=... instead")
    dir=""
    which=""
    for flag in flags:
        if flag.find("-b=")==0: dir=flag[3:]
        if flag.find("-a=")==0: which=flag[3:]
    return dir,which

def getch():
    """ get single character from terminal """
    import termios
    import sys, tty
    def _getch():
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    return _getch

def floatRange(Range):
    if Range.find("[")==-1 or Range.find("]")==-1:
        print("not a valid float range: ",Range,"specify as in, e.g. [1.e-4,100]")
        sys.exit(1)

    rang=Range[Range.find("[")+1:Range.find("]")].split(",")
    if len(rang)!=2:
        print("not a valid float range: ",Range,"specify as in, e.g. [1.e-4,100]")
        sys.exit(1)
    return float(rang[0]),float(rang[1])


def splitOutsideBrackets(str,sep,left,right):
    """ split string at sep, except where within left-right pairs of brackets"""
    if len(left)!=len(right):
        print("number of left brackets does not match number of right brackets")
        print(" left:",left)
        print("right:",right)

    # replace separators between brackets by PLACEHOLDER
    ph="PLACEHOLDER"
    ph=ph.replace(sep,"") # modify to avoid accidental agreement

    loc=str.find(sep)
    while loc!=-1:
        for k in range(len(left)):

            if str[loc:].count(left[k])!=str[loc:].count(right[k]):
                str=str[:loc]+str[loc:].replace(sep,ph)
                break;
        loc=str.find(sep,loc+1)

    # split and re-insert separators
    sp=str.split(sep)
    for k in range(len(sp)): sp[k]=sp[k].replace(ph,sep)

    return sp

def rangeInSquareBrackets(Range):
    beg=0
    end=-1
    if Range.find('[')==-1: return beg,end

    if Range.find('[')==-1 or Range.find(']')==-1 or Range.find(':')==-1:
        print("not a valid range: ",Range)
        sys.exit(1)

    r=Range[Range.rfind('[')+1:].strip()
    r=r[:r.rfind(']')].strip()
    if r[0:r.find(':')]!="": beg=int(r[:r.find(':')])
    if r[r.find(':')+1:]!="":end=int(r[r.find(':')+1:])

    return beg,end

def expandRange(colRange):
    """
    expand a range string
    example: "1,2,6,8-10,3" -> [1,2,6,8,9,10,3]
    """
    jobs=[]
    subr=colRange.split(",")
    for s in subr:
        rang=s.split("-")
        if len(rang)==1: jobs.append(int(rang[0]))
        else:
            end=10000
            if rang[1]!="": end=int(int(rang[1])+1)
            jobs=jobs+ list(range(int(rang[0]),end))
    return jobs

def findFirstLarger(Range,Val):
    for i in range(len(Range)):
        if Range[i]>Val: return i
    return len(Range)

class DataFile:
    """
    structured access to a data file
    specify files and columns
    examples:
        myDataFile... read all columns, plot column 0 against all other, or, if 2d cols 0,1 against the rest
        myDataFile[3,1:7,8,4]...cols 3,1,7,8,4 from myDataFile, interprete 3,1 as x- and y-axes of a 2d plots 7,8,4
                                note: file[7][2,3:4] ... last square bracket is column info
    """
    def __init__(self,FileCols):
        self.head=None
        self.name=FileCols
        if self.name.find("[")!=-1: self.name=FileCols[:FileCols.rfind("[")]
        openFile = open(self.name)
        self.kind="gnuplot"
        self.pars=RunParameters(FileCols[:FileCols.rfind("/")])

        # collect the header lines
        line=str(openFile.readline())
        self.head=[]
        while line[0]=="#":
            self.head.append(line)
            line=str(openFile.readline())

        # find how columns are separated: ',' or whitespace
        if self.firstDataLine(openFile).find(',')!=-1: self.colSep=","
        else:                                          self.colSep=""

        # detect multiple datasets separated by blank lines
        self.isMulti=False
        linePrev=self.firstDataLine(openFile)
        while linePrev!="":
            line=str(openFile.readline())
            if linePrev.strip()=="" and line.strip()!="":
                self.isMulti=True
                break
            linePrev=line

        # determine columns to be read
        line=str(self.firstDataLine(openFile))
        self.axCols=[0,]
        if self.isMulti: self.axCols=[0,1]
        if FileCols.find('[')==-1:
            # nothing selected - read all columns
            self.datCols=list(range(len(self.axCols),len(self.rowSplit(line))))
        else:
            # columns given --- expand ranges and get axes (if given)
            colstr=FileCols[FileCols.rfind('[')+1:FileCols.rfind(']')]
            self.datCols=expandRange(colstr.split(':')[-1])
            if 0 in self.datCols or self.isMulti and 1 in self.datCols:
                sys.exit("data columns start at 1 for 1d and at 2 for 2d")
            if colstr.find(':')!=-1:
                self.axCols=expandRange(colstr.split(':')[0])

        # extract colmn names; if no headers, just number
        colNames=[]
        if len(self.head)>0:
            colNames=self.head[-1][1:].split("=")[-1].split()

        # axis names
        self.axisName=["x","y","z"]
        nAxNam=0
        if len(self.head)>0:
            self.axisName=self.head[-1][1:].split(":")[0].split(",")

            # if "=", first 1 or 2 columns are for axes
            if self.head[-1].find("=")!=-1: nAxNam=len(self.axCols)
            else:                           nAxNam=0

        # number of column headers does not match number of columns - name by column number
        if(len(colNames)+nAxNam)!=len(self.rowSplit(self.firstDataLine(openFile))):
            colNames=list(range(len(self.rowSplit(line))))

        # strip brackets from column names
        for i in range(len(colNames)):
            n=str(colNames[i]).strip()
            if n.find("(") ==0:        n=n[1:]
            if n.rfind(")")==len(n)-1: n=n[:-1]
            colNames[i]=n

        # set up a dictionary for axis and data columns
        self.cols={}
        for c in self.axCols:  self.cols[c]=["ax"+str(c),[]]
        for c in self.datCols: self.cols[c]=[colNames[c-nAxNam],[]]

        # get axes and data columns
        self.datCols.extend(self.axCols)
        line=self.firstDataLine(openFile)
        while line!="":
            if line.strip()!="":
                data=self.rowSplit(line)
                for c in self.datCols:
                    self.cols[c][1].append(float(data[c]))
            line=openFile.readline()
        self.datCols=self.datCols[:-len(self.axCols)]

        # done, close file
        openFile.close()

    def xAxis(self):
        return np.array(self.cols[self.axCols[0]][1])

    def yAxis(self):
        return np.array(self.cols[self.axCols[1]][1])

    def column(self,col):
        return np.array(self.cols[col][1])

    def str(self):
        s=self.name+"\n"
        for key in self.cols:
            s+=str(key)+":"+str(self.cols[key][0])+"["+str(len(self.cols[key][1]))+"] "
        return s

    def colName(self,col):
        if not col in self.cols:
            print("Column",col,"not read",self.cols)
        return self.cols[col][0]


    def rowSplit(self,line):
        line=str(line)
        if self.colSep=="": return line.strip().split()
        else: return line.strip().split(self.colSep)

    def firstDataLine(self,file):
        file.seek(0)
        for n in range(len(self.head)): line=file.readline()
        return file.readline()

    def ranges(self,col):
        args,flags=argsAndFlags()

        ratio=1.e-5
        if "-linY" in flags: ratio=None
        vMin,vMax=getRange(flags,self.column(col),"-vrange",ratio)
        xMin,xMax=getRange(flags,self.xAxis(),"-xrange")
        if self.isMulti:
           yMin,yMax=getRange(flags,self.yAxis(),"-yrange")
        else:
           yMin,yMax=vMin,vMax
        return vMin,vMax,xMin,xMax,yMin,yMax




class RunParameters:
    """
    access to input parameters from a given tRecX "linp" file
    (the linp-file echos all inputs as actually read)
    """

    def __init__(self,RunDir):
        self.run=RunDir
        if os.path.exists(RunDir+'/linp-extract'):
            self.linp=open(RunDir+'/linp-extract','r')
        elif os.path.exists(RunDir+'/linp'):
            self.linp=open(RunDir+'/linp','r')
        else:
            self.linp=None;

    def hasParameters(self):
        return self.linp!=None

    def item(self,Name,I):
        """ string value of input item """
        if not self.hasParameters() or len(self.allItems(Name)[0])==0: return None
        return self.allItems(Name)[0][I]

    def short(self,Name,I):
        """ short name of input item """
        return self.allItems(Name)[1][I]


    def floatItem(self,Name,I):
        """ float value of input item """
        return float(self.allItems(Name)[0][I])

    def allItems(self,Name):
        """ all input items matching Name """
        if self.linp==None:
            print("no linp-file found, cannot display paramters")
            return
        self.linp.seek(0) # position to beginning of file
        items=[]
        short=[]
        line=str(self.linp.readline())
        while line!="":
            if line.find(Name)!=-1:
                val = line.split("=")
                if len(val)>1: items.append(val[1].rstrip()) # remove possible trailing whitspace
                else:          items.append("")
                if len(val)>2: short.append(val[2].rstrip()) # remove possible trailing whitspace
                else:          short.append(Name)
            line=str(self.linp.readline())
        return items,short


