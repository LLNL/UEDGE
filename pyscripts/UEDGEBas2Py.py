#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:06:22 2020
@author: jguterl
"""
import os
from uedge import *
import re
class UEDGEBas2Py(object):
    def __init__(self,BasFileName=None,PyFileName=None,Verbose=False,Doc=Doc):
        self.BasFileName=BasFileName
        self.Verbose=Verbose
        self.ActionKeyWord=['read','integer','real','character','real8']
        self.Doc=Doc
        self.MathSymbols='+-/*'
        if self.BasFileName is not None:
            if PyFileName is None:
                self.PyFileName=os.path.splitext(self.BasFileName)[0] + '.py'
            else:
                self.PyFileName=PyFileName
            if Verbose:
                print('## Converting UEDGE bas input file into py input file:')
                print('## {} -> {}'.format(self.BasFileName,self.PyFileName))
        
            self.ConvertBasFile()    
            
    def ConvertBasFile(self):
        with open(self.BasFileName,'r') as self.f:
            with open(self.PyFileName,'w') as self.fw:
                #Lines = f.readlines() 
                self.Line=next(self.f,None)
                while self.Line:
                    Out=self.LineParser(self.Line) 
                    if Out is not None:
                        for L in Out:
                            self.fw.write(L+'\n')
                    self.Line=next(self.f,None)  
                    
    def LineParser(self,Line):
            Out=[]
            if self.Verbose:
                print('Parsing: {}'.format(Line.rstrip()))
            # First we ignore all comments in line
            Comment=''
            Line=Line.strip()
            if Line.startswith('#'):
                Out.append(Line)
                return Out
            #Gather cut lines
            while Line.endswith('\\'):
               Lnext=next(self.f)
               self.Line=Lnext
               Line=Line.split('\\')[0]+Lnext.rstrip()
               
            if Line.count('#')>0:
                Comment=' #'+Line.split('#')[1]
                Line=Line.split('#')[0]
            
            # Second we look for several commands in one line 
            if Line.count(';')>0:
               Line=Line.split(';')
               Out=[]
               for L in Line:
                   Out.extend(self.LineParser(L)) 
                   print(Out)
               return Out       
            
            #Third Let do the special keywords
            if len(Line.split())>1:
                L=Line.strip().split()[0].strip()
                
                if L in 'read':
                    Arg=Line.split()[1].strip()
                    if len(Arg)>0:
                       Out.append('ReadInput("'+Arg+'")')
                       
                    else:
                        Out.append('##'+Line+Comment)
                    return Out
                    
                if L in 'real' or L in 'integer':
                    Out.append(Line.split()[1].strip()+Comment)
                    return Out
                
                if L in 'call':
                     Out.append(Line.split()[1].strip()+Comment)
                     
                     return Out
                
            
            if Line.count('=')==0: # not an assignment
                
                OutVar=self.ParseVariable(Line,ForceCheck=True)
                
                if OutVar is None:
                    Out.append('##'+Line+Comment) 
                else:
                    Out.append(OutVar+Comment)
            
            else:
                
                Var=Line.split('=')[0].strip()
                Arg=Line.split('=')[1].strip()
                Varout=self.ParseVariable(Var,ForceCheck=True)
                Argout=self.SubExpo(Arg)
                if Varout is not None:
                    Out.append(Varout+'='+Argout+Comment)
            
            return Out    
            
    def SubExpo(self,Arg):
         Out=Arg
         Out=Out.replace('D+','e+')
         Out=Out.replace('D-','e-')
         Out=Out.replace('d+','e+')
         Out=Out=Out.replace('d-','e-')
         Out=Out=Out.replace('.D','.e')
         Out=Out.replace('.d','.e')
         return Out
     
    def GetPackageMethod(self,VarName):
        Out=None
        for pkg in self.Doc.ListPkg:
            D={}
            exec('pkgobj='+pkg,globals(),D)
            obj=D['pkgobj']
            if hasattr(obj,VarName.strip()) and callable(getattr(obj,VarName.strip())):
                Out=pkg+'.'+VarName+'()'
        return Out
    
    def GetPackageVariable(self,VarName):
        Info=self.Doc.GetVarInfo(VarName,exact=True)
        if len(Info)>1:
            raise ValueError('More than one variable name found for {}'.format(VarName))
        if len(Info)<1:
           # if not a variable of a package, check if this is a package method      
            Out=None
        else:    
            Out=Info[0]['Package']+'.'+VarName
        return Out
    
    def ParseVariable(self,Var,ForceCheck=True):
        # Does Var has dimension?
        Out=None
        #print('#Var:',Var) 
        if Var.isnumeric():
            Out=Var
        else:    
            Out=self.ParseVariableName(Var.strip().split('(')[0].strip())
            
            if Out is None and not ForceCheck:
                Out=Var.strip().split('(')[0].strip()
            if Out is None and ForceCheck:
                Out=None
                return Out
            
            if '(' in Var:
                il=Var.find('(')+1
                ir=Var.rfind(')')
                Dim=Var[il:ir].strip()
                Dout=self.ParseVariableDim(Dim)
                #print(Dim,'->',Dout)
                Out=Out+'['+Dout+']'
        print('#:',Out)        
        return Out        
                
    def ParseVariableName(self,VarName):
        Out=self.GetPackageVariable(VarName)
        if Out is None:
            Out=self.GetPackageMethod(VarName)
        return Out
    
    def ParseMathSymbol(self,Str):
        Out=[]
        W=''
        
        for c in Str:
            if c in self.MathSymbols:
                Out.append(W)
                Out.append(c)
                W=''
            else:
                W=W+c
        Out.append(W)    
        return Out   
    
    def ParseVariableDim(self,Dim):
        Dout=[]
        idl=self.findall(Dim,'(')
        idr=self.findall(Dim,')')
        idc=self.findall(Dim,',')
        Ind=[-1]
        if (len(idl)!=len(idr)):
            raise ValueError('Cannot parse expression with incomplete bracket:',Dim)
        Ind.extend([i for i in idc if all([False for (il,ir) in zip(idl,idr) if i>=il and i<=ir])])
        DimArgs = [Dim[i+1:j] for i,j in zip(Ind, Ind[1:]+[None])]

        #print('DimArgs:',DimArgs)
        for D in DimArgs:
            #print('D:',D)
            d=D.strip()
            if d!='':
                rout=[]
                IsFirst=True
                for r in d.split(':'):
                    
                    rr=self.ParseMathSymbol(r)
                    print('rr:',rr)
                    outrr=[]
                    for rrr in rr:
                        if rrr in self.MathSymbols:
                            outrr.append(rrr)
                        else:
                            O=self.ParseVariable(rrr,ForceCheck=False)
                            outrr.extend(O)
                    routrr=outrr
                    if IsFirst:
                        if len(routrr)==1 and routrr[0].isnumeric(): 
                            routrr=[str(int(routrr[0])-1)]
                        else:
                            routrr.extend('-1')
                    IsFirst=False        
                    rout.append(''.join(routrr))
                Dout.append(':'.join(rout))
            else:
                Dout.append(':')
        Dout=','.join(Dout) 
        return Dout
    
    def findall(self,string,substring):
        """
        Function: Returning all the index of substring in a string
        Arguments: String and the search string
        Return:Returning a list
        """
        length = len(substring)
        c=0
        indexes = []
        while c < len(string):
            if string[c:c+length] == substring:
                indexes.append(c)
            c=c+1
        return indexes
# if __name__== "__main__":
#    UEDGEBas2Py('/home/jguterl/Dropbox/python/UEDGEInputDir/InputW.bas','/home/jguterl/Dropbox/python/UEDGEInputDir/InputW.py',Verbose=True)
#    UEDGEBas2Py('/home/jguterl/Dropbox/python/UEDGEInputDir/plate2.iter-feat.bas','/home/jguterl/Dropbox/python/UEDGEInputDir/plate2.iter-feat.py',Verbose=True)        
        
    