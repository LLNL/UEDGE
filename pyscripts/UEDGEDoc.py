#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 02:17:25 2020
@author: jguterl
"""

from colorama import Fore, Back, Style


import colorama
from colorama import Fore, Style,Back


import uedge
# from uedge import *
from . import UEDGEToolBox
class UedgeDoc(object):
    
    #---------------------------------------------------------------------------
    def __init__(self,PackageList=[],Verbose=False,Debug=True):
            self.ListPkg=UEDGEToolBox.GetListPackage()
            self.DocPkg={}
            self.DocGrp={}
            self.DocCom={}
            self._OnlyVar=False
            
            print('# Setup documentation of UEDGE')
            self.SetupDoc()
    #---------------------------------------------------------------------------
    def ListGrp(self,Str=None,OnlyVar=True):
        if Str is None:
            ListG=list(self.DocGrp.keys())
            ListG.sort()
            for G in ListG:
                self.PrintGrpInfo(G)
        else:
            print('Group in package:',Str)
            ListG=[g for g in list(self.DocGrp.keys()) if g['Package']==Str]
            ListG.sort()
            for G in ListG:
                self.PrintGrpInfo(G)
            
    #---------------------------------------------------------------------------
    def ListPkg(self,OnlyVar=True):
        ListP=list(self.DocPkg.keys())
        ListP.sort()
        for P in ListP:
            self.PrintPkgInfo(P)
                
    #---------------------------------------------------------------------------
    def ListVar(self,Str,InStr='',OnlyVar=True):
        ListV=[]
        if Doc.DocPkg.get(Str) is not None:
            ListV=[v for v in Doc.DocPkg.get(Str).keys() if v.startwith(InStr)]
            ListV.sort()
            for V in ListV:
                Doc.PrintVarInfo(Doc.DocPkg.get(Str)[V],OnlyVar) 
        if Doc.DocGrp.get(Str) is not None:
            ListV=[v for v in Doc.DocGrp.get(Str).keys() if v.startwith(InStr)]
            ListV.sort()
            for V in ListV:
                Doc.PrintVarInfo(Doc.DocGrp.get(Str)[V],OnlyVar) 
            
    #---------------------------------------------------------------------------  
    def ParseInfoVar(self,Str)->dict:
        d={}
        Fields=['Package','Group','Attributes','Dimension','Type','Address','Pyadress','Unit','Comment']
        for F in Fields:
            if F!='Comment':
                idx=Str.find(F+':')
                if idx>-1:
                    S=Str[idx+len(F+':'):]
                    idxS=S.find('\n')
                    S=S[:idxS].strip()
                else:
                    S=None
            else:
                idx=Str.find(F+':')
                if idx>-1: 
                    S=Str[idx+len(F+':'):].strip()
                else:
                    S=None
            d[F]=S
        return d        
    
    #---------------------------------------------------------------------------  
    def PrintGrpInfo(self,G):
        widthG=10
        S='{colorG}{:<{widthG}}{reset}'.format(G,colorG=Back.CYAN,reset=Style.RESET_ALL,widthG=widthG)
        print(S)
    #---------------------------------------------------------------------------  
    def PrintPkgInfo(self,P):
        widthP=5
        S='{colorG}{:<{widthP}}{reset}'.format(P,colorP=Back.RED,reset=Style.RESET_ALL,widthP=widthP)
        print(S)
    #---------------------------------------------------------------------------     
    def PrintVarInfo(self,Dic,OnlyVar=False):
        widthV=15
        widthG=10
        widthU=10
        if OnlyVar or self._OnlyVar:
            S='{colorV}{:<{widthV}}\
                {reset}'.format(Dic.get('Variable'),colorV=Back.GREEN,reset=Style.RESET_ALL,widthV=widthV)
            print(S)
        else:
                
            S='{colorV}{:<{widthV}}{reset} : {color}{:<}{reset}\
                / {colorG}{:<{widthG}}{reset}  [{:<}] : {} : {colorD}{:.10}{reset}\
                    '.format(Dic.get('Variable'),Dic.get('Package'),Dic.get('Group'),Dic.get('Unit'),\
                    Dic.get('Dimension'),str(Dic.get('Default')),widthV=widthV,colorV=Back.GREEN\
                        ,color=Back.RED,reset=Style.RESET_ALL,colorG=Back.CYAN,widthG=widthG,widthU=widthU,colorD=Back.MAGENTA)
            print(S)
    	#NLines=Dic.get('Comment').count('\n')
            for S in Dic.get('Comment').split('\n'):
                print('{:<{widthV}} : {:<}'.format(' ',S,widthV=widthV))
                print()
                
    #---------------------------------------------------------------------------  
    def SearchDoc(self,Str,OnlyVar=False):
       for pkg in self.ListPkg:
           Dic=self.DocPkg[pkg]
           for k in Dic.keys():
                   if Str in Dic[k]['Doc']: 
                       self.PrintVarInfo(Dic[k],OnlyVar)
                    
    #---------------------------------------------------------------------------    
    def Search(self,Str,exact=False,OnlyVar=False):
        for pkg in self.ListPkg:
            Dic=self.DocPkg[pkg]
            for k in Dic.keys():
                if exact:
                    if Str==k:
                        
                        self.PrintVarInfo(Dic[k],OnlyVar)
                else:
                    if Str in k: 
                        self.PrintVarInfo(Dic[k],OnlyVar) 
    #---------------------------------------------------------------------------                    
    def GetVarInfo(self,Str,exact=True):
        L=[]
        for pkg in self.ListPkg:
            Dic=self.DocPkg[pkg]
            for k in Dic.keys():
                if exact:
                    if Str==k:
                     L.append(Dic[k])
                else:
                    if Str in k: 
                     L.append(Dic[k])
        return L
    #---------------------------------------------------------------------------                     
    def SetupDoc(self):
        for pkg in self.ListPkg:
            loc={}
            self.DocPkg[pkg]={}
            exec('VarList=uedge.'+pkg +'py.'+ pkg+'.varlist()',globals(),loc)
            VarList=loc['VarList']
            for V in VarList:
                exec('VarStr=uedge.'+pkg +'py.'+ pkg+'.listvar("'+V+'")',globals(),loc)
                VarStr=loc['VarStr']
                
                if V in self.DocPkg[pkg].keys():
                    print('Variable is already defined:',self.DocVar[pkg][V])
                    print('V:',V,'pkg:',pkg)
                    raise ValueError('Variable is already defined')

                try :
                    exec('Val=' + pkg+'.'+V,globals(),loc)
                except: 
                   loc['Val']=None
                Val=loc['Val']

                
                self.DocPkg[pkg][V]=self.ParseInfoVar(VarStr)
                self.DocPkg[pkg][V]['Variable']=V
                if Val is None:
                    Val=''
                self.DocPkg[pkg][V]['Default']=Val
                self.DocPkg[pkg][V]['Doc']=VarStr
                #self.DocDict[pkg][V]
        for pkg in self.ListPkg:
            for V,D in (self.DocPkg[pkg].items()):
                if 'Group' in D.keys(): 
                    if not D['Group'] in self.DocGrp.keys():
                        self.DocGrp[D['Group']]={}
                    self.DocGrp[D['Group']][V]=D
                    self.DocGrp[D['Group']]['Package']=pkg
    #---------------------------------------------------------------------------
    def ToggleVarInfo(self):
        self._OnlyVar=not self._OnlyVar
    
    def HelpDoc(self):
        print('****** UEDGE Documentation ***** ')
        print('*** Documentation command ***')
        print('- ListGrp()')
        print('- ListPkg()')
        print('- ListVar(): list variable in package or group')
        print('- Search()')
        print('- SearchDoc()')                     
#---------------------------------------------------------------------------                
#---------------------------------------------------------------------------
                    
Doc=UedgeDoc(Verbose=False,Debug=False)
 
#---------------------------------------------------------------------------
def ListGrp(Str=None,OnlyVar=True):
    Doc.ListGrp(Str,OnlyVar)
    
#---------------------------------------------------------------------------
def ListPkg(Str,OnlyVar=True):
    Doc.ListPkg()
    
#---------------------------------------------------------------------------
def ListVar(Str,InStr='',OnlyVar=True):
    Doc.ListVar(Str,InStr,OnlyVar)
    
#---------------------------------------------------------------------------
def Search(Str,exact=False,OnlyVar=False):
    Doc.Search(Str,exact,OnlyVar)
#--------------------------------------------------------------------------
        
def SearchDoc(Str,OnlyVar=False):
    Doc.SearchDoc(Str,OnlyVar)

def HelpDoc():
    Doc.HelpDoc()