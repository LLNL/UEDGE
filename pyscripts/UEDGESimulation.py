#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:52:13 2020

@author: jguterl
"""
import uedge
import matplotlib.pyplot as plt
import numpy as np
import os, sys, math, string, re
from pathlib import Path
from . import UEDGEToolBox
from uedge import *
from colorama import Fore, Back, Style

class UEDGESimulation(object):
       
    def __init__(self,  *args, **kwargs):
        # Import Uedge packages as attribute of the class instance 
        # WARNING: this is not a deep copy!!!
        self.ListPkg=UEDGEToolBox.GetListPackage()
        print('Loaded Packages:',self.ListPkg)
        for pkg in self.ListPkg:
            exec('self.' + pkg + '=' + pkg,globals(),locals())
        
            
    def ReadInput(self,FileName:str,Verbose:bool=True):
        for pkg in self.ListPkg:
            exec('from uedge import '+pkg)
        print('### Looking for input file:',FileName)
        if not os.path.exists(FileName):
            print('### Cannot find {}'.format(os.path.abspath(FileName)))
            FileName_=os.path.join(Settings.InputDir,FileName)
            if not os.path.exists(FileName_):
                raise IOError('Cannot find the requested file in InputDir either:{} '.format(FileName))
            else:
                print('### File found in InputDir:{}'.format(FileName_))
        else: 
            FileName_=os.path.abspath(FileName)
            print('### File found :{}'.format(FileName_))
            
        print('### Loading {} '.format(FileName_))    
        
        f=open(FileName_,'r')
        lines=f.read()
        f.close()
        Lines=lines.splitlines()
        count=1
        for L in Lines:
            if not L.strip().startswith('#'):
                if Verbose:
                    print('{} : {}'.format(count,L))
                exec(L)
            count=count+1
            
        
         
            # except Exception as e:
            #     print(repr(e))
            #     print('Could not import:', os.path.join(Settings.InputDir,FileName)) 
            
            
    def ReadDivertorPlateFile(self,FileName,Verbose=False):
    
        Dic={'rplate1':[],'zplate1':[],'rplate2':[],'zplate2':[],'Comment':'No comment','FileName':FileName,'InnerDivPlateFile':None,
        'OuterDivPlateFile':None}
        
        
        try:
            with  open(os.path.join(FileName)) as f:
                exec(f.read(),Dic)
        except:
            raise ValueError('Cannot Read File:'+FileName)
        
        
        
        if Dic['InnerDivPlateFile'] is not None:
             InnerData=np.loadtxt(Dic['InnerDivPlateFile'])
             if Verbose: print(InnerData)   
             Dic['rplate1']=list(InnerData[:,0])
             Dic['zplate1']=list(InnerData[:,1])
        
        if Dic['OuterDivPlateFile'] is not None:
            OuterData=np.loadtxt(Dic['OuterDivPlateFile'])  
            Dic['rplate2']=list(OuterData[:,0])
            Dic['zplate2']=list(OuterData[:,1])
        # Check if data are correct
        # (len(rplate2)<2 or len(zplate2)<2 or len(rplate1)<2 or len(zplate1)<2):
         #       raise ValueError('Wrong size of plates coordinates')
        return Dic
    
    # outer is 2, inner is 1


    def SetDivertorPlates(self,FileName):
        Dic=ReadDivertorPlateFile(self,FileName)
        if len(Dic['rplate1'])<2 or len(Dic['rplate1'])!=len(Dic['zplate1']):
            raise ValueError('wrong dimension of coordinates of plate #1')
        if len(Dic['rplate2'])<2 or len(Dic['rplate2'])!=len(Dic['zplate2']):
            raise ValueError('wrong dimension of coordinates of plate #2')
        print('FileName:'+Dic['FileName'])
        print('Comment:'+Dic['Comment'])
        self.grd.nplate1=len(Dic['rplate1'])
        self.grd.nplate2=len(Dic['rplate2'])
        self.grd.gchange("Mmod")
        self.grd.rplate2=Dic['rplate2']
        self.grd.zplate2=Dic['zplate2']
        self.grd.rplate2=Dic['rplate1']
        self.grd.zplate1=Dic['zplate1']
    
    def PlotDivertorPlates(self,FileName=None):
        if type(FileName)==str:
            Dic=ReadDivertorPlateFile(FileName)
            if len(Dic['rplate1'])<2 or len(Dic['rplate1'])!=len(Dic['zplate1']):
                print('rplate1=',Dic['rplate1'])
                print('zplate1=',Dic['zplate1'])
                raise ValueError('wrong dimension of coordinates of plate #1')
            if len(Dic['rplate2'])<2 or len(Dic['rplate2'])!=len(Dic['zplate2']):
                raise ValueError('wrong dimension of coordinates of plate #2')
        
            r1=Dic['rplate1']
            r2=Dic['rplate2']
            z1=Dic['zplate1']
            z2=Dic['zplate2']
            plt.suptitle=Dic['Comment']
            plt.xlabel('R[m]')
            plt.ylabel('Z[m]')
            plt.axis('equal')
            plt.grid(True)
        else:
            r1=grd.rplate1
            r2=grd.rplate2
            z1=grd.zplate1
            z2=grd.zplate2
        
        plt.plot(r1,z1,color='r')
        plt.plot(r2,z2,'g')
    def PrintInfo(self,Str,color=Back.CYAN):
        print("*---------------------------------------------------------*")
        print("{color}{}{reset}".format(Str,color=color,reset=Style.RESET_ALL))
        print("*---------------------------------------------------------*")
              
    def InitRun(self):
            if (bbb.iterm == 1 and bbb.ijactot>1):
               self.PrintInfo("Initial successful time-step exists",Back.GREEN)
               return
            else:
               self.PrintInfo("Need to take initial step with Jacobian; trying to do here",Back.CYAN)
               bbb.icntnunk = 0
               bbb.exmain()
               
            if (bbb.iterm != 1):
                self.PrintInfo("Error: converge an initial time-step first",Back.RED)
                bbb.exmain_aborted=1
            else:
               self.PrintInfo("First initial time-step has converged",Back.GREEN) 
               return
    
    def Restart(self,dt=None,mul_dt=3):
        bbb.iterm=1
        bbb.dtreal/=mult_dt
        self.RunTime(dt,mul_dt)
              
    def GetMinTimeStep(self):
        #evaluate yldot
        # bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
        # dt=np.zeros(bbb.neq) 
        # for i in range(bbb.neq):
        #     if bbb.iseqalg[i]==1:
        #         print('eq al:',i)
        #         dt[i]=1e20
        #     else:
        #         if bbb.yldot[i]==0.0:
        #             dt[i]=1e20
        #         else:    
        #             dt[i]=abs(bbb.yl[i]/(bbb.yldot[i]*bbb.sfscal[0:bbb.neq]))
        # return dt.min()
        return -1
    def RunTime(self,dt=None,tstop=10,Imax=500,Jmax=5,ftol_dt=1e-10,itermx=7,rlx=0.9,incpset=7,method_dt=0,mult_dt_fwd=3.4,mult_dt_bwd=3):
        # this allow run restart 
        bbb.rlx=rlx
        bbb.ftol_dt=ftol_dt
        bbb.itermx=itermx
        bbb.incpset=incpset
        if dt is not None:
            bbb.dtreal=dt
        bbb.exmain_aborted=0    
        while bbb.exmain_aborted==0:
            self.InitRun()   
            if bbb.exmain_aborted==1:
                break
            self.PrintInfo('----Starting Main Loop ----')
            for imain in range(Imax):
                bbb.icntnunk = 0
                
                bbb.ylodt = bbb.yl #this is use in set_dt which provides option to select optimal timestep
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq]*bbb.sfscal[0:bbb.neq])**2))
                bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
                
               
                
                if (bbb.initjac == 0): 
                    bbb.newgeo=0
                self.PrintInfo('Number time-step changes = {} New time-step = {colorT}{:.4E}{reset}'.format(imain, bbb.dtreal,reset=Style.RESET_ALL,colorT=Back.MAGENTA)) 
                self.PrintInfo('Suggested timestep:{}'.format(self.GetMinTimeStep()))
                bbb.exmain() # take a single step at the present bbb.dtreal
                if bbb.exmain_aborted==1:
                    break
                
                if (bbb.iterm == 1):
                    bbb.ylodt = bbb.yl
                    bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                    fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                    if (bbb.dt_tot>=bbb.t_stop):
                        self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                        break
                    bbb.icntnunk = 1
                    bbb.isdtsfscal = 0
                    for ii2 in range( 1, bbb.ii2max+1): #take ii2max steps at the present time-step
                        if bbb.exmain_aborted==1:
                            break
                        bbb.itermx = itermx
                        bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
                        bbb.exmain()
                        if bbb.exmain_aborted==1:
                            break
                        if bbb.iterm == 1:
                            bbb.ylodt = bbb.yl
                            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                            fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                            bbb.dt_tot += bbb.dtreal
                        else:
                            break
                        if (bbb.dt_tot>=bbb.t_stop):
                            self.PrintInfo('SUCCESS: dt_tot >= t_stop')
                            bbb.exmain_aborted=1
                            break
                        
                if bbb.exmain_aborted==1:
                            break
                        
                if (bbb.iterm != 1):    #print bad eqn, cut dtreal by 3, set irev flag
                    #print('\n\n\n exmain_aborted:',bbb.exmain_aborted,'\n\n\n')
                    self.Itrouble()
                    
                    if (bbb.dtreal < bbb.dt_kill):
                        self.PrintInfo('FAILURE: time-step < dt_kill',Back.RED)
                        break
                    self.PrintInfo('Converg. fails for bbb.dtreal; reduce time-step by 3',Back.RED) 
                    bbb.dtreal /= mult_dt_bwd
                    
                    
                if (bbb.iterm == 1):
                    bbb.dtreal *= mult_dt_fwd
                    
                bbb.iterm = 1
                
    def Itrouble(self):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod,argmax
        from uedge import bbb
        # Set scaling factor
        scalfac = bbb.sfscal
        if (bbb.svrpkg[0].decode('UTF-8').strip() != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.
    
        # Find the fortran index of the troublemaking equation
        itrouble=argmax(abs(bbb.yldot[:bbb.neq]))+1
        print("** Fortran index of trouble making equation is:")
        print(itrouble)
    
        # Print equation information
        print("** Number of equations solved per cell:")
        print("numvar = {}".format(bbb.numvar))
        print(" ")
        iv_t = mod(itrouble-1,bbb.numvar) + 1 # Use basis indexing for equation number
        print("** Troublemaker equation is:")
        # Verbose troublemaker equation
        if abs(bbb.idxte-itrouble).min()==0:
            print('Electron energy equation: iv_t={}'.format(iv_t))           
        elif abs(bbb.idxti-itrouble).min()==0:
            print('Ion energy equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxphi-itrouble).min()==0:
            print('Potential equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxu-itrouble).min()==0:
            for species in range(bbb.idxu.shape[2]):
                if abs(bbb.idxu[:,:,species]-itrouble).min()==0:
                    print('Ion momentum equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxn-itrouble).min()==0:
            for species in range(bbb.idxn.shape[2]):
                if abs(bbb.idxn[:,:,species]-itrouble).min()==0:
                    print('Ion density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxg-itrouble).min()==0:
            for species in range(bbb.idxg.shape[2]):
                if abs(bbb.idxg[:,:,species]-itrouble).min()==0:
                    print('Gas density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxtg-itrouble).min()==0:
            for species in range(bbb.idxtg.shape[2]):
                if abs(bbb.idxtg[:,:,species]-itrouble).min()==0:
                    print('Gas temperature equation of species {}: iv_t={}'.format(species, iv_t))   
        # Display additional information about troublemaker cell
        print(" ")
        print("** Troublemaker cell (ix,iy) is:")
        print(bbb.igyl[itrouble-1,])
        print(" ")
        print("** Timestep for troublemaker equation:")
        print(bbb.dtuse[itrouble-1])
        print(" ")
        print("** yl for troublemaker equation:")
        print(bbb.yl[itrouble-1])
        print(" ")
        
    def WhichEq(self,itrouble):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod,argmax
        from uedge import bbb
        # Set scaling factor
        scalfac = bbb.sfscal
        if (bbb.svrpkg[0].decode('UTF-8').strip() != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.
    
        # Find the fortran index of the troublemaking equation
        print("** Fortran index of trouble making equation is:")
        print(itrouble)
    
        # Print equation information
        print("** Number of equations solved per cell:")
        print("numvar = {}".format(bbb.numvar))
        print(" ")
        iv_t = mod(itrouble-1,bbb.numvar) + 1 # Use basis indexing for equation number
        print("** Troublemaker equation is:")
        # Verbose troublemaker equation
        if abs(bbb.idxte-itrouble).min()==0:
            print('Electron energy equation: iv_t={}'.format(iv_t))           
        elif abs(bbb.idxti-itrouble).min()==0:
            print('Ion energy equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxphi-itrouble).min()==0:
            print('Potential equation: iv_t={}'.format(iv_t))   
        elif abs(bbb.idxu-itrouble).min()==0:
            for species in range(bbb.idxu.shape[2]):
                if abs(bbb.idxu[:,:,species]-itrouble).min()==0:
                    print('Ion momentum equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxn-itrouble).min()==0:
            for species in range(bbb.idxn.shape[2]):
                if abs(bbb.idxn[:,:,species]-itrouble).min()==0:
                    print('Ion density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxg-itrouble).min()==0:
            for species in range(bbb.idxg.shape[2]):
                if abs(bbb.idxg[:,:,species]-itrouble).min()==0:
                    print('Gas density equation of species {}: iv_t={}'.format(species, iv_t))   
        elif abs(bbb.idxtg-itrouble).min()==0:
            for species in range(bbb.idxtg.shape[2]):
                if abs(bbb.idxtg[:,:,species]-itrouble).min()==0:
                    print('Gas temperature equation of species {}: iv_t={}'.format(species, iv_t))   
        # Display additional information about troublemaker cell
        print(" ")
        print("** Troublemaker cell (ix,iy) is:")
        print(bbb.igyl[itrouble-1,])
        print(" ")
        print("** Timestep for troublemaker equation:")
        print(bbb.dtuse[itrouble-1])
        print(" ")
        print("** yl for troublemaker equation:")
        print(bbb.yl[itrouble-1])
        print(" ")
Sim=UEDGESimulation()

def RunTime(*arg,**kwargs):
    Sim.RunTime(*arg,**kwargs)

def ReadInput(FileName,Verbose=False):
    Sim.ReadInput(FileName,Verbose=Verbose)
    
def ir():
    Sim.InitRun()
    
                   
        
        
     
