#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 22:09:43 2020

@author: jguterl
"""

import sys,os
import easygui
import platform,inspect
import configparser
import getpass

import numpy as np
class UEDGESettings():
    def __init__(self):
        
        self.Platform=self.GetPlatform()
        self.config={'UEDGE':{'UserName':None,'RunDir':None,'SaveDir':None,'InputDir':None}}
        self.config=self.CheckUEDGEConfig()
        if self.config is None:
            print('Using default settings. Run InitConfig() to set the configuration...')
            self.UserName=getpass.getuser()
            self.RunDir=os.getcwd()
            self.SaveDir=os.getcwd()
            self.InputDir=os.getcwd() 
        else:
            print(' # Loading UEDGE settings from {}'.format(self.ConfigFileName))
            self.UserName=self.config['UEDGE'].get('UserName')
            self.RunDir=self.config['UEDGE'].get('RunDir')
            self.SaveDir=self.config['UEDGE'].get('SaveDir')
            self.InputDir=self.config['UEDGE'].get('InputDir') 

    def Config(self):
        '''Print current UEDGE configuration'''
        print('******** UEDGE configuration ********')
        print('Username:',self.UserName)
        print('RunDir:',self.RunDir)
        print('SaveDir:',self.SaveDir)
        print('InputDir:',self.InputDir)
    def GetPlatform(self):
        PF={}
        AllFunctions = inspect.getmembers(platform, inspect.isfunction)
        for (n,f) in AllFunctions:
            if not n.startswith('_'):
                try: 
                    PF[n]=f()
                except:
                    pass
        return PF
    
    def CheckUEDGEConfig(self):
        self.ConfigFileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
        if not os.path.exists(self.ConfigFileName):
            print ('Cannot find the config file ' + self.ConfigFileName +'....')
            return None
        config = configparser.ConfigParser()
        try:
            config.read(self.ConfigFileName)
            return config
        except: 
            print('Cannot parse the config file:',self.ConfigFileName)
            return None
    def CreateStamp(self): 
            from uedge import bbb
            Stamp={}
            Stamp['time'] = time.time()
            Stamp['ctime'] = time.ctime()
            Stamp['User']=self.UserName
            Stamp['Version'] = bbb.uedge_ver
            Stamp['PlatForm'] = self.Platform 
        
def CreateUEDGESettingsFile():
        config = configparser.ConfigParser()
        GetInput=True
        while GetInput:
            UserName=input('Enter the username:')
            
            RunDir=easygui.diropenbox(title='Select RunDir for pyUEDGE',default='.')
            SaveDir=easygui.diropenbox(title='Select SaveDir for pyUEDGE',default='.')
            InputDir=easygui.diropenbox(title='Select InputDir for pyUEDGE',default='.')
            if RunDir is None:
                RunDir=os.getcwd()
            if SaveDir is None:
                RunDir=os.getcwd()
            if InputDir is None:
                InputDir=os.getcwd()
            if UserName is None:
                UserName='unknown'
            RunDir=os.path.os.path.abspath(RunDir)
            SaveDir=os.path.os.path.abspath(SaveDir)
            InputDir=os.path.os.path.abspath(InputDir)
            print(' # Settings:')
            print('  - UserName:',UserName)
            print('  - RunDir:',RunDir)
            print('  - SaveDir:',SaveDir)
            print('  - InputDir:',InputDir)
            if QueryYesNo('Are these settings correct?'):
               config['UEDGE']={}
               config['UEDGE']['UserName']=UserName
               config['UEDGE']['RunDir']=RunDir
               config['UEDGE']['SaveDir']=SaveDir
               config['UEDGE']['InputDir']=InputDir
               ConfigFileName=os.path.join(os.path.expanduser("~"),'.UedgeSettings')
               with open(ConfigFileName, 'w') as f:
                   config.write(f)
               print(' # Creation of the config file:'+ConfigFileName)
               GetInput=False
               return
            elif not QueryYesNo('Enter the settings again?'):
               print(' # The config file .UEDGEInfo has not been created in the home folder...') 
               GetInput=False
               return 





def QueryYesNo(question, default="no"):
        """Ask a yes/no question via input() and return their answer.
    
        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).
    
        The "answer" return value is True for "yes" or False for "no".
        """
        valid = {"yes": True, "y": True, "ye": True,
                 "no": False, "n": False}
        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)
    
        while True:
            sys.stdout.write(question + prompt)
            choice = input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")  
                
Settings=UEDGESettings()
S=Settings

def CdRunDir():
    os.chdir(Settings.RunDir)
    print('Current working directory:',os.getcwd())
 
def CdSaveDir():
    os.chdir(Settings.SaveDir)
    print('Current working directory:',os.getcwd())
    
def CdInputDir():
    os.chdir(Settings.InputDir)
    print('Current working directory:',os.getcwd())
def Config():
    Settings.Config() 
    
def InitConfig():
    CreateUEDGESettingsFile()  
# start into the run folder
CdRunDir()    