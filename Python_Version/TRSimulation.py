""" Python Class To Simulate the FBG Time Response

Copyright (C) Gilmar Pereira & Ben Frey
2015 DTU Wind Energy

Author: Gilmar Pereira & Ben Frey
Email Gilmar: gfpe@dtu.dk; gilmar_fp@outlook.com
Email Ben: ben.frey@stthomas.edu; freynben@gmail.com
Last revision: 05-31-2020

***License***:

This file is part of FBG_SiMul.

FBG_SiMul is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FBG_SiMul is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>
"""

#Packages
import numpy as np
import scipy as sp
import sympy
import matplotlib.pyplot as plt
import math
import cmath, random

class TRSimulation(object):
    def __init__(self,InputList,NumberFBG,FBGLength,Tolerance,SkipRow,FBGPosition,InputUnits):
        """ Initialized the classe
        This will load all the input files with the stress and strain.
        
        It will be assigned to each FBG sensor, an average strain,a maximum
        and minumum strain, and average stress.

        Parameters:
        ----------
        InputList:List (Mandatory)
                 List containing a path for every input files
        NumberFBG: int (Mandatory)
                Number of FBG per Optical fibre
        FBG length:float (Mandatory)
                Length of the Gratting                          
        Tolerance: float (Mandatory)
                Tolerance in the FBG length 
        SkipRow:Int (Mandatory)
                How any rows will be skip to load the file
        FBGPosition:List (Mandatory)
                Position of each FBG from the beggin of the Path
        InputUnits:Int (Mandatory)
                0: meters; 1: mm
        
        In this Class everything will be converted to (mm)
        """
        self.InputList=InputList
        self.NumberFBG=NumberFBG
        self.FBGLength=FBGLength
        self.Tolerance=Tolerance
        self.SkipRow=SkipRow
        self.FBGPosition=FBGPosition
        self.InputUnits=InputUnits
        
        """---------------------------------------------------------------------
        Creates an empty FBG array, where it will be assigned to each FBG sensor,
        an average strain,a maximum and minumum strain, and average stress.
        FBGArrayTR['xx']['yy'][zz]
        xx- Number of the FBG-1
        yy- variable 'x','AV-LE11','AV-LE22','AV-LE33',
            'Max-LE11','Max-LE22','Max-LE33',
            'Min-LE11','Min-LE22','Min-LE33',
            'Av-S11','AV-S22','AV-S33'
        zz - Increment number """
        
        #Naming the Array
        self.FBGArrayTR={}
        for b in np.arange(1,self.NumberFBG+1):
            self.FBGArrayTR['FBG'+str(b)]={}
            self.FBGArrayTR['FBG'+str(b)]['Increment']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-LE11']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-LE22']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-LE33']=[]
            self.FBGArrayTR['FBG'+str(b)]['Max-LE11']=[]
            self.FBGArrayTR['FBG'+str(b)]['Max-LE22']=[]
            self.FBGArrayTR['FBG'+str(b)]['Max-LE33']=[]
            self.FBGArrayTR['FBG'+str(b)]['Min-LE11']=[]
            self.FBGArrayTR['FBG'+str(b)]['Min-LE22']=[]
            self.FBGArrayTR['FBG'+str(b)]['Min-LE33']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-S11']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-S22']=[]
            self.FBGArrayTR['FBG'+str(b)]['AV-S33']=[]        
        
        #Converting the tolerance and the FBG length to the mm if Input unitas are in meters
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength/1000.0
            self.Tolerance=self.Tolerance/1000.0        

        #Loading the File
        names =('x','LE11','LE22','LE33','S11','S22','S33')#File Format
        formats=('f8','f8','f8','f8','f8','f8','f8',)
        dtypes = {'names' : names,'formats': formats}
        
        
        for i in range (0,len(self.InputList)):
            self.RawData=np.loadtxt(self.InputList[i], dtype=dtypes,skiprows=self.SkipRow)

            #Sorting the raw data for each FBG according with the FBG position
            for b in np.arange(0,self.NumberFBG):
                #Here it will be calculated the average strain,a maximum and minumum strain, and average stress
                LE11=[] #This is a temp variable
                LE22=[]
                LE33=[]
                S11=[]
                S22=[]
                S33=[]

                for f in np.arange(0,len(self.RawData)):
                    #Check the lines inside the FBG length+tolerance
                    if self.RawData['x'][f]>self.FBGPosition[b]-self.Tolerance and self.RawData['x'][f]<self.FBGPosition[b]+self.FBGLength+self.Tolerance:
                        #Convert to SI(mm,MPA)
                        if self.InputUnits==0:
                            LE11.append(self.RawData['LE11'][f])
                            LE22.append(self.RawData['LE22'][f])
                            LE33.append(self.RawData['LE33'][f])
                            S11.append(self.RawData['S11'][f]/(10**6))
                            S22.append(self.RawData['S22'][f]/(10**6))
                            S33.append(self.RawData['S33'][f]/(10**6))
                       
                        else:
                            LE11.append(self.RawData['LE11'][f])
                            LE22.append(self.RawData['LE22'][f])
                            LE33.append(self.RawData['LE33'][f])
                            S11.append(self.RawData['S11'][f])
                            S22.append(self.RawData['S22'][f])
                            S33.append(self.RawData['S33'][f])

                #Calculating average strain,a maximum and minumum strain, and average stress
                self.FBGArrayTR['FBG'+str(b+1)]['Increment'].append(i)
                #Average Strain
                self.FBGArrayTR['FBG'+str(b+1)]['AV-LE11'].append(np.mean(LE11))
                self.FBGArrayTR['FBG'+str(b+1)]['AV-LE22'].append(np.mean(LE22))
                self.FBGArrayTR['FBG'+str(b+1)]['AV-LE33'].append(np.mean(LE33))
                #Maximum Strain
                self.FBGArrayTR['FBG'+str(b+1)]['Max-LE11'].append(np.max(LE11))
                self.FBGArrayTR['FBG'+str(b+1)]['Max-LE22'].append(np.max(LE22))
                self.FBGArrayTR['FBG'+str(b+1)]['Max-LE33'].append(np.max(LE33))
                #Maximum Strain
                self.FBGArrayTR['FBG'+str(b+1)]['Min-LE11'].append(np.min(LE11))
                self.FBGArrayTR['FBG'+str(b+1)]['Min-LE22'].append(np.min(LE22))
                self.FBGArrayTR['FBG'+str(b+1)]['Min-LE33'].append(np.min(LE33))              
                #Average Strain
                self.FBGArrayTR['FBG'+str(b+1)]['AV-S11'].append(np.mean(S11))
                self.FBGArrayTR['FBG'+str(b+1)]['AV-S22'].append(np.mean(S22))
                self.FBGArrayTR['FBG'+str(b+1)]['AV-S33'].append(np.mean(S33))
                            
        #Converting tolerance and FBG length back to mm
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength*1000.0
            self.Tolerance=self.Tolerance*1000.0
            


    def Calculate(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient,FBGDirection):
        """ 
        Function to Simulate the FBG signal
        Input:
        
        FBGOriginalWavel: List(Mandatory)
                        List with the FBG array original wavelength
        PhotoElasticParam: float (Mandatory)
                        Photo-Elastic parameter
        InitialRefractiveIndex: float (Mandatory)
                        Initial effective refractive index (neff)                        
        FBGDirection: int (mandatory)
                    0-xx 1-yy 2-zz
        DirectionalRefractiveP11: float (Mandatory)
                    Directional Photo-elastic coeffic
        DirectionalRefractiveP12: float (Mandatory)
                    Directional Photo-elastic coefficient      
        YoungsModule: float (Mandatory)
                    Optical Fibre Young's Module
        PoissionsCoefficient: float (Mandatory)
                    Optical Fibre Poisson Ration
        """        
        self.FBGOriginalWavel=FBGOriginalWavel
        self.PhotoElasticParam=PhotoElasticParam
        self.InitialRefractiveIndex=InitialRefractiveIndex
        self.FBGDirection=FBGDirection  
        self.DirectionalRefractiveP11=DirectionalRefractiveP11
        self.DirectionalRefractiveP12=DirectionalRefractiveP12
        self.YoungsModule=YoungsModule/10**6 #Change to MPA
        self.PoissionsCoefficient=PoissionsCoefficient
        
        #Directions
        if self.FBGDirection==0: #FBG longitudinal direction (xx)
            AVStrainDirection='AV-LE11'
            MaxStrainDirection='Max-LE11'
            MinStrainDirection='Min-LE11'
            AVStresszzDirection='AV-S33'
            AVStressyyDirection='AV-S22'

        if self.FBGDirection==1: #FBG longitudinal direction (yy)
            AVStrainDirection='AV-LE22'
            MaxStrainDirection='Max-LE22'
            MinStrainDirection='Min-LE22'
            AVStresszzDirection='AV-S33'
            AVStressyyDirection='AV-S11'
            
        if self.FBGDirection==2: #FBG longitudinal direction (yy)
            AVStrainDirection='AV-LE33'
            MaxStrainDirection='Max-LE33'
            MinStrainDirection='Min-LE33'
            AVStresszzDirection='AV-S11'
            AVStressyyDirection='AV-S22'

        #fixed component: Transverse stress"""
        fce=(((1+self.PoissionsCoefficient)*self.DirectionalRefractiveP12-(1+self.PoissionsCoefficient)*self.DirectionalRefractiveP11)*self.InitialRefractiveIndex**3)/self.YoungsModule
        
        #The sensor output will be presented as FBGTimeResponse
        self.FBGTimeResponse={}
        for b in np.arange(1,self.NumberFBG+1):
            self.FBGTimeResponse['FBG'+str(b)]={}
            self.FBGTimeResponse['FBG'+str(b)]['Increment']=[]
            self.FBGTimeResponse['FBG'+str(b)]['WaveShift']=[]
            self.FBGTimeResponse['FBG'+str(b)]['WaveWidth']=[]
        
        for b in np.arange(0,self.NumberFBG): #Cycle each FBG
            for i in np.arange(0,len(self.FBGArrayTR['FBG'+str(b+1)]['Increment'])):
                #Increment Numbering
                self.FBGTimeResponse['FBG'+str(b+1)]['Increment'].append(str(i))
                
                #--------Wavelength Shift Calculation----------------------
                WavlShift=(1.0-self.PhotoElasticParam)*self.FBGOriginalWavel[b]*self.FBGArrayTR['FBG'+str(b+1)][AVStrainDirection][i]
                self.FBGTimeResponse['FBG'+str(b+1)]['WaveShift'].append(WavlShift)
                
                #---Variation of the peak width: Non-uniform strain contribution----
                #Grating period caused by maximum strain
                graperiodmax=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*self.FBGArrayTR['FBG'+str(b+1)][MaxStrainDirection][i])
                #Grating period caused by minimum strain
                graperiodmin=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*self.FBGArrayTR['FBG'+str(b+1)][MinStrainDirection][i])
                #Peak width variation
                PeakWV1=2*self.InitialRefractiveIndex*(graperiodmax-graperiodmin)
                
                #---Variation of the peak width: Transverse stress contribution----
                PeakWV2=abs(self.FBGArrayTR['FBG'+str(b+1)][AVStresszzDirection][i]-self.FBGArrayTR['FBG'+str(b+1)][AVStressyyDirection][i])*(self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex))*fce
                
                #Add Both contributions
                PeakWV=PeakWV1+PeakWV2
                self.FBGTimeResponse['FBG'+str(b+1)]['WaveWidth'].append(PeakWV)
                
