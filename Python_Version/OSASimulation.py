""" Python Class To Simulate the FBG Spectrum (Specific Step)

Copyright (C) Gilmar Pereira & Ben Frey

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

class OSASimulation(object):
    def __init__(self,filename,NumberFBG,FBGLength,Tolerance,SkipRow,FBGPosition,InputUnits):
        """ Initialized the classe
        This will load the file with the stress and strain, and create a variable
        per FBG with the correspondent Strain and Stress
        
        Parameters:
        ----------
        Filename:string (Mandatory)
                 Path to the file
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
        self.filename=str(filename)
        self.NumberFBG=NumberFBG
        self.FBGLength=FBGLength
        self.Tolerance=Tolerance
        self.SkipRow=SkipRow
        self.FBGPosition=FBGPosition
        self.InputUnits=InputUnits
        
        
        #Loading the File
        names =('x','LE11','LE22','LE33','S11','S22','S33')
        formats=('f8','f8','f8','f8','f8','f8','f8',)
        dtypes = {'names' : names,'formats': formats}
        self.RawData=np.loadtxt(self.filename, dtype=dtypes,skiprows=self.SkipRow)
        
        """------------------------------------------------------------------------------
        Creating FBG array, and give them the correspondent Data
        FBGArray['xx']['yy'][zz]
        xx- Number of the FBG-1
        yy- variable 'x','LE11','LE22','LE33','S11','S22','S33'
        zz line(number of elements per grating)
        """
        #Naming the Array
        self.FBGArray={}
        for b in np.arange(1,self.NumberFBG+1):
            self.FBGArray['FBG'+str(b)]={}
            self.FBGArray['FBG'+str(b)]['x']=[]
            self.FBGArray['FBG'+str(b)]['LE11']=[]
            self.FBGArray['FBG'+str(b)]['LE22']=[]
            self.FBGArray['FBG'+str(b)]['LE33']=[]
            self.FBGArray['FBG'+str(b)]['S11']=[]
            self.FBGArray['FBG'+str(b)]['S22']=[]
            self.FBGArray['FBG'+str(b)]['S33']=[]
            
        #Converting tolerance and FBG length to the SI if in meters
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength/1000.0
            self.Tolerance=self.Tolerance/1000.0

        #Sorting Data for each FBG
        for b in np.arange(0,self.NumberFBG):
            for f in np.arange(0,len(self.RawData)):
                #Check the lines inside the FBG length+tolerance
                if self.RawData['x'][f]>self.FBGPosition[b]-self.Tolerance and self.RawData['x'][f]<self.FBGPosition[b]+self.FBGLength+self.Tolerance:
                    #Convert to SI(mm,MPA)
                    if self.InputUnits==0:
                        self.FBGArray['FBG'+str(b+1)]['x'].append(self.RawData['x'][f]*1000)
                        self.FBGArray['FBG'+str(b+1)]['LE11'].append(self.RawData['LE11'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE22'].append(self.RawData['LE22'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE33'].append(self.RawData['LE33'][f])
                        self.FBGArray['FBG'+str(b+1)]['S11'].append(self.RawData['S11'][f]/(10**6))
                        self.FBGArray['FBG'+str(b+1)]['S22'].append(self.RawData['S22'][f]/(10**6))
                        self.FBGArray['FBG'+str(b+1)]['S33'].append(self.RawData['S33'][f]/(10**6))
                   
                    else:                        
                        self.FBGArray['FBG'+str(b+1)]['x'].append(self.RawData['x'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE11'].append(self.RawData['LE11'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE22'].append(self.RawData['LE22'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE33'].append(self.RawData['LE33'][f])
                        self.FBGArray['FBG'+str(b+1)]['S11'].append(self.RawData['S11'][f])
                        self.FBGArray['FBG'+str(b+1)]['S22'].append(self.RawData['S22'][f])
                        self.FBGArray['FBG'+str(b+1)]['S33'].append(self.RawData['S33'][f])
                        
                #Converting tolerance and FBG length bacl to mm
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength*1000.0
            self.Tolerance=self.Tolerance*1000.0

    def UndeformedFBG(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution):
        """ 
        Function to Simulate the Undeformed FBG signal
        Input:
        
        FBGOriginalWavel: List(Mandatory)
                        List with the FBG array original wavelength
        PhotoElasticParam: float (Mandatory)
                        Photo-Elastic parameter
        InitialRefractiveIndex: float (Mandatory)
                        Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
                        Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
                        Fringe Visibility (Fv)
        MinBandWidth: float (Mandatory)
                        Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
        SimulationResolution: float (Mandatory)
                        Simulation resolution- Wavelength increment  
        """        
        self.FBGOriginalWavel=FBGOriginalWavel
        self.PhotoElasticParam=PhotoElasticParam
        self.InitialRefractiveIndex=InitialRefractiveIndex
        self.MeanChangeRefractiveIndex=MeanChangeRefractiveIndex
        self.FringeVisibility=FringeVisibility
        self.MinBandWidth=MinBandWidth
        self.MaxBandWidth=MaxBandWidth
        self.SimulationResolution=SimulationResolution
        
        #Array with the Original FBG period
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):
            self.APFBG.append(self.FBGOriginalWavel[i]*(1-self.MeanChangeRefractiveIndex/self.InitialRefractiveIndex)/(2.0*self.InitialRefractiveIndex))
            
        #Empty Original Reflec spectrum
        self.OReflect={}
        self.OReflect['wavelength']=[]
        self.OReflect['reflec']=[]
        
        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):    
            # Wavelength cycle (Here the simulation resolution is used)
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix
                M=20.0 #Sections the gratting is divided-- Transfer Matrix
                #FBG increment size (nm)
                deltz=(self.FBGLength*(10.0**6))/M
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*self.InitialRefractiveIndex*((1.0/l)-(1.0/(2.0*self.InitialRefractiveIndex*self.APFBG[i])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                self.OReflect['wavelength'].append(l)
                self.OReflect['reflec'].append(REF) 


    def UniformStrain(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,FBGDirection):
        """ 
        Function to Simulate the FGB Spectrum Considering only Uniform Strain
        Contribution
        
        The uniform strain is calculating by the average strain along the
        Optical fiber
        
        Input:
        
        FBGOriginalWavel: List(Mandatory)
                        List with the FBG array original wavelength
        PhotoElasticParam: float (Mandatory)
                        Photo-Elastic parameter
        InitialRefractiveIndex: float (Mandatory)
                        Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
                        Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
                        Fringe Visibility (Fv)
        MinBandWidth: float (Mandatory)
                        Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
        SimulationResolution: float (Mandatory)
                        Simulation resolution- Wavelength increment 
        FBGDirection: int (mandatory)
                    0-xx 1-yy 2-zz
        """        
        self.FBGOriginalWavel=FBGOriginalWavel
        self.PhotoElasticParam=PhotoElasticParam
        self.InitialRefractiveIndex=InitialRefractiveIndex
        self.MeanChangeRefractiveIndex=MeanChangeRefractiveIndex
        self.FringeVisibility=FringeVisibility
        self.MinBandWidth=MinBandWidth
        self.MaxBandWidth=MaxBandWidth
        self.SimulationResolution=SimulationResolution
        self.FBGDirection=FBGDirection
        
        #Array with the FBG period
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):  
            if self.FBGDirection==0: #FBG longitudinal direction (xx)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE11'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

            if self.FBGDirection==1: #FBG longitudinal direction (yy)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE22'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

            if self.FBGDirection==2: #FBG longitudinal direction (zz)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE33'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

        #Empty Reflec spectrum -Uniform Strain
        self.USReflect={}
        self.USReflect['wavelength']=[]
        self.USReflect['reflec']=[]
        
        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):    
            # Wavelength cycle (Here the simulation resolution is used)
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix 
                M=len(self.FBGArray['FBG'+str(i+1)]['x']) #Sections the gratting is divided-- Transfer Matrix
                #FBG increment size (nm)
                deltz=(self.FBGLength*(10.0**6))/M
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*self.InitialRefractiveIndex*((1.0/l)-(1.0/(2.0*self.InitialRefractiveIndex*self.APFBG[i])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                self.USReflect['wavelength'].append(l)
                self.USReflect['reflec'].append(REF) 



    def NonUniformStrain(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,FBGDirection):
        """ 
        Function to Simulate the FGB Spectrum Considering Non-Uniform Strain
        effects
        
        Input:
        
        FBGOriginalWavel: List(Mandatory)
                        List with the FBG array original wavelength
        PhotoElasticParam: float (Mandatory)
                        Photo-Elastic parameter
        InitialRefractiveIndex: float (Mandatory)
                        Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
                        Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
                        Fringe Visibility (Fv)
        MinBandWidth: float (Mandatory)
                        Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
        SimulationResolution: float (Mandatory)
                        Simulation resolution- Wavelength increment 
        FBGDirection: int (mandatory)
                    0-xx 1-yy 2-zz
        """        
        self.FBGOriginalWavel=FBGOriginalWavel
        self.PhotoElasticParam=PhotoElasticParam
        self.InitialRefractiveIndex=InitialRefractiveIndex
        self.MeanChangeRefractiveIndex=MeanChangeRefractiveIndex
        self.FringeVisibility=FringeVisibility
        self.MinBandWidth=MinBandWidth
        self.MaxBandWidth=MaxBandWidth
        self.SimulationResolution=SimulationResolution
        self.FBGDirection=FBGDirection
        
        #Array with the original FBG period
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):
            self.APFBG.append(self.FBGOriginalWavel[i]/(2.0*self.InitialRefractiveIndex))        

        #Empty Reflec spectrum -Uniform Strain
        self.NUSReflect={}
        self.NUSReflect['wavelength']=[]
        self.NUSReflect['reflec']=[]
        
        
        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):    
            #Sections the gratting is divided-- Transfer Matrix
            M=len(self.FBGArray['FBG'+str(i+1)]['x'])
            #FBG increment size (nm)
            deltz=(self.FBGLength*(10.0**6))/M
            
            #----Built the grating period chenged by the Non-uniform strain---
            FBGperiod=[]
            for j in np.arange(0,M):
                if self.FBGDirection==0: #FBG longitudinal direction (xx)
                    FBGperiod.append(self.APFBG[i]*(1+(1-self.PhotoElasticParam)*self.FBGArray['FBG'+str(i+1)]['LE11'][j]))
                if self.FBGDirection==1: #FBG longitudinal direction (yy)
                    FBGperiod.append(self.APFBG[i]*(1+(1-self.PhotoElasticParam)*self.FBGArray['FBG'+str(i+1)]['LE22'][j]))                    
                if self.FBGDirection==2: #FBG longitudinal direction (zz)
                    FBGperiod.append(self.APFBG[i]*(1+(1-self.PhotoElasticParam)*self.FBGArray['FBG'+str(i+1)]['LE33'][j]))

            # Wavelength cycle (Here the simulation resolution is used)
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix 
                #Cycle-Each FBG increment
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*self.InitialRefractiveIndex*((1.0/l)-(1.0/(2.0*self.InitialRefractiveIndex*FBGperiod[z])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                #Add to the Reflection file
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                self.NUSReflect['wavelength'].append(l)
                self.NUSReflect['reflec'].append(REF)
 


    def TransverseStrain(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,FBGDirection,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient):
        """ 
        Function to Simulate the FGB Spectrum Considering Uniform longitudinal strain
        and transverse Stress
        
        Input:
        
        FBGOriginalWavel: List(Mandatory)
                        List with the FBG array original wavelength
        PhotoElasticParam: float (Mandatory)
                        Photo-Elastic parameter
        InitialRefractiveIndex: float (Mandatory)
                        Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
                        Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
                        Fringe Visibility (Fv)
        MinBandWidth: float (Mandatory)
                        Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
        SimulationResolution: float (Mandatory)
                        Simulation resolution- Wavelength increment 
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
        self.MeanChangeRefractiveIndex=MeanChangeRefractiveIndex
        self.FringeVisibility=FringeVisibility
        self.MinBandWidth=MinBandWidth
        self.MaxBandWidth=MaxBandWidth
        self.SimulationResolution=SimulationResolution
        self.FBGDirection=FBGDirection  
        self.DirectionalRefractiveP11=DirectionalRefractiveP11
        self.DirectionalRefractiveP12=DirectionalRefractiveP12
        self.YoungsModule=YoungsModule/10**6 #Change to MPA
        self.PoissionsCoefficient=PoissionsCoefficient

      
        #Array with the FBG period- Longitudinal Strain
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):  
            if self.FBGDirection==0: #FBG longitudinal direction (xx)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE11'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

            if self.FBGDirection==1: #FBG longitudinal direction (yy)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE22'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

            if self.FBGDirection==2: #FBG longitudinal direction (zz)
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE33'] )#Strain average
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*TempMeanFBG) #weavelength at uniform strain
                self.APFBG.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period at strain state

        #Empty Reflec spectrum- Two waves
        self.TSYReflect={}#Y wave
        self.TSYReflect['wavelength']=[]
        self.TSYReflect['reflec']=[]      
      
        self.TSZReflect={}#Z wave
        self.TSZReflect['wavelength']=[]
        self.TSZReflect['reflec']=[]      

        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):    
            #Sections the gratting is divided-- Transfer Matrix
            M=len(self.FBGArray['FBG'+str(i+1)]['x'])
            #FBG increment size (nm)
            deltz=(self.FBGLength*(10.0**6))/M
                      
            #Build a List with the change in the refractive index DeltaNEffY
            self.Dneffy=[]
            self.Dneffz=[]
            for j in np.arange(0,M):
                if self.FBGDirection==0: #FBG longitudinal direction (xx)
                    DirecX='S11'
                    DirecY='S22'
                    DirecZ='S33'
                    self.Dneffy.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecY][j]
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecZ][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
                                      
                    
                    self.Dneffz.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecZ][j]
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecY][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
      
                if self.FBGDirection==1: #FBG longitudinal direction (yy)
                    DirecX='S22'
                    DirecY='S11' #need to be negative
                    DirecZ='S33'
                    self.Dneffy.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*(-self.FBGArray['FBG'+str(i+1)][DirecY][j])
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecZ][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
                                      
                    
                    self.Dneffz.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecZ][j]
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(-self.FBGArray['FBG'+str(i+1)][DirecY][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
                    
                if self.FBGDirection==2: #FBG longitudinal direction (zz)
                    DirecX='S33'
                    DirecY='S22'
                    DirecZ='S11'#Negative
                    self.Dneffy.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecY][j]
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(-self.FBGArray['FBG'+str(i+1)][DirecZ][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
                                      
                    
                    self.Dneffz.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissionsCoefficient*self.DirectionalRefractiveP12)*(-self.FBGArray['FBG'+str(i+1)][DirecZ][j])
                    +((1-self.PoissionsCoefficient)*self.DirectionalRefractiveP12-self.PoissionsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecY][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))

           # Wavelength cycle (Here the simulation resolution is used)
            #YWave
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix 
                #Cycle-Each FBG increment
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*(self.InitialRefractiveIndex+self.Dneffy[z])*((1.0/l)-(1.0/(2.0*(self.InitialRefractiveIndex+self.Dneffy[z])*self.APFBG[i])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                #Add to the Reflection file- YWave
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                self.TSYReflect['wavelength'].append(l) #Output File
                self.TSYReflect['reflec'].append(REF) #Output File

            #ZWave
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix 
                #Cycle-Each FBG increment
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*(self.InitialRefractiveIndex+self.Dneffz[z])*((1.0/l)-(1.0/(2.0*(self.InitialRefractiveIndex+self.Dneffz[z])*self.APFBG[i])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                #Add to the Reflection file- YWave
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                self.TSZReflect['wavelength'].append(l) #Output File
                self.TSZReflect['reflec'].append(REF)  #Output File

    def FBGOutputSum(self,FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,FBGDirection,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient):
        """
          Function to create a table with the wavelength shift and the peak spliting per sensor
          
          Similar calculation to the TR simulation

          Output:
            self.FBGOutSum['FBG'+str(b)]={}
            self.FBGOutSum['FBG'+str(b)]['WaveShift']=[]
            self.FBGOutSum['FBG'+str(b)]['WaveWidth']=[]
          
        """
        self.FBGOriginalWavel=FBGOriginalWavel
        self.PhotoElasticParam=PhotoElasticParam
        self.InitialRefractiveIndex=InitialRefractiveIndex
        self.FBGDirection=FBGDirection  
        self.DirectionalRefractiveP11=DirectionalRefractiveP11
        self.DirectionalRefractiveP12=DirectionalRefractiveP12
        self.YoungsModule=YoungsModule/10**6 #Change to MPA
        self.PoissionsCoefficient=PoissionsCoefficient
        

        #Calculating average strain,a maximum and minumum strain, and average stress  per FBG
        FBGmaxmin={}        
        for b in np.arange(1,self.NumberFBG+1):
            FBGmaxmin['FBG'+str(b)]={}
            FBGmaxmin['FBG'+str(b)]['AV-LE11']=np.mean(self.FBGArray['FBG'+str(b)]['LE11'])
            FBGmaxmin['FBG'+str(b)]['AV-LE22']=np.mean(self.FBGArray['FBG'+str(b)]['LE22'])
            FBGmaxmin['FBG'+str(b)]['AV-LE33']=np.mean(self.FBGArray['FBG'+str(b)]['LE33'])
            FBGmaxmin['FBG'+str(b)]['Max-LE11']=np.max(self.FBGArray['FBG'+str(b)]['LE11'])
            FBGmaxmin['FBG'+str(b)]['Max-LE22']=np.max(self.FBGArray['FBG'+str(b)]['LE22'])
            FBGmaxmin['FBG'+str(b)]['Max-LE33']=np.max(self.FBGArray['FBG'+str(b)]['LE33'])
            FBGmaxmin['FBG'+str(b)]['Min-LE11']=np.min(self.FBGArray['FBG'+str(b)]['LE11'])
            FBGmaxmin['FBG'+str(b)]['Min-LE22']=np.min(self.FBGArray['FBG'+str(b)]['LE22'])
            FBGmaxmin['FBG'+str(b)]['Min-LE33']=np.min(self.FBGArray['FBG'+str(b)]['LE33'])
            FBGmaxmin['FBG'+str(b)]['AV-S11']=np.mean(self.FBGArray['FBG'+str(b)]['S11'])
            FBGmaxmin['FBG'+str(b)]['AV-S22']=np.mean(self.FBGArray['FBG'+str(b)]['S22'])
            FBGmaxmin['FBG'+str(b)]['AV-S33']=np.mean(self.FBGArray['FBG'+str(b)]['S33'])           
        
        

        #Directions selection
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

        # Dictionary to store the data 
        self.FBGOutSum={}
        for b in np.arange(1,self.NumberFBG+1):
            self.FBGOutSum['FBG'+str(b)]={}
            self.FBGOutSum['FBG'+str(b)]['WaveShift']=[]
            self.FBGOutSum['FBG'+str(b)]['WaveWidth']=[]

        for b in np.arange(0,self.NumberFBG): #Cycle each FBG
  
            #--------Wavelength Shift Calculation----------------------
            WavlShift=(1.0-self.PhotoElasticParam)*self.FBGOriginalWavel[b]*FBGmaxmin['FBG'+str(b+1)][AVStrainDirection]
            self.FBGOutSum['FBG'+str(b+1)]['WaveShift'].append(WavlShift)

            #---Variation of the peak width: Non-uniform strain contribution----
            #Grating period caused by maximum strain
            graperiodmax=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][MaxStrainDirection])
            #Grating period caused by minimum strain
            graperiodmin=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][MinStrainDirection])
            #Peak width variation
            PeakWV1=2*self.InitialRefractiveIndex*(graperiodmax-graperiodmin)
            
            #---Variation of the peak width: Transverse stress contribution----
            PeakWV2=abs(FBGmaxmin['FBG'+str(b+1)][AVStresszzDirection]-FBGmaxmin['FBG'+str(b+1)][AVStressyyDirection])*(self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex))*fce
            
            #Add Both contributions
            PeakWV=PeakWV1+PeakWV2
            self.FBGOutSum['FBG'+str(b+1)]['WaveWidth'].append(PeakWV)
                
