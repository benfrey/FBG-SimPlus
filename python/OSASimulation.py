"""
Simulatation of reflected FBG spectrum using coupled-mode theory

Copyright (C) 2021 Ben Frey

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
Email: freynben@gmail.com
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
        names =('x','LE11','LE22','LE33','S11','S22','S33','T')
        formats=('f8','f8','f8','f8','f8','f8','f8','f8')
        dtypes = {'names' : names,'formats': formats}
        self.RawData=np.loadtxt(self.filename, dtype=dtypes,skiprows=self.SkipRow)
        print("Printing head of RawData from input file...")
        print(self.RawData[0:5])

        """------------------------------------------------------------------------------
        Creating FBG array, and give them the correspondent Data
        FBGArray['xx']['yy'][zz]
        xx- Number of the FBG-1
        yy- variable 'x','LE11','LE22','LE33','S11','S22','S33','T'
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
            self.FBGArray['FBG'+str(b)]['T']=[]

        #Converting tolerance and FBG length to the SI if in meters
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength/1000.0
            self.Tolerance=self.Tolerance/1000.0

        #Sorting Data for each FBG
        for b in np.arange(0,self.NumberFBG):
            for f in np.arange(0,len(self.RawData)):
                #Check the lines inside the FBG length+tolerance
                #if self.RawData['x'][f]>self.FBGPosition[b]-self.Tolerance and self.RawData['x'][f]<self.FBGPosition[b]+self.FBGLength+self.Tolerance:
                    #Convert to SI, previosly (mm,MPA)
                    if self.InputUnits==0:
                        #print("Conversion to non-SI units triggered!")
                        self.FBGArray['FBG'+str(b+1)]['x'].append(self.RawData['x'][f]*1000)
                        self.FBGArray['FBG'+str(b+1)]['LE11'].append(self.RawData['LE11'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE22'].append(self.RawData['LE22'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE33'].append(self.RawData['LE33'][f])
                        self.FBGArray['FBG'+str(b+1)]['S11'].append(self.RawData['S11'][f]/(10**6))
                        self.FBGArray['FBG'+str(b+1)]['S22'].append(self.RawData['S22'][f]/(10**6))
                        self.FBGArray['FBG'+str(b+1)]['S33'].append(self.RawData['S33'][f]/(10**6))
                        self.FBGArray['FBG'+str(b+1)]['T'].append(self.RawData['T'][f])

                    else:
                        #print("Using standard units!")
                        self.FBGArray['FBG'+str(b+1)]['x'].append(self.RawData['x'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE11'].append(self.RawData['LE11'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE22'].append(self.RawData['LE22'][f])
                        self.FBGArray['FBG'+str(b+1)]['LE33'].append(self.RawData['LE33'][f])
                        self.FBGArray['FBG'+str(b+1)]['S11'].append(self.RawData['S11'][f])
                        self.FBGArray['FBG'+str(b+1)]['S22'].append(self.RawData['S22'][f])
                        self.FBGArray['FBG'+str(b+1)]['S33'].append(self.RawData['S33'][f])
                        self.FBGArray['FBG'+str(b+1)]['T'].append(self.RawData['T'][f])

        #Converting tolerance and FBG length back to mm
        if self.InputUnits==0:
            self.FBGLength=self.FBGLength*1000.0
            self.Tolerance=self.Tolerance*1000.0

    def UndeformedFBG(self,SimulationResolution,MinBandWidth,MaxBandWidth,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,PoissonsCoefficient,FBGOriginalWavel):
        """
        
        Parameters
        ----------
        SimulationResolution: float (Mandatory)
            Simulation resolution- Wavelength increment
        MinBandWidth: float (Mandatory)
            Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
            Light Max. Bandwidth
        InitialRefractiveIndex: float (Mandatory)
            Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
            Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
            Fringe Visibility (FV)
        DirectionalRefractiveP11 : float (Mandatory)
            Pockel’s normal photoelastic constant
        DirectionalRefractiveP12 : float (Mandatory)
            Pockel’s shear photoelastic constant
        PoissonsCoefficient : float (Mandatory)
            Poisson's Coefficient of fiber
        FBGOriginalWavel: array (Mandatory)
            Original wavelengths of FBG regions

        Returns
        -------
        None.

        """
        self.SimulationResolution = SimulationResolution
        self.MinBandWidth = MinBandWidth
        self.MaxBandWidth = MaxBandWidth
        self.InitialRefractiveIndex = InitialRefractiveIndex
        self.MeanChangeRefractiveIndex = MeanChangeRefractiveIndex
        self.FringeVisibility = FringeVisibility
        self.DirectionalRefractiveP11 = DirectionalRefractiveP11
        self.DirectionalRefractiveP12 = DirectionalRefractiveP12
        self.PoissonsCoefficient = PoissonsCoefficient
        self.FBGOriginalWavel = FBGOriginalWavel

        #Array with the Original FBG period
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):
            self.APFBG.append(self.FBGOriginalWavel[i]/(2.0*self.InitialRefractiveIndex))

        #Empty Original Reflec spectrum
        self.OReflect={}
        self.OReflect['wavelength']=[]
        self.OReflect['reflec']=[]

        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):
            # Wavelength cycle (Here the simulation resolution is used)
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                #print("---------- WAVELENGTH: "+str(l))
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix
                M=20.0 #Sections the gratting is divided-- Transfer Matrix
                #FBG increment size (nm)
                deltz=(self.FBGLength*(10.0**6))/M
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*self.InitialRefractiveIndex*((1.0/l)-(1.0/(2.0*self.InitialRefractiveIndex*self.APFBG[i])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
                    #print("Sigma: "+str(sig))
                    #print("Lamb: "+str((2.0*self.InitialRefractiveIndex*self.APFBG[i])))
                    #Kaa Function
                    kaa=math.pi*self.FringeVisibility*self.MeanChangeRefractiveIndex/l
                    #print("Kaa: "+str(kaa))
                    #Gamma Function
                    gammab=cmath.sqrt(kaa**2.0-sig**2.0)
                    #print("Gammab: "+str(gammab))

                    #Transfer Matrix Calculation
                    f11=complex(cmath.cosh(gammab*deltz),-(sig/gammab)*cmath.sinh(gammab*deltz))
                    f22=complex(cmath.cosh(gammab*deltz),(sig/gammab)*cmath.sinh(gammab*deltz))
                    f12=complex(0,-(kaa/gammab)*cmath.sinh(gammab*deltz))
                    f21=complex(0,+(kaa/gammab)*cmath.sinh(gammab*deltz))
                    #print("TM: "+str(gammab))

                    #Assemble Transfer Matrix
                    f1=np.dot(f1,np.matrix([[f11, f12],[ f21, f22]]))
                    #print("Transfer Matrix "+str(z)+": "+str(f1))
                PO=f1[0,0]
                NO=f1[1,0]
                REF=abs(NO/PO)**2
                #print("Reflectivity "+str(REF))

                self.OReflect['wavelength'].append(l)
                self.OReflect['reflec'].append(REF)

        print(self.OReflect['wavelength'][0:15])
        print(self.OReflect['reflec'][0:15])

    def DeformedFBG(self,SimulationResolution,MinBandWidth,MaxBandWidth,AmbientTemperature,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissonsCoefficient,ThermoOptic,StrainType,StressType,EmulateTemperature,FiberThermalExpansionCoefficient,HostThermalExpansionCoefficient,FBGOriginalWavel):
        """

        Parameters
        ----------

        SimulationResolution: float (Mandatory)
            Simulation resolution- Wavelength increment
        MinBandWidth: float (Mandatory)
            Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
            Light Max. Bandwidth
        AmbientTemperature : float (Mandatory)
            Base in which our reference temperature is set.
        InitialRefractiveIndex: float (Mandatory)
            Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
            Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
            Fringe Visibility (FV)
        DirectionalRefractiveP11 : float (Mandatory)
            Pockel’s normal photoelastic constant
        DirectionalRefractiveP12 : float (Mandatory)
            Pockel’s shear photoelastic constant
        YoungsModule : float (Mandatory)
            Young's modulus of fiber
        PoissonsCoefficient : float (Mandatory)
            Poisson's Coefficient of fiber
        ThermoOptic : float (Mandatory)
            Thermo optic coefficient of fiber
        StrainType : int (Mandatory)
            0 for none, 1 for uniform, 2 for non-uniform
        StressType : int (Mandatory)
            0 for none, 1 for included
        EmulateTemperature : float (Mandatory)
            Theoretical emulated temperature of fiber
        FiberThermalExpansionCoefficient : float (Mandatory)
            Thermal expansion coefficient of fiber material
        HostThermalExpansionCoefficient : float (Mandatory)
            Thermal expansion coefficient of host material
        FBGOriginalWavel: array (Mandatory)
            Original wavelengths of FBG regions            

        Returns
        -------
        None.

        """
        
        self.SimulationResolution = SimulationResolution
        self.MinBandWidth = MinBandWidth
        self.MaxBandWidth = MaxBandWidth
        self.AmbientTemperature = AmbientTemperature
        self.InitialRefractiveIndex = InitialRefractiveIndex
        self.MeanChangeRefractiveIndex = MeanChangeRefractiveIndex
        self.FringeVisibility = FringeVisibility
        self.DirectionalRefractiveP11 = DirectionalRefractiveP11
        self.DirectionalRefractiveP12 = DirectionalRefractiveP12
        self.YoungsModule = YoungsModule
        self.PoissonsCoefficient = PoissonsCoefficient
        self.ThermoOptic = ThermoOptic
        self.EmulateTemperature = EmulateTemperature
        self.FiberThermalExpansionCoefficient = FiberThermalExpansionCoefficient
        self.HostThermalExpansionCoefficient = HostThermalExpansionCoefficient
        self.FBGOriginalWavel = FBGOriginalWavel

        #Calculate photoelastic coef from directional coefs.
        self.PhotoElasticParam = (self.InitialRefractiveIndex**2/2)*(self.DirectionalRefractiveP12-self.PoissonsCoefficient*(self.DirectionalRefractiveP11-self.DirectionalRefractiveP12))

        #Determine if we need to emulate a theoretical temperature value
        if self.EmulateTemperature != -1.0:
            for i in np.arange(0,self.NumberFBG):
                for j in range(len(self.FBGArray['FBG'+str(i+1)]['T'])):
                    self.FBGArray['FBG'+str(i+1)]['T'][j] = self.EmulateTemperature
                print(self.FBGArray['FBG'+str(i+1)]['T'][0:5])
            
        #Array with the original FBG periods
        self.APFBG=[]
        for i in np.arange(0,self.NumberFBG):
            self.APFBG.append(self.FBGOriginalWavel[i]/(2.0*self.InitialRefractiveIndex))

        #Empty Reflec spectrum- Two waves (individual contributions to account for any transverse stress)
        self.YReflect={}#Y wave
        self.YReflect['wavelength']=[]
        self.YReflect['reflec']=[]

        self.ZReflect={}#Z wave
        self.ZReflect['wavelength']=[]
        self.ZReflect['reflec']=[]
        
        #Empty Reflec spectrum- Composite wave of Y and Z contributions
        self.DReflect={}#Composite wave
        self.DReflect['wavelength']=[]
        self.DReflect['reflec']=[]       

        #Cycle all the FBG sensors
        for i in np.arange(0,self.NumberFBG):
            #Sections the gratting is divided-- Transfer Matrix
            M=len(self.FBGArray['FBG'+str(i+1)]['x'])
            #FBG increment size (nm)
            deltz=(self.FBGLength*(10.0**6))/M

            #--- STRAIN ---
            
            FBGperiod=[]
            if StrainType == 0:
                #--- Case of "no longitudinal strain" ---
                for j in np.arange(0,M):
                    FBGperiod.append(self.APFBG[i])
            
            elif StrainType == 1:
                #--- Case of "uniform longitudinal strain" ---
                StrainMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['LE11']) #Strain average over FBG region
                TempMeanFBG=np.mean(self.FBGArray['FBG'+str(i+1)]['T']) #Temperature average over FBG region
                TempnewWavelength=self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*StrainMeanFBG+(self.FiberThermalExpansionCoefficient+(1-self.PhotoElasticParam)*(self.HostThermalExpansionCoefficient-self.FiberThermalExpansionCoefficient)+self.ThermoOptic)*(TempMeanFBG-self.AmbientTemperature)) #weavelength at uniform strain and temperature
                for j in np.arange(0,M):
                    FBGperiod.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period 

            else:
                #--- Case of "non-uniform longitudinal strain" build the grating period changed ---
                for j in np.arange(0,M):
                    TempnewWavelength = self.FBGOriginalWavel[i]*(1+(1-self.PhotoElasticParam)*self.FBGArray['FBG'+str(i+1)]['LE11'][j]+(self.FiberThermalExpansionCoefficient+(1-self.PhotoElasticParam)*(self.HostThermalExpansionCoefficient-self.FiberThermalExpansionCoefficient)+self.ThermoOptic)*(self.FBGArray['FBG'+str(i+1)]['T'][j]-self.AmbientTemperature)) #weavelength at nonuniform strain and temperature
                    #print(self.FBGArray['FBG'+str(i+1)]['T'][0]-self.AmbientTemperature)
                    FBGperiod.append(TempnewWavelength/(2.0*self.InitialRefractiveIndex)) # Grating period 

            #--- STRESS ---
            
            #Build a List with the change in the refractive index DeltaNEffY
            self.Dneffy=[]
            self.Dneffz=[]
            
            if StressType == 0:
                #--- Case of "no transverse stress" ---
                for j in np.arange(0,M):
                     self.Dneffy.append(0)
                     self.Dneffz.append(0)
            else:
                #--- Case of "included transverse stress" ---
                for j in np.arange(0,M):
                    DirecX='S11'
                    DirecY='S22'
                    DirecZ='S33'
                    self.Dneffy.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissonsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecY][j]
                    +((1-self.PoissonsCoefficient)*self.DirectionalRefractiveP12-self.PoissonsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecZ][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))

                    self.Dneffz.append(-(self.InitialRefractiveIndex**3.0)/(2*self.YoungsModule)
                    *((self.DirectionalRefractiveP11-2*self.PoissonsCoefficient*self.DirectionalRefractiveP12)*self.FBGArray['FBG'+str(i+1)][DirecZ][j]
                    +((1-self.PoissonsCoefficient)*self.DirectionalRefractiveP12-self.PoissonsCoefficient*self.DirectionalRefractiveP11)
                    *(self.FBGArray['FBG'+str(i+1)][DirecY][j]+self.FBGArray['FBG'+str(i+1)][DirecX][j])))
            
            #--- SIMULATION --- 
            
            #YWave
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution): # Wavelength cycle (Here the simulation resolution is used)
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix
                #Cycle-Each FBG increment
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*(self.InitialRefractiveIndex+self.Dneffy[z])*((1.0/l)-(1.0/(2.0*(self.InitialRefractiveIndex+self.Dneffy[z])*FBGperiod[z])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
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
                self.YReflect['wavelength'].append(l) #Output File
                self.YReflect['reflec'].append(REF) #Output File

            #ZWave
            for l in np.arange(self.MinBandWidth,self.MaxBandWidth,self.SimulationResolution):
                f1 = np.matrix('1 0; 0 1') #empty transfer matrix
                #Cycle-Each FBG increment
                for z in np.arange(0,M):
                    #Sigma Function
                    sig=2.0*math.pi*(self.InitialRefractiveIndex+self.Dneffz[z])*((1.0/l)-(1.0/(2.0*(self.InitialRefractiveIndex+self.Dneffz[z])*FBGperiod[z])))+(2*math.pi*self.MeanChangeRefractiveIndex/l)
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
                self.ZReflect['wavelength'].append(l) #Output File
                self.ZReflect['reflec'].append(REF)  #Output File

        #Halve the amplitude of each of the y and z waves, then sum to find composite wave.
        self.DReflect['wavelength'] = self.YReflect['wavelength']
        self.YReflect['reflec'] = np.divide(self.YReflect['reflec'],2.0)
        self.ZReflect['reflec'] = np.divide(self.ZReflect['reflec'],2.0)
        self.DReflect['reflec'] = np.add(self.YReflect['reflec'],self.ZReflect['reflec'])

        #print(self.DReflect['wavelength'][0:15])
        #print(self.YReflect['reflec'][0:15])

    def FBGOutputSum(self,AmbientTemperature,InitialRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissonsCoefficient,ThermoOptic,StrainType,StressType,EmulateTemperature,FiberThermalExpansionCoefficient,HostThermalExpansionCoefficient,FBGOriginalWavel):
        """
        Function to create a table with the wavelength shift and the peak spliting per sensor

        Similar calculation to the TR simulation

        Parameters
        ----------

        SimulationResolution: float (Mandatory)
            Simulation resolution- Wavelength increment
        MinBandWidth: float (Mandatory)
            Light Min. Bandwidth
        MaxBandWidth: float (Mandatory)
            Light Max. Bandwidth
        AmbientTemperature : float (Mandatory)
            Base in which our reference temperature is set.
        InitialRefractiveIndex: float (Mandatory)
            Initial effective refractive index (neff)
        MeanChangeRefractiveIndex: float (Mandatory)
            Mean induced change in the refractive index (dneff)
        FringeVisibility: float (Mandatory)
            Fringe Visibility (FV)
        DirectionalRefractiveP11 : float (Mandatory)
            Pockel’s normal photoelastic constant
        DirectionalRefractiveP12 : float (Mandatory)
            Pockel’s shear photoelastic constant
        YoungsModule : float (Mandatory)
            Young's modulus of fiber
        PoissonsCoefficient : float (Mandatory)
            Poisson's Coefficient of fiber
        ThermoOptic : float (Mandatory)
            Thermo optic coefficient of fiber
        StrainType : int (Mandatory)
            0 for none, 1 for uniform, 2 for non-uniform
        StressType : int (Mandatory)
            0 for none, 1 for included
        EmulateTemperature : float (Mandatory)
            Theoretical emulated temperature of fiber
        FiberThermalExpansionCoefficient : float (Mandatory)
            Thermal expansion coefficient of core and cladding material
        HostThermalExpansionCoefficient : float (Mandatory)
            Thermal expansion coefficient of host material
        FBGOriginalWavel: array (Mandatory)
            Original wavelengths of FBG regions            

        Returns
        -------
        None.

        Output
        -------
        self.FBGOutSum['FBG'+str(b)]={}
        self.FBGOutSum['FBG'+str(b)]['WaveShift']=[]
        self.FBGOutSum['FBG'+str(b)]['WaveWidth']=[]
        
        """
        self.AmbientTemperature = AmbientTemperature
        self.InitialRefractiveIndex = InitialRefractiveIndex
        self.DirectionalRefractiveP11 = DirectionalRefractiveP11
        self.DirectionalRefractiveP12 = DirectionalRefractiveP12
        self.YoungsModule = YoungsModule
        self.PoissonsCoefficient = PoissonsCoefficient
        self.ThermoOptic = ThermoOptic
        self.EmulateTemperature = EmulateTemperature
        self.FiberThermalExpansionCoefficient = FiberThermalExpansionCoefficient
        self.HostThermalExpansionCoefficient = HostThermalExpansionCoefficient
        self.FBGOriginalWavel = FBGOriginalWavel
        self.FBGDirection = 0
        
        #Calculate photoelastic coef from directional coefs.
        self.PhotoElasticParam = (self.InitialRefractiveIndex**2/2)*(self.DirectionalRefractiveP12-self.PoissonsCoefficient*(self.DirectionalRefractiveP11-self.DirectionalRefractiveP12))

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
            FBGmaxmin['FBG'+str(b)]['AV-T']=np.mean(self.FBGArray['FBG'+str(b)]['T'])

        #Directions selection
        if self.FBGDirection==0: #FBG longitudinal direction (xx)
            AVStrainDirection='AV-LE11'
            MaxStrainDirection='Max-LE11'
            MinStrainDirection='Min-LE11'
            AVStresszzDirection='AV-S33'
            AVStressyyDirection='AV-S22'
            AVTemperature='AV-T'
        
        #Determine if we need to emulate a theoretical temperature value
        if self.EmulateTemperature != -1.0:
            for i in np.arange(0,self.NumberFBG):
                FBGmaxmin['FBG'+str(b)][AVTemperature]=self.EmulateTemperature
                    
        #fixed component: Transverse stress"""
        fce=(((1+self.PoissonsCoefficient)*self.DirectionalRefractiveP12-(1+self.PoissonsCoefficient)*self.DirectionalRefractiveP11)*self.InitialRefractiveIndex**3)/self.YoungsModule

        # Dictionary to store the data
        self.FBGOutSum={}
        for b in np.arange(1,self.NumberFBG+1):
            self.FBGOutSum['FBG'+str(b)]={}
            self.FBGOutSum['FBG'+str(b)]['WaveShift']=[]
            self.FBGOutSum['FBG'+str(b)]['WaveWidth']=[]

        for b in np.arange(0,self.NumberFBG): #Cycle each FBG

        
            #--- STRAIN ---
            if StrainType == 0:
                #--- Case of "no longitudinal strain" ---
                
                #--------Wavelength Shift Calculation----------------------
                WavlShift=self.FBGOriginalWavel[b]*((self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature)) #weavelength due to temp
                self.FBGOutSum['FBG'+str(b+1)]['WaveShift'].append(WavlShift)   
                
                self.graperiodmax = self.graperiodmin = self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature))
                PeakWV1=0
                
            elif StrainType == 1:
                #--- Case of "uniform longitudinal strain" 
                
                #--------Wavelength Shift Calculation----------------------
                WavlShift=self.FBGOriginalWavel[b]*((1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][AVStrainDirection]+(self.FiberThermalExpansionCoefficient+(1-self.PhotoElasticParam)*(self.HostThermalExpansionCoefficient-self.FiberThermalExpansionCoefficient)+self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature)) #weavelength at uniform strain and temperature
                self.FBGOutSum['FBG'+str(b+1)]['WaveShift'].append(WavlShift)   
                
                self.graperiodmax = self.graperiodmin = self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][AVStrainDirection]+(self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature))
                PeakWV1=0
                
            else:
                #--- Case of "non-uniform longitudinal strain"---
         
                #--------Wavelength Shift Calculation----------------------
                WavlShift=self.FBGOriginalWavel[b]*((1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][AVStrainDirection]+(self.FiberThermalExpansionCoefficient+(1-self.PhotoElasticParam)*(self.HostThermalExpansionCoefficient-self.FiberThermalExpansionCoefficient)+self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature)) #weavelength at uniform strain and temperature
                self.FBGOutSum['FBG'+str(b+1)]['WaveShift'].append(WavlShift)   
                
                #--- Variation of the peak width: Non-uniform strain contribution and temperature ----
                #Grating period caused by maximum strain
                self.graperiodmax=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][MaxStrainDirection]+(self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature))
                #Grating period caused by minimum strain
                self.graperiodmin=self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex)*(1+(1-self.PhotoElasticParam)*FBGmaxmin['FBG'+str(b+1)][MinStrainDirection]+(self.ThermoOptic)*(FBGmaxmin['FBG'+str(b+1)][AVTemperature]-self.AmbientTemperature))
                #Peak width variation
                PeakWV1=2*self.InitialRefractiveIndex*(self.graperiodmax-self.graperiodmin)


            #--- STRESS ---
            if StressType == 0:
                #--- Case of "no transverse stress" ---
                PeakWV2=0
            else:
                #---Variation of the peak width: Transverse stress contribution----
                PeakWV2=abs(FBGmaxmin['FBG'+str(b+1)][AVStresszzDirection]-FBGmaxmin['FBG'+str(b+1)][AVStressyyDirection])*(self.FBGOriginalWavel[b]/(2*self.InitialRefractiveIndex))*fce

            #Add Both contributions
            PeakWV=PeakWV1+PeakWV2
            self.FBGOutSum['FBG'+str(b+1)]['WaveWidth'].append(PeakWV)
            
            print("FBG "+str(b)+ " | Grating period max: "+str(self.graperiodmax)+" | Grating period min: "+str(self.graperiodmin)+" | PeakWV1: "+str(PeakWV1)+" | PeakWV2: "+str(PeakWV2))
