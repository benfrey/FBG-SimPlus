""" Python Class To Extract Stress and Strain Along a path in ABAQUS.

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
import os
import sys
import time
import numpy as np
import fnmatch
import subprocess
import shutil

class ExtractAbaqus(object):
    def __init__(self, Odbpath,Writefolder,PathCoordinates,RotateAxis,DisplacementStep,SpecificDisplacement,UserDefinedCoordinate,VectorNewCoordinate):
        """ Initialized the classe using the OdbPath

        Parameters:
        ----------
        Odbpath: string (Mandatory)
                 Path to Abaqus ODB file
        Writefolder: string (Mandatory)
                Path where the output will be writen
        PathCoordinates: Diction. (Mandatory)
                Dict. with the paths coordinates                          
        RotateAxis: Intg (Mandatory)
                0- Default; 1- User-Defined; 2- Vector with new coordinates
        DisplacementStep:Intg (Mandatory)
                0-all displacement; 1- Specific displacement step
        SpecificDisplacement:Intg (Mandatory)
                Number of the specific step time to extract the result
        UserDefinedCoordinate: String (Mandatory)
                String with the name of the user defined coordinatesystem
        VectorNewCoordinate: Dict (Mandatory)
                Vector with new coordinate system direction
                VectorNewCoordinate['0'],VectorNewCoordinate['1'],VectorNewCoordinate['2'],
        """
        self.Odbpath = Odbpath
        self.Writefolder=Writefolder.replace("\\","/",99)
        self.PathCoordinates=PathCoordinates
        self.RotateAxis=RotateAxis
        self.DisplacementStep=DisplacementStep
        self.SpecificDisplacement=SpecificDisplacement
        self.UserDefinedCoordinate=UserDefinedCoordinate
        self.VectorNewCoordinate=VectorNewCoordinate

        
    def createPyFile(self):
        """
        This function creates the Python file that is submited to Abaqus.
        
        Output file Format (.txt):
        Path Position(distance); LE11(Strain); LE22(Strain); LE33(Strain); s11(Stress);s22(Stress);s33(Stress) 
        """
        os.chdir(self.Writefolder) #Change to dir. to output folder
        if os.path.isdir("temp")==False: #Creates a temp directory
            os.makedirs("temp")
        AbaqusFile=open("temp\AbaqusSubmitFile.py","w+")
        #-------Fixed Text (Texted required for abaqus to run the file---------
        AbaqusFile.write("from abaqus import * \nfrom abaqusConstants import *\n")
        AbaqusFile.write("session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=30.0,height=30.0)\n")
        AbaqusFile.write("session.viewports['Viewport: 1'].makeCurrent()\n")
        AbaqusFile.write("session.viewports['Viewport: 1'].maximize()\n")
        AbaqusFile.write("from viewerModules import *\n")
        AbaqusFile.write("from driverUtils import executeOnCaeStartup\n")
        AbaqusFile.write("executeOnCaeStartup()\n")
                
        #-----------The ODB Path is inserted---------------------------------
        #Because Abaqus needs the path with rigth slash /, the next line will replace it
        NewPath=self.Odbpath.replace("\\","/",99)
        AbaqusFile.write("o2 = session.openOdb(name='"+str(NewPath)+"')\n")
        AbaqusFile.write("session.viewports['Viewport: 1'].setValues(displayedObject=o2) \n \n")
        #-----------Coordinate system Change----------------------------------
        #self.RotateAxis=0 do nothing- this is the default
        if self.RotateAxis == 1: #1- User-Defined
            AbaqusFile.write("session.viewports['Viewport: 1'].setValues(displayedObject=o2)\n")
            AbaqusFile.write("dtm = session.odbs['"+str(NewPath)+"'].rootAssembly.datumCsyses['"+str(self.UserDefinedCoordinate)+"']\n")
            AbaqusFile.write("session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, rigidTransformPrimary=True, rigidTransformDeformed=True, datumCsys=dtm)\n")
        if self.RotateAxis == 2: #2- Vector with new coordinates
            AbaqusFile.write("session.viewports['Viewport: 1'].setValues(displayedObject=o2)\n")
            AbaqusFile.write("odb = session.odbs['"+str(NewPath)+"']\n")
            AbaqusFile.write("scratchOdb = session.ScratchOdb(odb)\n")
            AbaqusFile.write("scratchOdb.rootAssembly.DatumCsysByThreePoints(name='Newaxis',coordSysType=CARTESIAN, origin="+str(self.VectorNewCoordinate['0'])+", point1="+str(self.VectorNewCoordinate['1'])+",point2="+str(self.VectorNewCoordinate['2'])+")\n")
            AbaqusFile.write("dtm = session.scratchOdbs['"+str(NewPath)+"'].rootAssembly.datumCsyses['Newaxis']\n")
            AbaqusFile.write("session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, rigidTransformPrimary=True,rigidTransformDeformed=True, datumCsys=dtm)\n")
            
        #-------------------Create Path---------------------------------------
        for i in np.arange(0,len(self.PathCoordinates)):
            PathName='Path: '+str(i+1)
            AbaqusFile.write("session.Path(name='"+str(PathName)+"', type=POINT_LIST, expression=(")
            for j in np.arange(0,len(self.PathCoordinates[str(i+1)])):
                AbaqusFile.write("("+str(self.PathCoordinates[str(i+1)][j][0])+','+str(self.PathCoordinates[str(i+1)][j][1])+','+str(self.PathCoordinates[str(i+1)][j][2])+'),')
            AbaqusFile.write("))\n \n")
                
        #Abaqus will create a varible (numberframes) with the number of steps in the model
        AbaqusFile.write("Stepname=str(o2.steps.keys()) \n")
        AbaqusFile.write("Stepname=Stepname.replace(\"['\",\"\") \n")
        AbaqusFile.write("Stepname=Stepname.replace(\"']\",\"\") \n")
        AbaqusFile.write("numberframes=len(o2.steps[Stepname].frames) \n \n")
        
        #Extract in a specific Displacement
        if self.DisplacementStep==1:
            #Select the specific Displacement
            AbaqusFile.write("session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame="+str(self.SpecificDisplacement)+") \n")
            #Create Var. with the setep Number
            AbaqusFile.write("DisplacementStep="+str(self.SpecificDisplacement)+'\n')
            AbaqusFile.write("for i in range(0,"+str(len(self.PathCoordinates))+"):\n") #Extract for all the paths
            #S11
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S11')) \n")
            AbaqusFile.write("\tpth = session.paths['Path: '+str(i+1)]\n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='s11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #S22
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S22')) \n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='s22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #S33
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S33')) \n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='s33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE1
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE11')) \n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='LE11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE2
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE22')) \n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='LE22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE3
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE33')) \n")
            AbaqusFile.write("\tsession.XYDataFromPath(name='LE33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            
            #Write in file
            AbaqusFile.write("\tx0 = session.xyDataObjects['LE11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\tx1 = session.xyDataObjects['LE22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\tx2 = session.xyDataObjects['LE33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\tx3 = session.xyDataObjects['s11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\tx4 = session.xyDataObjects['s22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\tx5 = session.xyDataObjects['s33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            
            AbaqusFile.write("\tsession.xyReportOptions.setValues(interpolation=ON) \n")
            AbaqusFile.write("\tsession.writeXYReport(fileName='"+str(self.Writefolder)+"/'+'Path'+str(i+1)+'_Step_'+str(DisplacementStep).zfill(4)+'.txt',appendMode=OFF, xyData=(x0, x1, x2, x3, x4, x5)) \n")
        #Extract for all the Displacements increments 
        if self.DisplacementStep==0: 
            #Cycle all the increments
            AbaqusFile.write("for l in range(0,numberframes): \n")
            AbaqusFile.write("\tsession.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=l) \n")
            #Create Var. with the setep Number
            AbaqusFile.write("\tDisplacementStep=str(l) \n")
            
            AbaqusFile.write("\tfor i in range(0,"+str(len(self.PathCoordinates))+"):\n") #Extract for all the paths
            #S11
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S11')) \n")
            AbaqusFile.write("\t\tpth = session.paths['Path: '+str(i+1)]\n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='s11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #S22
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S22')) \n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='s22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #S33
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S33')) \n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='s33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE1
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE11')) \n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='LE11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE2
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE22')) \n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='LE22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            #LE3
            AbaqusFile.write("\t\tsession.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'LE33')) \n")
            AbaqusFile.write("\t\tsession.XYDataFromPath(name='LE33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep), path=pth, includeIntersections=True, pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED,labelType=TRUE_DISTANCE) \n")
            
            #Write in file
            AbaqusFile.write("\t\tx0 = session.xyDataObjects['LE11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\t\tx1 = session.xyDataObjects['LE22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\t\tx2 = session.xyDataObjects['LE33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\t\tx3 = session.xyDataObjects['s11_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\t\tx4 = session.xyDataObjects['s22_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            AbaqusFile.write("\t\tx5 = session.xyDataObjects['s33_'+'Path'+str(i+1)+'_Step:'+str(DisplacementStep)] \n")
            
            AbaqusFile.write("\t\tsession.xyReportOptions.setValues(interpolation=ON) \n")
            AbaqusFile.write("\t\tsession.writeXYReport(fileName='"+str(self.Writefolder)+"/'+'Path'+str(i+1)+'_Step_'+str(DisplacementStep).zfill(4)+'.txt',appendMode=OFF, xyData=(x0, x1, x2, x3, x4, x5)) \n")
    
        AbaqusFile.close() #Close the file
        
         
    def submitAbaqus(self):  
          """
          This function submit the python file that was generated before to ABAQUS.
          Options: With or Without GUI
          """
          #NoGUI
          #process = subprocess.call('abaqus cae noGUI=temp\AbaqusSubmitFile.py',shell=True)
          #WithGUI
          process = subprocess.call('abaqus cae script=temp\AbaqusSubmitFile.py',shell=True)
          try:
              os.remove('abaqus.rpy')
          except:
              pass
          try:
              os.remove('abaqus_acis.log')
          except:
              pass
          try:
              shutil.rmtree(str(self.Writefolder)+'/temp')
          except:
              pass
        