# -*- coding: utf-8 -*-
""" This Python Class contains all the functions needed to build the GUI
interface, and to test all the INPUT parameters.

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
from __future__ import print_function
import numpy as np
import os.path as path
import scipy as sp
import sympy
import matplotlib.pyplot as plt
import time
import webbrowser
import sys

#Other Classes used
import GUI.MyPlotMainWindowUI
from OSASimulation import *
from TRSimulation import *
from QtGuiLoader import QtMainWindowLoader, QtWidgetLoader, QtDialogLoader
#Class to connect the UI to the Python Script
from PyQt5 import QtCore, QtGui, QtWidgets

#----------------Main Window Class---------------------------------------------
class MyPlotMainWindow(QtMainWindowLoader):
    """
    This is the Main Window UI Class.

    _Init_- It converts the QTdesign file in a Python Script using the class
    QTGUILoader. It connects all the actions in the UI with the Python Code.
    Also it creates the internal variables that will be submitted later to the
    other developed tools.
    """
    def __init__(self):
        """
        Initiation of the Main Window
        """
        self.module = GUI.MyPlotMainWindowUI
        try:self.ui = self.module.Ui_Form() # Enables autocompletion
        except: pass
        QtMainWindowLoader.__init__(self, self.module)

        #Start the MessageBoard used to plot errors
        self.ui.MessageBoard.insertPlainText('>>')

        #Empty Dict. with path coordinates (Extract Stress/strain)
        self.PathCoordinates={}
        self.PathNumber=0

        #Empty list with FBG position and FBG array original wavelength-OSA
        self.FBGPosition=[]
        self.FBGOriginalWavel=[]
        #Empty list with FBG position and FBG array original wavelength-TR
        self.FBGPositionTR=[]
        self.FBGOriginalWavelTR=[]

    """--------------------------Tool Bar Buttons----------------------"""
    #Push Button:Wedbsite
    def actionCopyrigth(self):
        msgbox=QtWidgets.QMessageBox()
        msgbox.setIcon(1)
        msgbox.setWindowTitle('About the FBG_SiMul/Copyrigth')
        msgbox.setWindowIcon(QtGui.QIcon("DESIGN_Resource/plot.png"))
        msgbox.setTextFormat(1)
        msgbox.setText('The FBG_SiMul is a software that simulates the reflected signal from a fibre Bragg grating array along a predefined path in a Finite Element Model. <br> <br>\
        This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.<br> <br>\
        If you found this software useful in any kind, please cite it together with article DOI.<br> <br>\
        The author would like to acknowledge:<br>\
        -The seventh Framework Programme (FP7) for funding the project MareWint (Project reference: 309395) as Marie-Curie Initial Training;<br>\
        -My supervisors, Lars P. Mikklesen and M. MCGugan, for all the support in the software development;<br>\
        -Denmark Technical University, Department of Wind Energy as my host institution (<a href="http://www.vindenergi.dtu.dk/english">DTU Wind Webpage</a>);')
        msgbox.exec_()

    def actionAbouttheauthor(self):
        msgbox=QtWidgets.QMessageBox()
        msgbox.setIcon(1)
        msgbox.setWindowTitle('About the author')
        msgbox.setWindowIcon(QtGui.QIcon("DESIGN_Resource/author.png"))
        msgbox.setTextFormat(1)
        msgbox.setText('Hi, <br> <br>\
        First, thanks for your interest in this software, my name is Gilmar F. Pereira and I\'m the FBG_SiMul main developer.<br><br>\
        I created the FBG_Simul while I was a PhD student at Denmark Technical University, as part of my thesis topic. My main motivation was to create a simulation tool to study the implementation of optical fibres in different structures; and ensure this tool could be used by any user even without a deep knowledge of optical fibres. <br><br>\
        If you want to know more about the theory behind the code please take a look in some of my articles. You can find my list of publications here: <a href="http://tinyurl.com/nmxehwd">DTU Orbit</a><br><br>\
        Please check out my Linkdin profil: <a href="https://dk.linkedin.com/in/gilmarfp">Linkdin</a> <br><br>\
        Contact: gfpe@dtu.dk; Gilmar_fp@outlook.com <br><br>\
        All the best and have a good FBG simulation,<br>\
        Gilmar P.<br>\
        Edited by Ben Frey in 2020')
        msgbox.exec_()

    def actionSoftwareDoc(self):
        webbrowser.open('Software_Documentation.pdf')

    def actionExitProgram(self):
        self.terminate()

    """--------------------------Tab Extract Stress/Strain-----------------
    In here, it is presented all the functions used in the Tab: Extract Stress/strain
    to extract a file which contains the stress and strain along an optical
    fibre path. This tool was developed for 2D and 3D models in ABAQUS."""

    #Push Button: Select Abaqus Path;
    def actionSelectAbaqusPath(self):
        selectedfile = str(QtWidgets.QFileDialog.getOpenFileName(self, "Select Abaqus ODB file (.odb)",'*.odb'))[0]
        self.ui.AbaqusODBFile.setText(selectedfile)

    #-------------------Path Coordinates Section------------------------------
    #Push Button: Add new Path
    def actionAdd_Path_Coordinates(self):
        vector, ok = QtWidgets.QInputDialog.getText(self, 'Insert the path coordinates','Coordinates Format: (x0,y0,z0);(x1,y1,z1);...,')
        if ok:
            #------In here, the format of the path will be checked.----------------
            #Check if it is empty
            if vector =='':
                QtWidgets.QMessageBox.warning(self, "Error",'No coordinates were inserted.')
                return
             #Check if the format is correct
            else:
                tempvector=vector
                tempvector=tempvector.split(';',99) # Splits the vector in points- Divided By ;
                self.PathCoordinates[str(self.PathNumber+1)]=[] #Creates an empty list to store the correct values
                for i in np.arange(0,len(tempvector)): #Cycle to all the points inserted
                    icoordinate=tempvector[i]
                    icoordinate=icoordinate.replace(' ','')#Remove all the White spaces
                    icoordinate=icoordinate.replace('(','')#Remove parenthesis(.
                    icoordinate=icoordinate.replace(')','')#Remove parenthesis).
                    try:
                        icoordList=map(float,icoordinate.split(','))#Convert the string in a list of float numbers
                    except:
                         QtWidgets.QMessageBox.warning(self, "Error",'The vector inserted is not valid. Please check the format of the coordinates.')
                         return
                    #Check if the list have 3 coordinate points (x,y and z)
                    if np.size(icoordList)!=3:
                        QtWidgets.QMessageBox.warning(self, "Error",'The vector inserted is not valid. Please check the format of the coordinates.')
                        return
                    else:
                        #Now that the input passed the test, let's add it to the dictonary: PathCoordinates
                        self.PathCoordinates[str(self.PathNumber+1)].append(icoordList)

            self.PathNumber+=1 #Add 1 to the path number counter
            #Write in the Plain Text Edit the created path coordinate
            self.ui.PathCoordinatesInputList.insertPlainText('Path No.'+str(self.PathNumber)+'\n'+ str(self.PathCoordinates[str(self.PathNumber)])+'\n')

    #Push Button: Remove path
    def actionRemovePath(self):
       #Remove the entry from the dictonary: PathCoordinates
       map(self.PathCoordinates.pop,str(self.PathNumber))
       self.PathNumber-=1 # Subtract 1 to the path number counter

       # Remove the entry from the Plain Text Edit
       self.ui.PathCoordinatesInputList.setPlainText('')
       tempPathNumber=1
       for o in np.arange(0,self.PathNumber):
           self.ui.PathCoordinatesInputList.insertPlainText('Path No.'+str(tempPathNumber)+'\n'+ str(self.PathCoordinates[str(tempPathNumber)])+'\n')
           tempPathNumber+=1

    #-------------------Radio Button Controllers ------------------------------
    def actionToggleDefaultCoordinate(self):
        self.ui.radioButton_UserDefinedCoordinate.setChecked(False)
        self.ui.radioButton_VectorNewCoordinate.setChecked(False)

    def actionToggleUserDefinedCoordinate(self):
        self.ui.radioButton_VectorNewCoordinate.setChecked(False)
        self.ui.radioButton_DefaultCoordinate.setChecked(False)

    def actionToggleVectorNewCoordinate(self):
        self.ui.radioButton_DefaultCoordinate.setChecked(False)
        self.ui.radioButton_UserDefinedCoordinate.setChecked(False)

    def actionToggleAllDisplacement(self):
        self.ui.radioButton_SpecificDisplacement.setChecked(False)

    def actionToggleSpecificDisplacement(self):
        self.ui.radioButton_AllDisplacement.setChecked(False)

    #---------------------Submit to Abaqus Section-----------------------------
    """In this section, the input data is checked before submit to Abaqus"""
    #Push Button: Select the path where the output will be saved.
    def actionSelectOutputPath(self):
        selectdir = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Abaqus ODB file (.odb)"))
        self.ui.lineEdit_OutputPath.setText(selectdir)

    #Push Button: Submit job to Abaqus
    def actionSubmitToAbaqus(self):
        #------Check if input is correct and all the options are selected------
        if self.ui.AbaqusODBFile.text()=='': # Check if Abaqus path is empty
            self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please select an Abaqus file (.odb). \n')
            return

        #Radio Button User-defined: Check if Coordinate name was inserted
        if self.ui.radioButton_UserDefinedCoordinate.isChecked() and self.ui.lineEdit_UserDefinedCoordinate.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please insert the name of the user-defined coordinate system. \n')
            return

        #Radio Button Vector new coordinates: Check the format of the vector
        if self.ui.radioButton_VectorNewCoordinate.isChecked():
            if self.ui.lineEdit_VectorNewCoordinate.text()=='': # Check if it is empty
                self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please insert a vector with the new coordinate system direction. \n')
                return
            else: # Check the format
                tempvector=self.ui.lineEdit_VectorNewCoordinate.text()
                tempvector=tempvector.split(';',3)
                if len(tempvector)!= 3.0: #Check if the vector have 3 points.
                    self.ui.MessageBoard.insertPlainText('>>ERROR!!: --Wrong format-- The vector with the new coordinate system direction doesn\'t match the required format. \n')
                    return
                #Split in three Coordinates (Origin, x-direction, y-direction)
                OCoordinate=tempvector[0]
                XCoordinate=tempvector[1]
                YCoordinate=tempvector[2]
                #Remove all white spaces
                OCoordinate=OCoordinate.replace(' ','')
                XCoordinate=XCoordinate.replace(' ','')
                YCoordinate=YCoordinate.replace(' ','')
                #Remove all the parentesis
                OCoordinate=OCoordinate.replace('(','')
                XCoordinate=XCoordinate.replace('(','')
                YCoordinate=YCoordinate.replace('(','')
                OCoordinate=OCoordinate.replace(')','')
                XCoordinate=XCoordinate.replace(')','')
                YCoordinate=YCoordinate.replace(')','')
                #From string to list of floats
                Olist=map(float,OCoordinate.split(','))
                Xlist=map(float,XCoordinate.split(','))
                Ylist=map(float,YCoordinate.split(','))

                #Check if the list is correct: each coordinate have 3 points x,y,z
                if np.size(Olist)!=3 or np.size(Xlist)!=3 or np.size(Ylist)!=3:
                    self.ui.MessageBoard.insertPlainText('>>ERROR!!: --Wrong format-- The vector with the new coordinate system direction doesn\'t match the required format.\n')
                    return

                #Create a global variable_dictionary with the "correct" vector for the new coordinate system direction
                self.NewCoordinateSystem={}
                self.NewCoordinateSystem['0']=Olist
                self.NewCoordinateSystem['1']=Xlist
                self.NewCoordinateSystem['2']=Ylist
          #-----------Check the path coordinate-------------------------------
        if self.ui.PathCoordinatesInputList.toPlainText()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please insert a Path.\n')
            return

          #-----------Check the specific displacement step---------------------
        if self.ui.radioButton_SpecificDisplacement.isChecked():
            if self.ui.lineEdit_SpecificDisplacement.text()=='':
                self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please insert the specific displacement step. \n')
                return
            #Check if the input is an integer
            if str(self.ui.lineEdit_SpecificDisplacement.text()).isdigit()==False:
                self.ui.MessageBoard.insertPlainText('>>ERROR!!: Specific displacement step invalid format. Please insert an integer number.')
                return

          #-----------Check the output path----------------------------------------
        if self.ui.lineEdit_OutputPath.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR!!: Please select a folder to write the output file. \n')
            return
            """---------------------------------------------------------------
            Create and Submit the Input file to Abaqus
            """
        #Generate data that will be submitted to the class Extract_Abaqus_Stress
        Odbpath=str(self.ui.AbaqusODBFile.text())

        Writefolder=str(self.ui.lineEdit_OutputPath.text())

        #PathCoordinates is self.PathCoordinates
        PathCoordinates=self.PathCoordinates
        #RotateAxis: 0- Default; 1- User-Defined; 2- Vector with new coordinates
        if self.ui.radioButton_DefaultCoordinate.isChecked():
            RotateAxis=0
            UserDefinedCoordinate=None
            VectorNewCoordinate=None
        if self.ui.radioButton_UserDefinedCoordinate.isChecked():
            RotateAxis=1
            VectorNewCoordinate=None
            UserDefinedCoordinate=str(self.ui.lineEdit_UserDefinedCoordinate.text())
        if self.ui.radioButton_VectorNewCoordinate.isChecked():
            RotateAxis=2
            UserDefinedCoordinate=None
            VectorNewCoordinate=self.NewCoordinateSystem

        # DisplacementStep: 0-all displacement; 1- Specific displacement step
        if self.ui.radioButton_AllDisplacement.isChecked():
           DisplacementStep=0
           SpecificDisplacement=None
        if self.ui.radioButton_SpecificDisplacement.isChecked():
            DisplacementStep=1
            SpecificDisplacement=str(self.ui.lineEdit_SpecificDisplacement.text())

       #------------------Generate Abaqus Python File--------------------------
        from ExtractAbaqusStress import ExtractAbaqus#Class to extract stress and strain along a path in ABAQUS

        sendabaqus=ExtractAbaqus(Odbpath,Writefolder,PathCoordinates,RotateAxis,DisplacementStep,SpecificDisplacement,UserDefinedCoordinate,VectorNewCoordinate)
        sendabaqus.createPyFile()
        self.ui.MessageBoard.insertPlainText('>> Input file successful generated. \n')

        #Send to Abaqus
        sendabaqus.submitAbaqus()
        self.ui.MessageBoard.insertPlainText('>> Stress and Strain along path successful created. \n')

    """Help pushButton-------------------------------------------------------"""
    #Help Button in Extract Stress/Strain Tab
    def actionHelpStressStrainExtract(self):
        QtWidgets.QMessageBox.information(self, "Tool to extract the stress/strain along a pre-defined path",
        "Tool to extract the stress and strain along a pre-defined path in a finite element model.  The output file/s generated in this tab will be used later to simulate the FBG spectrum.\
        \n\nThis tool have the following options:\n\
        \n -Rotate the axis coordinate system, to match the direction of the model with the directions of the optical fibre.\
        \n -Create multiple paths, to extract multiple fibre lines.\
        \n -Extract the stress/strain data from specific or all step increments.\
        \n\nFor more information check the user-manual: Section 4.")

    def actionHelpRotateAxisCoordinate(self):
        QtWidgets.QMessageBox.information(self, "Help - Rotate Axis",
        "In the station (2) Rotate Coordinate System, the user have the possibility to rotate the coordinate axis. This allows the user to match the model axis direction with the optical fibre direction.\n \
         \nThe user have the following options:\n\n\
         -Default: uses the default default model coordinate system.\n\
         -User-defined coordinate sytem: uses a user-defined coordinate system, created during themodel development; Input: name of the user-defined coordinate system.\n\
         -Vector with new coordinate system directions: uses a vector to defined the new coordinate sytem direction. Input format: (x0,y0,y0);(x1,y1,y1);(x2,y2,y2); 0- origin; 1- x direction; 2- y direction.\
         \n\nFor more information check the user-manual: Section 4.2.")

    def actionHelpPathCoordinates(self):
        QtWidgets.QMessageBox.information(self, "Help - Optical Fibre Path Coordinates",
        "In the station (3) Optical Fibre Path Coordinates, the user is asked to insert the optical fibre path coordinates. This path represents the virtual location of the optical fibre line in the FEM model. Each path should have a minimum of two points, but multiple paths can be inserted.\n \
        \nInput format: (x0,y0,y0);(x1,y1,y1);...(x...,y...,y...);\
         \n\nFor more information check the user-manual: Section 4.3.")

    def actionHelpDisplacementStep(self):
        QtWidgets.QMessageBox.information(self, "Help - Time increment (Step)",
        "In station (4) Time increment, the user can select between two time increment options to extract the stress/strain: All time increments or specific time increment.\n \
        \nThe software will name the output files accordingly with the increment number.\
         \n\nFor more information check the user-manual: Section 4.4.")

    """-------------------Signal Simulation Section (OSA)-------------------"""
    #Push Button: Load the .txt files with stress and strain along a fiber
    def actionLoadFolderOSA(self):
        selectedfile = QtWidgets.QFileDialog.getOpenFileNames(self, "Select file(s) with stress and strain along the path (.txt).",'*.txt')[0]
        self.ui.listWidget_stressStrainFiles.addItems(selectedfile)
        #self.ui.listWidget_stressStrainFiles.sortItems(0) #Sort the items

    #Push Button: Remove loaded files
    def actionRemoveFileStressStrain(self):
        SelectedItem=self.ui.listWidget_stressStrainFiles.currentRow()
        self.ui.listWidget_stressStrainFiles.takeItem(SelectedItem)

    #Push Button: Add FBG array position
    def actionAddFBGPosition(self):
        #Empty List with FBG array position
        self.FBGPosition=[]
        #Clear the EditBox
        self.ui.textEdit_FBGPosition.setText('')
        #Number of FBG per fibre
        numberFBG=self.ui.spinBox_NumberFBG.value()
        for i in range(0,numberFBG):
            position, ok = QtWidgets.QInputDialog.getText(self, 'FBG:'+str(i+1)+' position','Insert the FBG:'+str(i+1)+' position. (mm or m from the begin of the line)')
            if ok:
                #------In here, the format of the path is checked.----------------
                #Check if it is empty
                if position =='':
                    QtWidgets.QMessageBox.warning(self, "Error",'No value inserted.')
                    return
                 #Try to convert to float
                else:
                    try:position=float(position)
                    except:
                        QtWidgets.QMessageBox.warning(self, "Error",'Invalid format.')
                        return
                    if self.FBGPosition:
                        if position <= max(self.FBGPosition):
                            QtWidgets.QMessageBox.warning(self, "Error",'The FBG position should be larger than the previous inserted.')
                            self.FBGPosition=[] #Clear the list with FBG array position
                            self.ui.textEdit_FBGPosition.setText('') #Clear the EditBox
                            return
                    self.FBGPosition.append(position)#save in a local variable
                    #Write in the EditBox
                    self.ui.textEdit_FBGPosition.insertPlainText('FBG '+str(i+1)+': '+str(position)+' (mm or m) \n')
            else:
                return
    #Push Button: Clear the FBG array position
    def actionClearFBGPosition(self):
        #Clear the list with FBG array position
        self.FBGPosition=[]
        #Clear the EditBox
        self.ui.textEdit_FBGPosition.setText('')

    #Push Button: Auto- Add FBG original wavelength automatically based on the number of FBG sensors per fibre
    def actionAutoOSA(self):
        #Empty List with FBG array original wavelength (nm)
        self.FBGOriginalWavel=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength.setText('')
        #Number of FBG per fibre
        numberFBG=self.ui.spinBox_NumberFBG.value()
        #minimum bandwidth wavelength
        minwav=float(self.ui.lineEdit_MinBandWidth.text())
        #mmaximum bandwidth wavelength
        maxwav=float(self.ui.lineEdit_MaxBandWidth.text())
        #Distance between FBG
        FBGincrem=(maxwav-minwav)/(numberFBG+1)
        for i in range(0,numberFBG):
            self.FBGOriginalWavel.append(minwav+FBGincrem*(i+1))
            #Write in the EditBox
            self.ui.textEdit_FBGOriginalWavelength.insertPlainText('FBG '+str(i+1)+': '+str(minwav+FBGincrem*(i+1))+' (nm) \n')

    #Push Button: Add FBG original wavelength
    def actionAddFBGOriginalWavelength(self):
        #Empty List with FBG array original wavelength (nm)
        self.FBGOriginalWavel=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength.setText('')
        #Number of FBG per fibre
        numberFBG=self.ui.spinBox_NumberFBG.value()
        for i in range(0,numberFBG):
            wavelength, ok = QtWidgets.QInputDialog.getText(self, 'FBG:'+str(i+1)+' original wavelength','Insert the FBG:'+str(i+1)+' original wavelength. (in nm)')
            if ok:
                #------In here, the format of the path is checked.----------------
                #Check if it is empty
                if wavelength =='':
                    QtWidgets.QMessageBox.warning(self, "Error",'No value inserted.')
                    return
                 #Try to convert to float
                else:
                    try:wavelength=float(wavelength)
                    except:
                        QtWidgets.QMessageBox.warning(self, "Error",'Invalid format.')
                        return
                    if self.FBGOriginalWavel:
                        if wavelength <= max(self.FBGOriginalWavel):
                            QtWidgets.QMessageBox.warning(self, "Error",'The FBG original wavelength should be larger than the previous inserted.')
                            self.FBGOriginalWavel=[] #Clear the list with FBG original wavelength position
                            self.ui.textEdit_FBGOriginalWavelength.setText('') #Clear the EditBox
                            return
                    self.FBGOriginalWavel.append(wavelength)#save in a local variable
                    #Write in the EditBox
                    self.ui.textEdit_FBGOriginalWavelength.insertPlainText('FBG '+str(i+1)+': '+str(wavelength)+' (nm) \n')
            else:
                return

    #Push Button: Clear the FBG Original Wavelength
    def actionClearFBGOriginalWavelength(self):
        #Clear the list with FBG array position
        self.FBGOriginalWavel=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength.setText('')


    #-------------------Radio Button Type of Simulation -----------------------
    def actionToggleTypeSimulationLS(self):
        self.ui.radioButton_TypeSimulationLSNUS.setChecked(False)
        self.ui.radioButton_TypeSimulationLSC.setChecked(False)
        self.ui.radioButton_TypeSimulationT.setChecked(False)
        self.ui.radioButton_TypeSimulationLSCT.setChecked(False)

    def actionToggleTypeSimulationLSNUS(self):
        self.ui.radioButton_TypeSimulationLS.setChecked(False)
        self.ui.radioButton_TypeSimulationLSC.setChecked(False)
        self.ui.radioButton_TypeSimulationT.setChecked(False)
        self.ui.radioButton_TypeSimulationLSCT.setChecked(False)

    def actionToggleTypeSimulationLSC(self):
        self.ui.radioButton_TypeSimulationLSNUS.setChecked(False)
        self.ui.radioButton_TypeSimulationLS.setChecked(False)
        self.ui.radioButton_TypeSimulationT.setChecked(False)
        self.ui.radioButton_TypeSimulationLSCT.setChecked(False)

    def actionToggleTypeSimulationT(self):
        self.ui.radioButton_TypeSimulationLSNUS.setChecked(False)
        self.ui.radioButton_TypeSimulationLS.setChecked(False)
        self.ui.radioButton_TypeSimulationLSC.setChecked(False)
        self.ui.radioButton_TypeSimulationLSCT.setChecked(False)

    def actionToggleTypeSimulationLSCT(self):
        self.ui.radioButton_TypeSimulationLSNUS.setChecked(False)
        self.ui.radioButton_TypeSimulationLS.setChecked(False)
        self.ui.radioButton_TypeSimulationLSC.setChecked(False)
        self.ui.radioButton_TypeSimulationT.setChecked(False)

    #-------------------Radio Button Input Units -----------------------
    def actionToggleInputUnitsSIM(self):
        self.ui.radioButton_InputUnitsSIMM.setChecked(False)

    def actionToggleInputUnitsSIMM(self):
        self.ui.radioButton_InputUnitsSIM.setChecked(False)

    #---------------------Generate the Spectrum Section-----------------------------
    """In this section, the input data is checked before generate the spectrum"""
    #Push Button: Generate
    def actionOsaGenerate(self):
        #Progress Bar
        self.ui.progressBar.setValue(0)
        #Check if input file is generated
        SelectedInput=self.ui.listWidget_stressStrainFiles.currentRow()
        if SelectedInput==-1:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (1)!!: Please Select an Input file. \n')
            return
        #Check if all the Optical Fibre Parameters were inserted and are float
        if self.ui.lineEdit_PhotoElasticParam.text()=='' \
            or self.ui.lineEdit_InitialRefractiveIndex.text()=='' \
            or self.ui.lineEdit_MeanChangeRefractiveIndex.text()=='' \
            or self.ui.lineEdit_FringeVisibility.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP11.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP12.text()=='' \
            or self.ui.lineEdit_YoungsModule.text()=='' \
            or self.ui.lineEdit_PoissionsCoefficient.text()=='' \
            or self.ui.lineEdit_MinBandWidth.text()=='' \
            or self.ui.lineEdit_MaxBandWidth.text()==''\
            or self.ui.lineEdit_ThermalExpansion.text()=='' \
            or self.ui.lineEdit_ThermoOptic.text()=='' \
            or self.ui.lineEdit_Temperature.text()=='':
                self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert all optical fibre parameters. \n')
                self.ui.progressBar.setValue(0)
                return
        try:
            float(self.ui.lineEdit_PhotoElasticParam.text())
            float(self.ui.lineEdit_InitialRefractiveIndex.text())
            float(self.ui.lineEdit_MeanChangeRefractiveIndex.text())
            float(self.ui.lineEdit_FringeVisibility.text())
            float(self.ui.lineEdit_DirectionalRefractiveP11.text())
            float(self.ui.lineEdit_DirectionalRefractiveP12.text())
            float(self.ui.lineEdit_YoungsModule.text())
            float(self.ui.lineEdit_PoissionsCoefficient.text())
            float(self.ui.lineEdit_MinBandWidth.text())
            float(self.ui.lineEdit_MaxBandWidth.text())
            float(self.ui.lineEdit_ThermalExpansion.text())
            float(self.ui.lineEdit_ThermoOptic.text())
            float(self.ui.lineEdit_Temperature.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Invalid format!! --Optical fibre parameters. \n')
            self.ui.progressBar.setValue(0)
            return

        #Progress Bar
        self.ui.progressBar.setValue(10)

        #Check if tolerance and FBG length were inserted
        if self.ui.lineEdit_FBGlength.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Please insert the FBG length. \n')
            self.ui.progressBar.setValue(0)
            return
        try: float(self.ui.lineEdit_FBGlength.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Invalid format!! --FBG length. \n')
            self.ui.progressBar.setValue(0)
            return

        if self.ui.lineEdit_tolerance.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Please insert the tolerance. \n')
            return
        try: float(self.ui.lineEdit_tolerance.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Invalid format!! --Tolerance. \n')
            self.ui.progressBar.setValue(0)
            return

        #Check if FBG position was inserted
        if self.FBGPosition==[]:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Please insert the FBG position. \n')
            self.ui.progressBar.setValue(0)
            return

        #Check if FBG original wavelength was inserted
        if self.FBGOriginalWavel==[]:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (4)!!: Please insert the FBG array original wavelength. \n')
            self.ui.progressBar.setValue(0)
            return

        #Check if simulation resolution was inserted
        if self.ui.lineEdit_SimulationResolution.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR in (5)!!: Please insert the simulation resolution. \n')
            self.ui.progressBar.setValue(0)
            return
        try: float(self.ui.lineEdit_SimulationResolution.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (2)!!: Invalid format!! --Simulation resolution. \n')
            self.ui.progressBar.setValue(0)
            return
        """ Now that the input data is corrected lets prepare it to submit
         to the OSASimulation Class"""
        self.ui.progressBar.setValue(30)#Progress Bar

        #First call the function OSASimulation.LoadFile- To generate the FBG Data that contains the stress and strain per FBG
        filename=self.ui.listWidget_stressStrainFiles.item(SelectedInput).text().replace("\\","/",99) # string Name of the input file
        NumberFBG=self.ui.spinBox_NumberFBG.value() #Int Number FBG
        FBGlength=float(self.ui.lineEdit_FBGlength.text())
        Tolerance=float(self.ui.lineEdit_tolerance.text())
        SkipRow=self.ui.spinBox_SkipRow.value()
        #Input Units
        if self.ui.radioButton_InputUnitsSIM.isChecked():
            self.InputUnits=0
        else:
            self.InputUnits=1

        #Start the OSASimulation function
        self.Osa=OSASimulation(filename,NumberFBG,FBGlength,Tolerance,SkipRow,self.FBGPosition,self.InputUnits)

        #Load other parameters
        PhotoElasticParam=float(self.ui.lineEdit_PhotoElasticParam.text())
        InitialRefractiveIndex=float(self.ui.lineEdit_InitialRefractiveIndex.text())
        MeanChangeRefractiveIndex=float(self.ui.lineEdit_MeanChangeRefractiveIndex.text())
        FringeVisibility=float(self.ui.lineEdit_FringeVisibility.text())
        DirectionalRefractiveP11=float(self.ui.lineEdit_DirectionalRefractiveP11.text())
        DirectionalRefractiveP12=float(self.ui.lineEdit_DirectionalRefractiveP12.text())
        YoungsModule=float(self.ui.lineEdit_YoungsModule.text())
        PoissionsCoefficient=float(self.ui.lineEdit_PoissionsCoefficient.text())
        MinBandWidth=float(self.ui.lineEdit_MinBandWidth.text())
        MaxBandWidth=float(self.ui.lineEdit_MaxBandWidth.text())
        ThermalExpansion=float(self.ui.lineEdit_ThermalExpansion.text())
        ThermoOptic=float(self.ui.lineEdit_ThermoOptic.text())
        Temperature=float(self.ui.lineEdit_Temperature.text())
        SimulationResolution=float(self.ui.lineEdit_SimulationResolution.text())

        self.ui.progressBar.setValue(50)#Progress Bar

        #Generate the undeformed FBG signal
        if self.ui.checkBox_UndeformedSignal.isChecked():
            self.Osa.UndeformedFBG(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution)

        self.ui.progressBar.setValue(75)#Progress Bar

        #FBG sensor direction 0-xx 1-yy 2-zz
        self.FBGDirection=self.ui.comboBox_FBGDirection.currentIndex()

        #Generate the FBG signal from the uniform Strain Contribution
        if self.ui.radioButton_TypeSimulationLS.isChecked():
            self.Osa.UniformStrain(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,self.FBGDirection)

        #Generate the FBG signal: Non-Uniform Strain Contribution
        if self.ui.radioButton_TypeSimulationLSNUS.isChecked():
            self.Osa.NonUniformStrain(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,self.FBGDirection)

        #Generate the FBG signal: Longitudinal Strain + Transverse Strain
        if self.ui.radioButton_TypeSimulationLSC.isChecked():
            self.Osa.TransverseStrain(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,self.FBGDirection,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient)

        #Generate the FBG signal: Temperature
        if self.ui.radioButton_TypeSimulationT.isChecked():
            self.Osa.Temperature(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,self.FBGDirection,ThermalExpansion,ThermoOptic,Temperature)

        #Generate the FBG signal: Longitudinal Strain + Transverse Strain + Temperature
        if self.ui.radioButton_TypeSimulationLSCT.isChecked():
            self.Osa.TransverseStrainTemperature(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,MinBandWidth,MaxBandWidth,SimulationResolution,self.FBGDirection,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient,ThermalExpansion,ThermoOptic,Temperature)

        #Generates the summaryzed Data
        self.Osa.FBGOutputSum(self.FBGOriginalWavel,PhotoElasticParam,InitialRefractiveIndex,self.FBGDirection,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissionsCoefficient,ThermalExpansion,ThermoOptic,Temperature)

        self.ui.progressBar.setValue(100) #Progress Bar
        #Message
        self.ui.MessageBoard.insertPlainText(">> The FBG reflected spectrum was successfully simulated. \n ")

    #Push Button: Load/Plot
    def actionPlotOSA(self):
        #Check if data was generated
        try:PlotWindowOSA(None,False,self.Osa).start()
        except:
            self.ui.MessageBoard.insertPlainText(">> Error!! Please simulated the FBG spectrum before ploting. \n ")

    #Push Button: Save as file
    def actionSaveAsFileOSA(self):
        fsize=None
        #Size of the file
        try:fsize=len(self.Osa.OReflect['wavelength'])
        except:pass
        try:fsize=len(self.Osa.USReflect['wavelength'])
        except:pass
        try:fsize=len(self.Osa.NUSReflect['wavelength'])
        except:pass
        try:fsize=len(self.Osa.TSYReflect['wavelength'])
        except:pass

        if fsize:
            saveFile = str(QtWidgets.QFileDialog.getSaveFileName(self, "Save FBG Spectrum Plot as a file.",'*.txt'))[0]
            if saveFile!='':
                savefile = open(saveFile,"w")
                #Write the Head
                savefile.write("File name: "+str(self.Osa.filename)+"\n") #file name
                savefile.write("Type of Simulation: ") #Type of simulation
                if self.ui.radioButton_TypeSimulationLS.isChecked():
                    savefile.write("Longitudinal Strain (Uniform Strain) \n")
                if self.ui.radioButton_TypeSimulationLSNUS.isChecked():
                    savefile.write("Longitudinal Strain (Non-Uniform Strain) \n")
                if self.ui.radioButton_TypeSimulationLSC.isChecked():
                    savefile.write("Longitudinal strain (US) and Transverse stress \n")
                if self.ui.radioButton_TypeSimulationT.isChecked():
                    savefile.write("Temperature \n")
                if self.ui.radioButton_TypeSimulationLSCT.isChecked():
                    savefile.write("Longitudinal strain (US) and Transverse stress and Temperature \n")
                savefile.write("\n ------FBG array configuration------ \n")
                savefile.write("Number of FBG sensors: "+str(self.Osa.NumberFBG)+" \n")
                savefile.write("FBG length (mm): "+str(self.Osa.FBGLength)+" \n")
                savefile.write("Tolerance (mm): "+str(self.Osa.Tolerance)+" \n")
                savefile.write("FBG position (mm): "+str(self.Osa.FBGPosition)+" \n")
                savefile.write("FBG original wavelength (nm): "+str(self.Osa.FBGOriginalWavel)+" \n \n")

                #Write FBG summary
                savefile.write('FBG number | Original Wavelength | Wavelength shift | Width variation \n')
                for b in np.arange(0,self.Osa.NumberFBG):
                    savefile.write('FBG: '+str(b)+' | '+str(self.Osa.FBGOriginalWavel[b])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveShift'][0])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveWidth'][0])+'\n')
                savefile.write('\n')

                if self.ui.radioButton_TypeSimulationLS.isChecked(): #Longitudinal Strain
                    if self.ui.checkBox_UndeformedSignal.isChecked():
                        savefile.write("Wavelength (nm) \t Reflectivity:Undeformed Shape  \t Reflectivity: Longitudinal Strain (Uniform Strain) \n")
                    else:
                        savefile.write("Wavelength (nm) \t Reflectivity: Longitudinal Strain (Uniform Strain) \n")

                if self.ui.radioButton_TypeSimulationLSNUS.isChecked():#Longitudinal Strain (Non-uniform)
                    if self.ui.checkBox_UndeformedSignal.isChecked():
                        savefile.write("Wavelength (nm) \t Reflectivity:Undeformed Shape  \t Reflectivity: Longitudinal Strain (Non-Uniform Strain) \n")
                    else:
                        savefile.write("Wavelength (nm) \t Reflectivity: Longitudinal Strain (Non-Uniform Strain) \n")

                if self.ui.radioButton_TypeSimulationLSC.isChecked():# Longitudinal strain + transverse strain
                    if self.ui.checkBox_UndeformedSignal.isChecked():
                        savefile.write("Wavelength (nm) \t Reflectivity:Undeformed Shape  \t Reflectivity: Longitudinal Strain + Transverse Stress (Y Wave) \t Reflectivity: Longitudinal Strain + Transverse Stress (Z Wave) \n")
                    else:
                        savefile.write("Wavelength (nm)\t Reflectivity: Longitudinal Strain + Transverse Stress (Y Wave) \t Reflectivity: Longitudinal Strain + Transverse Stress (Z Wave) \n")

                # Write the Data
                for b in np.arange(0,fsize):
                    if self.ui.radioButton_TypeSimulationLS.isChecked(): #Longitudinal Strain
                        if self.ui.checkBox_UndeformedSignal.isChecked():
                            savefile.write('%5f \t %5f \t %5f \n' % (self.Osa.USReflect['wavelength'][b],self.Osa.OReflect['reflec'][b],self.Osa.USReflect['reflec'][b]))
                        else:
                            savefile.write('%5f \t %5f \n' % (self.Osa.USReflect['wavelength'][b],self.Osa.USReflect['reflec'][b]))

                    if self.ui.radioButton_TypeSimulationLSNUS.isChecked():#Longitudinal Strain (Non-uniform)
                        if self.ui.checkBox_UndeformedSignal.isChecked():
                            savefile.write('%5f \t %5f \t %5f \n' % (self.Osa.NUSReflect['wavelength'][b],self.Osa.OReflect['reflec'][b],self.Osa.NUSReflect['reflec'][b]))
                        else:
                            savefile.write('%5f \t %5f \n' % (self.Osa.NUSReflect['wavelength'][b],self.Osa.NUSReflect['reflec'][b]))

                    if self.ui.radioButton_TypeSimulationLSC.isChecked():# Longitudinal strain + transverse strain
                        if self.ui.checkBox_UndeformedSignal.isChecked():
                            savefile.write('%5f \t %5f \t %5f \t %5f \n' % (self.Osa.TSYReflect['wavelength'][b],self.Osa.OReflect['reflec'][b],self.Osa.TSYReflect['reflec'][b],self.Osa.TSZReflect['reflec'][b]))
                        else:
                            savefile.write('%5f \t %5f \t %5f \n' % (self.Osa.TSYReflect['wavelength'][b],self.Osa.TSYReflect['reflec'][b],self.Osa.TSZReflect['reflec'][b]))
                savefile.close()

        else:
            self.ui.MessageBoard.insertPlainText(">> Error!! Please simulated the FBG spectrum before saving as file. \n ")
            return

    """Help pushButton-------------------------------------------------------"""
    def actionAboutFBGSpectrumSimulation(self):
        QtWidgets.QMessageBox.information(self, "About FBG Spectrum Simulation (Specific Step)",
        "In this tab, the user can simulate the FBGreflected spectrumby using the stress/strain file generated in the previous tab. The output from this tab is a plot or a .txt file of the FBG reflected signal(s) at a specific time increment.\
        \n\nTool options:\n\
        \n-Type of simulation;\
        \n-User-defined optical fibre parameters;\
        \n-Multiple FBGs per fibre;\
        \n-User-defined FBG array configuration;\
        \n\nFor more information check the user-manual: Section 5.")


    def actionHelpFileFormatOSA(self):
        QtWidgets.QMessageBox.information(self, "Help - File Format OSA",
        "The input must have the following format:\
        \n\n-8 columns separated by tab;\
        \n\n-Variables per column: -1st- FBG path length/distance; -2nd-  Strain in direction 11; -3rd- Strain in direction 22; -4th- Strain in direction 33; -5th- Stress in direction 11; -6th- Stress in direction 22; -7th- Stress in direction 33;\
        \n\nFor more information check the user-manual: Section 4.0")

    def actionHelpFBGPositionOSA(self):
        QtWidgets.QMessageBox.information(self, "Help - FBG Position (mm or m)",
        "Here, the user defines the FBGs position along the path. By pressing the Add button the user is prompted to input the length from the beginning of the line to the beginning of that particular FBG grating. This distance is defined by the first column of the input file and it should be inserted with the input units (mm or m) defined in station (1).\
        \n\nFor more information check the user-manual: Section 5.4")

    def actionHelpFBGWavelengthOSA(self):
        QtWidgets.QMessageBox.information(self, "Help - FBG Original Wavelength OSA",
        "Here, the user should insert the FBG array wavelength in its unstrained state. The values must lie within the light bandwidth specified in station (3) (by default:between 1500 to 1600nm).\
        The user can either press the Auto button to distribute the original wavelength along the available bandwidth with equal spacing, or may choose to input each FBG original wavelength manually by pressing Manual button.\
        \n\nFor more information check the user-manual: Section 5.4")

    def actionHelpTypeSimulationOSA(self):
        QtWidgets.QMessageBox.information(self, "Help - Type of Simulation OSA",
        "In station (2) Type of simulation (contribution), the user chooses the type of simulation to be performed. \
        \n\n Type of simulations:\
        \n-Longitudinal Strain (Uniform Strain): In this option only uniform strain acting in the longitudinal direction of the FBG sensor is considered;\
        \n-Longitudinal Strain (Uniform+Non-UniformStrain): In this option the effect of Non-uniform strain acting in the longitudinal direction of the FBG sensor is considered;\
        \n-Longitudinal Strain (Uniform Strain) and Transverse Stress: In this option the effect of uniformstrain along the longitudinal direction of the FBG sensor and the transverse stress acting perpendicularly to the Sensor are considered;\
        \n-Temperature: In this option only temperature effects are considered;\
        \n-Long. Strain (Uni), Trans. Stress, Temperature: In this option the effect of uniformstrain along the longitudinal direction of the FBG sensor and the transverse stress acting perpendicularly and temperature effects are considered;\
        \n\nThe simulation resolution correspond to the light bandwidth discretization, i.e. amount of points that the reflected spectrum will be calculated.\
        \n\nFor more information check the user-manual: Section 5.2")

        """--------------------------Tab FBG Signal Variation-----------------
    In here, it is presented all the functions in the Tab: FBG Signal Variation
    used to simulate the FBG response along multiple time increments

    The sensor response is Wavelength shift: givin by the uniform longitudinal
    strain; and peak width variation givin by the non-uniform strain and transverse
    stress contribution
    """

    #Push Button: Add .txt files  with stress strain along fiber
    def actionLoadFolderTR(self):
        selectedfile = QtWidgets.QFileDialog.getOpenFileNames(self, "Select file(s) with stress and strain along the path (.txt).",'*.txt')[0]
        self.ui.listWidget_stressStrainFiles_TV.addItems(selectedfile)
        #self.ui.listWidget_stressStrainFiles.sortItems(0) #Sort the items

    #Push Button: Sort the Input file
    def actionSortFiles(self):
        self.ui.listWidget_stressStrainFiles_TV.sortItems(0)

    #Push Button: Remove loaded files
    def actionRemoveFileStressStrainTR(self):
        SelectedItem=self.ui.listWidget_stressStrainFiles_TV.currentRow()
        self.ui.listWidget_stressStrainFiles_TV.takeItem(SelectedItem)

    #Push Button: Clear all files
    def actionClearAllTR(self):
        self.ui.listWidget_stressStrainFiles_TV.clear()


#-------------------Fibre Bragg Grating Array Configuration-------------------
    #Push Button: Add FBG array position
    def actionAddFBGPositionTR(self):
        #Empty List with FBG array position (TR- Time Response tab)
        self.FBGPositionTR=[]
        #Clear the EditBox
        self.ui.textEdit_FBGPosition_TR.setText('')
        #Number of FBG per fibre
        numberFBGTR=self.ui.spinBox_NumberFBG_TR.value()
        for i in range(0,numberFBGTR):
            position, ok = QtWidgets.QInputDialog.getText(self, 'FBG:'+str(i+1)+' position','Insert the FBG:'+str(i+1)+' position. (mm or m from the begin of the line)')
            if ok:
                #------Here the format of the path will be checked.----------------
                #Check if it is empty
                if position =='':
                    QtWidgets.QMessageBox.warning(self, "Error",'No value inserted.')
                    return
                 #Try to convert to float
                else:
                    try:position=float(position)
                    except:
                        QtWidgets.QMessageBox.warning(self, "Error",'Invalid format.')
                        return
                    if self.FBGPositionTR:
                        if position <= max(self.FBGPositionTR):
                            QtWidgets.QMessageBox.warning(self, "Error",'The FBG position should be larger than the previous inserted.')
                            self.FBGPositionTR=[] #Clear the list with FBG array position
                            self.ui.textEdit_FBGPosition_TR.setText('') #Clear the EditBox
                            return
                    self.FBGPositionTR.append(position)#save in a local variable
                    #Write in the EditBox
                    self.ui.textEdit_FBGPosition_TR.insertPlainText('FBG '+str(i+1)+': '+str(position)+' (mm or m) \n')
            else:
                return
    #Push Button: Clear the FBG array position
    def actionClearFBGPositionTR(self):
        #Clear the list with FBG array position
        self.FBGPositionTR=[]
        #Clear the EditBox
        self.ui.textEdit_FBGPosition_TR.setText('')

    #Push Button: Auto- Add FBG original wavelength automatically based on the number of FBG sensors per fibre
    def actionAutoTR(self):
        #Empty List with FBG array original wavelength (nm)
        self.FBGOriginalWavelTR=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength_TR.setText('')
        #Number of FBG per fibre
        numberFBG=self.ui.spinBox_NumberFBG_TR.value()
        #minimum bandwidth wavelength
        minwav=1500.00
        #mmaximum bandwidth wavelength
        maxwav=1600.00
        #Distance between FBG
        FBGincrem=(maxwav-minwav)/(numberFBG+1)
        for i in range(0,numberFBG):
            self.FBGOriginalWavelTR.append(minwav+FBGincrem*(i+1))
            #Write in the EditBox
            self.ui.textEdit_FBGOriginalWavelength_TR.insertPlainText('FBG '+str(i+1)+': '+str(minwav+FBGincrem*(i+1))+' (nm) \n')

    #Push Button: Add FBG original wavelength
    def actionAddFBGOriginalWavelengthTR(self):
        #Empty List with FBG array original wavelength (nm)
        self.FBGOriginalWavelTR=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength_TR.setText('')
        #Number of FBG per fibre
        numberFBGTR=self.ui.spinBox_NumberFBG_TR.value()
        for i in range(0,numberFBGTR):
            wavelength, ok = QtWidgets.QInputDialog.getText(self, 'FBG:'+str(i+1)+' original wavelength','Insert the FBG:'+str(i+1)+' original wavelength. (in nm)')
            if ok:
                #------Here the format of the path will be checked.----------------
                #Check if it is empty
                if wavelength =='':
                    QtWidgets.QMessageBox.warning(self, "Error",'No value inserted.')
                    return
                 #Try to convert to float
                else:
                    try:wavelength=float(wavelength)
                    except:
                        QtWidgets.QMessageBox.warning(self, "Error",'Invalid format.')
                        return
                    if self.FBGOriginalWavelTR:
                        if wavelength <= max(self.FBGOriginalWavelTR):
                            QtWidgets.QMessageBox.warning(self, "Error",'The FBG original wavelength should be larger than the previous inserted.')
                            self.FBGOriginalWavelTR=[] #Clear the list with FBG original wavelength position
                            self.ui.textEdit_FBGOriginalWavelength_TR.setText('') #Clear the EditBox
                            return
                    self.FBGOriginalWavelTR.append(wavelength)#save in a local variable
                    #Write in the EditBox
                    self.ui.textEdit_FBGOriginalWavelength_TR.insertPlainText('FBG '+str(i+1)+': '+str(wavelength)+' (nm) \n')
            else:
                return

    #Push Button: Clear the FBG Original Wavelength
    def actionClearFBGOriginalWavelengthTR(self):
        #Clear the list with FBG array position
        self.FBGOriginalWavelTR=[]
        #Clear the EditBox
        self.ui.textEdit_FBGOriginalWavelength_TR.setText('')

    #-------------------Radio Button Input Units -----------------------
    def actionToggleInputUnitsSIMTR(self):
        self.ui.radioButton_InputUnitsSIMM_TR.setChecked(False)

    def actionToggleInputUnitsSIMMTR(self):
        self.ui.radioButton_InputUnitsSIM_TR.setChecked(False)



    #---------------------Generate the Spectrum Section-----------------------------
    """
    In this section, the input data will be check before generate the spectrum
    """
    #Push Button: Generate
    def actionTimeResponseGenerate(self):
        #Progress Bar
        self.ui.progressBar_TR.setValue(0)

        #Check if input files were inserted
        InputFileItemNumber=self.ui.listWidget_stressStrainFiles_TV.count()
        if InputFileItemNumber==0:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (1)!!: Please insert the input file(s). \n')
            return

        #Check if all the Optical Fibre Parameters were inserted and are float
        if self.ui.lineEdit_PhotoElasticParam_TR.text()=='' \
            or self.ui.lineEdit_InitialRefractiveIndex_TR.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP11_TR.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP12_TR.text()=='' \
            or self.ui.lineEdit_YoungsModule_TR.text()=='' \
            or self.ui.lineEdit_PoissionsCoefficient_TR.text()=='':
                self.ui.MessageBoard.insertPlainText('>>ERROR in (2)!!: Please insert all optical fibre parameters. \n')
                self.ui.progressBar_TR.setValue(0)
                return

        try:
            float(self.ui.lineEdit_PhotoElasticParam_TR.text())
            float(self.ui.lineEdit_InitialRefractiveIndex_TR.text())
            float(self.ui.lineEdit_DirectionalRefractiveP11_TR.text())
            float(self.ui.lineEdit_DirectionalRefractiveP12_TR.text())
            float(self.ui.lineEdit_YoungsModule_TR.text())
            float(self.ui.lineEdit_PoissionsCoefficient_TR.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (2)!!: Invalid format!! --Optical fibre parameters. \n')
            self.ui.progressBar_TR.setValue(0)
            return

        #Progress Bar
        self.ui.progressBar_TR.setValue(10)


        #Check if tolerance and FBG length were inserted
        if self.ui.lineEdit_FBGlength_TR.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert the FBG length. \n')
            self.ui.progressBar_TR.setValue(0)
            return
        try: float(self.ui.lineEdit_FBGlength_TR.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Invalid format!! --FBG length. \n')
            self.ui.progressBar_TR.setValue(0)
            return

        if self.ui.lineEdit_tolerance_TR.text()=='':
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert the tolerance. \n')
            return
        try: float(self.ui.lineEdit_tolerance_TR.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Invalid format!! --Tolerance. \n')
            self.ui.progressBar_TR.setValue(0)
            return

        #Check if FBG position was inserted
        if self.FBGPositionTR==[]:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert the FBG position. \n')
            self.ui.progressBar_TR.setValue(0)
            return

        #Check if FBG original wavelength was inserted
        if self.FBGOriginalWavelTR==[]:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert the FBG array original wavelength. \n')
            self.ui.progressBar_TR.setValue(0)
            return

        """ Now that the input data is corrected lets prepare to submit
        it to the TRSimulation Class"""
        self.ui.progressBar_TR.setValue(30)#Progress Bar

        #Generate list with all the input files
        InputList=[]
        for i in range(0,InputFileItemNumber):
            InputList.append(str(self.ui.listWidget_stressStrainFiles_TV.item(i).text().replace("\\","/",99)))
        #Get the other variables
        NumberFBGTR=self.ui.spinBox_NumberFBG_TR.value() #Int Number FBG
        FBGlengthTR=float(self.ui.lineEdit_FBGlength_TR.text())
        ToleranceTR=float(self.ui.lineEdit_tolerance_TR.text())
        SkipRowTR=self.ui.spinBox_SkipRow_TR.value()
        PhotoElasticParamTR=float(self.ui.lineEdit_PhotoElasticParam_TR.text())
        InitialRefractiveIndexTR=float(self.ui.lineEdit_InitialRefractiveIndex_TR.text())
        DirectionalRefractiveP11TR=float(self.ui.lineEdit_DirectionalRefractiveP11_TR.text())
        DirectionalRefractiveP12TR=float(self.ui.lineEdit_DirectionalRefractiveP12_TR.text())
        YoungsModuleTR=float(self.ui.lineEdit_YoungsModule_TR.text())
        PoissionsCoefficientTR=float(self.ui.lineEdit_PoissionsCoefficient_TR.text())
        self.FBGDirectionTR=self.ui.comboBox_FBGDirection_TR.currentIndex()

        #Input Units
        if self.ui.radioButton_InputUnitsSIM_TR.isChecked():
            self.InputUnitsTR=0
        else:
            self.InputUnitsTR=1

        #Start the Time Response Simulation function
        self.TR=TRSimulation(InputList,NumberFBGTR,FBGlengthTR,ToleranceTR,SkipRowTR,self.FBGPositionTR,self.InputUnitsTR)

        self.ui.progressBar_TR.setValue(50)#Progress Bar

        #Calculate The response
        self.TR.Calculate(self.FBGOriginalWavelTR,PhotoElasticParamTR,InitialRefractiveIndexTR,DirectionalRefractiveP11TR,DirectionalRefractiveP12TR,YoungsModuleTR,PoissionsCoefficientTR,self.FBGDirectionTR)

        self.ui.progressBar_TR.setValue(100)#Progress Bar
        #Message
        self.ui.MessageBoard.insertPlainText(">> The FBG time response was successfully simulated. \n ")

    #Push Button: Load/Plot
    def actionPlotTR(self):
        #Check if data was generated
        try:PlotWindowTR(None,False,self.TR).start()
        except:
            self.ui.MessageBoard.insertPlainText(">> Error!! Please simulated the FBG time response ploting. \n ")

    #Push Button: Save as file
    def actionSaveAsFileTR(self):
        #Check if data was generated
        fsize=None
        #Size of the file
        try:fsize=len(self.TR.FBGTimeResponse['FBG1']['WaveShift'])
        except:pass

        if fsize:
            saveFile = str(QtWidgets.QFileDialog.getSaveFileName(self, "Save FBG time response as a file.",'*.txt'))[0]
            if saveFile!='':
                savefile = open(saveFile,"w")
                #Write the Head
                savefile.write("Number of input files (increment number):  "+str(len(self.TR.InputList))+"\n") #file name
                savefile.write("\n ------FBG array configuration------ \n")
                savefile.write("Number of FBG sensors: "+str(self.TR.NumberFBG)+" \n")
                savefile.write("FBG length (mm): "+str(self.TR.FBGLength)+" \n")
                savefile.write("Tolerance (mm): "+str(self.TR.Tolerance)+" \n")
                savefile.write("FBG position (mm): "+str(self.TR.FBGPosition)+" \n")
                savefile.write("FBG original wavelength (nm): "+str(self.TR.FBGOriginalWavel)+" \n \n")

                #Write variables names
                savefile.write("Increment number \t")
                for c in np.arange(0,self.TR.NumberFBG):
                     savefile.write("FBG "+str(c+1)+"-Wavelength shift (nm)\t")
                for c in np.arange(0,self.TR.NumberFBG):
                     savefile.write("FBG "+str(c+1)+"-Peak width variation (nm)\t")
                savefile.write("\n")

                #Write the Data
                for b in np.arange(0,fsize): #Cycle each line
                    savefile.write('%d \t' % (int(self.TR.FBGTimeResponse['FBG1']['Increment'][b]))) #write the increment number
                    #Cycle each FBG
                    for c in np.arange(0,self.TR.NumberFBG):
                        savefile.write('%5f \t' % (self.TR.FBGTimeResponse['FBG'+str(c+1)]['WaveShift'][b]))
                    for c in np.arange(0,self.TR.NumberFBG):
                        savefile.write('%5f \t' % (self.TR.FBGTimeResponse['FBG'+str(c+1)]['WaveWidth'][b]))
                    savefile.write("\n")
               #close file
                savefile.close()
        else:
            self.ui.MessageBoard.insertPlainText(">> Error!! Please simulated the FBG spectrum before saving as file. \n ")
            return


    """Help pushButton-------------------------------------------------------"""
    def actionAboutTimeResponse(self):
        QtWidgets.QMessageBox.information(self, "About FBG Signal Variation (Time Response)",
        "In this tab, the user can simulate the FBG response using multiple input files (time response). Similar to Tab 3- FBG Spectrum Simulation, this tool calculates the reflected signal of the FBG based on the stress/strain; however multiple input files are selected, and the output is presented as wavelength shift variation and peak width variation.\
        \n\nTool options:\n\
        \n-User-defined optical fibre parameters;\
        \n-Multiple FBGs per fibre;\
        \n-User-defined FBG array configuration;\
        \n\nFor more information check the user-manual: Section 6.")


"""----------------------------------------------------------------------------
-------------------------------------------------------------------------------
---------------------------PLOT CLass------------------------------------------
-------------------------------------------------------------------------------
"""

"""----------------------------Class to plot: TR tab------------------------
"""
import GUI.PlotWindow_TR
class PlotWindowTR(QtDialogLoader):
    def __init__(self, parent, modal,TR):
        #Local Variables
        self.TR=TR
        self.linewidth=0.5
        self.LineColorUndeformed='grey'
        self.LineColorSimulated='blue'
        #This will start the dialog
        self.ui_module = GUI.PlotWindow_TR
        try: self.ui = self.ui_module.Ui_Form()  #enable autocomplete
        except: pass
        QtDialogLoader.__init__(self, self.ui_module, parent, modal)

        #Start the Canvas for Time Response Plot
        self.TRplot = TRCanvas()

        self.ui.gridLayout_TR_Plot.addWidget(self.TRplot)#Insert The canvas

        #Set the Axes names
        self.TRplot.axes0.set_xlabel('Increments ',fontsize=10)
        self.TRplot.axes0.set_ylabel('Wavelength shift (nm)',fontsize=14)
        self.TRplot.axes1.set_xlabel('Increments ',fontsize=10)
        self.TRplot.axes1.set_ylabel('Peak width variation (nm)',fontsize=14)

        #Set maxIncrements
        self.ui.doubleSpinBox_MaxIncrement.setValue(len(self.TR.FBGTimeResponse['FBG1']['Increment'])-1)

        #Set max and Min Wavelength shift and peak width
        #Calculate maximum
        minwav=0.0
        maxwav=0.0
        minwid=0.0
        maxwid=0.0
        for b in np.arange(0,self.TR.NumberFBG):
            tempminwav=np.min(self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveShift'])
            tempmaxwav=np.max(self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveShift'])
            tempminwid=np.min(self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveWidth'])
            tempmaxwid=np.max(self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveWidth'])
            if tempminwav<minwav:minwav=tempminwav
            if tempmaxwav>maxwav:maxwav=tempmaxwav
            if tempminwid<minwid:minwid=tempminwid
            if tempmaxwid>maxwid:maxwid=tempmaxwid

        self.ui.spinBox_WSMin.setValue(minwav)
        self.ui.spinBox_WSMax.setValue(maxwav)
        self.ui.spinBox_WDMin.setValue(minwid)
        self.ui.spinBox_WDMax.setValue(maxwid)


        #Connect widget signals to actionUpdatePlot
        self.ui.doubleSpinBox_MaxIncrement.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_MinIncrement.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_WSMin.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_WSMax.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_WDMin.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_WDMax.valueChanged.connect(self.actionUpdatePlot)
        self.ui.checkBox_legend.clicked.connect(self.actionUpdatePlot)
        self.ui.checkBox_Grid.clicked.connect(self.actionUpdatePlot)
        self.ui.doubleSpinBox_LineWidth.valueChanged.connect(self.actionUpdatePlot)

        #Plot the TR
        #self.TRplot.axes0.hold(True)    #---BEN_EDIT, check later
        #self.TRplot.axes1.hold(True)    #---BEN_EDIT, check later
        self.Plot()

    def Plot(self):
        #Cycle all the FBG sensors
        for b in np.arange(0,self.TR.NumberFBG):
            #plot Wavelength Shift
            self.TRplot.axes0.plot(self.TR.FBGTimeResponse['FBG'+str(b+1)]['Increment'],self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveShift'],linewidth=self.linewidth,label='FBG'+str(b+1))
            #plot peak width variation
            self.TRplot.axes1.plot(self.TR.FBGTimeResponse['FBG'+str(b+1)]['Increment'],self.TR.FBGTimeResponse['FBG'+str(b+1)]['WaveWidth'],linewidth=self.linewidth)

    def actionUpdatePlot(self):
        self.ui.gridLayout_TR_Plot.removeWidget(self.TRplot)
        self.TRplot = TRCanvas()
        self.ui.gridLayout_TR_Plot.addWidget(self.TRplot)
        #self.TRplot.axes0.hold(True)    #---BEN_EDIT, check later
        #self.TRplot.axes1.hold(True)    #---BEN_EDIT, check later

       #Set the Axes names
        self.TRplot.axes0.set_xlabel('Increments ',fontsize=10)
        self.TRplot.axes0.set_ylabel('Wavelength shift (nm)',fontsize=14)
        self.TRplot.axes1.set_xlabel('Increments ',fontsize=10)
        self.TRplot.axes1.set_ylabel('Peak width variation (nm)',fontsize=14)


        #Set xlimits (Increments)
        self.TRplot.axes0.set_xlim([self.ui.spinBox_MinIncrement.value(), self.ui.doubleSpinBox_MaxIncrement.value()])
        self.TRplot.axes1.set_xlim([self.ui.spinBox_MinIncrement.value(), self.ui.doubleSpinBox_MaxIncrement.value()])
        #Set ylimits
        self.TRplot.axes0.set_ylim([self.ui.spinBox_WSMin.value(), self.ui.spinBox_WSMax.value()])
        self.TRplot.axes1.set_ylim([self.ui.spinBox_WDMin.value(), self.ui.spinBox_WDMax.value()])
        #Set LineWidth
        self.linewidth=self.ui.doubleSpinBox_LineWidth.value()

        #Plot the Time Response
        #self.TRplot.axes0.hold(True)    #---BEN_EDIT, check later
        #self.TRplot.axes1.hold(True)    #---BEN_EDIT, check later
        self.Plot()
        #Legend
        if self.ui.checkBox_legend.isChecked():
            self.TRplot.axes0.legend(fontsize=10, loc="best")
         #Grid
        if self.ui.checkBox_Grid.isChecked():
            self.TRplot.axes0.grid()
            self.TRplot.axes1.grid()

    def actionSavePicture(self):
        savePictureFile = str(QtWidgets.QFileDialog.getSaveFileName(self, "Save FBG Spectrum Plot",'*.png'))[0]
        if savePictureFile!='':
            self.TRplot.figure.savefig(str(savePictureFile))


"""----------------------------Class to plot: OSA tab------------------------
"""
import GUI.PlotWindow_OSA
class PlotWindowOSA(QtDialogLoader):
    def __init__(self, parent, modal,Osa):
        #Local Variables
        self.Osa=Osa
        self.linewidth=0.5
        self.LineColorUndeformed='grey'
        self.LineColorSimulated='blue'
        #This will start the dialog
        self.ui_module = GUI.PlotWindow_OSA
        try: self.ui = self.ui_module.Ui_Form()  #enable autocomplete
        except: pass
        QtDialogLoader.__init__(self, self.ui_module, parent, modal)

        #Start the Canvas for OSA Plot
        self.osaplot = OSACanvas()
        self.ui.gridLayout_OSA_Plot.addWidget(self.osaplot)#Insert The canvas
        #Set the Axes names
        self.osaplot.axes.set_xlabel('Wavelength (nm)',fontsize=14)
        self.osaplot.axes.set_ylabel('Reflectivity (R)',fontsize=14)
        #Connect widget signals to actionUpdatePlot
        self.ui.spinBox_YlimMin.valueChanged.connect(self.actionUpdatePlot)
        self.ui.doubleSpinBox_YlimMax.valueChanged.connect(self.actionUpdatePlot)
        self.ui.spinBox_XlimMin.valueChanged.connect(self.actionUpdatePlot)
        self.ui.doubleSpinBox_XlimMax.valueChanged.connect(self.actionUpdatePlot)
        self.ui.checkBox_legend.clicked.connect(self.actionUpdatePlot)
        self.ui.checkBox_Grid.clicked.connect(self.actionUpdatePlot)
        self.ui.doubleSpinBox_LineWidth.valueChanged.connect(self.actionUpdatePlot)
        self.ui.comboBox_LineColorUndeformed.currentIndexChanged.connect(self.actionUpdatePlot)
        self.ui.comboBox_LineColorSimulated.currentIndexChanged.connect(self.actionUpdatePlot)
        #Plot the Spectrum
        #self.osaplot.axes.hold(True)  #---BEN_EDIT, check later
        self.Plot()

        #FBG Output Table
        self.ui.textEdit_FBGoutput.insertPlainText('FBG number | Original Wavelength | Wavelength shift | Width variation \n')
        for b in np.arange(0,self.Osa.NumberFBG):
            self.ui.textEdit_FBGoutput.insertPlainText('FBG: '+str(b)+' | '+str(self.Osa.FBGOriginalWavel[b])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveShift'][0])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveWidth'][0])+'\n')

    def Plot(self):
        #Plot the Undeformed shape
        try:self.osaplot.axes.plot(self.Osa.OReflect['wavelength'],self.Osa.OReflect['reflec'],color=self.LineColorUndeformed,linewidth=self.linewidth, label="Undeformed FBG Spectrum")
        except:pass
        #Plot the Uniform Strain contribution
        try:self.osaplot.axes.plot(self.Osa.USReflect['wavelength'],self.Osa.USReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Longitudinal Strain (Uniform Strain)")
        except:pass
        #Plot the Non-Uniform Strain contribution
        try:self.osaplot.axes.plot(self.Osa.NUSReflect['wavelength'],self.Osa.NUSReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Longitudinal Strain (Non-Uniform Strain)")
        except:pass
        #Plot the trasnverse Stress Contribution
        try:
            self.osaplot.axes.plot(self.Osa.TSYReflect['wavelength'],self.Osa.TSYReflect['reflec'],color="red",linewidth=self.linewidth, label="Longitudinal Strain and Transverse Stress (Y Wave)")
            self.osaplot.axes.plot(self.Osa.TSZReflect['wavelength'],self.Osa.TSZReflect['reflec'],color="blue",linewidth=self.linewidth, label="Longitudinal Strain and Transverse Stress (Z Wave)")
        except:pass
        #Plot the Temperature contribution
        try:self.osaplot.axes.plot(self.Osa.TReflect['wavelength'],self.Osa.TReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Temperature Contribution")
        except:pass
        #Plot the trasnverse Stress Contribution and Temperature
        try:
            self.osaplot.axes.plot(self.Osa.TSTYReflect['wavelength'],self.Osa.TSTYReflect['reflec'],color="red",linewidth=self.linewidth, label="Longitudinal Strain, Transverse Stress, and Temperature (Y Wave)")
            self.osaplot.axes.plot(self.Osa.TSTZReflect['wavelength'],self.Osa.TSTZReflect['reflec'],color="blue",linewidth=self.linewidth, label="Longitudinal Strain, Transverse Stress, and Temperature (Z Wave)")
        except:pass

    def actionUpdatePlot(self):
        self.ui.gridLayout_OSA_Plot.removeWidget(self.osaplot)
        self.osaplot = OSACanvas()
        self.ui.gridLayout_OSA_Plot.addWidget(self.osaplot)
        #self.osaplot.axes.hold(True)  #---BEN_EDIT, check later
        #Set the Axes names
        self.osaplot.axes.set_xlabel('Wavelength (nm)',fontsize=14)
        self.osaplot.axes.set_ylabel('Reflectivity (R)',fontsize=14)
        #Set xlimits
        self.osaplot.axes.set_xlim([self.ui.spinBox_XlimMin.value(), self.ui.doubleSpinBox_XlimMax.value()])
        #Set ylimits
        self.osaplot.axes.set_ylim([self.ui.spinBox_YlimMin.value(), self.ui.doubleSpinBox_YlimMax.value()])
        #Set LineWidth
        self.linewidth=self.ui.doubleSpinBox_LineWidth.value()
        #Set Line Colour
        self.LineColorUndeformed=str(self.ui.comboBox_LineColorUndeformed.currentText())
        self.LineColorSimulated=str(self.ui.comboBox_LineColorSimulated.currentText())
        #Plot the Spectrum
        #self.osaplot.axes.hold(True)     #---BEN_EDIT, check later
        self.Plot()
        #Legend
        if self.ui.checkBox_legend.isChecked():
            self.osaplot.axes.legend(fontsize=10, loc="best")
         #Grid
        if self.ui.checkBox_Grid.isChecked():
            self.osaplot.axes.grid()

    def actionSavePicture(self):
        savePictureFile = str(QtWidgets.QFileDialog.getSaveFileName(self, "Save FBG Spectrum Plot",'*.png'))[0]
        if savePictureFile!='':
            self.osaplot.figure.savefig(str(savePictureFile))

"""
------------------Class for the OSA Canvas------------------------------------
"""
from PyQt5.QtWidgets import QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['font.size'] = 9

class TRCanvas(Canvas):
    """
    MatplotlibWidget inherits PyQt5.QtWidgets.QWidget
    and matplotlib.backend_bases.FigureCanvasBase
    -------
    parent (None): parent widget
    title (''): figure title
    xlabel (''): X-axis label
    ylabel (''): Y-axis label
    xlim (None): X-axis limits ([min, max])
    ylim (None): Y-axis limits ([min, max])
    xscale ('linear'): X-axis scale
    yscale ('linear'): Y-axis scale
    width (4): width in inches
    height (3): height in inches
    dpi (100): resolution in dpi
    hold (False): if False, figure will be cleared each time plot is called

    Widget attributes:
    -----------------
    figure: instance of matplotlib.figure.Figure
    axes: figure axes
    """
    def __init__(self, parent=None, title='', xlabel='', ylabel='',
                 xlim=None, ylim=None, xscale='linear', yscale='linear',
                 width=10, height=3, dpi=100, hold=False):
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        #Top Figure
        self.axes0 = self.figure.add_subplot(211)
        #self.axes0.set_title('Wavelength shift')
        self.axes0.set_xlabel(xlabel)
        self.axes0.set_ylabel(ylabel)
        #Bot Figure
        self.axes1 = self.figure.add_subplot(212)
        #self.axes1.set_title('Peak width variation')
        self.axes1.set_xlabel(xlabel)
        self.axes1.set_ylabel(ylabel)

        if xscale is not None:
            self.axes0.set_xscale(xscale)
            self.axes1.set_xscale(xscale)
        if yscale is not None:
            self.axes0.set_yscale(yscale)
            self.axes1.set_yscale(yscale)
        if xlim is not None:
            self.axes0.set_xlim(*xlim)
            self.axes1.set_xlim(*xlim)
        if ylim is not None:
            self.axes0.set_ylim(*ylim)
            self.axes1.set_ylim(*ylim)
        #self.axes0.hold(hold)    #---BEN_EDIT, check later
        #self.axes1.hold(hold)    #---BEN_EDIT, check later

        Canvas.__init__(self, self.figure)
        self.setParent(parent)

        Canvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        Canvas.updateGeometry(self)

class OSACanvas(Canvas):
    """
    The same as TRCanvas
    """
    def __init__(self, parent=None, title='', xlabel='', ylabel='',
                 xlim=None, ylim=None, xscale='linear', yscale='linear',
                 width=10, height=3, dpi=100, hold=False):
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_title('FBG Reflected Spectrum Simulation')
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        if xscale is not None:
            self.axes.set_xscale(xscale)
        if yscale is not None:
            self.axes.set_yscale(yscale)
        if xlim is not None:
            self.axes.set_xlim(*xlim)
        if ylim is not None:
            self.axes.set_ylim(*ylim)
        #self.axes.hold(hold)    #---BEN_EDIT, check later

        Canvas.__init__(self, self.figure)
        self.setParent(parent)

        Canvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        Canvas.updateGeometry(self)
