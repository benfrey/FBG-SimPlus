""" 
Build program GUIs and check input fields for valid format

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
from __future__ import print_function
import numpy as np
import scipy as sp
import sympy
import matplotlib.pyplot as plt
import time
import webbrowser
import sys

#Other Classes used
import GUI.MyPlotMainWindowUI
from OSASimulation import *
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

    """--------------------------Tool Bar Buttons----------------------"""
    #Push Button:Wedbsite
    def actionCopyrigth(self):
        msgbox=QtWidgets.QMessageBox()
        msgbox.setIcon(1)
        msgbox.setWindowTitle('About the Software')
        msgbox.setWindowIcon(QtGui.QIcon("GUI/resource/help.png"))
        msgbox.setTextFormat(1)
        msgbox.setText("FBG-SimPlus is a software that simulates the reflected spectrum of a FBG array placed virtually along a predefined path in a Finite Element Model.<br> <br>\
        This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.<br> <br>\
        If you found this software useful in any kind, please cite it together with article DOI from our 2021 letter as well as Dr. Pereira's previous article(s).")
        msgbox.exec_()

    def actionAbouttheauthor(self):
        msgbox=QtWidgets.QMessageBox()
        msgbox.setIcon(1)
        msgbox.setWindowTitle('About the author')
        msgbox.setWindowIcon(QtGui.QIcon("GUI/resource/author.png"))
        msgbox.setTextFormat(1)
        msgbox.setText("Hi, <br> <br>\
        My name is Ben Frey and I am an undergraduate student at the University of St. Thomas in St. Paul, MN.</br>\
        Work on FBG-SimPlus was conducted in order to add further functionality to Dr. Gilmar Pereira's work on FBG_SiMul.</br>\
        The main addition to this software includes the ability to model the dependence of temperature on the reflected</br>\
        spectrum produced by the FBG Array. The code has also been heavily modified to allow users greater flexibility</br>\
        in viewing contribution components such as strain, stress, and temperature. Please visit the software</br>\
        documentation before program usage. Thanks!")
        msgbox.exec_()

    def actionSoftwareDoc(self):
        webbrowser.open('https://github.com/benfrey/FBG-SimPlus/blob/master/documentation.pdf')

    def actionExitProgram(self):
        self.terminate()

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

    #Push file format help button
    def actionHelpFileFormatOSA(self):
        msgbox=QtWidgets.QMessageBox()
        msgbox.setIcon(1)
        msgbox.setWindowTitle('File Format Help')
        msgbox.setWindowIcon(QtGui.QIcon("GUI/resources/document.png"))
        msgbox.setTextFormat(1)
        msgbox.setText('The expected input format for each column is given by the following:<br>\
        1.) Distance from start of path in [m] or [mm]<br>\
        2.) True strain (logarithmic strain) in longitudinal direction (x direction)<br>\
        3.) True strain (logarithmic strain) in y direction<br>\
        4.) True strain (logarithmic strain) in z direction<br>\
        5.) True stress (normal from stress tensor) in longitudinal direction (x direction) in [Pa]<br>\
        6.) True stress (normal from stress tensor) in y direction in [Pa]<br>\
        7.) True stress (normal from stress tensor) in z direction in [Pa]<br>\
        8.) Temperature at point in [K]<br><br>\
        For more information, please see the provided documentation.')
        msgbox.exec_()

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
    
    #-------------------Simulation Contribution-----------------------
    def actionToggleUniform(self):
        self.ui.radioButton_TypeNonUniform.setChecked(False)
     
    def actionToggleNonUniform(self):
        self.ui.radioButton_TypeUniform.setChecked(False)
        
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
        #Check if contribution settings were correct
        if self.ui.lineEdit_EmulateTemperature.text()=='' \
            or self.ui.lineEdit_HostThermalExpansion.text()=='':
                self.ui.MessageBoard.insertPlainText('>>ERROR in (2)!!: Please insert all contribution parameters. \n')
                self.ui.progressBar.setValue(0)
                return
        try:
            float(self.ui.lineEdit_EmulateTemperature.text())
            float(self.ui.lineEdit_HostThermalExpansion.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (2)!!: Invalid format!! --Contribution parameters. \n')
            self.ui.progressBar.setValue(0)
            return    
        
        #Check if all the Simulation Parameters were inserted and are float
        if self.ui.lineEdit_InitialRefractiveIndex.text()=='' \
            or self.ui.lineEdit_MeanChangeRefractiveIndex.text()=='' \
            or self.ui.lineEdit_FringeVisibility.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP11.text()=='' \
            or self.ui.lineEdit_DirectionalRefractiveP12.text()=='' \
            or self.ui.lineEdit_YoungsModule.text()=='' \
            or self.ui.lineEdit_PoissionsCoefficient.text()=='' \
            or self.ui.lineEdit_MinBandWidth.text()=='' \
            or self.ui.lineEdit_MaxBandWidth.text()==''\
            or self.ui.lineEdit_ThermoOptic.text()=='' \
            or self.ui.lineEdit_AmbientTemperature.text()=='':
                self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert all optical fibre parameters. \n')
                self.ui.progressBar.setValue(0)
                return
        try:
            float(self.ui.lineEdit_InitialRefractiveIndex.text())
            float(self.ui.lineEdit_MeanChangeRefractiveIndex.text())
            float(self.ui.lineEdit_FringeVisibility.text())
            float(self.ui.lineEdit_DirectionalRefractiveP11.text())
            float(self.ui.lineEdit_DirectionalRefractiveP12.text())
            float(self.ui.lineEdit_YoungsModule.text())
            float(self.ui.lineEdit_PoissionsCoefficient.text())
            float(self.ui.lineEdit_MinBandWidth.text())
            float(self.ui.lineEdit_MaxBandWidth.text())
            float(self.ui.lineEdit_ThermoOptic.text())
            float(self.ui.lineEdit_AmbientTemperature.text())
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
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Please insert the simulation resolution. \n')
            self.ui.progressBar.setValue(0)
            return
        try: float(self.ui.lineEdit_SimulationResolution.text())
        except:
            self.ui.MessageBoard.insertPlainText('>>ERROR in (3)!!: Invalid format!! --Simulation resolution. \n')
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
        InitialRefractiveIndex=float(self.ui.lineEdit_InitialRefractiveIndex.text())
        MeanChangeRefractiveIndex=float(self.ui.lineEdit_MeanChangeRefractiveIndex.text())
        FringeVisibility=float(self.ui.lineEdit_FringeVisibility.text())
        DirectionalRefractiveP11=float(self.ui.lineEdit_DirectionalRefractiveP11.text())
        DirectionalRefractiveP12=float(self.ui.lineEdit_DirectionalRefractiveP12.text())
        YoungsModule=float(self.ui.lineEdit_YoungsModule.text())
        PoissonsCoefficient=float(self.ui.lineEdit_PoissionsCoefficient.text())
        MinBandWidth=float(self.ui.lineEdit_MinBandWidth.text())
        MaxBandWidth=float(self.ui.lineEdit_MaxBandWidth.text())
        ThermoOptic=float(self.ui.lineEdit_ThermoOptic.text())
        AmbientTemperature=float(self.ui.lineEdit_AmbientTemperature.text())
        SimulationResolution=float(self.ui.lineEdit_SimulationResolution.text())
        self.EmulateTemperature=float(self.ui.lineEdit_EmulateTemperature.text())
        self.HostThermalExpansionCoefficient=float(self.ui.lineEdit_HostThermalExpansion.text())

        self.ui.progressBar.setValue(50)#Progress Bar

        #Generate the undeformed FBG signal
        if self.ui.checkBox_UndeformedSignal.isChecked():
            self.Osa.UndeformedFBG(SimulationResolution,MinBandWidth,MaxBandWidth,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,PoissonsCoefficient,self.FBGOriginalWavel)

        self.ui.progressBar.setValue(75)#Progress Bar

        #FBG sensor direction 0-xx 1-yy 2-zz
        self.FBGDirection=0

        #Determine contributions to add
        
        #--- STRAIN ---
        if self.ui.checkBox_TypeStrain.isChecked():
            if self.ui.radioButton_TypeUniform.isChecked():
                self.StrainType = 1
            if self.ui.radioButton_TypeNonUniform.isChecked():
                self.StrainType = 2
        else:
            self.StrainType = 0

        #--- STRESS ---
        if self.ui.checkBox_TypeTransverseStress.isChecked():
            self.StressType = 1
        else:
            self.StressType = 0
            
        #--- EMULATE TEMPERATURE ---
        if self.ui.checkBox_TypeEmulateTemperature.isChecked() == False:
             self.HostThermalExpansionCoefficient = 0
             self.EmulateTemperature = -1.0
             
        #Generate the FBG signal
        self.Osa.DeformedFBG(SimulationResolution,MinBandWidth,MaxBandWidth,AmbientTemperature,InitialRefractiveIndex,MeanChangeRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissonsCoefficient,ThermoOptic,self.StrainType,self.StressType,self.EmulateTemperature,self.HostThermalExpansionCoefficient,self.FBGOriginalWavel)

        #Generates the summarized Data
        self.Osa.FBGOutputSum(AmbientTemperature,InitialRefractiveIndex,FringeVisibility,DirectionalRefractiveP11,DirectionalRefractiveP12,YoungsModule,PoissonsCoefficient,ThermoOptic,self.StrainType,self.StressType,self.EmulateTemperature,self.HostThermalExpansionCoefficient,self.FBGOriginalWavel)

        self.ui.progressBar.setValue(100) #Progress Bar
        #Message
        self.ui.MessageBoard.insertPlainText(">> The FBG reflected spectrum was successfully simulated. \n ")

    #Push Button: Load/Plot
    def actionPlotOSA(self):
        #Check if data was generated
        try:PlotWindowOSA(None,False,self.Osa,self.StressType).start()
        except:
            self.ui.MessageBoard.insertPlainText(">> Error!! Please simulated the FBG spectrum before ploting. \n ")

    #Push Button: Save as file
    def actionSaveAsFileOSA(self):
        fsize=None
        #Size of the file
        #try:fsize=len(self.Osa.OReflect['wavelength'])
        #except:pass
        #try:fsize=len(self.Osa.USReflect['wavelength'])
        #except:pass
        #try:fsize=len(self.Osa.NUSReflect['wavelength'])
        #except:pass
        try:fsize=len(self.Osa.YReflect['wavelength'])
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

                """
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
                """

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
    
"""----------------------------------------------------------------------------
-------------------------------------------------------------------------------
---------------------------PLOT CLass------------------------------------------
-------------------------------------------------------------------------------
"""

"""----------------------------Class to plot: OSA tab------------------------
"""
import GUI.PlotWindow_OSA
class PlotWindowOSA(QtDialogLoader):
    def __init__(self, parent, modal,Osa,StressType):
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
        self.ui.checkBox_ViewSplitWaves.clicked.connect(self.actionUpdatePlot)
        self.ui.checkBox_legend.clicked.connect(self.actionUpdatePlot)
        self.ui.checkBox_Grid.clicked.connect(self.actionUpdatePlot)
        self.ui.doubleSpinBox_LineWidth.valueChanged.connect(self.actionUpdatePlot)
        self.ui.comboBox_LineColorUndeformed.currentIndexChanged.connect(self.actionUpdatePlot)
        self.ui.comboBox_LineColorSimulated.currentIndexChanged.connect(self.actionUpdatePlot)
        #Plot the Spectrum
        #self.osaplot.axes.hold(True)  #---BEN_EDIT, check later
        self.Plot()
                
        #FBG Output Table
        self.ui.textEdit_FBGoutput.insertPlainText('FBG # | Orig. Wavelength [nm] | Shift [nm] | Width Variation [nm]\n')
        for b in np.arange(0,self.Osa.NumberFBG):
            self.ui.textEdit_FBGoutput.insertPlainText('FBG: '+str(b)+' | '+str(self.Osa.FBGOriginalWavel[b])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveShift'][0])+' | '+str(self.Osa.FBGOutSum['FBG'+str(b+1)]['WaveWidth'][0])+'\n')

        print("here")

        #Enable button to view split waves due to Birefringence
        if StressType == 1:
            print("here2")
            self.ui.checkBox_ViewSplitWaves.setEnabled(True)

        print("Plot has been displayed!")

    def Plot(self):
        #Plot the Undeformed shape
        try:self.osaplot.axes.plot(self.Osa.OReflect['wavelength'],self.Osa.OReflect['reflec'],color=self.LineColorUndeformed,linewidth=self.linewidth, label="Undeformed FBG Spectrum")
        except:pass
 
        #Plot the deformed shape
        try:self.osaplot.axes.plot(self.Osa.DReflect['wavelength'],self.Osa.DReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Deformed FBG Spectrum")
        except:pass
        
        if self.ui.checkBox_ViewSplitWaves.isChecked():
            #Plot the individual contributions
            try:self.osaplot.axes.plot(self.Osa.YReflect['wavelength'],self.Osa.YReflect['reflec'],color='red',linewidth=self.linewidth, label="Y-Wave Contribution")
            except:pass
            try:self.osaplot.axes.plot(self.Osa.ZReflect['wavelength'],self.Osa.ZReflect['reflec'],color='red',linewidth=self.linewidth, label="Z-Wave Contribution")
            except:pass
        #Plot the individual contributions
        #try:self.osaplot.axes.plot(self.Osa.YReflect['wavelength'],self.Osa.YReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Y-Wave Contribution")
        #except:pass
        #try:self.osaplot.axes.plot(self.Osa.ZReflect['wavelength'],self.Osa.ZReflect['reflec'],color=self.LineColorSimulated,linewidth=self.linewidth, label="Z-Wave Contribution")
        #except:pass
                
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
            self.osaplot.axes.legend(fontsize=10, loc="upper right")
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
