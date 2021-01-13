"""
Classes for loading a QWidget designed in QT Designer as MainWindow, Dialog or Widget
Examples of how to use can be found in UseQtGuiLoader.py

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

from __future__ import print_function
from PyQt5 import QtCore, QtGui, QtWidgets
import os
import sys
import time
import importlib

class QtGuiLoader(object):
    def compile_ui(self, ui_module):
        basename = ui_module.__name__.replace(".",os.path.sep)
        ui_file = basename + ".ui"
        py_file = basename + ".py"
        if os.path.exists(ui_file):
            if not os.path.exists(py_file) or \
            os.path.getmtime(ui_file) > os.path.getmtime(py_file) or \
            os.path.getsize(py_file)==0:
                print("compile %s > %s" % (ui_file, py_file))
                # quick and dirty hack
                if os.name=='posix':
                    # someone has to make that more general
                    pyuic_path='/usr/lib/python2.7/dist-packages/PyQt5/uic/pyuic.py'
                else:
                    # windows
                    pyuic_path = os.path.join(os.path.dirname(sys.executable), 'Lib/site-packages/PyQt5/uic/pyuic.py')
                os.system("%s %s %s > %s" % (sys.executable, pyuic_path, ui_file, py_file))
        importlib.reload(ui_module)
        
    def connect_actions(self,action_receiver=None):
        for name, action in [(n,a) for n,a in vars(self.ui).items() if isinstance(a,QtWidgets.QAction)]:
            if action_receiver is None:
                action_receiver = self
            if hasattr(action_receiver, name):
                action.triggered.connect(getattr(action_receiver,name))
            elif action.receivers(QtCore.SIGNAL("triggered()"))==0:
                raise Warning("Action %s not connected. Method with name '%s' not found"%(action.text(), name))

    def setupUI(self, widget):
        self.ui.setupUi(widget)
        root_widgets = [w for w in widget.children() if w.__class__.__name__ == "QWidget"]
        if len(root_widgets) == 0 or widget.layout() is not None:
            self.ui_widget = self
        else:
            self.ui_widget = root_widgets[-1]
            g = QtWidgets.QGridLayout()
            if isinstance(self, QtWidgetLoader):
                g.setContentsMargins(0, 0, 0, 0)
                g.setSpacing(0)
            widget.setLayout(g)
            g.addWidget(self.ui_widget)

class QtGuiApplication():
    def __init__(self,ui_module):
        self.ui_module = ui_module
        self.app_filename = os.path.basename(sys.argv[0])
        self.app_name = os.path.splitext(self.app_filename)[0]
        if QtWidgets.QApplication.startingUp():
            self.app = QtWidgets.QApplication(sys.argv)
        self.compile_ui(self.ui_module)
        
    def save_settings(self):
        settings = QtCore.QSettings("QtGuiApplication", "%s_%s" % (self.app_name, self.__class__.__name__))
        settings.setValue(self.ui_module.__name__ +"/geometry", self.saveGeometry())
        
    def load_settings(self):
        settings = QtCore.QSettings("QtGuiApplication", "%s_%s" % (self.app_name, self.__class__.__name__))
        geometry = settings.value(self.ui_module.__name__ + "/geometry")
        try:
            geometry = geometry.toByteArray()
        except:
            pass  # Fails in PySide
        if geometry:
            self.restoreGeometry(geometry)
    
    def save_setting(self,key,value):
        settings = QtCore.QSettings("QtGuiApplication", "%s_%s" % (self.app_name, self.__class__.__name__))
        settings.setValue(self.ui_module.__name__ +"/" +key, value)
    
    def load_setting(self,key,default_value=None):
        settings = QtCore.QSettings("QtGuiApplication", "%s_%s" % (self.app_name, self.__class__.__name__))
        return settings.value(self.ui_module.__name__ + "/"+key,default_value)
    
class QtMainWindowLoader(QtGuiLoader,QtGuiApplication, QtWidgets.QMainWindow):
    def __init__(self, ui_module, parent=None, connect_actions=True):
        QtGuiApplication.__init__(self, ui_module)
        QtWidgets.QMainWindow.__init__(self, parent)
        if "Ui_Form" in dir(ui_module):
            print ("Ui_Form")
            self.ui = ui_module.Ui_Form()
            centralWidget = QtWidgets.QWidget(self)
            self.setCentralWidget(centralWidget)
            self.setupUI(centralWidget)        
        elif "Ui_MainWindow" in dir(ui_module):
            print ("Ui_MainWindow")
            self.ui = ui_module.Ui_MainWindow()
            self.ui.setupUi(self)
        if connect_actions: 
            print ("Connect actions")
            self.connect_actions()

    def start(self):
        self.load_settings()
        self.show()
        if hasattr(self, "app"):
            self.app.exec_()
            
    def terminate(self):
        QtWidgets.QApplication.quit()
        QtCore.QCoreApplication.quit()  
        sys.exit()
        
    def closeEvent(self, *args, **kwargs):
        self.save_settings()
        #Enable paste of clipboard after termination
        clipboard = QtWidgets.QApplication.clipboard()
        event = QtCore.QEvent(QtCore.QEvent.Clipboard)
        QtWidgets.QApplication.sendEvent(clipboard, event)
        return QtWidgets.QMainWindow.closeEvent(self, *args, **kwargs)

class QtDialogLoader(QtGuiLoader,QtGuiApplication, QtWidgets.QDialog):    
    def __init__(self, ui_module, parent, modal=True,connect_actions=True):
        QtGuiApplication.__init__(self, ui_module)
        QtWidgets.QDialog.__init__(self, parent)
        self.modal = modal
        self.setModal(modal)
        self.ui = ui_module.Ui_Form()
        self.setupUI(self)        
        if connect_actions: 
            self.connect_actions()
        
    def start(self):
        self.load_settings()
        self.show()
        self.raise_()
        if hasattr(self, "app"):
            return self.app.exec_()
        elif self.modal:
            return self.exec_()
            
    def hideEvent(self, *args, **kwargs):
        self.save_settings()
        return QtWidgets.QDialog.hideEvent(self, *args, **kwargs)
       
class QtWidgetLoader(QtGuiLoader,QtWidgets.QWidget):
    def __init__(self, ui_module,action_receiver=None,parent=None,connect_actions=True):
        if "ui_module" not in vars(self):
            QtWidgets.QWidget.__init__(self, parent)
            self.ui_module = ui_module
            self.compile_ui(ui_module)
            self.ui = ui_module.Ui_Form()
            try:
                self.setupUI(self)
            except:
                self.compile_ui(ui_module, True)
                self.ui = ui_module.Ui_Form()
                self.setupUI(self)
            if connect_actions: 
                self.connect_actions(action_receiver)


