#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-

# Creating the GUI


# Reading the data, processing, displaying them
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=18
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QPushButton, QMessageBox,QToolBar,QAction,QStatusBar,
                             QHBoxLayout, QVBoxLayout, QApplication, QListWidget,QSplitter)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

import os, sys
import numpy as np
old_settings = np.seterr(all='ignore')  # Avoid warning messages

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class PositionCanvas(MplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        t = np.arange(0.0, 3.0, 0.01)
        s = np.sin(2*np.pi*t)
        self.axes.plot(t, s)


class FluxCanvas(MplCanvas):
    """A canvas that updates itself every 3 seconds with a new plot."""

    #def __init__(self, *args, **kwargs):
    #    MyMplCanvas.__init__(self, *args, **kwargs)
    #    timer = QtCore.QTimer(self)
    #    timer.timeout.connect(self.update_figure)
    #    timer.start(3000)

    def compute_initial_figure(self):
        self.axes.set_xlim([0,20*20]) # Set at least 20 observations
        self.axes.set_ylabel('Flux [V/s]')
        self.axes.plot([0,0],[np.nan,np.nan],'.b')

    def updateFigure(self,nod,fn,spec):
        # get number of grating positions
        ng = (np.shape(spec))[1]
        coverage = 16+0.5*ng
        sp16 = np.arange(16)
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        #l = [random.randint(0, 10) for i in range(4)]
        labels = []
        for i in np.arange(len(fn)):
            if nod[i] == 'A':
                c = 'blue'
            else:
                c = 'red'
            s = spec[i]
            for j in np.arange(ng):
                self.axes.plot(sp16+i*coverage+j*0.5, s[j,:], color=c)
            self.axes.axvline(i*coverage, color='gray', linewidth=0.2)
        self.axes.axvline(len(fn)*coverage, color='gray', linewidth=0.2)
        labels = np.array(fn, dtype='str')    
        plt.xticks((np.arange(len(fn))+0.5)*coverage, labels, rotation='vertical',fontsize=15)
        self.axes.yaxis.grid(True)
        self.axes.autoscale(enable=True,axis='y')
        self.draw()

class myListWidget(QListWidget):

    def Clicked(self,item):
#      QMessageBox.information(self, "ListWidget", "You clicked: "+item.text())
        # 2000 means message erased after 2 seconds
        #self.parent().parent().statusBar().showMessage("Clicked item "+item.text(),2000)
        mw = self.parent().parent()
        mw.sb.showMessage("You selected the FileGroupID: "+item.text(),2000)
#        self.parent().parent().lf.hide()
        mw.lf.setVisible(False)
        # Trigger event related to item list ....
        print "item is ",item.text()
        mw.addObs(item.text())



class ApplicationWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        #        self.setWindowTitle("application main window")
        path0 = sys.path[0]

        # Start exploring directory
        from fifitools import exploreDirectory
        cwd = os.getcwd()
        self.files, self.start, self.fgid = exploreDirectory(cwd+"/")
        # Compile the list of unique File group IDs
        self.fgidList = list(set(self.fgid))
        # This step should be repeated in case new files appear, to update the plots

        
        # Start list of observations
        self.fileNames = []
        self.obs = []
        
        # Menu
        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&About', self.about)

        # Define main widget
        self.main_widget = QWidget(self)

        # Define sub widgets
        self.pc = PositionCanvas(self.main_widget, width=3, height=4, dpi=100)
        self.fc = FluxCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.mpl_toolbar = NavigationToolbar(self.pc, self)
        self.mpl_toolbar2 = NavigationToolbar(self.fc, self)

        # Actions
        exitAction = QAction(QIcon(path0+'/icons/exit24.png'), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(qApp.quit)
        hideAction = QAction(QIcon(path0+'/icons/list.png'), 'List of File Group IDs', self)
        hideAction.setShortcut('Ctrl+L')
        hideAction.triggered.connect(self.changeVisibility)

        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.addAction(exitAction)
        self.tb.addAction(hideAction)

        # Widget list
        self.lf = myListWidget()
        for item in self.fgidList:
            self.lf.addItem(item); 
        self.lf.setWindowTitle('File Group ID')
        self.lf.itemClicked.connect(self.lf.Clicked)
        
        self.sb = QStatusBar()
        

        # Layout
        mainLayout = QHBoxLayout(self.main_widget)
        splitter1 = QSplitter(Qt.Horizontal)

        radecLayout = QVBoxLayout()
        radecWidget = QWidget()
        radecWidget.setLayout(radecLayout)
        fluxWidget = QWidget()
        fluxLayout = QVBoxLayout()
        fluxWidget.setLayout(fluxLayout)

        radecLayout.addWidget(self.pc)
        radecLayout.addWidget(self.mpl_toolbar)
        radecLayout.addWidget(self.tb)
        fluxLayout.addWidget(self.fc)
        fluxLayout.addWidget(self.mpl_toolbar2)
        fluxLayout.addWidget(self.sb)

        splitter1.addWidget(radecWidget)
        splitter1.addWidget(fluxWidget)

        mainLayout.addWidget(splitter1)
        mainLayout.addWidget(self.lf)
        self.lf.hide()

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)


        self.sb.showMessage("Welcome to FIFI Monitor!", 10000)

    def changeVisibility(self):
        state = self.lf.isVisible()
        self.lf.setVisible(not state)

        
    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        file=open("copyright.txt","r")
        message=file.read()
        QtWidgets.QMessageBox.about(self, "About", message)

    def addObs(self, fileGroupId):
        from fifitools import readData, multiSlopes, Obs
        mask = self.fgid == fileGroupId
        selFiles = self.files[mask]
        selFileNames = [os.path.splitext(os.path.basename(f))[0] for f in selFiles]
        selFileNames = np.array(selFileNames)
        print "Files selected are "
        print selFileNames
        # Check if file names already appear in previous list, otherwise append them and read/process/display relative data
        for infile in selFileNames:
            if infile not in self.fileNames:
                print "Reading file: ", infile
                self.fileNames.append(infile)
                aor, hk, gratpos, flux = readData(infile+".fits")
                spectra = multiSlopes(flux)
                spectrum=np.nanmedian(spectra,axis=2)
                detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum = aor
                obsdate, coords, offset, angle, za, altitude, wv = hk
                self.obs.append(Obs(spectrum,coords,offset,angle,altitude,za,wv,nodbeam,filegpid,filenum,gratpos))
                # Create lists of spectra, nod, and file number
                fn,nod,spectra = zip(*((o.n,o.nod,o.spec)  for o in self.obs if o.fgid == fileGroupId))
                # Display the data
                self.fc.updateFigure(nod,fn,spectra)

""" Main code """
        
if __name__ == '__main__':
    qApp = QApplication(sys.argv)

    screen_resolution = qApp.desktop().screenGeometry()
    width, height = screen_resolution.width(), screen_resolution.height()
    aw = ApplicationWindow()
    aw.setGeometry(100, 100, width*0.9, width*0.35)
    progname = 'FIFI Monitor'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    sys.exit(qApp.exec_())
