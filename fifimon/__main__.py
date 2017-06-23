#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-

# Creating the GUI


# Reading the data, processing, displaying them
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
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
from PyQt5.QtCore import Qt, QThread, QTimer


import os, sys
import numpy as np
old_settings = np.seterr(all='ignore')  # Avoid warning messages

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
#        self.axes = self.fig.add_axes([0.1, 0.1, 0.8, 0.8])

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.compute_initial_figure()

    def compute_initial_figure(self):
        pass


class PositionCanvas(MplCanvas):
    """Simple canvas with a sine plot."""
    
    def compute_initial_figure(self,ra=None,dec=None):
        from astropy.wcs import WCS
        self.w = WCS(naxis=2)
        self.w.wcs.ctype=["RA---TAN","DEC--TAN"]
        self.w.wcs.crpix=[1,1]
        self.w.wcs.cdelt = np.array([1,1])
        self.w.wcs.pc[0]=[1,0]
        self.w.wcs.pc[1]=[-0,1]
        #self.w.wcs.naxis1=1000
        #self.w.wcs.naxis2=1000
        self.w.wcs.cunit=['deg','deg']
        if ra == None:
            ra = 0
        if dec == None:
            dec = 0
        self.fig.clear()
        self.w.wcs.crval=[ra,dec]
        #        self.axes = self.fig.add_subplot(111, projection=self.w)
        #        #        self.axes.set_xlim([0,20*20]) # Set at least 20 observations

        self.axes = self.fig.add_subplot(111, projection=self.w)
        self.t = self.axes.get_transform(self.w)
        self.axes.axis('equal')
        self.axes.coords[0].set_major_formatter('hh:mm:ss')
        self.axes.plot(np.nan,np.nan)

    def updateFigure(self, nod, ra, dec, dx, dy, angle):
        n = len(nod)
        i = n-1
        ra=np.asarray(ra)
        dec=np.asarray(dec)
        dx=np.asarray(dx)
        dy=np.asarray(dy)
        x = ra+dx-self.w.wcs.crval[0]
        y = dec+dy-self.w.wcs.crval[1]
        colors = ['blue' if a =='A' else 'red' for a in nod]
        self.axes.autoscale(enable=True, axis='both')
        self.axes.scatter(x,y,facecolors='none',edgecolors='none',marker=(4,0,angle[0]),s=0,transform=self.t)
        self.axes.clear()
        # Get new limits
        x0,y0 = self.axes.transAxes.transform(([0,0]))
        x1,y1 = self.axes.transAxes.transform(([1,1]))
        xl0,xl1 = self.axes.get_xlim()
        yl0,yl1 = self.axes.get_ylim()
        x2p = (x1-x0)/(xl0-xl1)
        #y2p = (y1-y0)/(yl1-yl0)


        bbox = self.axes.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        print 'height is: ',height
        print 'delta y is: ',yl1-yl0
        y2p = height*72/(yl1-yl0)*self.axes.get_data_ratio()  # from data to points (heights in inch, 72 points per inch)

        print 'pixels ',self.fig.dpi*height,yl1-yl0
        # Adjust marker size 1 arcmin for red, 1/2 arcmin for blue
        radiant = np.pi/180.
        theta = (angle[0]+45.)*radiant
        fact = np.sin(theta)+np.cos(theta)
        size = (0.5/60.*y2p*fact)**2
        deltay = (yl1-yl0)*3600.
        
        size = (30./deltay * height*72 *self.axes.get_data_ratio()*fact)**2
        print "size ",size
        # size should be in points^2, a points is 1/72 of an inch
        self.axes.scatter(x,y,facecolors='none',edgecolors=colors,marker=(4,0,45+angle[0]),s=size,transform=self.t)
        self.axes.scatter(x[i],y[i],facecolors='none',edgecolors='green',marker=(4,0,45+angle[0]),s=size,transform=self.t)


        # Invert x axis
        xlim = self.axes.get_xlim()
        if xlim[0] < xlim[1]:
            self.axes.set_xlim([xlim[1],xlim[0]])

        self.draw()

class FluxCanvas(MplCanvas):
    """A canvas that updates itself every 3 seconds with a new plot."""

    #def __init__(self, *args, **kwargs):
    #    MyMplCanvas.__init__(self, *args, **kwargs)
    #    timer = QTimer(self)
    #    timer.timeout.connect(self.update_figure)
    #    timer.start(3000)

    def compute_initial_figure(self):
        self.axes = self.fig.add_subplot(111)
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
        n = len(nod)

        i = n-1
        if nod[i] == 'A':
            c = 'blue'
        else:
            c = 'red'
        s = spec[i]
        for j in np.arange(ng):
            self.axes.plot(sp16+i*coverage+j*0.5, s[j,:], color=c)
        self.axes.axvline(i*coverage, color='gray', linewidth=0.2)            
        self.axes.axvline(len(fn)*coverage, color='gray', linewidth=0.2)
        self.axes.yaxis.grid(True)
        self.axes.autoscale(enable=True,axis='y')
        if i*coverage > 400:
            self.axes.autoscale(enable=True, axis='x')
        labels = np.array(fn, dtype='str')    
        #plt.xticks((np.arange(len(fn))+0.5)*coverage, labels, rotation='vertical',fontsize=15)
        self.axes.set_xticks((np.arange(len(fn))+0.5)*coverage)
        self.axes.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)

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
    '''
    Notes about threading pyqt using also multiprocessing
    https://stackoverflow.com/questions/15675043/multiprocessing-and-gui-updating-qprocess-or-multiprocessing
    
    Other examples:
    https://stackoverflow.com/questions/15698251/multiprocessing-gui-schemas-to-combat-the-not-responding-blocking
    '''

            
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
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
                                 Qt.CTRL + Qt.Key_Q)
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
        exitAction.triggered.connect(self.close)
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

        # Periodical update
        timer = QTimer(self)
        timer.timeout.connect(self.update_fifimon)
        timer.start(3000)

        

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

    def addObs(self, fileGroupId = None):
        from fifitools import readData, multiSlopes, Obs
        import os
        import psutil
        firstRun = 0
        if fileGroupId != None:
            self.fileGroupId = fileGroupId
            firstRun = 1

            
        #print "selected file group id is: ", self.fileGroupId
        mask = self.fgid == self.fileGroupId
        selFiles = self.files[mask]
        selFileNames = [os.path.splitext(os.path.basename(f))[0] for f in selFiles]
        selFileNames = np.array(selFileNames)
        #print "Files selected are "
        #print selFileNames
        # Check if file names already appear in previous list, otherwise append them and read/process/display relative data

        if firstRun:
            # Grab RA-Dec of first file and restart the WCS
            aor, hk, gratpos, flux = readData(selFileNames[0]+".fits")
            obsdate, coords, offset, angle, za, altitude, wv = hk
            self.pc.compute_initial_figure(coords[0],coords[1])


        
        for infile in selFileNames:
            if infile not in self.fileNames:
                print "Reading file: ", infile
                self.fileNames.append(infile)
                aor, hk, gratpos, flux = readData(infile+".fits")
                spectra = multiSlopes(flux)
                spectrum = np.nanmedian(spectra,axis=2)
                detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum = aor
                obsdate, coords, offset, angle, za, altitude, wv = hk
                self.obs.append(Obs(spectrum,coords,offset,angle,altitude,za,wv,nodbeam,filegpid,filenum,gratpos))
                # Create lists of spectra, nod, positions, and file numbers
                fn,nod,ra,dec,x,y,angle,spectra = zip(*((o.n,o.nod,o.ra,o.dec,o.x,o.y,o.angle,o.spec)  for o in self.obs if o.fgid == self.fileGroupId))
                # Display the data
                self.fc.updateFigure(nod,fn,spectra)
                self.pc.updateFigure(nod,ra,dec,x,y,angle)
                QApplication.processEvents()
                pid = os.getpid()
                py = psutil.Process(pid)
                #memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
                #print('memory use:', memoryUse)

    def update_fifimon(self):
        '''
        This module is called automatically every n seconds (n=10, typically) to look
        for new files in the directory, update the lists and, eventually, the plot
        '''
        from fifitools import exploreDirectory
        cwd = os.getcwd()
        files, start, fgid = exploreDirectory(cwd+"/")
        # Check if new fileGroupID or files added
        for f,fg in zip(files,fgid):
            if f not in self.files:
                self.files = np.append(self.files,f)
                self.fgid = np.append(self.fgid, fg)
                # update plot (if same fileGroupID)
                self.addObs()
            if fg not in self.fgidList:
                # update fileGroupID list    
                self.fgidList.append(fg)
                # update widget list
                self.lf.addItem(fg); 


                
""" Main code """
        
def main():
    app = QApplication(sys.argv)
    screen_resolution = app.desktop().screenGeometry()
    width, height = screen_resolution.width(), screen_resolution.height()
    aw = ApplicationWindow()
    aw.setGeometry(100, 100, width*0.9, width*0.35)
    progname = 'FIFI Monitor'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()

# Ensure that the app is created once 
if __name__ == '__main__':
    main()
