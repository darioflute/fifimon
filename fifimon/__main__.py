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
from matplotlib.patches import Rectangle

# Make sure that we are using QT5
#matplotlib.use('Qt5Agg')
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
#old_settings = np.seterr(all='ignore')  # Avoid warning messages
import warnings
#with warnings.catch_warnings():
#    warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
# To avoid excessive warning messages
warnings.filterwarnings('ignore')

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
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        self.cidPress = self.fig.canvas.mpl_connect('scroll_event', self.onWheel)
    

    def onWheel(self, event):
        eb = event.button
        curr_xlim = self.axes.get_xlim()
        curr_ylim = self.axes.get_ylim()
        curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
        curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        self.axes.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.axes.set_ylim([curr_y0-new_height,curr_y0+new_height])
        # Adjust size of array
        
        self.draw()

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
        self.rects = []    
        self.w.wcs.crval=[ra,dec]
        #        self.axes = self.fig.add_subplot(111, projection=self.w)
        #        #        self.axes.set_xlim([0,20*20]) # Set at least 20 observations

        self.axes = self.fig.add_subplot(111, projection=self.w)
        self.t = self.axes.get_transform(self.w)
        self.axes.axis('equal')
        self.axes.coords[0].set_major_formatter('hh:mm:ss')
        self.axes.plot(np.nan,np.nan)

    def updateFigure(self, nod, fn, ra, dec, dx, dy, angle, infile):
        n = len(nod)
        # Find index of infile
        i = fn.index(int(infile[:5]))
        
        ra=np.asarray(ra)
        dec=np.asarray(dec)
        dx=np.asarray(dx)
        dy=np.asarray(dy)
        x = ra+dx-self.w.wcs.crval[0]
        y = dec+dy-self.w.wcs.crval[1]
        colors = ['blue' if a =='A' else 'red' for a in nod]
        
        #self.axes.scatter(x,y,facecolors='none',edgecolors='none',marker=(4,0,angle[0]),s=0,transform=self.t)
        #self.axes.clear()
        # Get new limits
        #x0,y0 = self.axes.transAxes.transform(([0,0]))
        #x1,y1 = self.axes.transAxes.transform(([1,1]))
        #xl0,xl1 = self.axes.get_xlim()
        #yl0,yl1 = self.axes.get_ylim()
        #bbox = self.axes.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
        #width, height = bbox.width, bbox.height
        # The size of the marker is the side of a square circumscribing the rotated square
        #radiant = np.pi/180.
        #theta = (angle[0]+45.)*radiant
        #fact = np.sin(theta)+np.cos(theta)
        #deltay = (yl1-yl0)*3600.
        ## Adjust marker size 1 arcmin for red, 30 arcsec for blue
        ## size should be in points^2, a points is 1/72 of an inch
        #size = (30./deltay * height*72 *self.axes.get_data_ratio()*fact)**2
        #self.arrays = self.axes.scatter(x,y,facecolors='none',edgecolors=colors,marker=(4,0,45+angle[0]),s=size,transform=self.t)
        ##        if i == (len(fn)-1):
        ##            self.axes.scatter(x[i],y[i],facecolors=colors[i],edgecolors=colors[i],marker=(4,0,45+angle[0]),s=size,transform=self.t)
            
        # Add rectangle patch
        side = 30./3600.
        theta = -angle[i]*np.pi/180.
        dx = side*0.5*(np.cos(theta)-np.sin(theta))
        dy = side*0.5*(np.sin(theta)+np.cos(theta))
        rect = Rectangle((x[i] - dx, y[i] - dy), side,side,angle=-angle[i],fc='none',ec=colors[i])
        patch = self.axes.add_patch(rect)
        # collect the patches to modify them later
        self.rects.append(patch)
        self.axes.autoscale(enable=True, axis='both')
        # Invert x axis
        xlim = self.axes.get_xlim()
        if xlim[0] < xlim[1]:
            self.axes.set_xlim([xlim[1],xlim[0]])
        # Update figure
        self.draw()

class FluxCanvas(MplCanvas):
    """A canvas that updates itself every 3 seconds with a new plot."""

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        #    timer = QTimer(self)
        #    timer.timeout.connect(self.update_figure)
        #    timer.start(3000)
        # Mouse wheel zooming
        self.cidPress = self.fig.canvas.mpl_connect('scroll_event', self.onWheel)
        #self.picPress = self.fig.canvas.mpl_connect('pick_event', self.onPick)
        self.picPress = self.fig.canvas.mpl_connect('button_press_event', self.onPick)
        self.cidmotion = self.fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
        self.picRelease = self.fig.canvas.mpl_connect('button_release_event', self.onRelease)

        # Factor controlling the width of the figure
        self.w = 10
        self.dragged = None

        
    def onWheel(self, event):
        eb = event.button
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        self.w *= factor

        modifiers = QApplication.keyboardModifiers()
        if modifiers == QtCore.Qt.ShiftModifier:
            curr_ylim = self.axes.get_ylim()
            curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
            new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
            self.axes.set_ylim([curr_y0-new_height,curr_y0+new_height])
        else:
            curr_xlim = self.axes.get_xlim()
            curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
            new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
            self.axes.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.draw()

    def onPick(self, event):
        if event.button == 2:    
            self.dragged = event
            self.pick_pos = (event.xdata, event.ydata)
            print "pick position: ", self.pick_pos

    def onMotion(self, event):
        if self.dragged is not None and self.pick_pos[0] is not None:
            #old_pos = self.dragged.get_position()
            new_pos = (event.xdata,event.ydata)
            deltax = new_pos[0] - self.pick_pos[0]
            deltay = new_pos[1] - self.pick_pos[1]
            curr_xlim = self.axes.get_xlim()
            curr_ylim = self.axes.get_ylim()
            self.axes.set_xlim(curr_xlim-deltax)
            self.axes.set_ylim(curr_ylim-deltay)
            self.pick_pos = new_pos
            # Not clear how to use blit
            #            self.draw()
            #self.fig.canvas.blit(self.axes.bbox)
            #self.update()
            self.draw_idle()
        return True

    def onRelease(self, event):
        if self.dragged is not None:
            #deltax = event.xdata - self.pick_pos[0]
            #deltay = event.ydata - self.pick_pos[1]
            #print "delta ", deltax,deltay
            #curr_xlim = self.axes.get_xlim()
            #curr_ylim = self.axes.get_ylim()
            #self.axes.set_xlim(curr_xlim-deltax)
            #self.axes.set_ylim(curr_ylim-deltay)
            self.dragged = None
            self.draw()
        return True



    def compute_initial_figure(self,fileGroupId=None):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.set_xlim([0,20*20]) # Set at least 20 observations
        self.axes.set_ylabel('Flux [V/s]')
        self.axes.plot([0,0],[np.nan,np.nan],'.b')
        if fileGroupId is not None:
            self.fig.suptitle(fileGroupId)


    def updateFigure(self,nod,fn,spec,infile):
        # get number of grating positions
        ng = (np.shape(spec))[1]
        coverage = 16+0.5*ng
        sp16 = np.arange(16)
        n = len(nod)
        #print fn
        i = fn.index(int(infile[:5]))
        colors = ['blue' if a =='A' else 'red' for a in nod]
        if nod[i] == 'A':
            c = 'blue'
        else:
            c = 'red'

        s = spec[i]
        for j in np.arange(ng):
            self.axes.plot(sp16+i*coverage+j*0.5, s[j,:], color=colors[i])

        self.axes.axvline(i*coverage, color='gray', linewidth=0.2)            
        self.axes.axvline(len(fn)*coverage, color='gray', linewidth=0.2)
        self.axes.yaxis.grid(True)
        labels = np.array(fn, dtype='str')    
        self.axes.set_xticks((np.arange(len(fn))+0.5)*coverage)
        self.axes.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation
        self.axes.set_xlim([coverage*(i-1.5*self.w),coverage*(i+0.5*self.w)]) # Set at least 20 observations
        self.axes.autoscale(enable=True,axis='y')
        self.draw()

    def updateFigureMedians(self,nod,fn,spec,infile):
        # get number of grating positions
        ng = (np.shape(spec))[1]
        coverage = 16*ng
        sp16 = np.arange(16)
        n = len(nod)
        i = fn.index(int(infile[:5]))
        colors = ['blue' if a =='A' else 'red' for a in nod]
        if nod[i] == 'A':
            c = 'blue'
        else:
            c = 'red'

        s = spec[i]
        x = []
        y = []
        for j in np.arange(ng):
            n0 = int((ng-j)*0.5)
            n1 = int(16+ng*0.5-n0)
            x.append(i*ng+j)
            y.append(np.median(s[j,n0:n1]))
            #            self.axes.plot(sp16+i*coverage+j*16, s[j,:], color=colors[i])

        self.axes.plot(x,y,color=colors[i] )

        self.axes.axvline(i*ng, color='gray', linewidth=0.2)            
        self.axes.axvline(len(fn)*ng, color='gray', linewidth=0.2)
        self.axes.yaxis.grid(True)
        labels = np.array(fn, dtype='str')    
        self.axes.set_xticks((np.arange(len(fn))+0.5)*ng)
        #self.axes.set_xticks((np.arange(len(fn))+0.5)*coverage)
        self.axes.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation
        self.axes.set_xlim([ng*(i-1.5*self.w),ng*(i+0.5*self.w)]) # Set at least 20 observations
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


# We should probably add a threaded class to add observations.
# In this way, when we are adding observations, the GUI is not blocked and the plot is continuously updated.        
# check example at: https://nikolak.com/pyqt-threading-tutorial/
# class AddObs(QThread):
#     sec_signal = pyqtSignal(str)
#     def __init__(self, parent=None):
#         super(Logger, self).__init__(parent)
#         self.current_time = 0
#         self.go = True

#     def __del__(self):
#         self.wait()
        
#     def run(self):
#         #this is a special fxn that's called with the start() fxn
#         while self.go:
#             time.sleep(1)
#             self.sec_signal.emit(str(self.current_time))
#             self.current_time += 1
    


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

        # Background color
        self.setStyleSheet("""
        QMainWindow {
        background-color: 'lemon chiffon';
        }
        QMenu {
        background-color: 'light goldenrod yellow';
        color: 'black';
        }
        QMenuBar {
        background-color: 'light goldenrod yellow';
        }
        """)
        # Set window background color
#        self.setAutoFillBackground(True)
#        p = self.palette()
#        p.setColor(self.backgroundRole(), Qt.white)
#        self.setPalette(p)
        
        # Start exploring directory
        from fifitools import exploreDirectory
        cwd = os.getcwd()
        self.files, self.start, self.fgid = exploreDirectory(cwd+"/")
        # Compile the list of unique File group IDs
        self.fgidList = list(set(self.fgid))
        
        # Start list of observations
        self.fileGroupId = None
        self.fileNames = []
        self.obs = []
        
        # Menu
        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit, Qt.CTRL + Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&About', self.about)

        # Define main widget
        self.main_widget = QWidget(self)

        # Define sub widgets
        self.pc = PositionCanvas(self.main_widget, width=5, height=4, dpi=100)
        self.fc = FluxCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.mpl_toolbar = NavigationToolbar(self.pc, self)
        self.mpl_toolbar.pan('on')
        self.mpl_toolbar2 = NavigationToolbar(self.fc, self)
        self.mpl_toolbar2.pan('on')

        # Actions
        exitAction = QAction(QIcon(path0+'/icons/exit.png'), 'Exit', self)
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
        toolWidget = QWidget()
        toolLayout = QHBoxLayout()
        toolWidget.setLayout(toolLayout)
        
        radecLayout.addWidget(self.pc)
        radecLayout.addWidget(self.mpl_toolbar)
        fluxLayout.addWidget(self.fc)

        #fluxLayout.addWidget(self.mpl_toolbar2)
        toolLayout.addWidget(self.mpl_toolbar2)
        toolLayout.addWidget(self.tb)
        #radecLayout.addWidget(self.tb)
        fluxLayout.addWidget(toolWidget)
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
        path0 = sys.path[0]
        file=open(path0+"/copyright.txt","r")
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
        if self.fileGroupId == None:
            return
            
        #print "selected file group id is: ", self.fileGroupId
        #print "current group id is: ", self.fileGroupId    
        mask = self.fgid == self.fileGroupId
        selFiles = self.files[mask]
        selFileNames = [os.path.splitext(os.path.basename(f))[0] for f in selFiles]
        selFileNames = np.array(selFileNames)
        #print "Files selected are "
        #print selFileNames
        # Check if file names already appear in previous list, otherwise append them and read/process/display relative data

        if firstRun:
            #print "This is the first run "
            # Grab RA-Dec of first file and restart the WCS
            aor, hk, gratpos, flux = readData(selFileNames[0]+".fits")
            obsdate, coords, offset, angle, za, altitude, wv = hk
            self.pc.compute_initial_figure(coords[0],coords[1])
            # Clear the plot
            self.fc.compute_initial_figure(self.fileGroupId)

        
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
                # Give time to the system to update the GUI
                QApplication.processEvents()
            # Create lists of spectra, nod, positions, and file numbers
            fn,nod,ra,dec,x,y,angle,spectra = map(list, zip(*((o.n,o.nod,o.ra,o.dec,o.x,o.y,o.angle,o.spec)  for o in self.obs if o.fgid == self.fileGroupId)))
            # Display the data
            #print "Adding to the plot ",infile
            self.fc.updateFigure(nod,fn,spectra,infile)
            self.pc.updateFigure(nod,fn,ra,dec,x,y,angle,infile)
            # Used memory check
            #pid = os.getpid()
            #py = psutil.Process(pid)
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
        try:
            if not self.files:
                print 'There were no files !'
                self.files=[]
                self.fgid=[]
        except:
            pass

        for f,fg in zip(files,fgid):
            if f not in self.files:
                #print "updating file list ..."
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
