#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-

# Creating the GUI


# Reading the data, processing, displaying them
#import matplotlib
#import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1
from matplotlib.patches import Rectangle
from matplotlib import collections  as mc

# Make sure that we are using QT5
#matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QMessageBox,QToolBar,QAction,QStatusBar,
                             QHBoxLayout, QVBoxLayout, QApplication, QListWidget,QSplitter,QMenu)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QTimer


import os, sys
import numpy as np
#old_settings = np.seterr(all='ignore')  # Avoid warning messages
import warnings
#with warnings.catch_warnings():
#    warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
# To avoid excessive warning messages
warnings.filterwarnings('ignore')

#from multiprocessing import Process

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
        self.offsets = []
        self.filename = []
        self.w.wcs.crval=[ra,dec]
        #        self.axes = self.fig.add_subplot(111, projection=self.w)
        #        #        self.axes.set_xlim([0,20*20]) # Set at least 20 observations

        self.axes = self.fig.add_subplot(111, projection=self.w)
        self.t = self.axes.get_transform(self.w)
        self.axes.axis('equal')
        self.axes.coords[0].set_major_formatter('hh:mm:ss')
        self.axes.plot(np.nan,np.nan)

        
    def updateFigure(self, nod, fn, ra, dec, dx, dy, angle, infile):
        #n = len(nod)
        # Find index of infile
        i = fn.index(int(infile[:5]))
        
        ra=np.asarray(ra)
        dec=np.asarray(dec)
        dx=np.asarray(dx)
        dy=np.asarray(dy)
        x = ra+dx-self.w.wcs.crval[0]
        y = dec+dy-self.w.wcs.crval[1]
        colors = ['blue' if a =='A' else 'red' for a in nod]
        self.offsets.append((dx[i],dy[i]))
        self.filename.append(infile)
        
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
        self.axes.add_patch(rect)
        # Add patch on current position
        try:
            self.greenpatch.remove()
        except:
            pass
        greenrect = Rectangle((x[i] - dx, y[i] - dy), side,side,angle=-angle[i],fc='#C1FFC1',ec='none',alpha=0.5)
        self.greenpatch = self.axes.add_patch(greenrect)

        
        # collect the patches to modify them later
        self.rects.append(rect)
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
            curr_ylim = self.axes2.get_ylim()
            curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
            new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
            self.axes2.set_ylim([curr_y0-new_height,curr_y0+new_height])
        else:
            curr_xlim = self.axes2.get_xlim()
            curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
            new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
            self.axes2.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.draw()

    def onPick(self, event):
        if event.button == 1:
            x = event.xdata
            # compute the position number, return if it is not defined
            try:
                n = int(x // self.coverage)
            except:
                return
            if n >= 0 and n < self.n:
                # draw rectangle on flux plot
                self.gr1.remove()
                self.gr2.remove()
                self.gr1 = self.axes1.axvspan(n*self.coverage, (n+1)*self.coverage, alpha=0.5, color='#E7FFE7')
                self.gr2 = self.axes2.axvspan(n*self.coverage, (n+1)*self.coverage, alpha=0.5, color='#E7FFE7')
                # update position 
                pc = self.parent().parent().parent().parent().pc
                rect = pc.rects[n]
                pc.greenpatch.set_xy(rect.get_xy())
                pc.draw()
                # write offset in the status bar
                offset = pc.offsets[n]
                fname  = pc.filename[n]
                mw = self.parent().parent().parent().parent()
                mw.sb.showMessage("File: "+fname+" --  Offset ("+"{:.0f}".format(offset[0]*3600.)+","+"{:.0f}".format(offset[1]*3600.)+") ",3000)

        elif event.button == 2:    
            self.dragged = event
            self.pick_pos = (event.xdata, event.ydata)
            print "pick position: ", self.pick_pos
        elif event.button == 3:
            # Call popup menu to show housekeeping values (altitude, zenith angle, water vapor)
            pass
        else:
            pass

            
    def contextMenuEvent(self, event):
        drawZA = QAction('Draw zenith angle', self)
        drawZA.triggered.connect(self.drawZA)
        drawAlt = QAction('Draw altitude', self)
        drawAlt.triggered.connect(self.drawAlt)
        drawWV = QAction('Draw water vapor', self)
        drawWV.triggered.connect(self.drawWV)
        menu = QMenu(self)
        menu.addAction(drawZA)
        menu.addAction(drawAlt)
        menu.addAction(drawWV)
        menu.exec_(event.globalPos())


    def drawZA(self):
        print "Draw zenith angle"    
        self.displayZA ^= True
        self.zaLayer.set_visible(self.displayZA)
        self.draw()

    def drawAlt(self):
        self.displayAlt ^= True
        self.altLayer.set_visible(self.displayAlt)
        self.draw()
        print "Draw altitude"    

    def drawWV(self):
        self.displayWV ^= True
        self.wvLayer.set_visible(self.displayWV)
        self.draw()
        print "Draw water vapor"    
        
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
        # Clear figure    
        self.fig.clear()
        # Initialize display
        self.displayZA = True
        self.displayAlt = True
        self.displayWV = True
        # Create box
        self.axes1 = self.fig.add_subplot(211)
        self.axes1.set_xlim([0,20*20]) # Set at least 20 observations
        self.axes1.set_ylim([25,75])
        self.axes1.set_ylabel('Zenith angle [degs]')
        self.axes1b = self.axes1.twinx()
        self.axes1b.set_ylim([36000,45000])
        self.axes1b.get_yaxis().set_tick_params(labelright='on',right='on')            
        self.axes1b.get_yaxis().set_tick_params(which='both', direction='out',colors='green')
        self.axes1b.yaxis.set_label_coords(-0.13,0.5)
        self.axes1b.set_ylabel('Altitude [ft]',color='green')
        self.axes1c = self.axes1.twinx()
        self.axes1c.set_ylim([1,100])
        self.axes1c.tick_params(labelright='on',right='on',direction='in',pad=-20,colors='orange')
        self.axes1c.yaxis.set_label_coords(-0.11,0.5)
        self.axes1c.set_ylabel('Water vapor [$\mu$m]',color='orange')
        self.axes2 = self.fig.add_subplot(212,sharex=self.axes1)
        self.axes2.set_ylabel('Flux [V/s]')
        self.axes2.plot([0,0],[np.nan,np.nan],'.b')
        self.fig.subplots_adjust(hspace=0., bottom=0.1)
        if fileGroupId is not None:
            self.fig.suptitle(fileGroupId)


    def updateFigure(self,nod,fn,spec,infile,za,alti,wv):
        # get number of grating positions
        #print "shape spec is: ",np.shape(spec)
        ng = (np.shape(spec))[1]
        self.coverage = 16+0.5*ng
        sp16 = np.arange(16)
        self.n = len(nod)
        #print fn
        i = fn.index(int(infile[:5]))
        colors = ['blue' if a =='A' else 'red' for a in nod]
        #if nod[i] == 'A':
        #    c = 'blue'
        #else:
        #    c = 'red'

        s = spec[i]
        #print "ng is: ",ng
        for j in np.arange(ng):
            self.axes2.plot(sp16+i*self.coverage+j*0.5, s[j,:], color=colors[i])

        # clear axes 1
        self.axes1.clear()    
        self.axes1b.clear()    
        self.axes1c.clear()
        self.axes1.set_ylabel('Zenith angle [degs]')
        self.axes1b.set_ylabel('Altitude [ft]',color='green')
        self.axes1c.set_ylabel('Water vapor [$\mu$m]',color='orange')
            
        self.axes2.axvline(i*self.coverage, color='gray', linewidth=0.2)            
        self.axes2.axvline(len(fn)*self.coverage, color='gray', linewidth=0.2)
        for ii in np.arange(i+1):
            self.axes1.axvline(ii*self.coverage, color='gray', linewidth=0.2)            
        self.axes1.axvline(len(fn)*self.coverage, color='gray', linewidth=0.2)
        self.axes2.yaxis.grid(True)
        labels = np.array(fn, dtype='str')    
        self.axes2.set_xticks((np.arange(len(fn))+0.5)*self.coverage)
        self.axes2.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation
        self.axes2.set_xlim([self.coverage*(i-1.5*self.w),self.coverage*(i+0.5*self.w)]) # Set at least 20 observations
        self.axes2.autoscale(enable=True,axis='y')

        # Shade background of last position
        try:
            # Try removing previous shaded region
            self.gr2.remove()
        except:
            pass
        self.gr1 = self.axes1.axvspan(i*self.coverage, (i+1)*self.coverage, alpha=0.5, color='#E7FFE7')
        self.gr2 = self.axes2.axvspan(i*self.coverage, (i+1)*self.coverage, alpha=0.5, color='#E7FFE7')
        
        # Display curves
        x1 = np.arange(i+1)*self.coverage
        x2 = np.arange(1,i+2)*self.coverage
        zalines = [[(xs,y[0]),(xe,y[1])] for xs,xe,y in zip(x1,x2,za)]
        altlines = [[(xs,y[0]),(xe,y[1])] for xs,xe,y in zip(x1,x2,alti)]
        wvlines = [[(xs,y[0]),(xe,y[1])] for xs,xe,y in zip(x1,x2,wv)]
        zamin = 32
        zamax = 67
        cza = ['black' if (y[0] < zamax and y[1] < zamax and y[0] > zamin and y[1] > zamin) else 'red' for y in za]
        self.zaLayer = mc.LineCollection(zalines, colors=cza, linewidths=1)
        self.altLayer = mc.LineCollection(altlines, colors='green', linewidths=1)
        self.wvLayer = mc.LineCollection(wvlines, colors='orange', linewidths=1)

        # Redraw
        self.axes1.add_collection(self.zaLayer)
        self.axes1b.add_collection(self.altLayer)
        self.axes1c.add_collection(self.wvLayer)
        
        # Hide/show curves
        self.altLayer.set_visible(self.displayAlt)
        self.zaLayer.set_visible(self.displayZA)
        self.wvLayer.set_visible(self.displayWV)

        self.draw()

    def updateFigureMedians(self,nod,fn,spec,infile):
        # get number of grating positions
        ng = (np.shape(spec))[1]
        self.coverage = 16*ng
        #sp16 = np.arange(16)
        #n = len(nod)
        i = fn.index(int(infile[:5]))
        colors = ['blue' if a =='A' else 'red' for a in nod]
        #if nod[i] == 'A':
        #    c = 'blue'
        #else:
        #    c = 'red'

        s = spec[i]
        x = []
        y = []
        for j in np.arange(ng):
            n0 = int((ng-j)*0.5)
            n1 = int(16+ng*0.5-n0)
            x.append(i*ng+j)
            y.append(np.median(s[j,n0:n1]))

        self.axes2.plot(x,y,color=colors[i] )

        self.axes2.axvline(i*ng, color='gray', linewidth=0.2)            
        self.axes2.axvline(len(fn)*ng, color='gray', linewidth=0.2)
        self.axes2.yaxis.grid(True)
        labels = np.array(fn, dtype='str')    
        self.axes2.set_xticks((np.arange(len(fn))+0.5)*ng)
        self.axes2.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation
        self.axes2.set_xlim([ng*(i-1.5*self.w),ng*(i+0.5*self.w)]) # Set at least 20 observations
        self.axes2.autoscale(enable=True,axis='y')
        self.draw()

class myListWidget(QListWidget):

    def Clicked(self,item):
        mw = self.parent().parent()
        # 2000 means message erased after 2 seconds
        mw.sb.showMessage("You selected the FileGroupID: "+item.text(),2000)
        mw.lf.setVisible(False)
        # Trigger event related to item list ....
        mw.addObs(item.text(),False)


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
        # Get the path of the package
        path0, file0 = os.path.split(__file__)

        # Background color (FFF7C0 is buttermilk, DCAE1D is honey, F2D388 is butter)
        # Colors from https://designschool.canva.com/blog/website-color-schemes/
        self.setStyleSheet("""
        QMainWindow {
        background-color: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #FFF6BA, stop: 1 #F2D388);
        }
        QMenu {
        background-color: #D8AB4E;
        background: '#FFF6BA';
        color: 'black';
        }
        QMenuBar {
        background-color: #F2D388;
        background: #F2D388;
        color: 'black';
        }
        QStatusBar {
        background-color: #FFF6BA;
        }
        QToolBar#tb1 {
        background-color: #FFF6BA;
        }
        QToolBar#tb2 {
        background-color: #FFF6BA;
        }
        QToolBar#tb {
        background-color: #FFF6BA;
        }
        QIcon {
        background: 'transparent';
        }
        """)
        
        # Start exploring directory
        from fifitools import exploreDirectory
        cwd = os.getcwd()
        self.files, self.start, self.fgid = exploreDirectory(cwd+"/")
        # Compile the list of unique File group IDs
        self.fgidList = list(set(self.fgid))
        
        # Start list of observations
        self.fileGroupId = None

        # Variables of the session
        self.fileNames = []
        self.obs = []
        try:
            self.loadData()
        except:
            pass

        # Menu
        # self.file_menu = QMenu('&File', self)
        # self.file_menu.addAction('&Quit', self.fileQuit, Qt.CTRL + Qt.Key_Q)

        # self.help_menu = QMenu('&Help', self)
        # self.help_menu.addAction('&About', self.about)

        # #menubar = self.menuBar()
        # menubar = QMenuBar()
        # self.setMenuBar(menubar)
        # # Request for the MacOSX - do not use native menubar
        # #menubar.setNativeMenuBar(False)
        # menubar.addMenu(self.file_menu)
        # menubar.addSeparator()
        # menubar.addMenu(self.help_menu)

        # Menu
        self.file_menu = self.menuBar().addMenu('&File')
        self.quit_program = QAction('Quit',self,shortcut='Ctrl+q',triggered=self.fileQuit)
        self.file_menu.addAction(self.quit_program)
        
        self.file_menu = self.menuBar().addMenu('&Help')
        self.about_code = QAction('About',self,shortcut='Ctrl+h',triggered=self.about)
        self.file_menu.addAction(self.about_code)

        self.menuBar().setNativeMenuBar(False)


        # Define main widget
        self.main_widget = QWidget(self)

        # Define sub widgets
        self.pc = PositionCanvas(self.main_widget, width=5, height=4, dpi=100)
        self.fc = FluxCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.mpl_toolbar = NavigationToolbar(self.pc, self)
        self.mpl_toolbar.pan('on')
        self.mpl_toolbar.setObjectName('tb1')
        self.mpl_toolbar2 = NavigationToolbar(self.fc, self)
        self.mpl_toolbar2.pan('on')
        self.mpl_toolbar2.setObjectName('tb2')

        # Actions
        exitAction = QAction(QIcon(path0+'/icons/exit.png'), 'Exiting', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(self.closeEvent)
        exitAction.setMenuRole(QAction.NoRole)
        hideAction = QAction(QIcon(path0+'/icons/list.png'), 'List of File Group IDs', self)
        hideAction.setShortcut('Ctrl+L')
        hideAction.triggered.connect(self.changeVisibility)

        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.addAction(hideAction)
        self.tb.addAction(exitAction)
        self.tb.setObjectName('tb')
        
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
        radecLayout.addWidget(self.sb)
        fluxLayout.addWidget(self.fc)

        #fluxLayout.addWidget(self.mpl_toolbar2)
        toolLayout.addWidget(self.mpl_toolbar2)
        toolLayout.addWidget(self.tb)
        #radecLayout.addWidget(self.tb)
        fluxLayout.addWidget(toolWidget)

        splitter1.addWidget(radecWidget)
        splitter1.addWidget(fluxWidget)

        mainLayout.addWidget(splitter1)
        mainLayout.addWidget(self.lf)
        self.lf.hide()

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)


        self.sb.showMessage("Welcome to FIFI Monitor!", 10000)

        # Periodical update
        #self.reduction=False
        timer = QTimer(self)
        timer.timeout.connect(self.update_fifimon)
        timer.start(3000)

    def changeVisibility(self):
        state = self.lf.isVisible()
        self.lf.setVisible(not state)

    def saveData(self):
        import json, io
        data = {}
        for fn,o in zip(self.fileNames,self.obs):
            data[fn] = {
                'n': o.n,
                'ra': o.ra,
                'dec': o.dec,
                'x': o.x,
                'y': o.y,
                'angle': o.angle,
                'alt': o.alt,
                'za': o.za,
                'gp': o.gp.tolist(),
                'wv': o.wv,
                'fgid': o.fgid,
                'nod': o.nod,
                'spec': o.spec.tolist()
            }
        with io.open('fifimon.json', mode='w') as f:
            str_= json.dumps(data,indent=2,sort_keys=True,
                             separators=(',',': '), ensure_ascii=False)
            f.write(unicode(str_))
        pass    

    def loadData(self):
        import json
        from collections import OrderedDict
        from fifitools import Obs
        print 'loading previous reduced data'
        # Read the json file as ordered dictionary to preserve the
        # order the objects were saved
        with open('fifimon.json') as f:
            data = json.load(f, object_pairs_hook=OrderedDict)

        for key,value in data.iteritems():
            self.fileNames.append(key)
            spec=np.array(value['spec'])
            x=value['x']
            y=value['y']
            ra=value['ra']
            dec=value['dec']
            angle=value['angle']
            alt=value['alt']
            za=value['za']
            wv=value['wv']
            nod=value['nod']
            fgid=value['fgid']
            n=value['n']
            gp=np.array(value['gp'])
            self.obs.append(Obs(spec,(ra,dec),(x,y),angle,(alt[0],alt[1]),(za[0],za[1]),(wv[0],wv[1]),nod,fgid,n,gp))
            

        # key = []
        # for k,v in data.iteritems():
        #     key.append(k)
        # key.sort()
            
        # for k in key:
        #     value = data[k]
        #     print k, np.shape(spec)

        #print "Added ", len(key)," objects"
        
    def fileQuit(self):
        self.saveData()    
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        # Get path of the package
        path0,file0 = os.path.split(__file__)
        file=open(path0+"/copyright.txt","r")
        message=file.read()
        QMessageBox.about(self, "About", message)

    def addObs(self, fileGroupId = None, newfiles=False):
        from fifitools import readData, multiSlopes, Obs
        import os
        #import psutil
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
            try:
                aor, hk, gratpos, flux = readData(selFileNames[0]+".fits")
                obsdate, coords, offset, angle, za, altitude, wv = hk
            except:
                print "Failed to read file ",selFileNames[0]
            # Clear the plot
            # Grab RA-Dec of first file and restart the WCS
            self.pc.compute_initial_figure(coords[0],coords[1])
            self.fc.compute_initial_figure(self.fileGroupId)

            
        for infile in selFileNames:
            #print "infile", infile    
            if infile not in self.fileNames:
                print "Reading file: ", infile
                try:
                    self.fileNames.append(infile)
                    aor, hk, gratpos, flux = readData(infile+".fits")
                    print "flux shape is: ",np.shape(flux)
                    spectra = multiSlopes(flux)
                    #print "shape of spectrum ", np.shape(spectra)
                    spectrum = np.nanmedian(spectra,axis=2)
                    detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum = aor
                    obsdate, coords, offset, angle, za, altitude, wv = hk
                    #print "shape of spectrum ", np.shape(spectrum)
                    self.obs.append(Obs(spectrum,coords,offset,angle,altitude,za,wv,nodbeam,filegpid,filenum,gratpos))
                    # Give time to the system to update the GUI
                    QApplication.processEvents()
                    self.update_figures(infile)
                except:
                    print "Problems with file: ", infile
                    #newfiles=True
            else:
                if not newfiles:
                    #print "call update figure"
                    self.update_figures(infile)
            # Used memory check
            #pid = os.getpid()
            #py = psutil.Process(pid)
            #memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
            #print('memory use:', memoryUse)
        # When exiting reset reduction
        #self.reduction = False

    def update_figures(self,infile):
        # Create lists of spectra, nod, positions, and file numbers
        fn,nod,ra,dec,x,y,angle,spectra,za,alti,wv = map(list, zip(*((o.n,o.nod,o.ra,o.dec,o.x,o.y,o.angle,o.spec,o.za,o.alt,o.wv)
                                                                     for o in self.obs if o.fgid == self.fileGroupId)))
        # Display the data
        print "Adding to the plot ",infile
        print np.shape(spectra)
        self.fc.updateFigure(nod,fn,spectra,infile,za,alti,wv)
        self.pc.updateFigure(nod,fn,ra,dec,x,y,angle,infile)
                
    
    def update_fifimon(self):
        '''
        This module is called automatically every n seconds (n=10, typically) to look
        for new files in the directory, update the lists and, eventually, the plot
        '''

        # Check if running reductions
#        if self.reduction == False:
#            self.reduction == True
#        else:
#            return

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
                print "updating file list ..."
                self.files = np.append(self.files,f)
                self.fgid = np.append(self.fgid, fg)
                # update plot (if same fileGroupID)
                # Check this part. It works if reading from fifimon.json and adding new files in the directory.
                # But does not work if Timer finds new files in the directory ... and restart replotting all the list of files slowing everything down
                self.addObs(None,True)
            if fg not in self.fgidList:
                # update fileGroupID list    
                self.fgidList.append(fg)
                # update widget list
                self.lf.addItem(fg)
                # pop-up window alerting and asking if processing the new data
                choice = QMessageBox.question(self,
                                              'New set of data appeared !','Do you want to show them ?',
                                              QMessageBox.Yes | QMessageBox.No, QMessageBox.No )
                if choice == QMessageBox.Yes:
                    # Write on status
                    self.sb.showMessage("You selected the FileGroupID: "+fg,2000)
                    # Eventually hide list
                    self.lf.setVisible(False)
                    # Start new plot
                    self.addObs(fg,False)
                else:
                    pass    


                
""" Main code """
        
def main():
    app = QApplication(sys.argv)
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    aw = ApplicationWindow()
    aw.setGeometry(100, 100, width*0.9, width*0.35)
    progname = 'FIFI Monitor'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()

# Ensure that the app is created once 
if __name__ == '__main__':
#    p = Process(target=main)
#    p.start()
#    p.join()
    main()
