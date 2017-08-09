#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-
# encoding=utf8
import sys
#reload(sys)
#sys.setdefaultencoding("utf-8")
# Creating the GUI


# Reading the data, processing, displaying them
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection

# Make sure that we are using QT5
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QMessageBox,QToolBar,QAction,QStatusBar,
                             QHBoxLayout, QVBoxLayout, QApplication, QListWidget,QSplitter,QMenu)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QTimer, QThread, pyqtSignal, QObject


import os
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
    """Simple canvas for positions in the sky"""
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
        self.draw_idle()
        
    def compute_initial_figure(self,ra=None,dec=None):
        from astropy.wcs import WCS
        self.w = WCS(naxis=2)
        self.w.wcs.ctype=["RA---TAN","DEC--TAN"]
        self.w.wcs.crpix=[1,1]
        self.w.wcs.cdelt = np.array([1,1])
        self.w.wcs.pc[0]=[1,0]
        self.w.wcs.pc[1]=[-0,1]
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

        self.axes = self.fig.add_subplot(111, projection=self.w)
        self.t = self.axes.get_transform(self.w)
        self.axes.axis('equal')
        self.axes.coords[0].set_major_formatter('hh:mm:ss')
        self.axes.plot(np.nan,np.nan)
        #plt.ion()
        #plt.show(False)

        
    def updateFigure(self, nod, ra, dec, dx, dy, angle, infile):

        x = ra+dx-self.w.wcs.crval[0]
        y = dec+dy-self.w.wcs.crval[1]
        color = 'blue' if nod == 'A' else 'red'
        self.offsets.append((dx,dy))
        self.filename.append(infile)
        
        # Add rectangle patch
        side = 30./3600.
        theta = -angle*np.pi/180.
        dx = side*0.5*(np.cos(theta)-np.sin(theta))
        dy = side*0.5*(np.sin(theta)+np.cos(theta))
        rect = Rectangle((x - dx, y - dy), side,side,angle=-angle,fc='none',ec=color)
        self.axes.add_patch(rect)
        # Add patch on current position
        try:
            self.greenpatch.remove()
        except:
            pass
        greenrect = Rectangle((x - dx, y - dy), side,side,angle=-angle,fc='#C1FFC1',ec='none',alpha=0.6)
        self.greenpatch = self.axes.add_patch(greenrect)
        
        # collect the patches to modify them later
        self.rects.append(rect)
        self.axes.autoscale(enable=True, axis='both')
        # Invert x axis
        xlim = self.axes.get_xlim()
        if xlim[0] < xlim[1]:
            self.axes.set_xlim([xlim[1],xlim[0]])
        # Update figure
        self.draw_idle()

class FluxCanvas(MplCanvas):
    """A canvas that show spectra and some housekeeping"""

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
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
        # ShiftModifier does not work with MAC OSX ! I have to use control (the Apple icon)
        if modifiers == Qt.ControlModifier:
            curr_ylim = self.axes2.get_ylim()
            curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
            new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
            self.axes2.set_ylim([curr_y0-new_height,curr_y0+new_height])
        else:
            curr_xlim = self.axes2.get_xlim()
            curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
            new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
            self.axes2.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.draw_idle()

    def onPick(self, event):
        if event.button == 1:
            x = event.xdata
            # Search the label position within the mininum distance
            try:
                pl = np.array(self.labpos)
                n = np.argmin(np.abs(pl-x))
                # draw rectangle on flux plot
                self.gr1.remove()
                self.gr2.remove()
                self.gr1 = self.axes1.axvspan(pl[n]-0.5*self.coverage[n], pl[n]+0.5*self.coverage[n], alpha=0.6, color='#E7FFE7')
                self.gr2 = self.axes2.axvspan(pl[n]-0.5*self.coverage[n], pl[n]+0.5*self.coverage[n], alpha=0.6, color='#E7FFE7')
                # update position 
                pc = self.parent().parent().parent().parent().pc
                rect = pc.rects[n]
                pc.greenpatch.set_xy(rect.get_xy())
                pc.draw_idle()
                # write offset in the status bar
                offset = pc.offsets[n]
                fname  = pc.filename[n]
                mw = self.parent().parent().parent().parent()
                mw.sb.showMessage("File: "+fname+" --  Offset ("+"{:.0f}".format(offset[0]*3600.)+","+"{:.0f}".format(offset[1]*3600.)+") ",3000)
            except:
                print "No data displayed"

                
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
        self.draw_idle()

    def drawAlt(self):
        self.displayAlt ^= True
        self.altLayer.set_visible(self.displayAlt)
        self.draw_idle()
        print "Draw altitude"    

    def drawWV(self):
        self.displayWV ^= True
        self.wvLayer.set_visible(self.displayWV)
        self.draw_idle()
        print "Draw water vapor"    
        
    def onMotion(self, event):
        if self.dragged is not None and self.pick_pos[0] is not None:
            #old_pos = self.dragged.get_position()
            new_pos = (event.xdata,event.ydata)
            deltax = new_pos[0] - self.pick_pos[0]
            deltay = new_pos[1] - self.pick_pos[1]
            curr_xlim = self.axes1.get_xlim()
            curr_ylim = self.axes1.get_ylim()
            self.axes1.set_xlim(curr_xlim-deltax)
            self.axes1.set_ylim(curr_ylim-deltay)
            self.pick_pos = new_pos
            # Draw only when idle
            self.draw_idle()
        return True

    def onRelease(self, event):
        if self.dragged is not None:
            self.dragged = None
            self.draw_idle()
        return True



    def compute_initial_figure(self,fileGroupId=None):
        # Clear figure    
        self.fig.clear()
        # Initialize display
        self.displayZA = True
        self.displayAlt = True
        self.displayWV = True

        # Initialize variables
        self.zalines = []
        self.altlines =[]
        self.wvlines = []
        self.cza = []
        self.labels = []
        self.labpos = []
        self.coverage = []
        self.wvl = []
        self.wdict = {'0':0}
        self.zamin = 32
        self.zamax = 67
        # Create box
        self.axes1 = self.fig.add_subplot(211)
        self.axes1.set_xlim([0,20*20]) # Set at least 20 observations
        self.axes1.set_ylim([25,75])
        self.axes1.set_ylabel('Zenith angle [degs]')
        self.axes1.get_yaxis().set_tick_params(which='both', direction='in',colors='black', right='on',pad=-20)
        self.axes1.yaxis.set_label_coords(-0.07,0.5)
        self.axes1b = self.axes1.twinx()
        self.axes1b.set_ylim([36000,45000])
        self.axes1b.get_yaxis().set_tick_params(labelright='on',right='on')            
        self.axes1b.get_yaxis().set_tick_params(which='both', direction='out',colors='green')
        self.axes1b.yaxis.set_label_coords(-0.15,0.5)
        self.axes1b.set_ylabel('Altitude [ft]',color='green')
        self.axes1c = self.axes1.twinx()
        self.axes1c.set_ylim([1,100])
        self.axes1c.tick_params(labelright='on',right='on',direction='in',pad=-20,colors='orange')
        self.axes1c.yaxis.set_label_coords(-0.13,0.5)
        self.axes1c.set_ylabel('Water vapor [$\mu$m]',color='orange')
        self.axes1d = self.axes1.twinx()
        self.axes1d.get_yaxis().set_tick_params(which='both', direction='out',colors='blue')
        self.axes1d.get_yaxis().set_tick_params(labelleft='on', left='on',labelright='off',right='off')
        self.axes1d.yaxis.set_label_coords(-0.11,0.5)
        self.axes1d.set_ylabel('Wavelength [$\mu$m]',color='blue')    
        self.axes2 = self.fig.add_subplot(212,sharex=self.axes1)
        self.axes2.set_ylabel('Flux [V/s]')
        self.axes2.plot([0,0],[np.nan,np.nan],'.b')
        self.axes2.yaxis.grid(True)
        self.fig.subplots_adjust(hspace=0., bottom=0.1)
        if fileGroupId is not None:
            mw = self.parent().parent().parent().parent()
            self.fig.suptitle(fileGroupId+" ("+mw.channel+")")


    def updateFigure(self,nod,spec,infile,za,alti,wv,gp,order,obsdate,dichroic):
        from fifitools import waveCal        
        # get number of grating positions
        ng = (np.shape(spec))[0]
        start = sum(self.coverage)
        self.coverage.append(16+0.5*(ng-1))
        sp16 = np.arange(16)
        self.labels.append(int(infile[:5]))
        i = len(self.labels)
        color = 'blue' if nod == 'A' else 'red'

        lines  = []
        for j in np.arange(ng):
            x = sp16+start+j*0.5
            y = spec[j,:]
            lines.append(zip(x,y))
        lc = LineCollection(lines, colors=color, linewidths=1)
        self.lines = self.axes2.add_collection(lc)
            
        self.axes1.axvline(start, color='gray', linewidth=0.2)            
        self.axes1.axvline(start+self.coverage[i-1], color='gray', linewidth=0.2)
        self.axes2.axvline(start, color='gray', linewidth=0.2)            
        self.axes2.axvline(start+self.coverage[i-1], color='gray', linewidth=0.2)
        # Update labels
        labels = np.array(self.labels, dtype='str')
        self.labpos.append(sum(self.coverage)-0.5*self.coverage[i-1])
        self.axes2.set_xticks(self.labpos)
        self.axes2.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation (at least 20 observations)
        self.axes2.set_xlim([self.coverage[i-1]*(i-1.5*self.w),self.coverage[i-1]*(i+0.5*self.w)])
        self.axes2.autoscale(enable=True,axis='y')

        # Shade background of last position
        try:
            # Try removing previous shaded region
            self.gr1.remove()
            self.gr2.remove()
        except:
            pass
        self.gr1 = self.axes1.axvspan(start, start+self.coverage[i-1], alpha=0.6, color='#E7FFE7')
        self.gr2 = self.axes2.axvspan(start, start+self.coverage[i-1], alpha=0.6, color='#E7FFE7')
        
        # Display curves
        xs = start
        xe = start+self.coverage[i-1]
        self.zalines.append([(xs,za[0]),(xe,za[1])])
        self.altlines.append([(xs,alti[0]),(xe,alti[1])])
        self.wvlines.append([(xs,wv[0]),(xe,wv[1])])
        if (za[0] < self.zamax and za[1] < self.zamax and za[0] > self.zamin and za[1] > self.zamin):
            self.cza.append('black')
        else:
            self.cza.append('red')
            
        try:
            self.zaLayer.remove()
            self.altLayer.remove()
            self.wvLayer.remove()
        except:
            pass

        self.zaLayer = LineCollection(self.zalines, colors=self.cza, linewidths=1)
        self.altLayer = LineCollection(self.altlines, colors='green', linewidths=1)
        self.wvLayer = LineCollection(self.wvlines, colors='orange', linewidths=1)

        # Redraw
        self.axes1.add_collection(self.zaLayer)
        self.axes1b.add_collection(self.altLayer)
        self.axes1c.add_collection(self.wvLayer)
        
        # Hide/show curves
        self.altLayer.set_visible(self.displayAlt)
        self.zaLayer.set_visible(self.displayZA)
        self.wvLayer.set_visible(self.displayWV)

        # Add grating position in blue
        #self.axes1d.plot(xs+0.5*np.arange(ng), gp*.001,'.',color=color)
        # Compute wavelenghts and plot them
        mw = self.parent().parent().parent().parent()
        j=0
        for g in gp:
            if g in self.wdict:
                waves = self.wdict[g]
            else:
                l,lw = waveCal(gratpos=g, order=order, array=mw.channel,dichroic=dichroic,obsdate=obsdate)
                waves = l[12,:]
                self.wdict[g]=waves
            x = sp16+start+j*0.5
            self.axes1d.plot(x,waves,'.',color=color)
            self.wvl.append(waves)  # Conserve the central wavelengths
            j += 1
        self.axes1d.autoscale(enable=True,axis='y')        
        self.draw_idle()
        ww = np.array(self.wvl)
        mw.sb.showMessage("Wavelength: "+"{:.1f}".format(np.min(ww))+" - "+"{:.1f}".format(np.median(ww))+" - "+"{:.1f}".format(np.max(ww)),5000)

class myListWidget(QListWidget):

    def Clicked(self,item):
        mw = self.parent().parent()
        # 2000 means message erased after 2 seconds
        mw.sb.showMessage("You selected the FileGroupID: "+item.text(),2000)
        mw.lf.setVisible(False)
        # Trigger event related to item list ....
        mw.addObs(item.text(),False)

    
class UpdateObjects(QObject):
    from fifitools import Obs
    newObj = pyqtSignal([Obs])


class AddObsThread(QThread):

    updateObjects = UpdateObjects()
    updateFigures = pyqtSignal('QString')
    updateFilenames = pyqtSignal('QString')
    updateStatus = pyqtSignal('QString')

    def __init__(self, selFileNames, fileNames, processAll, parent=None):
        super(AddObsThread, self).__init__(parent)
        self.selFileNames = selFileNames
        self.fileNames = fileNames
        self.processAll = processAll
        
    def run(self):
        from fifitools import readData, multiSlopes, Obs
        from timeit import default_timer as timer
        for infile in self.selFileNames:
            if infile not in self.fileNames:
                try:
                    t1=timer()
                    aor, hk, gratpos, flux = readData(infile+".fits")
                    spectra = multiSlopes(flux)
                    spectrum = np.nanmedian(spectra,axis=2)
                    detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum = aor
                    obsdate, coords, offset, angle, za, altitude, wv = hk
                    obj = Obs(spectrum,coords,offset,angle,altitude,za,wv,nodbeam,filegpid,filenum,gratpos,detchan,order,obsdate,dichroic)
                    t2=timer()
                    print "Fitted: ", infile, " ",np.shape(spectrum), " in: ", t2-t1," s"
                    # Call this with a signal from thread
                    self.updateObjects.newObj.emit(obj)
                    self.updateFilenames.emit(infile)
                    self.updateFigures.emit(infile)
                except:
                    print "Problems with file: ", infile
            else:
                self.updateFigures.emit(infile)
        print "Done adding observations thread"
        if self.processAll:
            self.updateStatus.emit('next')
        # Disconnect signal at the end of the thread
        self.updateObjects.newObj.disconnect()


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
        # To work on MAC-OSX for QToolBar, one has to set the border to force the style on the system
        # https://bugreports.qt.io/browse/QTBUG-12717
        # Qt documentation: http://doc.qt.io/qt-5/stylesheet-examples.html
        self.setStyleSheet("""
        QMainWindow {
        background-color: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 LemonChiffon, stop: 1 #F2D388);
        }
        QMenu {
        background-color: #D8AB4E;
        background: '#FFF6BA';
        color: 'black';
        }
        QMenuBar {
        background-color: QLinearGradient(x1:0, y1:0, x2:0, y2:1, stop:0 #FFF6BA, stop:1 #F2D388);
        background: #F2D388;
        color: 'black';
        }
        QMenuBar::item {
        background: transparent;
        spacing: 3px;
        padding: 1px 4px;
        border-radius: 4px;
        }
        QMenuBar::item:selected { /* when selected using mouse or keyboard */
        background: #FFF6BA;
        }
        QMenuBar::item:pressed {
        background: #DCAE1D;
        }
        QStatusBar {
        background-color: #FFF6BA;
        border: 1px solid black;
        border-radius: 3px;
        }
        QToolBar#tb1, QToolBar#tb2, QToolBar#tb {
        background-color: transparent;
        border: 1px transparent;
        }
        QToolBar::separator{
        background-color: transparent;
        }
        QToolButton:pressed {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:hover {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:focused {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:checked {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolTip {
        border: 1px solid black;
        padding: 2px;
        border-radius: 3px;
        opacity: 200;
        background-color: LemonChiffon;
        }
        """)

        # Start exploring directory
        from fifitools import exploreDirectory
        cwd = os.getcwd()
        self.files, self.start, self.fgid, self.ch = exploreDirectory(cwd+"/")
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
        self.mpl_toolbar = NavigationToolbar(self.pc, self)
        self.mpl_toolbar.pan('on')
        self.mpl_toolbar.setObjectName('tb1')

        self.fc = FluxCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.mpl_toolbar2 = NavigationToolbar(self.fc, self)
        self.mpl_toolbar2.pan('on')
        self.mpl_toolbar2.setObjectName('tb2')

        # Actions
        exitAction = QAction(QIcon(path0+'/icons/exit.png'), 'Exit the program', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(self.closeEvent)
        exitAction.setMenuRole(QAction.NoRole)
        hideAction = QAction(QIcon(path0+'/icons/list.png'), 'List of File Group IDs', self)
        hideAction.setShortcut('Ctrl+L')
        hideAction.triggered.connect(self.changeVisibility)
        procAction = QAction(QIcon(path0+'/icons/gears.png'), 'Process all data', self)
        procAction.setShortcut('Ctrl+P')
        procAction.triggered.connect(self.processAll)

        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.addAction(hideAction)
        self.tb.addAction(procAction)
        self.tb.addAction(exitAction)
        self.tb.setObjectName('tb')
        
        # Widget list
        self.lf = myListWidget()
        for item in self.fgidList:
            self.lf.addItem(item); 
        self.lf.setWindowTitle('File Group ID')
        self.lf.itemClicked.connect(self.lf.Clicked)
        # Status Bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to FIFI Monitor!", 10000)
        
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

        toolLayout.addWidget(self.mpl_toolbar2)
        toolLayout.addWidget(self.tb)
        fluxLayout.addWidget(toolWidget)

        splitter1.addWidget(radecWidget)
        splitter1.addWidget(fluxWidget)

        mainLayout.addWidget(splitter1)
        mainLayout.addWidget(self.lf)
        self.lf.hide()

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        # Thread status
        self.threadStatus = 'done'    
        
        # Periodical update
        self.numberOfFiles = 0
        timer = QTimer(self)
        timer.timeout.connect(self.update_fifimon)
        timer.start(5000)

    def changeVisibility(self):
        state = self.lf.isVisible()
        self.lf.setVisible(not state)

    def processAll(self):
        ''' Action to process everything in the directory '''

        self.processItem = 0
        item = self.fgidList[self.processItem]
        self.sb.showMessage("Processing FileGroupID: "+item.tostring(),5000)
        self.addObs(item, False, True)

        
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
                'ch': o.ch,
                'order': o.order,
                'obsdate': o.obsdate,
                'dichroic': o.dichroic,
                'spec': o.spec.tolist()
            }
        with io.open('fifimon.json', mode='w') as f:
            str_= json.dumps(data,indent=2,sort_keys=True,
                             separators=(',',': '), ensure_ascii=False,encoding='utf-8')
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
            data = json.load(f, object_pairs_hook=OrderedDict, encoding='utf-8')

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
            ch=value['ch']
            # Damn observation date has a character in latin-1 !
            obsdate=value['obsdate'].decode('latin-1').encode('utf-8')
            order=value['order']
            dichroic=value['dichroic']
            print obsdate
            self.obs.append(Obs(spec,(ra,dec),(x,y),angle,(alt[0],alt[1]),(za[0],za[1]),(wv[0],wv[1]),nod,fgid,n,gp,ch,order,obsdate,dichroic))
        print "Loading completed "

        
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

    def addObs(self, fileGroupId = None, newfiles=False, processAll=False):
        from astropy.io import fits
        import os
        firstRun = 0
        if fileGroupId != None:
            self.fileGroupId = fileGroupId
            firstRun = 1
        if self.fileGroupId == None:
            return

        mask = self.fgid == self.fileGroupId
        # Select the channel
        channels = list(set(self.ch[mask]))
        self.channel = channels[0]
        # Check if there are files
        print "there are ",np.sum(mask)," observations"
        if np.sum(mask) == 0:
            print "There no files in this channel"
            return
        
        
        selFiles = self.files[mask]
        selFileNames = [os.path.splitext(os.path.basename(f))[0] for f in selFiles]
        selFileNames = np.array(selFileNames)
            
        if firstRun:
            try:
                # Grab RA-Dec of first file and restart the WCS
                hdulist = fits.open(selFileNames[0]+".fits")
                header = hdulist[0].header
                ra   = float(header['OBSLAM'])
                dec  = float(header['OBSBET'])
                hdulist.close()
                self.pc.compute_initial_figure(ra,dec)
                self.fc.compute_initial_figure(self.fileGroupId)
            except:
                print "Failed to read file ",selFileNames[0]
                return
        
            
        self.addObsThread = AddObsThread(selFileNames,self.fileNames,processAll)
        self.addObsThread.updateObjects.newObj.connect(self.update_objects)
        self.addObsThread.updateFigures.connect(self.update_figures)
        self.addObsThread.updateFilenames.connect(self.update_filenames)
        self.addObsThread.updateStatus.connect(self.update_status)
        self.addObsThread.start()
        
        # Used memory check
        #pid = os.getpid()
        #py = psutil.Process(pid)
        #memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
        #print('memory use:', memoryUse)

    def update_objects(self, obj):
        self.obs.append(obj)

    def update_status(self, status):
        ''' process next item is available '''
        self.processItem += 1
        if self.processItem < len(self.fgidList):
            item = self.fgidList[self.processItem]
            self.sb.showMessage("Processing FileGroupID: "+item.tostring(),5000)
            self.addObs(item, False, True)
            
    def update_filenames(self, infile):
        #print "updating filenames with file ", infile
        self.fileNames.append(infile)
      
    def update_figures(self, infile):
        from timeit import default_timer as timer

        # Select obs corresponding to infile
        n = self.fileNames.index(infile)
        o = self.obs[n]
        #print "updating figure with filename ", n
        nod,ra,dec,x,y,angle,spectra,za,alti,wv,gp,order,obsdate,dichroic = o.nod,o.ra,o.dec,o.x,o.y,o.angle,o.spec,o.za,o.alt,o.wv,o.gp,o.order,o.obsdate,o.dichroic
        t1=timer()
        self.fc.updateFigure(nod,spectra,infile,za,alti,wv,gp,order,obsdate,dichroic)
        t2=timer()
        self.pc.updateFigure(nod,ra,dec,x,y,angle,infile)
        t3=timer()
        print "Plotted ",infile," in ",t2-t1," and ",t3-t2," s"


        
    def update_fifimon(self):
        '''
        This module is called automatically every n seconds (n=10, typically) to look
        for new files in the directory, update the lists and, eventually, the plot
        '''

        from fifitools import exploreDirectory

        # Check if there are new files
        cwd = os.getcwd()
        files = os.listdir(cwd)
        if len(files) == self.numberOfFiles:
            return
        else:
            self.numberOfFiles = len(files)

        
        files, start, fgid, ch = exploreDirectory(cwd+"/")
        # Check if new fileGroupID or files added
        try:
            if not self.files:
                print 'There were no files !'
                self.files=[]
                self.fgid=[]
                self.ch=[]
        except:
            pass

            
        for f,fg in zip(files,fgid):
            if f not in self.files:
                print "updating file list ..."
                self.files = np.append(self.files,f)
                self.fgid = np.append(self.fgid, fg)
                self.ch = np.append(self.ch, ch)
                # update plot (if same fileGroupID)
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
    main()
