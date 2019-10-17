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
matplotlib.use('QT5Agg')
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

# import for multiprocessing
from fifimon.fifitools import readData, multiSlopes, Obs, applyFlats, waveCal
#from fifitools import readData, multiSlopes, Obs
from timeit import default_timer as timer


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
        self.fmode = False
    

    def onWheel(self, event):
        eb = event.button
        curr_xlim = self.axis.get_xlim()
        curr_ylim = self.axis.get_ylim()
        curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
        curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        self.axis.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.axis.set_ylim([curr_y0-new_height,curr_y0+new_height])
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

        self.axis = self.fig.add_subplot(111, projection=self.w)
        self.t = self.axis.get_transform(self.w)
        self.axis.axis('equal')
        self.axis.coords[0].set_major_formatter('hh:mm:ss')
        self.axis.plot(np.nan,np.nan)
        #plt.ion()
        #plt.show(False)

    def flightmode(self, flightmode):
        if flightmode:
            self.fig.suptitle('Flight mode')
            self.fmode = True
        else:
            self.fig.suptitle('Analysis mode')
            self.fmode = False
        self.draw_idle()
                
        
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
        self.axis.add_patch(rect)
        # Add patch on current position
        try:
            self.greenpatch.remove()
        except:
            pass
        greenrect = Rectangle((x - dx, y - dy), side,side,angle=-angle,fc='#C1FFC1',ec='none',alpha=0.6)
        self.greenpatch = self.axis.add_patch(greenrect)
        
        # collect the patches to modify them later
        self.rects.append(rect)
        self.axis.autoscale(enable=True, axis='both')
        # Invert x axis
        xlim = self.axis.get_xlim()
        if xlim[0] < xlim[1]:
            self.axis.set_xlim([xlim[1],xlim[0]])
        # Update figure
        if self.fmode:
            self.fig.suptitle('Fligth mode')
        else:
            self.fig.suptitle('Analysis mode')
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
            curr_ylim = self.axis2.get_ylim()
            curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
            new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
            self.axis2.set_ylim([curr_y0-new_height,curr_y0+new_height])
        else:
            curr_xlim = self.axis2.get_xlim()
            curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
            new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
            self.axis2.set_xlim([curr_x0-new_width,curr_x0+new_width])
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
                self.gr1 = self.axis1.axvspan(pl[n]-0.5*self.coverage[n], pl[n]+0.5*self.coverage[n], alpha=0.6, color='#E7FFE7')
                self.gr2 = self.axis2.axvspan(pl[n]-0.5*self.coverage[n], pl[n]+0.5*self.coverage[n], alpha=0.6, color='#E7FFE7')
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
                print ("No data displayed")

                
        elif event.button == 2:    
            self.dragged = event
            self.pick_pos = (event.xdata, event.ydata)
            print ("pick position: ", self.pick_pos)
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
        #print "Draw zenith angle"    
        self.displayZA ^= True
        self.zaLayer.set_visible(self.displayZA)
        self.draw_idle()

    def drawAlt(self):
        self.displayAlt ^= True
        self.altLayer.set_visible(self.displayAlt)
        self.draw_idle()
        #print "Draw altitude"    

    def drawWV(self):
        self.displayWV ^= True
        self.wvLayer.set_visible(self.displayWV)
        self.draw_idle()
        #print "Draw water vapor"    
        
    def onMotion(self, event):
        if self.dragged is not None and self.pick_pos[0] is not None:
            #old_pos = self.dragged.get_position()
            new_pos = (event.xdata,event.ydata)
            deltax = new_pos[0] - self.pick_pos[0]
            deltay = new_pos[1] - self.pick_pos[1]
            curr_xlim = self.axis1.get_xlim()
            curr_ylim = self.axis1.get_ylim()
            self.axis1.set_xlim(curr_xlim-deltax)
            self.axis1.set_ylim(curr_ylim-deltay)
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
        self.dwdict = {'0':0}
        self.zamin = 32
        self.zamax = 67
        # Create box
        
        self.grid1 = self.fig.add_gridspec(nrows=2, ncols=1, top=.95, bottom=0.4, hspace=0.0)
        self.axis1 = self.fig.add_subplot(self.grid1[0,0])
        self.axis2 = self.fig.add_subplot(self.grid1[1,0], sharex=self.axis1)
        #self.axis1 = self.fig.add_subplot(211)
        self.axis1.set_xlim([0,20*20]) # Set at least 20 observations
        self.axis1.set_ylim([25,75])
        self.axis1.set_ylabel('Zenith angle [degs]')
        self.axis1.get_yaxis().set_tick_params(which='both', direction='in',colors='black', right='on',pad=-20)
        self.axis1.yaxis.set_label_coords(-0.07,0.5)
        self.axis1b = self.axis1.twinx()
        self.axis1b.set_ylim([36000,45500])
        self.axis1b.get_yaxis().set_tick_params(labelright='on',right='on')            
        self.axis1b.get_yaxis().set_tick_params(which='both', direction='out',colors='green')
        self.axis1b.yaxis.set_label_coords(-0.16,0.5)
        self.axis1b.set_ylabel('Altitude [ft]',color='green')
        self.axis1c = self.axis1.twinx()
        self.axis1c.set_ylim([1,100])
        self.axis1c.tick_params(labelright='on',right='on',direction='in',pad=-20,colors='orange')
        self.axis1c.yaxis.set_label_coords(-0.14,0.5)
        self.axis1c.set_ylabel('Water vapor [$\mu$m]',color='orange')
        self.axis1d = self.axis1.twinx()
        self.axis1d.yaxis.tick_left()
        self.axis1d.get_yaxis().set_tick_params(which='both', direction='out',colors='blue')
        self.axis1d.yaxis.set_label_coords(-0.12,0.5)
        self.axis1d.set_ylabel('Wavelength [$\mu$m]',color='blue')    
        self.axis1d.set_ylim([50,200])
        #self.axis2 = self.fig.add_subplot(212,sharex=self.axis1)
        self.axis2.set_ylabel('Flux [V/s/Hz]')
        self.axis2.plot([0,0],[np.nan,np.nan],'.b')
        self.axis2.yaxis.grid(True)
        #self.fig.subplots_adjust(hspace=0., bottom=0.1)
        self.grid2 = self.fig.add_gridspec(nrows=1, ncols=2, top=0.35, bottom=0.05)
        self.axis3 = self.fig.add_subplot(self.grid2[0,0])
        self.axis4 = self.fig.add_subplot(self.grid2[0,1])
        self.axis3.set_xlabel('Wavelength [$\mu$m]')    
        self.axis4.set_xlabel('Wavelength [$\mu$m]')    
        self.axis3.set_ylabel('Sky flux [V/s/Hz]')
        self.axis4.set_ylabel('Chop diff flux [V/s/Hz]')
        
        if fileGroupId is not None:
            mw = self.parent().parent().parent().parent()
            self.fig.suptitle(fileGroupId+" ("+mw.channel+")")


    def updateFigure(self,nod,wave,spec,sky,infile,za,alti,wv,gp,order,obsdate,dichroic,ra,dec):
        #from fifimon.fifitools import waveCal        
        # get number of grating positions

        ng = (np.shape(spec))[0]
        start = sum(self.coverage)
        self.coverage.append(16+0.5*(ng-1))
        sp16 = np.arange(16)
        self.labels.append(int(infile[:5]))
        i = len(self.labels)
        color = 'blue' if nod == 'A' else 'red'

        lines  = []
        ly = []
        for j in np.arange(ng):
            x = sp16+start+j*0.5
            y = spec[j,:]
            ly.append(y)
            lines.append(list(zip(x,y)))
        lc = LineCollection(lines, colors=color, linewidths=1)
        self.lines = self.axis2.add_collection(lc)

        # Median spectrum for all the grating positions
        ly = np.array(ly)
        medline = np.nanmedian(ly,axis=0)
        self.axis2.plot(sp16+start, medline, color='green')
            
        self.axis1.axvline(start, color='gray', linewidth=0.2)            
        self.axis1.axvline(start+self.coverage[i-1], color='gray', linewidth=0.2)
        self.axis2.axvline(start, color='gray', linewidth=0.2)            
        self.axis2.axvline(start+self.coverage[i-1], color='gray', linewidth=0.2)
        # Update labels
        labels = np.array(self.labels, dtype='str')
        self.labpos.append(sum(self.coverage)-0.5*self.coverage[i-1])
        self.axis2.set_xticks(self.labpos)
        self.axis2.set_xticklabels(labels,rotation=90,ha='center',fontsize=10)
        # Set limits around last observation (at least 20 observations)
        self.axis2.set_xlim([self.coverage[i-1]*(i-1.5*self.w),self.coverage[i-1]*(i+0.5*self.w)])
        self.axis2.autoscale(enable=True,axis='y')

        # Shade background of last position
        try:
            # Try removing previous shaded region
            self.gr1.remove()
            self.gr2.remove()
        except:
            pass
        self.gr1 = self.axis1.axvspan(start, start+self.coverage[i-1], alpha=0.6, color='#E7FFE7')
        self.gr2 = self.axis2.axvspan(start, start+self.coverage[i-1], alpha=0.6, color='#E7FFE7')
        
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
        self.axis1.add_collection(self.zaLayer)
        self.axis1b.add_collection(self.altLayer)
        self.axis1c.add_collection(self.wvLayer)
        
        # Hide/show curves
        self.altLayer.set_visible(self.displayAlt)
        self.zaLayer.set_visible(self.displayZA)
        self.wvLayer.set_visible(self.displayWV)

        # plot wavelengths, sky and object fluxes
        mw = self.parent().parent().parent().parent()
        for l in self.axis3.lines:
            l.set_alpha(.2)
            l.set_linestyle(':')
        for l in self.axis4.lines:
            l.set_alpha(.2)
            l.set_linestyle(':')
        # To hide use linestyle('None')
        j=0
        for g, w, sk, sp in zip(gp, wave, sky, spec):
            x = sp16+start+j*0.5
            self.axis1d.plot(x, w,'.',color=color)
            self.axis3.plot(w, sk, color=color)
            self.axis4.plot(w, sp, color=color)
            self.wvl.append(w)  # Conserve the central wavelengths
            j += 1
            
        
        self.axis1d.autoscale(enable=True,axis='y')        
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
    #from fifimon.fifitools import Obs
    newObj = pyqtSignal([Obs])


class AddObsThread(QThread):

    updateObjects = UpdateObjects()
    updateFigures = pyqtSignal('QString')
    updateExclude = pyqtSignal('QString')
    updateFilenames = pyqtSignal('QString')
    updateStatus = pyqtSignal('QString')

    def __init__(self, selFileNames, fileNames, processAll, firstRun, parent=None):
        #super(AddObsThread, self).__init__(parent)
        super().__init__()
        self.selFileNames = selFileNames
        self.fileNames = fileNames
        self.processAll = processAll
        self.firstRun = firstRun
        
    def run(self):

        print ('files are: ',self.selFileNames)
        print ('previous files are: ', self.fileNames)
        c = 299792458.e+6 # um/s
        for infile in self.selFileNames:
            print (infile)
            if infile not in self.fileNames:
                try:
                    t1=timer()
                    aor, hk, gratpos, flux = readData(infile+".fits")
                    detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum = aor
                    obsdate, coords, offset, angle, za, altitude, wv = hk
                    #print("Data read from ",infile)
                    onspectra, offspectra = multiSlopes(flux)
                    #offspectra, onspectra = multiSlopes(flux)
                    #print("Slope fitted")
                    # Take the median over the space dimension of the detector
                    # flat
                    # We can apply flats here    
                    wave = []
                    dw  = []
                    for gp in gratpos:
                        l,lw = waveCal(gratpos=gp, order=order, array=detchan,dichroic=dichroic,obsdate=obsdate)
                        wave.append(np.transpose(l))
                        dw.append(np.transpose(lw))
                    # Compute flux
                    dw = np.array(dw)
                    wave = np.array(wave)
                    dnu = c/wave * dw/wave
                    print(np.shape(wave), np.shape(onspectra))
                    onspectra /= dnu
                    offspectra /= dnu
                    print('units V/s/Hz')
                    onspectra = applyFlats(wave, onspectra, detchan, order, dichroic, obsdate)    
                    offspectra = applyFlats(wave, offspectra, detchan, order, dichroic, obsdate)    
                    print('flats applied')
                    onspectrum = np.nanmedian(onspectra,axis=2)
                    offspectrum = np.nanmedian(offspectra,axis=2)
                    wave = np.nanmedian(wave, axis=2)
                    if nodbeam == 'A':
                        print('nodbeam is: ', nodbeam)
                        spectrum = onspectrum - offspectrum
                        skyspectrum = offspectrum
                    else:
                        spectrum = offspectrum - onspectrum
                        skyspectrum = onspectrum
                        
                    #obj = Obs(skyspectrum,coords,offset,angle,altitude,za,wv,nodbeam,filegpid,filenum,gratpos,detchan,order,obsdate,dichroic)
                    obj = Obs(wave,spectrum, skyspectrum, coords,offset,angle,altitude,
                              za,wv,nodbeam,filegpid,filenum,gratpos,detchan,order,obsdate,dichroic)
                    t2=timer()
                    print ("Fitted: ", infile, " ",np.shape(spectrum), " in: ", t2-t1," s")
                    # Call this with a signal from thread
                    self.updateObjects.newObj.emit(obj)
                    self.updateFilenames.emit(infile)
                    self.updateFigures.emit(infile)
                except:
                    print ("Problems with file: ", infile)
                    self.updateExclude.emit(infile)
            else:
                print ('updating figure')
                if self.firstRun:
                    self.updateFigures.emit(infile)
        if self.processAll:
            self.updateStatus.emit('next')
        else:
            self.updateStatus.emit('All processed')
        # Disconnect signal at the end of the thread
        self.updateObjects.newObj.disconnect()
        print ("Done adding observations thread")


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


        #os.register_at_fork()
        # Start exploring directory
        from fifimon.fifitools import exploreDirectory
        cwd = os.getcwd()
        self.files, self.start, self.fgid, self.ch = exploreDirectory(cwd+"/")
        print ("Directory scanned")
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
        self.flightmode = False
        self.excludedFiles = []
            
        # Menu
        self.file_menu = self.menuBar().addMenu('&File')
        self.quit_program = QAction('Quit',self,shortcut='Ctrl+q',triggered=self.fileQuit)
        self.file_menu.addAction(self.quit_program)
        
        self.help_menu = self.menuBar().addMenu('&Help')
        self.about_code = QAction('About',self,shortcut='Ctrl+h',triggered=self.about)
        self.help_menu.addAction(self.about_code)

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
        exitAction = QAction(QIcon(os.path.join(path0,'icons','exit.png')), 'Exit the program', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(self.closeEvent)
        exitAction.setMenuRole(QAction.NoRole)
        hideAction = QAction(QIcon(os.path.join(path0,'icons','list.png')), 'List of File Group IDs', self)
        hideAction.setShortcut('Ctrl+L')
        hideAction.triggered.connect(self.changeVisibility)
        procAction = QAction(QIcon(os.path.join(path0,'icons','gears.png')), 'Process all data', self)
        procAction.setShortcut('Ctrl+P')
        procAction.triggered.connect(self.processAll)
        modsaveAction = QAction(QIcon(os.path.join(path0,'icons','save.png')), 'Scale and update data', self)
        modsaveAction.setShortcut('Ctrl+M')
        modsaveAction.triggered.connect(self.modsaveData)
        flyAction = QAction(QIcon(os.path.join(path0,'icons','plane.png')), 'Switch between flight and QA mode', self)
        flyAction.setShortcut('Ctrl+F')
        flyAction.triggered.connect(self.flightMode)
        
        
        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.addAction(hideAction)
        self.tb.addAction(procAction)
        self.tb.addAction(modsaveAction)
        self.tb.addAction(flyAction)
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

        # Periodical update
        self.computeStatus = False
        cwd = os.getcwd()
        files = os.listdir(cwd)
        self.numberOfFiles = len(files)
        self.timer = QTimer(self)
        
    def changeVisibility(self):
        state = self.lf.isVisible()
        self.lf.setVisible(not state)

    def processAll(self):
        ''' Action to process everything in the directory '''

        self.processItem = 0
        item = self.fgidList[self.processItem]
        #self.sb.showMessage("Processing FileGroupID: "+item.tostring(),5000)
        self.sb.showMessage("Processing FileGroupID: "+item,5000)
        self.addObs(item, False, True)

    def flightMode(self):
        ''' switch between flight and QA mode'''    
        self.flightmode ^= True
        if self.flightmode:
            self.timer.timeout.connect(self.update_fifimon)
            self.timer.start(5000)
            # In the future it should have the possibility to reduce
            # blue and red channel at once (and display them)
        else:
            self.timer.timeout.disconnect()
        self.pc.flightmode(self.flightmode)


    def modsaveData(self):
        from astropy.io import fits


        buttonReply = QMessageBox.question(self, 'Rescale and save', "Do you want to rescale the data and update them ?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if buttonReply == QMessageBox.Yes:
            print('Yes clicked.')
            for fn,o in zip(self.fileNames,self.obs):
                hdulist = fits.open(fn+".fits")
                header = hdulist[0].header
                scidata = hdulist[1].data
                sciheader = hdulist[1].header
                hdulist.close()
                data = np.float32(scidata.DATA)+2**15
                ncycles = header['C_CYC_B']
                ngrat = header['G_PSUP_B']
                ndata = ncycles*4*32
                # Compute scaling factors
                spec = o.spec
                ly = []
                for j in np.arange(ngrat):
                    y = spec[j,:]
                    ly.append(y)
                ly = np.array(ly)
                medline = np.nanmedian(ly,axis=0)
                alpha = np.ones(ngrat)
                for j in np.arange(ngrat):
                    alpha[j] = np.nanmedian(medline/spec[j,:])
                # Apply to raw data
                for j in np.arange(ngrat):
                    data[j*ndata:(j+1)*ndata,:,:] *= alpha[j]
                # Save back to the file
                data = np.int16(data-2**15)
                scidata.DATA = data
                fits.update(fn+".fits", scidata, sciheader, 'FIFILS_rawdata')
                print(fn+' updated')
        else:
            print('No clicked.')
            pass

            
        
        
        
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
                'wave': o.wave.tolist(),
                'spec': o.spec.tolist(),
                'sky': o.sky.tolist()
            }
        with io.open('fifimon.json', mode='w') as f:
            str_= json.dumps(data,indent=2,sort_keys=True,separators=(',',': '), ensure_ascii=False)
            #,encoding='utf-8')
            #f.write(unicode(str_))
            #f.write(str(str_, 'utf-8'))
            f.write(str_)
        pass    

    def loadData(self):
        import json
        from collections import OrderedDict
        #from fifimon.fifitools import Obs
        print ('loading previous reduced data')
        # Read the json file as ordered dictionary to preserve the
        # order the objects were saved
        with open('fifimon.json') as f:
            #data = json.load(f, object_pairs_hook=OrderedDict, encoding='utf-8')
            data = json.load(f, object_pairs_hook=OrderedDict)

        for key,value in data.items():
            print (key)
            self.fileNames.append(key)
            spec=np.array(value['spec'])
            sky=np.array(value['sky'])
            wave=np.array(value['wave'])
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
            obsdate=value['obsdate']#.decode('latin-1').encode('utf-8')
            order=value['order']
            dichroic=value['dichroic']
            self.obs.append(Obs(wave,spec,sky,(ra,dec),(x,y),angle,(alt[0],alt[1]),(za[0],za[1]),(wv[0],wv[1]),nod,fgid,n,gp,ch,order,obsdate,dichroic))
        print ("Loading completed ")


    def saveExcludedFiles(self):
        ''' Save the excluded file list if some is found '''    
        if len(self.excludedFiles) == 0:
            return
        else:
            # Open file to collect problematic files after renaming previous file
            try:
                os.rename('fifimon.exclude','fifimon.exclude.old')
            except:
                pass
            self.excludefile = open('fifimon.exclude','w')
            # Order list of excluded files
            self.excludedFiles.sort()
            # Save them
            for f in self.excludedFiles:
                self.excludefile.write(f)
                self.excludefile.write('\n')
            # Close the file
            self.excludefile.close()
        
        
    def fileQuit(self):
        self.saveData()
        self.saveExcludedFiles()
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        from fifimon import __version__
        # Get path of the package
        path0,file0 = os.path.split(__file__)
        file=open(os.path.join(path0,"copyright.txt"),"r")
        message=file.read()
        message = 'FIFIMON - version: ' + __version__ + '\n' + message
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
        print('mask ',np.shape(mask), 'ch ',np.shape(self.ch))
        # Select the channel
        channels = list(set(self.ch[mask]))
        self.channel = channels[0]
        # Check if there are files
        print ("there are ", np.sum(mask)," observations")
        if np.sum(mask) == 0:
            print ("There no files in this channel")
            return
        
        
        selFiles = self.files[mask]
        selFileNames = [os.path.splitext(os.path.basename(f))[0] for f in selFiles]
        selFileNames = np.array(selFileNames)
            
        if self.computeStatus:
            return
        else:
            self.computeStatus = True

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
                print ("Failed to read file ", selFileNames[0])
                return
        
        # Starts thread to compute slopes
        self.addObsThread = AddObsThread(selFileNames, self.fileNames, processAll, firstRun)
        self.addObsThread.updateObjects.newObj.connect(self.update_objects)
        self.addObsThread.updateFigures.connect(self.update_figures)
        self.addObsThread.updateExclude.connect(self.update_exclude)
        self.addObsThread.updateFilenames.connect(self.update_filenames)
        self.addObsThread.updateStatus.connect(self.update_status)
        #self.addObsThread.finished.connect(self.addObsThread.quit) # Quit once finished
        print('starting new thread with files ', selFileNames)
        self.addObsThread.start()
        
        # Used memory check
        #pid = os.getpid()
        #py = psutil.Process(pid)
        #memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
        #print('memory use:', memoryUse)

    def update_objects(self, obj):
        self.obs.append(obj)

    def update_exclude(self, infile):
        '''We could save the list of files, order it, and save it at the end in one shot'''
        self.excludedFiles.append(infile)
        
    def update_status(self, status):
        ''' process next item is available '''
        if status == 'All processed':
            print('Thread completed')
            self.addObsThread.quit()
            self.addObsThread.wait()
            self.computeStatus = False
            return
        
        self.processItem += 1
        if self.processItem < len(self.fgidList):
            item = self.fgidList[self.processItem]
            #self.sb.showMessage("Processing FileGroupID: "+item.tostring(),5000)
            self.sb.showMessage("Processing FileGroupID: "+item,5000)
            self.addObs(item, False, True)
        self.addObsThread.quit()
         
    def update_filenames(self, infile):
        #print "updating filenames with file ", infile
        self.fileNames.append(infile)
      
    def update_figures(self, infile):
        #from timeit import default_timer as timer

        # Select obs corresponding to infile
        n = self.fileNames.index(infile)
        o = self.obs[n]
        #print "updating figure with filename ", n
        
        t1=timer()
        self.fc.updateFigure(o.nod,o.wave,o.spec,o.sky,infile,o.za,o.alt,o.wv,o.gp,o.order,o.obsdate,o.dichroic,o.ra,o.dec)
        t2=timer()
        self.pc.updateFigure(o.nod,o.ra,o.dec,o.x,o.y,o.angle,infile)
        t3=timer()
        print ("Plotted ",infile," in ",t2-t1," and ",t3-t2," s")


        
    def update_fifimon(self):
        '''
        This module is called automatically every n seconds (n=10, typically) to look
        for new files in the directory, update the lists and, eventually, the plot
        '''

        from fifimon.fifitools import exploreDirectory
        
        # Check if there are new files
        cwd = os.getcwd()
        files = os.listdir(cwd)
        if len(files) == self.numberOfFiles:
            return
        else:
            self.numberOfFiles = len(files)

        
        files, start, fgid, ch = exploreDirectory(cwd)
        print('files ', files)
        print('ch ', ch)
        print('fgid ' ,fgid)
        # Check if new fileGroupID or files added
        try:
            if not self.files:
                print ('There were no files !')
                self.files=[]
                self.fgid=[]
                self.ch=[]
        except:
            pass

        newfiles = 0
        for f, fg, ch_ in zip(files, fgid, ch):
            if f not in self.files:
                print ("updating file list ...")
                self.files = np.append(self.files, f)
                self.fgid = np.append(self.fgid, fg)
                self.ch = np.append(self.ch, ch_)
                newfiles += 1
                # update plot (if same fileGroupID)
                # self.addObs(None, True)
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
                    self.addObs(fg, False)
                else:
                    pass  
        if newfiles > 0:
            self.addObs(self.fileGroupId, True)

                
""" Main code """
        
def main():
    from fifimon import __version__
    app = QApplication(sys.argv)
    app.setApplicationVersion(__version__)
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    aw = ApplicationWindow()
    aw.setGeometry(100, 100, width*0.9, width*0.45)
    progname = 'FIFI Monitor'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()
    app.deleteLater() # to avoid weird QThread messages
