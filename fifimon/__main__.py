#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-

# Creating the GUI


# Reading the data, processing, displaying them
import matplotlib
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

import sys
import numpy as np

class MyMplCanvas(FigureCanvas):
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


class RaDecCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        t = np.arange(0.0, 3.0, 0.01)
        s = np.sin(2*np.pi*t)
        self.axes.plot(t, s)


class FluxCanvas(MyMplCanvas):
    """A canvas that updates itself every 3 seconds with a new plot."""

    #def __init__(self, *args, **kwargs):
    #    MyMplCanvas.__init__(self, *args, **kwargs)
    #    timer = QtCore.QTimer(self)
    #    timer.timeout.connect(self.update_figure)
    #    timer.start(3000)

    def compute_initial_figure(self):
        self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        l = [random.randint(0, 10) for i in range(4)]
        self.axes.cla()
        self.axes.plot([0, 1, 2, 3], l, 'r')
        self.draw()

class myListWidget(QListWidget):

    def Clicked(self,item):
#      QMessageBox.information(self, "ListWidget", "You clicked: "+item.text())
        # 2000 means message erased after 2 seconds
        #self.parent().parent().statusBar().showMessage("Clicked item "+item.text(),2000)
        self.parent().parent().sb.showMessage("Clicked item "+item.text(),2000)
#        self.parent().parent().listFiles.hide()
        self.parent().parent().listFiles.setVisible(False)


        

class ApplicationWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

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
        sc = RaDecCanvas(self.main_widget, width=3, height=4, dpi=100)
        dc = FluxCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.mpl_toolbar = NavigationToolbar(sc, self)
        self.mpl_toolbar2 = NavigationToolbar(dc, self)

        # Actions
        exitAction = QAction(QIcon('icons/exit24.png'), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(qApp.quit)
        hideAction = QAction(QIcon('icons/list.png'), 'Hide', self)
        hideAction.setShortcut('Ctrl+L')
        hideAction.triggered.connect(self.changeVisibility)

        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.addAction(exitAction)
        self.tb.addAction(hideAction)

        # Widget list
        self.listFiles = myListWidget()
        self.listFiles.addItem("Item 1"); 
        self.listFiles.addItem("Item 2");
        self.listFiles.addItem("Item 3");
        self.listFiles.addItem("Item 4");
        self.listFiles.setWindowTitle('PyQT QListwidget Demo')
        self.listFiles.itemClicked.connect(self.listFiles.Clicked)

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

        radecLayout.addWidget(sc)
        radecLayout.addWidget(self.mpl_toolbar)
        radecLayout.addWidget(self.tb)
        fluxLayout.addWidget(dc)
        fluxLayout.addWidget(self.mpl_toolbar2)
        fluxLayout.addWidget(self.sb)

        splitter1.addWidget(radecWidget)
        splitter1.addWidget(fluxWidget)

        mainLayout.addWidget(splitter1)
        mainLayout.addWidget(self.listFiles)
        self.listFiles.hide()

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)


        self.sb.showMessage("Welcome to FIFI Monitor!", 10000)

    def changeVisibility(self):
        state = self.listFiles.isVisible()
        self.listFiles.setVisible(not state)

        
    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        file=open("copyright.txt","r")
        message=file.read()
        QtWidgets.QMessageBox.about(self, "About", message)

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
