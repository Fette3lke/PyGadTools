#!/usr/bin/env python
from embeddedIpython import QIPythonWidget 
import pyqtgraph as pg
from PyQt4 import QtGui, QtCore
#from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
from pyqtgraph.widgets.MatplotlibWidget import MatplotlibWidget
import numpy as np
import pyqtgraph.opengl as gl
from astro import *
from sphvis import *
from optparse import OptionParser
import matplotlib



from pyqtgraph.dockarea import *

class myGLViewWidget(gl.GLViewWidget):
    def init(self, **kwargs):
        gl.GLViewWidget.__init__(self, **kwargs)
        self.reset()
                        
    def reset(self):
        self.opts['azimuth'] = -90
        self.opts['elevation'] = 00
#        self.opts['azimuth'] = 45
#        self.opts['elevation'] = 30 
 
        self.opts['center'] = QtGui.QVector3D(0.0, 0.0, 0.0)
        self.opts['fov'] = 60
        self.opts['distance'] = 30
        self.update()


class popupPlot(QtGui.QDialog):
    def __init__(self, parent=None):
        super(popupPlot, self).__init__()


class mainWindow(QtGui.QMainWindow):
    def __init__(self, load=None, loadsnapshot=None, options={}, **kwargs):
        super(mainWindow, self).__init__()
        self.setupUI()
        if load is not None:
            self.loaddata_gadget(load)
        if loadsnapshot is not None:
            self.loadsnapshot_gadget(loadsnapshot, basic=options.basic)        
        

    def setupUI(self):
        area = DockArea()
        self.setCentralWidget(area)
        self.resize(1000,800)
        self.setWindowTitle('Gadget Visualization')

        self.setupMenu()

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        d1 = Dock("Main Window", size=(800, 800))
        d2 = Dock("Console", size=(800,200))
        d3 = Dock("View", size=(200,800))
        d4 = Dock("Options", size=(200,800))
        d1.hideTitleBar()
        d2.hideTitleBar()
#        d3.hideTitleBar()

        area.addDock(d1, 'left')
        area.addDock(d3, 'right', d1)
        area.addDock(d4, 'top', d3)
        area.addDock(d2, 'bottom', d1)
        area.moveDock(d3, 'above', d4)

        ## Test ability to move docks programatically after they have been placed
        #area.moveDock(d1, 'top', d2)     ## move d4 to top edge of d2

        ## Add widgets into each dock
    #    self.sn = None
        self.filePath = '.'
        self.sn = snapshot()
        self.vis = sphvis()

        self.ViewWidget = myGLViewWidget()
        self.ViewWidget.addItem(self.vis)
        self.ViewWidget.reset()
        d1.addWidget(self.ViewWidget)

        ## setup embedded IPython terminal
        namespace={'pg': pg, 'np': np, 'vis': self.vis, 'sn': self.sn, 'vw': self.ViewWidget, 'main': self }
        #w2 = pg.console.ConsoleWidget(namespace=namespace)
        w2 = QIPythonWidget()
        w2.pushVariables(namespace)
        d2.addWidget(w2)


        ## setup control widget
        w3 = pg.LayoutWidget()
        w3.nextRow()
        d3.addWidget(w3)

        ##buttons
        buttons = [('Reset', self.ViewWidget.reset, 'r')]
        for text, action, shortcut in buttons:
            btn = QtGui.QPushButton(text)
            w3.addWidget(btn)
            w3.nextRow()
            if shortcut:
                btn.setShortcut(shortcut)
            btn.clicked.connect(action)

        ## spins
        spins = [
            ("Pointsize", pg.SpinBox(value=0.1, bounds=[0, 100]), self.sizechange),
            ("Pointsize", pg.SpinBox(value=0.1, bounds=[0, 100]), self.sizechange)
            ]        
        for text, spin, slot in spins:
            label = QtGui.QLabel("%16s" % text)
            w3.addWidget(label)
            w3.addWidget(spin)
            w3.nextRow()
            spin.sigValueChanged.connect(slot)
        label =  QtGui.QLabel("") 
        w3.addWidget(label)

    def sizechange(self, sb):
        self.vis.sizes[:] = sb.value()
        self.vis.update()

    def initpopup(self):
        self.popupPlot = popupPlot()
        vbox = QtGui.QVBoxLayout()
        self.popupPlot.setLayout(vbox)
        self.popupPlot.show()     
        self.mw = MatplotlibWidget()
        vbox.addWidget(self.mw)
        self.subplot = self.mw.getFigure().add_subplot(111)
        self.subplot.plot(np.arange(10))
        mw.draw()

    def setupMenu(self):
        openFile = QtGui.QAction('&Open File', self)        
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new Gadget File')
        openFile.triggered.connect(self.fileSelect)

        popupAction = QtGui.QAction('&PlotPopup', self)        
        popupAction.setStatusTip('Exit application')
        popupAction.triggered.connect(self.initpopup)

        exitAction = QtGui.QAction('&Exit', self)        
        exitAction.setShortcut('Esc')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(exitAction)

        self.toolbar = self.addToolBar('Exit')
        self.toolbar.addAction(openFile)
        self.toolbar.addAction(popupAction)
        self.toolbar.addAction(exitAction)


    def fileSelect(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.filePath)
        if fname:
            self.filePath =  QtCore.QFileInfo(fname).path()
            self.loadsnapshot_gadget(fname)

    def loadsnapshot_gadget(self, fname, **kwargs):
        self.sn.readgadget(fname, **kwargs)
        self.loaddata_gadget(self.sn)

    def loaddata_gadget(self, *args):
        self.vis.loaddata_gadget(*args)
        self.statusBar().showMessage('Snapshot loaded')



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        parser = OptionParser()
        parser.add_option("-b", "--basic", dest="basic", default=False,
                  help="only read basic properties")
        (options, args) = parser.parse_args()
        app = QtGui.QApplication([])
        if len(args):
            mainApp = mainWindow(loadsnapshot=args[0], options=options)
        else:
            mainApp = mainWindow(options=options)
        mainApp.show()
        app.exec_()
