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
import numpy.linalg
import OpenGL.GL

VERSION = "0.1"


from pyqtgraph.dockarea import *

class myGLViewWidget(gl.GLViewWidget):
    def __init__(self, **kwargs):
        gl.GLViewWidget.__init__(self, **kwargs)
        self.reset()        
        self.initview = np.array([[ 1.,  0.,  0.],
                                  [ 0.,  0.,  1.],
                                  [ 0., -1.,  0.]])

                        
    def reset(self):
        self.opts['azimuth'] = -90
        self.opts['elevation'] = 00
#        self.opts['azimuth'] = 45
#        self.opts['elevation'] = 30 
 
        self.opts['center'] = QtGui.QVector3D(0.0, 0.0, 0.0)
        self.opts['fov'] = 60
        self.opts['distance'] = 30        
        self.update()

    def rotateXZ(self):
        self.opts['azimuth'] = -90
        self.opts['elevation'] = 00
        self.update()

    def rotateXY(self):
        self.opts['azimuth'] = -90
        self.opts['elevation'] = 90
        self.update()

    def rotateYZ(self):
        self.opts['azimuth'] = 0
        self.opts['elevation'] = 00
        self.update()

    def currentView(self):
        return np.dot(OpenGL.GL.glGetFloatv(OpenGL.GL.GL_MODELVIEW_MATRIX)[:3, :3], self.initview)


class popupPlot(QtGui.QDialog):
    def __init__(self, parent=None):
        super(popupPlot, self).__init__()


class mainWindow(QtGui.QMainWindow):
    def __init__(self, load=None, loadsnapshot=None, options={}, **kwargs):
        super(mainWindow, self).__init__()
        self._setupUI()
        if load is not None:
            self.loaddata_gadget(load)
        if loadsnapshot is not None:
            self.loadsnapshot_gadget(loadsnapshot, basic=options.basic)        
        

    def _setupUI(self):
        area = DockArea()
        self.setCentralWidget(area)
        self.resize(1000,800)
        self.setWindowTitle('Gadget Visualization')

        self.setupMenu()

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        self.dock1 = Dock("Main Window", size=(800, 800))
        self.dock2 = Dock("Console", size=(800,200))
        self.dock3 = Dock("View", size=(200,800))
        self.dock4 = Dock("Options", size=(200,800))
        self.dock1.hideTitleBar()
        self.dock2.hideTitleBar()
#        dock3.hideTitleBar()

        area.addDock(self.dock1, 'left')
        area.addDock(self.dock3, 'right', self.dock1)
        area.addDock(self.dock4, 'top', self.dock3)
        area.addDock(self.dock2, 'bottom', self.dock1)
        area.moveDock(self.dock3, 'above', self.dock4)

        ## Test ability to move docks programatically after they have been placed
        #area.moveDock(dock1, 'top', dock2)     ## move dock4 to top edge of dock2

        ## Add widgets into each dock
    #    self.sn = None
        self.filePath = '.'
        self.sn = snapshot()
        self.vis = sphvis()

        self.ViewWidget = myGLViewWidget()
        self.ViewWidget.addItem(self.vis)
        self.dock1.addWidget(self.ViewWidget)

        ## setup embedded IPython terminal
        self._setupIpythonWidget()

        ## setup control widget
        self._setupControls()

    def _setupIpythonWidget(self):
        namespace={'pg': pg, 'np': np, 'vis': self.vis, 'sn': self.sn, 'vw': self.ViewWidget, 'main': self }
        #w2 = pg.console.ConsoleWidget(namespace=namespace)
        
        self.ipythonWidget = QIPythonWidget()
        self.ipythonWidget.pushVariables(namespace)
        self.dock2.addWidget(self.ipythonWidget)
        welcomeText = "\nPySPHViewer v%s\n" % (VERSION) \
            + "The data of the currently visible Snapshot can be accessed with th 'sn' object\n" \
            + "import numpy as np, pyqtgraph as pg\n" \
            + "\n"
        self.ipythonWidget.printText(welcomeText)

    def _setupControls(self):
        w3 = pg.LayoutWidget()
        self.dock3.addWidget(w3)

        label =  QtGui.QLabel("") 
        w3.addWidget(label)
        w3.nextRow()
        ##buttons
        buttons = [('Reset', self.ViewWidget.reset, 'r', -1),
                   ('Save Rotation', self.saveCurrentView, None, -1),
                   ('X-Z', self.ViewWidget.rotateXZ, 'c', None),
                   ('X-Y', self.ViewWidget.rotateXY, 'v', None),
                   ('Y-Z', self.ViewWidget.rotateYZ, 'b', 1),
                   ]
        for text, action, shortcut, cs in buttons:
            btn = QtGui.QPushButton(text)
            btn.setMinimumWidth(15)
            if cs:
                w3.addWidget(btn, colspan=cs)
                w3.nextRow()
            else:
                w3.addWidget(btn)
            if shortcut:
                btn.setShortcut(shortcut)
            btn.clicked.connect(action)
            
        w4 = pg.LayoutWidget()
        self.dock3.addWidget(w4)
        ## spins
        spins = [
            ("Pointsize", pg.SpinBox(value=1.0, bounds=[0, 100]), self.vis.sizeChange),
            ("Alpha", pg.SpinBox(value=0.5, bounds=[0., 1.]), self.vis.alphaChange)
            ]        
        for text, spin, slot in spins:
            label = QtGui.QLabel("%16s" % text)
            w4.addWidget(label)
            w4.addWidget(spin)
            w4.nextRow()
            spin.sigValueChanged.connect(slot)
        label =  QtGui.QLabel("") 
        w4.addWidget(label)

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
        openFile.triggered.connect(self.openGadget)

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
            return fname

    def openGadget(self):
        fname = self.fileSelect()
        if fname:
            self.loadsnapshot_gadget(fname)

    def saveCurrentView(self):
        fname = self.fileSelect()
        with (open(fname, 'w')) as f:
            self.ViewWidget.currentView().transpose().astype(np.float64).tofile(f)

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
