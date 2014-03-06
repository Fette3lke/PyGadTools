import pyqtgraph.opengl as gl
import numpy as np

class sphvis(gl.GLScatterPlotItem):
    def __init__(self, load=None, **kwargs):
        gl.GLScatterPlotItem.__init__(self, **kwargs)
        if not load == None:
            self.loaddata_gadget(load)
        else:
            self.setData(pos = np.empty((1,3)), pxMode=False)

    def loaddata_gadget(self, snapshot, lum = 0.5):
        self.startgas = snapshot.startgas
        self.starthalo = snapshot.starthalo
        self.startdisk = snapshot.startdisk
        self.startbulge = snapshot.startbulge
        self.startstars = snapshot.startstars
        self.startbh = snapshot.startbh
        self.endgas = snapshot.endgas
        self.endhalo = snapshot.endhalo
        self.enddisk = snapshot.enddisk
        self.endbulge = snapshot.endbulge
        self.endstars = snapshot.endstars
        self.endbh = snapshot.endbh

        self.npart = np.sum(snapshot.head.npart)
        self.pos = snapshot.pos
        self.colors= np.ones([self.npart, 4])
        self.sizes = np.ones(self.npart) * 0.1
        self.setData(pos = self.pos, color=self.colors, size=self.sizes, pxMode=False)
        self.lum = lum
        self.initColors()
        self.initSizes()

    def initColors(self):
        self.gasalpha = 1.
        self.haloalpha = 1.
        self.diskalpha = 1.
        self.bulgealpha = 1.
        self.staralpha = 1.
        
        self.gascol = [1., 0., 0.]
        self.halocol = [0., 0., 1.]
        self.diskcol = [0., 1., 0.]
        self.bulgecol = [0., 1., 0.]
        self.starcol = [1., 1., 0.]
        self.setColors()

    def initSizes(self):
        self.gassize = 0.1
        self.halosize = 0.1
        self.disksize = 0.1
        self.bulgesize = 0.1
        self.starsize = 0.1
        self.mastersize = 1.0
        self.setSizes()

    def setColors(self):                
        self.colors[self.startgas:self.endgas, :] = self.gascol + [self.gasalpha * self.lum ]
        self.colors[self.starthalo:self.endhalo, :] = self.halocol + [self.haloalpha * self.lum]
        self.colors[self.startdisk:self.enddisk, :] = self.diskcol + [self.haloalpha * self.lum]
        self.colors[self.startbulge:self.endbulge, :] = self.bulgecol + [self.bulgealpha * self.lum]
        self.colors[self.startstars:self.endstars, :] = self.starcol + [self.staralpha * self.lum]
        self.setData(color=self.colors)

    def setSizes(self):
        self.sizes[self.startgas:self.endgas] = self.gassize * self.mastersize
        self.sizes[self.starthalo:self.endhalo] = self.halosize * self.mastersize
        self.sizes[self.startdisk:self.enddisk] = self.disksize * self.mastersize
        self.sizes[self.startbulge:self.endbulge] = self.bulgesize * self.mastersize
        self.sizes[self.startstars:self.endstars] = self.starsize * self.mastersize
        self.setData(size=self.sizes)
        
    def sizeChange(self, sb):
        self.mastersize = sb.value()
        self.setSizes()
        self.update()

    def alphaChange(self, sb):
        self.lum = sb.value()
        self.setColors()
        self.update()
