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
        self.setColors()

    def setColors(self, gascol=[1.,0,0,0.5], halocol=[1., 1., 0., 0.5]):        
        lum = 0.5
        self.gasalpha = lum
        self.haloalpha = lum
        self.diskalpha = lum
        self.bulgealpha = lum
        self.staralpha = lum
        
#        self.gascol = [1., 0., 0., self.gasalpha]
        self.gascol = gascol
        self.halocol = [1., 1., 0., self.staralpha]
        self.diskcol = [1., 1., 0., self.staralpha]
        self.bulgecol = [1., 1., 0., self.staralpha]
        self.starcol = [1., 1., 0., self.staralpha]

        self.colors[self.startgas:self.endgas, :] = self.gascol
        self.colors[self.starthalo:self.endhalo, :] = self.halocol
        self.colors[self.startdisk:self.enddisk, :] = self.diskcol
        self.colors[self.startbulge:self.endbulge, :] = self.bulgecol
        self.colors[self.startstars:self.endstars, :] = self.starcol
        self.setData(color=self.colors)

