import pyqtgraph.opengl as gl
import numpy as np
import sys

class sphvis(gl.GLScatterPlotItem):
    def __init__(self, load=None, **kwargs):
        gl.GLScatterPlotItem.__init__(self, **kwargs)
        if not load == None:
            self.loaddata_gadget(load)
        else:
            self.setData(pos = np.empty((1,3)), pxMode=False)
        self.mastersize = 1.
        self.lum = 0.5
        self.step = kwargs.pop("step", 1)
        self.ntypes = kwargs.pop("ntypes", 6)
        self.typenames = ["gas", "halo", "disk", "bulge", "stars", "bh"]
        self.starttype = np.zeros(self.ntypes, dtype=np.int64)
        self.endtype = np.zeros(self.ntypes, dtype=np.int64)
        self.alpha = np.ones(self.ntypes) * kwargs.pop("alpha", 1)
        self.col = [[1., 1., 1]] * self.ntypes
        self.psize = [kwargs.pop("size", 0.1)] * self.ntypes        
        self.showTypes = set([])
        self.drawtype = []

    def loaddata_gadget(self, snapshot):
        self.starttype[0] = snapshot.startgas
        self.starttype[1] = snapshot.starthalo
        self.starttype[2] = snapshot.startdisk
        self.starttype[3] = snapshot.startbulge
        self.starttype[4] = snapshot.startstars
        self.starttype[5] = snapshot.startbh
        self.endtype[0] = snapshot.endgas
        self.endtype[1] = snapshot.endhalo
        self.endtype[2] = snapshot.enddisk
        self.endtype[3] = snapshot.endbulge
        self.endtype[4] = snapshot.endstars
        self.endtype[5] = snapshot.endbh

        self.npart = np.sum(snapshot.head.npart)
        self.posdata = snapshot.pos 
        self.colors= np.ones([self.npart, 4])
        self.sizes = np.ones(self.npart) * 0.1
        self.setData(pos = self.posdata[::self.step], color=self.colors[::self.step], size=self.sizes[::self.step], pxMode=False)
        self.initColors()
        self.showTypes = set(range(self.ntypes))
        self.setDrawTypes()
        self.Update()

    def loaddata_custom(self, **kwargs):
        self.setData(**kwargs)

    def Update(self):
        self.setPos()
        self.setColors()
        self.setSizes()
        self.update()

    def initColors(self):                
        self.col[0] = [1., 0., 0.]
        self.col[1] = [0., 0., 1.]
        self.col[2] = [0., 1., 0.]
        self.col[3] = [0., 1., 0.] 
        self.col[4] = [1., 1., 0.]
        self.col[5] = [1., .5, 0.]
        self.setColors()

    def setDrawTypes(self):
        self.drawtype = np.array([], dtype=np.int64)
        for t in self.showTypes:
            self.drawtype = np.concatenate((self.drawtype, np.arange(self.starttype[t], self.endtype[t])))        
        self.Update()

    def changeDrawTypes(self, state, ptype):
        self.showTypes ^= set([ptype])
        self.setDrawTypes()

    def setStep(self, step):
        if not step > 0:
            return
        self.step = step
        self.Update()

    def setPos(self):
        self.setData(pos = self.posdata[self.drawtype[::self.step]])

    def setColors(self):                
        for t in range(self.ntypes):
            self.colors[self.starttype[t]:self.endtype[t], :] = self.col[t] + [self.alpha[t] * self.lum ]
        self.setData(color=self.colors[self.drawtype[::self.step]])

    def setSizes(self):
        for t in range(self.ntypes):
            self.sizes[self.starttype[t]:self.endtype[t]] = self.psize[t] * self.mastersize
        self.setData(size=self.sizes[self.drawtype[::self.step]])
        
    def sizeChange(self, val):
        self.mastersize = (val/100. * 4)**2
        self.Update()

    def alphaChange(self, val):
        self.lum = val/100.
        self.Update()
        
    def stepChange(self, val):
        self.setStep(val)
        self.Update()



