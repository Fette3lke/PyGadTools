import numpy as np
import astro

def readParticleDump(filename, **kwargs):
    sn = astro.snapshot(**kwargs)
    file = open(filename,'rb') 
    npart = np.fromfile(file, dtype=np.int32, count=1)
    print npart
    sn.pos = np.fromfile(file, dtype=np.dtype((np.float64,3)), count=npart)
    sn.vel = np.fromfile(file, dtype=np.dtype((np.float64,3)), count=npart)
    sn.id  = np.fromfile(file, dtype=np.int32, count=npart)
    file.close()
    return sn
