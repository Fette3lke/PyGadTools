import astro
from . import readParticleDump 
from . import readSubfind

__all__ = ["readParticleDump", "readSubfind"]

readParticleDump = readParticleDump.readParticleDump

astro.snapshot.readSubfind = readSubfind.readSubfind
astro.snapshot.getSubHalo  = readSubfind.getSubHalo

