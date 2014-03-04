# routines for reading headers and data blocks from Gadget snapshot files
# usage e.g.:
#
# import readsnap as rs
# header = rs.snapshot_header("snap_063.0") # reads snapshot header
# print header.massarr
# mass = rs.read_block("snap_063","MASS",parttype=5) # reads mass for particles of type 5, using block names should work for both format 1 and 2 snapshots
# print "mass for", mass.size, "particles read"
# print mass[0:10]
#
# before using read_block, make sure that the description (and order if using format 1 snapshot files) of the data blocks
# is correct for your configuration of Gadget 
#
# for mutliple file snapshots give e.g. the filename "snap_063" rather than "snap_063.0" to read_block
# for snapshot_header the file number should be included, e.g."snap_063.0", as the headers of the files differ
#
# the returned data block is ordered by particle species even when read from a multiple file snapshot

import numpy as np
import os
import sys
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import copy

parsec_in_cm = 3.08567758e018
kpc_in_cm = 3.08567758e021
msun_in_g = 1.9891e33
G_in_cgs = 6.6742e-8
crit_dens = (2.78e-8) / kpc_in_cm**3 * (1e10*msun_in_g)          # in cgs * h**2 at z=0
GAS  = 0
HALO = 1
DISK = 2
BULGE= 3
STARS= 4
BH   = 5

HUB = 0.72
BFRAC = 0.044 / 0.26

def power_law(x, a, b):
  return a*x**b

def linear(x,a,b):
  return a*x+b

  
# ----- class for snapshot ----- 

class snapshot(object):
  def __init__(self, filename=None, unit_length_in_cm=kpc_in_cm, unit_mass_in_g=1.9891e43, **kwargs):
    self.ntypes = 6
    self.center = np.array([0., 0., 0.])
    self.unit_length_in_cm = unit_length_in_cm
    self.unit_mass_in_g = unit_mass_in_g
    if filename != None:
      self.readgadget(filename, **kwargs)
    
  def readgadget(self, filename, export=False, convert=True, basic=False, calcDist=False, loadrot=None):
    self.basic = basic
    self.head = snapshot_header(filename)
    self.head.npart = self.head.nall
    self.npart = self.head.npart
    self.ntot = sum(self.head.npart)
    self.type = np.zeros(self.ntot, dtype=np.int)    

    self.pos  = read_block(filename, "POS ")  
    self.vel  = read_block(filename, "VEL ")  
    self.ids  = read_block(filename, "ID  ")  
    self.mass = read_block(filename, "MASS") 
    
    if not basic > 1:
      if self.head.npart[0]:
        self.u    = read_block(filename, "U   ")
        self.rho  = read_block(filename, "RHO ")
        if self.head.cooling:
          self.ne   = read_block(filename, "NE  ")
          self.nh   = read_block(filename, "NH  ")
        self.hsml = read_block(filename, "HSML")
      if self.head.sfr:
        self.sfr  = read_block(filename, "SFR ")
      if self.head.age and self.head.npart[4]:
        self.age  = read_block(filename, "AGE ")
      if self.head.metals and not basic and self.head.npart[4]:
        self.let = read_block(filename, "LET ")
        self.imass = read_block(filename, "INIM")
        self.metstars = read_block(filename, "Z   ", parttype=4, csformat = 1)
        self.metgas = read_block(filename, "Z   ", parttype=0, csformat = 1)
        self.temp = read_block(filename, "CSTE", csformat = 1)
        self.pot = read_block(filename, "POT ", csformat = 1)

      if self.head.npart[5]:
        self.bhmass = read_block(filename, "BHMA", csformat = self.head.metals)
        self.bhmdot = read_block(filename, "BHMD", csformat = self.head.metals)
        self.acrb = read_block(filename, "ACRB", csformat = self.head.metals)


    self.setOffsets()
    
    if convert:
      self.convert()

    if calcDist:
      self.calcDistances()

    self.assign_names()
      
    if export:
      self.export()

    if not loadrot == None:
      self.loadrotfile(loadrot)

  def convert(self):
    self.mass /= self.head.hubble
    self.pos *= self.head.time / self.head.hubble
    self.vel *= np.sqrt(self.head.time)
#    self.physical = True

  def setOffsets(self):
    self.offset = {'pos' : [0, self.ntot], 'vel' : [0, self.ntot],
                   'ids' : [0, self.ntot], 'mass' : [0, self.ntot],
                   'type': [0, self.ntot]}
    if not self.basic > 1:
      if self.head.cooling:
        self.offset.update({'ne' : [0, self.npart[0]], 'nh' : [0, self.npart[0]]})
      self.offset.update({'u' : [0, self.npart[0]], 'rho' : [0, self.npart[0]], 'hsml' : [0, self.npart[0]]})
      if self.head.sfr:
        self.offset.update({'sfr' : [0, self.npart[0]]})
      if self.head.age and sum(self.head.npart[4:6]):
        self.offset.update({'age' : [sum(self.npart[0:4]), sum(self.npart[4:6])]})
      if self.head.metals and not self.basic and self.head.npart[4]:
        self.offset.update({'let' : [sum(self.npart[0:4]), self.npart[4]], 
                            'imass' : [sum(self.npart[0:4]), self.npart[4]],
                            'metstars' : [sum(self.npart[0:4]), self.npart[4]],
                            'metgas' : [0, self.npart[0]],
                            'temp' : [0, self.npart[0]],
                            'pot' : [0, self.ntot]
                            })
      if self.head.npart[5]:  
        self.offset.update({'bhmass' : [sum(self.npart[0:5]), self.npart[5]],
                            'bhmdot' : [sum(self.npart[0:5]), self.npart[5]],
                            'acrb' : [sum(self.npart[0:5]), self.npart[5]]
                            })

  def slice(self, indices, docopy=False, assign=True):
    """
    create a slice of the snapshot, containing the data of a subset of the particles given by their indices
    """    
    sn = snapshot()
    if len(indices) == self.ntot and type(indices[0]) == np.bool_:
      indices = np.arange(self.ntot)[indices]
    ind = indices
    for entry in self.__dict__:
      if entry in self.dontcopy:
        continue
      if entry in self.offset:
        ind, = np.where ((indices >= self.offset[entry][0]) & (indices < sum(self.offset[entry])))
        ind  = indices[ind] - self.offset[entry][0]
#        print entry,  self.offset[entry][0], ind
        if docopy:
          sn.__dict__[entry] = self.__dict__[entry][ind].copy()
        else:
          sn.__dict__[entry] = self.__dict__[entry][ind]
      else:
#        print 'deepcopy', entry
#        print type(self.__dict__[entry])
#        if not entry in self.parttypes:
        sn.__dict__[entry] = copy.deepcopy(self.__dict__[entry])
    for t in range(sn.ntypes):
      i, = np.where(sn.type == t)
      sn.head.npart[t] = len(i)
      sn.head.nall[t] = len(i)
    sn.npart = sn.head.npart
    sn.ntot = sum(sn.npart)
    sn.setOffsets()
    if assign:
      sn.assign_names()
    return sn

#  def gas(self):
#    return self.slice(np.arange(self.startgas,self.endgas))

  def _makevarname(self, str1, str2):
    varname =  str1 + str2
    self.dontcopy.add(varname)
    return varname
    
  def assign_names(self):
    parttypes = ['gas', 'halo', 'disk', 'bulge', 'stars', 'bh']
    self.parttypes = parttypes
    startind = 0
    itype = 0
    self.dontcopy = set()
    for i in xrange(len(parttypes)):
      endind = startind + self.head.npart[i]
      varname = self._makevarname('start', parttypes[i])
      self.__dict__[varname] = startind
      varname = self._makevarname('end', parttypes[i])
      self.__dict__[varname] = endind
      varname = self._makevarname('n', parttypes[i])
      self.__dict__[varname] = self.head.npart[i]
      varname = self._makevarname('x', parttypes[i])
      self.__dict__[varname] = self.pos[startind:endind, 0]
      varname = self._makevarname('y', parttypes[i])
      self.__dict__[varname] = self.pos[startind:endind, 1]
      varname = self._makevarname('z', parttypes[i])
      self.__dict__[varname] = self.pos[startind:endind, 2]
      varname = self._makevarname('vx', parttypes[i])
      self.__dict__[varname] = self.vel[startind:endind, 0]
      varname = self._makevarname('vy', parttypes[i])
      self.__dict__[varname] = self.vel[startind:endind, 1]
      varname = self._makevarname('vz', parttypes[i])
      self.__dict__[varname] = self.vel[startind:endind, 2]
      varname = self._makevarname('m', parttypes[i])
      self.__dict__[varname] = self.mass[startind:endind]
      varname = self._makevarname('id', parttypes[i])
      self.__dict__[varname] = self.ids[startind:endind]
      self.type[startind:endind] = itype
      itype += 1
#      self.__setattr__(parttypes[i], self.slice(np.arange(startind,endind), assign=False)) 
      self.dontcopy.add(parttypes[i])
      self.__dict__[parttypes[i]] = self.slice(np.arange(startind,endind), assign=False) 
      startind = endind
#
#      self.__setattr__(parttypes[i], self.test)


  def loadrotfile(self, fname):
    with open(fname, 'r') as f:
      self.rot = np.fromfile(f, dtype=np.dtype((np.float64, (3,3))), count=1)[0]
      self.rot = np.transpose(self.rot)
    self.rotate(self.rot.astype(np.float32))

  def rotate(self, matrix):
    posdum = np.dot( self.pos, matrix)
    self.pos[:] = posdum
    veldum = np.dot( self.vel, matrix)
    self.vel[:] = veldum
#    self.assig_names()

  def calcTemp(self):
      gamma_minus_1 = 5.0/3 -1
      boltzmann = 1.3806e-16
      protonm= 1.6726e-24
      hydr_frac= 0.76
      yhelium = ( 1 - hydr_frac ) / ( 4 * hydr_frac )
      mass_in_g= 1.989e43
      length_cm=3.085678e21
      vel_in_cm_per_s=1e5
      time_in_s= length_cm / vel_in_cm_per_s
      energy_cgs = mass_in_g * (length_cm**2) / (time_in_s**2);
      mu = (1 + 4 * yhelium) / (1 + yhelium + self.ne);          
      temp = gamma_minus_1 / boltzmann * self.u * protonm * mu;
      temp *= energy_cgs / mass_in_g;
      self.temp = temp

  def calcDistances(self):
    self.dist = np.sqrt(np.sum((self.pos-self.center)**2, dtype=np.float64, axis=1))
    self.sind = np.argsort(self.dist)

  def export(self):
    """
    Adds snapshot data to global namespace
    """
    for k in self.__dict__:
      sys.modules['__builtin__'].__dict__[k] = self.__dict__[k]

  def getSfr(self, redshifts=None, masses=None, bins=np.logspace(-1, 1, 50)):
    """
    returns Star formation rate (inferred from stellar ages) and bincenters (redshift)
    """
    if redshifts == None:
      redshifts = (a2z(self.stars.age))
    if masses == None:
      masses = (self.mstars)

    sf, bins = np.histogram(redshifts, weights=masses, bins=bins)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    time = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        time[i] = galage(bins[i]) - galage(bins[i+1])
    sfr = sf/time
    return sfr, bincenters

  def bins2d(self, z=None, bins=np.arange(10, dtype=np.float32), type=None):
    """
    Create projected bins for particles of type given by array (e.g. type = [0,4] for baryons)
    z          - cutoff distance in 3rd dimension

    return binind, bincenters, area, sigma
    binind     - binned indices
    bincenters - center of bins
    vol        - surface area of bins
    sigma      - 1/sqrt(number of particles in bin)
    """
    ind = []
    if z > 0:
      ind, = np.where (np.abs(self.pos[:,2]) < z)
    else:
      ind = np.arange(self.ntot)

    if not type == None:
      bytetype = np.sum(1<<np.array(type))
      i, = np.where( (1<<self.type[ind]) & bytetype)
      ind = ind[i]

    dist2d = np.sqrt(self.pos[ind,0]**2 + self.pos[ind,1]**2)
    low = bins[0]
    binind = []
    area = []
    sigma = []
    bincenters = []
    for high in bins[1:]:
      i, = np.where( (dist2d >= low) & (dist2d < high) )
      binind.append(ind[i])
      area = np.append(area, np.pi * (high**2 - low**2))
      if len(i) > 0:
        err = 1/np.sqrt(len(i))
      else:
        err = 1
      sigma = np.append(sigma, err)
      bincenters = np.append(bincenters, (low+high)/2.)
      low = high
    return binind, bincenters, area, sigma

  def bins3d(self, bins=np.arange(10, dtype=np.float32), type=None):
    """
    Create spherical bins for particles of type given by array (e.g. type = [0,4] for baryons)
    return binind, bincenters, vol, sigma
    binind     - binned indices
    bincenters - center of bins
    vol        - volume of bins
    sigma      - 1/sqrt(number of particles in bin)
    """
    if not hasattr(self, 'dist'):
      self.calcDistances()
    ind = []
    if not type == None:
      bytetype = np.sum(1<<np.array(type))
      ind, = np.where( (1<<self.type) & bytetype)
    else:
      ind = np.arange(self.ntot)

    low = bins[0]
    binind = []
    vol = []
    sigma = []
    bincenters = []
    for high in bins[1:]:
      i, = np.where( (self.dist[ind] >= low) & (self.dist[ind] < high) )
      binind.append(ind[i])
      vol.append( (4./3.) * np.pi * (high**3 - low**3) )
      if len(i) > 0:
        err = 1/np.sqrt(len(i))
      else:
        err = 1
      sigma = np.append(sigma, err)
      bincenters = np.append(bincenters, (low+high)/2.)
      low = high
    return binind, bincenters, vol, sigma

  def densityProfile(self, **kwargs):
    binind, bincenters, vol, sigma = self.bins3d(**kwargs)
    density = []
    for i in range(len(binind)):
      ind = binind[i]
      mass = np.sum(self.mass[ind], dtype=np.float64)
      density = np.append(density, mass/vol[i])
    return density, bincenters

  def surfaceDensity(self, **kwargs):
    binind, bincenters, area, sigma = self.bins2d(**kwargs)
    surfdens = []
    for i in range(len(binind)):
      ind = binind[i]
      mass = np.sum(self.mass[ind], dtype=np.float64)
      surfdens = np.append(surfdens, mass/area[i])
    return surfdens, bincenters

  def surfaceDensitySFR(self, fit=None, **kwargs):
    binind, bincenters, area, sigma = self.bins2d(type=[0], **kwargs)
#    for high in bins[1:]:
    sfr = []
    surfdens = []
    area *= self.unit_length_in_cm**2 / parsec_in_cm**2
    for i in range(len(binind)):
      ind = binind[i]
      mass = np.sum(self.mass[ind], dtype=np.float64) * self.unit_mass_in_g / msun_in_g
      sfr  = np.append(sfr, np.sum(self.sfr[ind], dtype=np.float64)/area[i])
      surfdens = np.append(surfdens, mass/area[i])

    i, = np.where (sfr > 0)
    surfdens = surfdens[i]
    sfr = sfr[i]
    bincenters=bincenters[i]
    sigma = sigma[i]
    if type(fit) == list:
#      fit[0], fit[1] = curve_fit(power_law, surfdens*msun_in_g/1e6, sfr)
      fit[0], fit[1] = curve_fit(linear, np.log10(surfdens), np.log10(sfr), sigma=sigma)
    return surfdens, sfr, bincenters

  def vcirc(self, bins=np.arange(0, 10, 0.1, dtype=np.float32), plot = True, **kwargs):
    """
    Calculate velocity curves for different particle types.
    """
    if not hasattr(self, 'dist'):
      self.calcDistances()
#    dist = np.sqrt(np.sum((self.pos-self.center)**2, axis=1))
    binmass = np.zeros([self.ntypes, len(bins)-1])
    vcirc = np.zeros([self.ntypes+1, len(bins)-1])
    bincenters = np.zeros(len(bins)-1)
    
    low = bins[0]
    i = 0

    for j in range(self.ntypes):
      ind, = np.where (self.type == j )
      binmass[j,:], binedges = np.histogram(self.dist[ind], bins=bins, weights=self.mass[ind])

    for high in bins[1:]:
        bincenters[i] = (low+high)/2.
        for j in range(self.ntypes):
            vcirc[j,i] = np.sqrt( (np.sum(binmass[j,0:(i+1)], dtype=np.float64)   / high ) *(G_in_cgs * self.unit_mass_in_g / self.unit_length_in_cm))  * 1e-5 # cm/s to km/s                      
        vcirc[self.ntypes, i] = np.sqrt( (np.sum(binmass[:,0:(i+1)], dtype=np.float64)   / high ) *(G_in_cgs * self.unit_mass_in_g / self.unit_length_in_cm))  * 1e-5             
        low = high
        i += 1

    if plot:
      plt.plot(bins[1:], vcirc[1,:], linewidth=2, **kwargs)
      plt.plot(bins[1:], vcirc[4,:], linewidth=2, **kwargs)
      plt.plot(bins[1:], vcirc[6,:], linewidth=2, **kwargs)
      plt.legend(["Halo", "Stars", "Total"], loc="lower right")
      plt.xlabel("$R$ [kpc]")
      plt.ylabel("$V_c$ [km/s]")

    return vcirc, bins[1:], binmass

  def r200(self):
    rvir = 0
    cdens = crit_dens * self.head.hubble**2 * (self.head.omega_l + self.head.omega_m * self.head.time**(-3) );
    cdens /= self.unit_mass_in_g / self.unit_length_in_cm**3
    if not hasattr(self, 'dist'):
      self.calcDistances()
    i = 0
    od = 201
    masstot = 0
    while ( ((od>200) or (i<5)) and (i < self.ntot) ):
      ind = self.sind[i]
      masstot += self.mass[ind]
      dens = masstot /( (self.dist[ind])**3 * (4./3.) * np.pi )
      od = dens / cdens
      i += 1
    ind = self.sind[i-1]
    masstot -= self.mass[ind]
    rvir = self.dist[ind]
    return rvir, masstot, od

  def test(self):
    print "1"

  def GasVelHist(self, **histargs):
    gasvel = np.sqrt( np.sum( self.gas.vel**2, axis=1 ) )
    return np.histogram(gasvel, **histargs)


class snapshot_header:
  def __init__(self, filename):

    if os.path.exists(filename):
      curfilename = filename
    elif os.path.exists(filename+".0"):
      curfilename = filename+".0"
    else:
      print "file not found:", filename
#      print "and:", curfilename
      sys.exit()
    filename = curfilename

    self.filename = filename  
    f = open(filename,'rb')    
    blocksize = np.fromfile(f,dtype=np.int32,count=1)
    if blocksize[0] == 8:
      swap = 0
      format = 2
    elif blocksize[0] == 256:
      swap = 0
      format = 1  
    else:
      blocksize.byteswap(True)
      if blocksize[0] == 8:
        swap = 1
        format = 2
      elif blocksize[0] == 256:
        swap = 1
        format = 1
      else:
        print "incorrect file format encountered when reading header of", filename
        sys.exit()
    
    self.format = format
    self.swap = swap
    
    if format==2:
      f.seek(16, os.SEEK_CUR)
    
    self.npart = np.fromfile(f,dtype=np.int32,count=6)
    self.massarr = np.fromfile(f,dtype=np.float64,count=6)
    self.time = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.redshift = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.sfr = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.feedback = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.nall = np.fromfile(f,dtype=np.int32,count=6)
    self.cooling = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.filenum = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.boxsize = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_m = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_l = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.hubble = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.age =  (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.metals =  (np.fromfile(f,dtype=np.int32,count=1))[0]
    
    if swap:
      self.npart.byteswap(True)
      self.massarr.byteswap(True)
      self.time = self.time.byteswap()
      self.redshift = self.redshift.byteswap()
      self.sfr = self.sfr.byteswap()
      self.feedback = self.feedback.byteswap()
      self.nall.byteswap(True)
      self.cooling = self.cooling.byteswap()
      self.filenum = self.filenum.byteswap()
      self.boxsize = self.boxsize.byteswap()
      self.omega_m = self.omega_m.byteswap()
      self.omega_l = self.omega_l.byteswap()
      self.hubble = self.hubble.byteswap()
      self.age =  self.age.byteswap()
      self.metals =  self.metals.byteswap()
     
    f.close()



    

# ----- find offset and size of data block ----- 

def find_block(filename, format, swap, block, block_num, only_list_blocks=False):
  if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()
            
  f = open(filename,'rb')
  f.seek(0, os.SEEK_END)
  filesize = f.tell()
  f.seek(0, os.SEEK_SET)
  
  found = False
  curblock_num = 1
  while ((not found) and (f.tell()<filesize)):
    if format==2:
      f.seek(4, os.SEEK_CUR)
      curblock = f.read(4)
      if (block == curblock):
        found = True
      f.seek(8, os.SEEK_CUR)  
    else:
      if curblock_num==block_num:
        found = True
        
    curblocksize = (np.fromfile(f,dtype=np.int32,count=1))[0]
    if swap:
      curblocksize = curblocksize.byteswap()
    
    # - print some debug info about found data blocks -
    #if format==2:
    #  print curblock, curblock_num, curblocksize
    #else:
    #  print curblock_num, curblocksize
    
    if only_list_blocks:
      print curblock_num,curblock,f.tell(),curblocksize
      found = False
    
    if found:
      blocksize = curblocksize
      offset = f.tell()
    else:
      f.seek(curblocksize, os.SEEK_CUR)
      blocksize_check = (np.fromfile(f,dtype=np.int32,count=1))[0]
      if swap: blocksize_check = blocksize_check.byteswap()
      if (curblocksize != blocksize_check):
        print "something wrong"
        sys.exit()
      curblock_num += 1
      
  f.close()
      
  if ((not found) and (not only_list_blocks)):
    print "Error: block not found %s" % block
    sys.exit()
    
  if (not only_list_blocks):
    return offset,blocksize
 
# ----- read data block -----
 
def read_block(filename, block, parttype=-1, physical_velocities=False, arepo=0, no_masses=False, verbose=False, csformat=0, shift=0):
  if (verbose):
	  print "reading block", block
  
  blockadd=0
  blocksub=0

  if os.path.exists(filename):
    curfilename = filename
  elif os.path.exists(filename+".0"):
    curfilename = filename+".0"
  else:
    print "file not found:", filename
    print "and:", curfilename
    sys.exit()

  head = snapshot_header(curfilename)
#  print type(head.npart)
#  print (head.npart != 0)
#  print (head.massarr == 0)
#  verbose = True
  massblocks = (head.npart != 0) & (head.massarr == 0)
  if not (massblocks).any():
    no_massses = True
  
  if arepo==0:
    if (verbose):	
	    print "Gadget format"
    blockadd=0
  if arepo==1:
    if (verbose):	
	    print "Arepo format"
    blockadd=1	
  if arepo==2:
    if (verbose):
	   print "Arepo extended format"
    blockadd=4	
  if no_masses==True:
    if (verbose):	
	    print "No mass block present"    
    blocksub=1
    
  if (shift):
    blocksub -= shift
		 
  if parttype not in [-1,0,1,2,3,4,5]:
    print "wrong parttype given"
    sys.exit()
  
  
  
  format = head.format
  swap = head.swap
  npart = head.npart
  massarr = head.massarr
  nall = head.nall
  filenum = head.filenum
  if filenum < 1:
    filenum = 1
  redshift = head.redshift
  time = head.time
  if head.npart[5]:
    bh_blockadd = 4
  else:
    bh_blockadd = 0
  del head
  
  # - description of data blocks -
  # add or change blocks as needed for your Gadget version
  data_for_type = np.zeros(6,bool) # should be set to "True" below for the species for which data is stored in the data block
  dt = np.float32 # data type of the data in the block
  if block=="POS ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 2
  elif block=="VEL ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 3
  elif block=="ID  ":
    data_for_type[:] = True
    dt = np.uint32
    block_num = 4
  elif block=="MASS":
    data_for_type[np.where(massarr==0)] = True
    block_num = 5
    if parttype>=0 and massarr[parttype]>0:   
      if (verbose):	
	      print "filling masses according to massarr"   
      return np.ones(nall[parttype],dtype=dt)*massarr[parttype]
  elif block=="U   ":
    data_for_type[0] = True
    block_num = 6-blocksub
  elif block=="RHO ":
    data_for_type[0] = True
    block_num = 7-blocksub
  elif block=="VOL ":
    data_for_type[0] = True
    block_num = 8-blocksub 
  elif block=="CMCE":
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 9-blocksub 
  elif block=="AREA":
    data_for_type[0] = True
    block_num = 10-blocksub
  elif block=="NFAC":
    data_for_type[0] = True
    dt = np.dtype(np.int32)	
    block_num = 11-blocksub
  elif block=="NE  ":
    data_for_type[0] = True
    block_num = 8+blockadd-blocksub
  elif block=="NH  ":
    data_for_type[0] = True
    block_num = 9+blockadd-blocksub
  elif block=="HSML":
    data_for_type[0] = True
    block_num = 10+blockadd-blocksub
  elif block=="SFR ":
    data_for_type[0] = True
    block_num = 11+blockadd-blocksub
  elif block=="AGE ":
    data_for_type[4:6] = True
    block_num = 12+blockadd-blocksub
  elif block=="LET ":
    data_for_type[4] = True
    block_num = 13+blockadd-blocksub
  elif block=="INIM":
    data_for_type[4] = True
    block_num = 14+blockadd-blocksub
  elif block=="Z   ":
    data_for_type[0] = True
    data_for_type[4] = True
    block_num = 13+blockadd-blocksub
    if csformat:
      dt = np.dtype((np.float32,12))
      block_num = 15+blockadd-blocksub
  elif block=="POT ":
    data_for_type[:] = True
    block_num = 13+blockadd-blocksub
    if csformat:
      block_num = 16+blockadd-blocksub+bh_blockadd
  elif block=="CSTE":
    data_for_type[0] = True
    block_num = 17+blockadd-blocksub+bh_blockadd
  elif block=="BHMA":
    data_for_type[5] = True
    block_num = 14+blockadd-blocksub
    if csformat:
      block_num += 2
  elif block=="BHMD":
    data_for_type[5] = True
    block_num = 15+blockadd-blocksub
    if csformat:
      block_num += 2
  elif block=="ACRB":
    data_for_type[5] = True
    block_num = 16+blockadd-blocksub
    if csformat:
      block_num += 2
  elif block=="COOR":
    data_for_type[0] = True
    block_num = -1
  else:
    print "Sorry! Block type", block, "not known!"
    sys.exit()
  # - end of block description -

  if (block_num < 0 and format==1):
    print "Sorry! Block number of", block, "not known! Unable to read this block from format 1 file!"
    sys.exit() 
    
  actual_data_for_type = np.copy(data_for_type)  
  if parttype >= 0:
    actual_data_for_type[:] = False
    actual_data_for_type[parttype] = True
    if data_for_type[parttype]==False:
      print "Error: no data for specified particle type", parttype, "in the block", block   
      sys.exit()
  elif block=="MASS":
    actual_data_for_type[:] = True
    
  allpartnum = np.int64(0)
  species_offset = np.zeros(6,np.int64)
  for j in range(6):
    species_offset[j] = allpartnum
    if actual_data_for_type[j]:
      allpartnum += nall[j]
    
  for i in range(filenum): # main loop over files
    if filenum>1:
      curfilename = filename+"."+str(i)
      
    if i>0:
      head = snapshot_header(curfilename)
      npart = head.npart  
      del head
      
    curpartnum = np.int32(0)
    cur_species_offset = np.zeros(6,np.int64)
    for j in range(6):
      cur_species_offset[j] = curpartnum
      if data_for_type[j]:
        curpartnum += npart[j]
    
    if parttype>=0:
      actual_curpartnum = npart[parttype]      
      add_offset = cur_species_offset[parttype] 
    else:
      actual_curpartnum = curpartnum
      add_offset = np.int32(0)

    curdat = None
    if curpartnum:
      offset,blocksize = find_block(curfilename,format,swap,block,block_num)

      if i==0: # fix data type for ID if long IDs are used
        if block=="ID  ":
          if blocksize == np.dtype(dt).itemsize*curpartnum * 2:
            dt = np.uint64 

      if np.dtype(dt).itemsize*curpartnum != blocksize:
        print "something wrong with blocksize! expected =",np.dtype(dt).itemsize*curpartnum,"actual =",blocksize, " in block ", block
        sys.exit()

      f = open(curfilename,'rb')
      f.seek(offset + add_offset*np.dtype(dt).itemsize, os.SEEK_CUR)  
      curdat = np.fromfile(f,dtype=dt,count=actual_curpartnum) # read data
      f.close()  
      if swap:
        curdat.byteswap(True)  

    if i==0:
      data = np.empty(allpartnum,dt)
    
    for j in range(6):
      if actual_data_for_type[j]:
        if block=="MASS" and massarr[j]>0: # add mass block for particles for which the mass is specified in the snapshot header
          data[species_offset[j]:species_offset[j]+npart[j]] = massarr[j]
        else:
          if parttype>=0:
            data[species_offset[j]:species_offset[j]+npart[j]] = curdat
          else:
            if not (block=="MASS" and npart[j] == 0):
              data[species_offset[j]:species_offset[j]+npart[j]] = curdat[cur_species_offset[j]:cur_species_offset[j]+npart[j]]
        species_offset[j] += npart[j]

    if not curdat == None:
      del curdat

  if physical_velocities and block=="VEL " and redshift!=0:
    data *= math.sqrt(time)

  return data
  
# ----- list all data blocks in a format 2 snapshot file -----

def list_format2_blocks(filename):
  if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()
  
  head = snapshot_header(filename)
  format = head.format
  swap = head.swap
  del head
  
  if (format != 2):
    print "not a format 2 snapshot file"
    sys.exit()
            
  print "#   BLOCK   OFFSET   SIZE"
  print "-------------------------"
  
  find_block(filename, format, swap, "XXXX", 0, only_list_blocks=True)
  
  print "-------------------------"


# ***************************************************************************************************************
#  Non-member functions
# ***************************************************************************************************************

def galage(z, omegam=0.26, omegal=0.74, h=0.72, lookback=False):
    """
    compute time between z (may be np.array) and the Big Bang, or look-back time
    """
    omegak = 1.0 - omegam - omegal;
    h0 = 3.24078e-18 * h;    

    if ((omegam < 0) or (omegam > 1)):
        return 0;

    if (omegam == 1):
        tsec = 2.0/(3.0*h0)/((1.0+z)**1.5);
    elif (omegak < 1.0e-10):
        tmp1 = 2.0/(3.0*np.sqrt(1.0-omegam)*h0);
        tmp2 = np.sqrt((1.0 - omegam)/omegam)/((1.0+z)**1.5);
        tmp3 = np.log(tmp2 + np.sqrt(1.0+tmp2*tmp2));
        tsec = tmp1 * tmp3;
    elif (omegal < 1.0e-10):
        tmp1 = 1.0/h0 * omegam/(2.0*((1.0-omegam)**1.5));
        tmp2 = 2.0*(1.0-omegam)/(omegam*(1.0+z))+1.0;
        tmp3 = 2.0 * \
            np.log(np.sqrt((tmp2+1.0)/2.0) \
            + np.sqrt((tmp2-1.0)/2.0));

        tsec = tmp1 * ( np.sqrt((tmp2**2)-1.0) - tmp3 );

    else:
        return 0;

    if lookback:
        thub = galage(0, omegam, omegal, h)
        return thub - (tsec / 3.155815e7)

    return tsec / 3.155815e7;


def a2z(a):
    return ((1 / a) - 1)


def z2a(z):    
    return (1 / (z + 1))

def conversionEfficiency(doplot = False):
  dark_mass = 10.**(np.arange(24)/4.+10.)
  # star_mass = fltarr(24)
  # star_mass = np.zeros(24)
  beta = 1.057
  gamma = 0.556
  M_1  = 10.**11.884
  m_M_0 = 0.02820
  #for i in range(24):
  #    star_mass[i] = 2* m_M_0* ( (dark_mass[i]/M_1)**(-beta) + (dark_mass[i]/M_1)**(gamma) )**(-1.) * dark_mass[i]
  star_mass = 2* m_M_0* ( (dark_mass/M_1)**(-beta) + (dark_mass/M_1)**(gamma) )**(-1.) * dark_mass
  # oplot,alog10(dark_mass),star_mass/(0.176*dark_mass), thick = 6, linestyle=2, color = 0
  if doplot:
    plt.plot(np.log10(dark_mass), star_mass/(BFRAC * dark_mass), 'k--',linewidth=2 )
  # z=1.
  beta = 0.91
  gamma = 0.43
  M_1  = 10.**11.98
  m_M_0 = 0.0142
  #for i in range(24):
  #    star_mass[i] = 2* m_M_0* ( (dark_mass[i]/M_1)^(-beta) + (dark_mass[i]/M_1)^(gamma) )^(-1.) * dark_mass[i]
  star_mass = 2* m_M_0* ( (dark_mass/M_1)**(-beta) + (dark_mass/M_1)**(gamma) )**(-1.) * dark_mass
  # oplot,alog10(dark_mass),star_mass/(0.176*dark_mass), thick = 6, linestyle=2, color = 200
  if doplot:
    plt.plot(np.log10(dark_mass), star_mass/(BFRAC * dark_mass), 'r--',linewidth=2 )
  # z=1.8
  beta = 1.53
  gamma = 0.41
  M_1  = 10.**12.28
  m_M_0 = 0.0116
  #for i in range(24):
  star_mass = 2* m_M_0* ( (dark_mass/M_1)**(-beta) + (dark_mass/M_1)**(gamma) )**(-1.) * dark_mass
  if doplot:
    plt.plot(np.log10(dark_mass), star_mass/(BFRAC * dark_mass), 'g--',linewidth=2 )
    plt.xlim([10, 14])
    plt.ylim([0,0.5])
    plt.legend(['z=0', 'z=1', 'z=2'], 'upper right')
