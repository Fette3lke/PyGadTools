import numpy as np
import astro
import matplotlib.pyplot as plt

def phaseDiagram(self, *args, **kwargs):
	alpha = kwargs.pop("alpha", 0.7)
	if len(args) == 0:
		symbol = 'b,'
	else:
		symbol = args[0]

	plt.plot(np.log(self.rho / (astro.protonmass / astro.hydrogen_fraction / self.unit_density_in_cgs)), np.log(self.temp), 
		symbol, alpha=alpha , **kwargs)
	plt.xlabel(r'$log(\rho) \; [cm^{-3}]$')
	plt.ylabel(r'$log(T) \; [K]$')
	return plt.gcf(), plt.gca()