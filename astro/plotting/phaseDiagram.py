import numpy as np
import astro
import matplotlib.pyplot as plt

def phaseDiagram(self, *args, **kwargs):
	alpha = kwargs.pop("alpha", 0.7)
	if len(args) == 0:
		symbol = 'b,'
	else:
		symbol = args[0]

	plt.plot(np.log10(self.rho / (astro.protonmass / astro.hydrogen_fraction / self.unit_density_in_cgs)), np.log10(self.temp), 
		symbol, alpha=alpha , **kwargs)
	plt.xlabel(r'$log(\rho) \; [cm^{-3}]$', size=20)
	plt.ylabel(r'$log(T) \; [K]$', size=20)
        

	return plt.gcf(), plt.gca()
