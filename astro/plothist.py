import numpy as np

def plothist(axes, data, *args, **kwargs):
    bins  = kwargs.pop('bins', 10)
    label = kwargs.pop('label', None)
    vals, binedges = np.histogram(data, bins=bins)
    binmids = (binedges[1:] + binedges[:-1]) / 2.
    vals = vals / (1. * sum(vals))    
    axes.hist(binmids, weights=vals, bins=binedges, histtype='step', linewidth=2, *args, **kwargs)
    return axes.hist(binmids, weights=vals, bins=binedges, label=label, *args, **kwargs)
