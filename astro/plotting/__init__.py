import astro
import phaseDiagram
import matplotlib

font = {
#    'family' : 'normal',
#    'weight' : 'bold',
    'size'   : 16
}

matplotlib.rc('font', **font)

astro.snapshot.phaseDiagram = phaseDiagram.phaseDiagram
