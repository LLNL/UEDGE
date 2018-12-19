#!/bin/sh

#from matplotlib import pylab
import matplotlib.pyplot as plt
import numpy
import tables
import sys

f=tables.openFile("2200_edge_1.h5")
#f=tables.openFile("2200_edge_10.h5")
vsf=tables.openFile("2200_coupled-lineout_9.vsh5")

# Electron temperature
tes=f.root.bbb.tesCore.read()[:,1:]
tesCoreAve=numpy.average(tes,0)
tesCoreStd=numpy.std(tes,0)
maxTeLims=tes.max()              #for figure limits
#tes=f.root.bbb.tesSol.read()
tes=vsf.root.outer_temperature_electron.read()
tesSolAve=tes
tesSolStd=0.

# Ion temperature
tis=f.root.bbb.tisCore.read()[:,1:]
tisCoreAve=numpy.average(tis,0)
tisCoreStd=numpy.std(tis,0)
maxTiLims=tis.max()
#tis=f.root.bbb.tisSol.read()
tis=vsf.root.outer_temperature_H2p1.read()
tisSolAve=tis
tisSolStd=0.

# Density
nis=f.root.bbb.nisCore.read()[:,1:,0]/2.e4
nisCoreAve=numpy.average(nis,0)
nisCoreStd=numpy.std(nis,0)
nis=vsf.root.outer_density.read()/2.e4
#nis=f.root.bbb.nisSol.read()[:,:,0]/2.e4
nisSolAve=nis
nisSolStd=0.

# Neutral gas
ngs=f.root.bbb.ngsCore.read()[:,1:,0]
ngsCoreAve=numpy.average(ngs,0)
ngsCoreStd=numpy.std(ngs,0)
ngs=f.root.bbb.ngsSol.read()
#ngs=f.root.bbb.ngsSol.read()
ngsSolAve=numpy.average(ngs,0)
ngsSolStd=numpy.std(ngs,0)

# Grid
gridCore=vsf.root.mpEdgeGrid.read()[1:]
gridEdge=vsf.root.mpOuterGrid.read()

maxTlim=numpy.maximum(maxTeLims,maxTiLims)*1.1

# Now switch to a more OO interface to exercise more features.
fig, axs = plt.subplots(nrows=2, ncols=1)
ax = axs[0]
ax.errorbar(gridCore, tesCoreAve, yerr=tesCoreStd, label="Electron temperature")
ax.errorbar(gridCore, tisCoreAve, yerr=tisCoreStd, label="Ion temperature")
#ax.errorbar(gridEdge, tesSolAve, yerr=tesSolStd, label="Electron temperature")
#ax.errorbar(gridEdge, tisSolAve, yerr=tisSolStd, label="Ion temperature")
ax.set_title('Average Temperatures')
ax.set_ylim(0,maxTlim)
ax.legend(loc='lower left',frameon=False)

ax = axs[1]
ax.errorbar(gridCore, nisCoreAve, yerr=nisCoreStd, label="Ion density/2.e4")
ax.errorbar(gridCore, ngsCoreAve, yerr=ngsCoreStd, label="Neutral gas density")
#ax.errorbar(gridEdge, nisSolAve, yerr=nisSolStd, label="Ion density")
#ax.errorbar(gridEdge, ngsSolAve, yerr=0., label="Neutral gas density")
ax.set_title('Average Densities')
ax.legend(loc='upper left',frameon=False)

fig.suptitle('Average quantities showing standard deviation.')

plt.show()
