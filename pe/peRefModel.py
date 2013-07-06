#!/usr/bin/env python
"""
Nearfield Lloyd-Mirror Pattern, Section 6.4.3, Eq. 6.109
Computational Ocean Acoustics, 2ed, Jensen, Kuperman, Porter, Schmidt

A single point source in a fluid halfspace, the boundary condition at
the surface (z=0) is pressure-release (p=0).
"""

from pylab import *

## Reference Half-space solution
# JKPS Eq. (6.109)
# All distances in meters.

f0 = 100.   # Hz
zs = 75.    # source depth
c = 1500.   # m/s
k = 2 * np.pi * f0 / c;

rgridfrac, zgridfrac = 0.25, 0.20
rstart, rend = 0.1, 500
zstart, zend = 0.1, 500
dr = rgridfrac * 2 * np.pi /k # grid at a fraction of a wavelength
dz = rgridfrac * 2 * np.pi /k; print "Grid Size:", dz

rpoints = (rend-rstart)/dr
zpoints = (zend-zstart)/dz

R = arange(rstart, rend, dr); print shape(R)
Z = arange(zstart, zend, dz); print shape(Z)

r, z = np.meshgrid(R,Z); print "Stop 1"

R1 = sqrt(r**2 + (z-zs)**2)
R2 = sqrt(r**2 + (z+zs)**2)

p = exp(1j*k*R1)/R1 - exp(1j*k*R2)/R2
p  = -10*log(abs(p))

# Contour levels chosen to best match figure in text
v = [15, 29, 35.7, 41, 46.5, 52, 57.5, 63.3,105]


contourf(r, z, p, v, cmap=cm.coolwarm_r)
contour(r, z, p, v, colors='#2F4F4F', linewidth=0.001)

xlabel('Depth (m)')
ylabel('Range (m)')
title('F = 100 Hz, SD = 75 m.     Eq. (6.109), Comp. Ocean Acous.',fontsize=10)
ylim(zend,zstart)
xlim(rstart,rend)
axes().set_aspect('equal')
show()

##savefig('lloyd-mirror.ps',format='PS',dpi=300)
##savefig('lloyd-mirror.jpg',format='PNG',dpi=300)
