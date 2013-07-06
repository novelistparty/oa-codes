#!/usr/bin/env python
"""
Nearfield Lloyd-Mirror Pattern, Section 6.4.3, Eq. 6.109
Computational Ocean Acoustics, 2ed, Jensen, Kuperman, Porter, Schmidt

A single point source in a fluid halfspace, the boundary condition at
the surface (z=0) is pressure-release (p=0).
"""

from pylab import *
import os, csv

## Reference Half-space solution
# JKPS Eq. (6.109)
# All distances in meters.

f0 = 100.   # Hz
zs = 75.    # source depth
c = 1500.   # m/s
k = 2 * np.pi * f0 / c;

rgridfrac, zgridfrac = 0.25, 0.25
rstart, rend = 0.01, 500
zstart, zend = 0.01, 500
dr = rgridfrac * 2 * np.pi /k # grid at a fraction of a wavelength
dz = rgridfrac * 2 * np.pi /k; print "Grid Size:", dz
print(k)
rpoints = (rend-rstart)/dr
zpoints = (zend-zstart)/dz

R = arange(rstart, rend, dr); print shape(R)
Z = arange(zstart, zend, dz); print shape(Z)

##R = np.linspace(rstart, rend, rpoints); print shape(R)
##Z = np.linspace(zstart, zend, zpoints); print shape(Z)

r, z = np.meshgrid(R,Z); print "Stop 1"

R1 = sqrt(r**2 + (z-zs)**2)
R2 = sqrt(r**2 + (z+zs)**2)

p = exp(1j*k*R1)/R1 - exp(1j*k*R2)/R2
p  = -10*log10(abs(p)**2)

# Contour levels chosen to best match figure in text
v = [5, 25, 31, 37, 43, 49, 55, 105]

############# Reference solution #########
figure("Reference Solution")
contourf(r, z, p, v, cmap=cm.coolwarm_r)
colorbar()
contour(r, z, p, v, colors='#2F4F4F', linewidth=0.001)

ylabel('Depth (m)')
xlabel('Range (m)')
title('F = 100 Hz, SD = 75 m. Exact Reference Solution.',fontsize=14)
ylim(zend,zstart)
xlim(rstart,rend)
axes().set_aspect('equal')

##savefig('lloyd-mirror-ref.pdf',format='PDF',dpi=300)
savefig('lloyd-mirror-ref.jpg',format='JPEG',dpi=300)

############ MATLAB Solution #####################
###refpath = os.getcwd()
###basepath = 'U:\classes\ECE576_compmethods\pe_model';
###os.chdir(basepath)

#psimat = loadtxt('psimat.txt', dtype=double)
#R = loadtxt('R.txt', dtype=double)
#Z = loadtxt('Z.txt', dtype=double)
## Change back to script directory
###os.chdir(refpath)

#figure("MATLAB Solution")
#contourf(R, Z, psimat, v, cmap=cm.coolwarm_r)
#colorbar(orientation='horizontal')
#contour(R, Z, psimat, v, colors='#2F4F4F', linewidth=0.001)

#ylabel('Depth (m)')
#xlabel('Range (m)')
#title('F = 100 Hz, SD = 75 m. Tappert Equation, Gaussian Source',fontsize=14)
#ylim(200,0)
#xlim(rstart,1500)
#axes().set_aspect('equal')

#savefig('lloyd-mirror-wrong-large.png',format='PNG',dpi=300)

################# Diff Solution #################
##diff1 = p - psimat[0:-1,:]
##print(shape(diff1))
##figure("Difference")
##contourf(r, z, diff1,  cmap=cm.coolwarm_r)
##colorbar()
##contour(r, z, diff1, colors='#2F4F4F', linewidth=0.001)
##
##ylabel('Depth (m)')
##xlabel('Range (m)')
##title('F = 100 Hz, SD = 75 m. Difference in Solutions',fontsize=14)
##ylim(500,zstart)
##xlim(rstart,500)
##axes().set_aspect('equal')
##show()
