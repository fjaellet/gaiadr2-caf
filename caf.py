# Compute dN/dt, using the formalism in Lamers et al. (2005, L05), 
# assuming a truncated power-law cluster IMF
#
# Author: M. Gieles, with additions by F. Anders (ICCUB)
# Last modified: 27.11.2020
#
# Reference: https://ui.adsabs.harvard.edu/abs/2020arXiv200601690A/abstract

from __future__ import division
import numpy
from pylab import log
from scipy.integrate import quad
from scipy import interpolate
import pylab as plt
from pylab import log10
from astropy.table import Table

def mu_ev(t):
    """
    Remaining mass fraction due to stellar evolution 
    (equations 2 & 3 and Table 1 in L05)
    """
    a_ev, b_ev, c_ev = 7, 0.255, -1.805
    q_ev  = 0 if t<12.5e6 else 10**((log10(t)-a_ev)**b_ev + c_ev)
    mu_ev = 1 - q_ev
    return mu_ev
    
def dNdm(M, t, Mlim, Mmax, gamma, t0, CFR, Nwise=True):
    """
    Evolved mass function at some time, 
    
    CFR_M = int M*phi dm = A ln(Mmax/Mlim)
    CFR_N = int phi dm = A [1/Mlim - Mmax] ~= A/Mlim
    """
    # Derive MF constant A from the CFR
    if Nwise:
        A = CFR * (1/Mlim - 1/Mmax)
    else:
        A = CFR/log(Mmax/Mlim)
    # Eq 15 in L05 and/or Eq (28) in G09 for alpha = -2 and without the exp part
    Delta_g = gamma*t/t0
    dNdm    = A * M**(gamma-1) * mu_ev(t) / (M**gamma + Delta_g)**(1+1/gamma)
    return dNdm
 
def caf_L05(CFR, lMmax, t0):
    """ 
    Generate age distribution (CAF) following L05.
    
    Input:
        CFR:   cluster formation rate [number of clusters per yr]
        lMmax: log Mmax [upper cutoff of the Cluster IMF]
        t0:    characteristic destruction time
    Output:
        t:     time [yr]
        dNdt:  CAF [number of clusters per yr]
    """
    # Assume CIMF is a -2 power-law powerlaw
    # Follow the evolve MF notation in Gieles 2009, MNRAS

    # Fixed model parameters:
    lMlim = 2 # minimum mass of mass function [Msun]
    ltmin, ltmax = 6, 10
    gamma = 0.62
    
    tmin, tmax, Mlim, Mmax = 10**ltmin, 10**ltmax, 10**lMlim, 10**lMmax

    # Define some arrays
    nt, nm = 300, 300
    lte = numpy.linspace(ltmin, ltmax, nt+1)
    lt = 0.5*(lte[0:-1] + lte[1:])
    t = 10**lt

    lMe = numpy.linspace(lMlim, lMmax, nm+1)
    lM = 0.5*(lMe[0:-1] + lMe[1:])
    M = 10**lM
    
    dNdt = numpy.zeros(nt)

    # Integrate over mass functions at each age
    for it in range(nt):
        # Note that the upper limit is decreasing because of stellar evolution
        dNdt[it] = quad(dNdm, Mlim, Mmax, args=(t[it], Mlim, Mmax*mu_ev(t[it]), gamma, t0, CFR))[0]
    return t, dNdt

def caf_L05burst2(CFR, lMmax, t0, burstamp, ltbegin, ltend):
    """ 
    Generate age distribution (CAF) following L05 & add a starburst
    """
    # Assume CIMF is a -2 power-law powerlaw
    # Follow the evolve MF notation in Gieles 2009, MNRAS

    # Fixed model parameters:
    lMlim = 2 # minimum mass of mass function [Msun]
    ltmin, ltmax = 6, 10
    gamma = 0.62
    
    tmin, tmax, Mlim, Mmax = 10**ltmin, 10**ltmax, 10**lMlim, 10**lMmax

    # Define some arrays
    nt, nm = 100, 100
    lte = numpy.linspace(ltmin, ltmax, nt+1)
    lt = 0.5*(lte[0:-1] + lte[1:])
    t = 10**lt

    lMe = numpy.linspace(lMlim, lMmax, nm+1)
    lM = 0.5*(lMe[0:-1] + lMe[1:])
    M = 10**lM
        
    dNdt = numpy.zeros(nt)
    # Model CFR as constant + Gaussian burst
    CFRt = CFR * numpy.ones(len(lt))
    CFRt[(lt>ltbegin) & (lt<ltend)] = CFR + burstamp
    
    # Integrate over mass functions at each age
    for it in range(nt):
        # Note that the upper limit is decreasing because of stellar evolution      
        dNdt[it] = quad(dNdm, Mlim, Mmax, args=(t[it], Mlim, Mmax*mu_ev(t[it]), gamma, t0, CFRt[it]))[0]
    return t, dNdt


if __name__ ==  "__main__":
    # Run an example
    t0 = 3.e6    # yr
    CFR = 5e-6   # Total cluster formation rate [/yr]
    lMmax = 4.5  # Maximum mass of mass function [Msun]


    plt.ion()
    plt.clf()
    plt.ylabel(r'$\log\,dN/dt$ [yr$^{-1}$]')
    plt.xlabel(r'$\log\,t$ [yr]')

    t, dNdt = caf_L05(CFR, lMmax, t0)
    plt.plot(log10(t), log10(dNdt),'b-')
