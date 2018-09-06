#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np

NANG = 91
NMXX = 150e3

def mie_bohren_huffman(x, refrel, nang=NANG):
    """
    Compute mie scattering based on Bohren and Huffman theory

    Input
    ------
        x        Size parameter = k*radius = 2pi/lambda * radius   
                  (lambda is the wavelength in the medium around the scatterers)
        refrel   Refraction index (n in complex form for example:  1.5+0.02*i;
        nang     Number of angles for S1 and S2 function in range from 0 to pi/2

    Output
    -------
        S1, S2   Funtion which correspond to the (complex) phase functions
        Qext     Extinction efficiency
        Qsca     Scattering efficiency 
        Qback    Backscatter efficiency
        gsca     Asymmetry parameter

    Note
    ----
    This file is converted from [mie.m](http://atol.ucsd.edu/scatlib/index.htm)
    Bohren and Huffman originally published the code in their book on light scattering.

    Source: http://scatterlib.googlecode.com/files/bhmie_herbert_kaiser_july2012.py
    """
    if (nang > 1000):
        raise AttributeError(f"Require NANG = {nang} <= 1000")

    if (nang < 2):
        raise AttributeError(
            f"Require NANG = {nang} > 1 in order to calculate scattering intensities")

    ang = .5 * np.pi / (nang-1)
    mu = np.cos(np.arange(0, nang, 1) * ang)

    # Series expansion terminated after NSTOP terms
    # Logarithmic derivatives calculated from NMX on down

    xstop = x + 4 * np.power(x, 1/3) + 2
    # xstop = x + 4 * np.power(x, 1/3) + 10  # Old form

    ymod = abs(x * refrel)

    nmx = np.fix( max(xstop, ymod) + 15 )

    # BTD experiment 91/1/15: add one more term to series and compare resu<s
    #      NMX = AMAX1(XSTOP, YMOD) + 16
    # test: compute 7001 wavelen > hs between .0001 and 1000 micron
    # for a = 1.0 micron SiC grain.  When NMX increased by 1, only a single
    # computed number changed (out of 4*7001) and it only changed by 1/8387
    # Conclusion: we are indeed retaining enough terms in series!

    if (nmx > NMXX):
        raise ValueError(f"nmx = {nmx} > NMXX = {NMXX} for |m|x = {ymod}")
    
    s1_1 = np.zeros(nang, dtype=np.complex128)
    s1_2 = np.zeros(nang, dtype=np.complex128)
    s2_1 = np.zeros(nang, dtype=np.complex128)
    s2_2 = np.zeros(nang, dtype=np.complex128)
    pi = np.zeros(nang, dtype=np.complex128)
    tau = np.zeros(nang, dtype=np.complex128)
    pi0 = np.zeros(nang, dtype=np.complex128)
    pi1 = np.ones(nang, dtype=np.complex128)

    # Logarithmic derivative D(J) calculated by downward recurrence
    # beginning with initial value (0,0) at J = NMX

    nn = int(nmx) - 1
    d = np.zeros(nn + 1, dtype=np.complex128)
    for n in range(0, nn):
        en = (nmx - n)/(x * refrel)
        d[nn-n-1] = en - 1 / (d[nn-n] + en)

    #*** Riccati-Bessel functions with real argument X
    #    calculated by upward recurrence

    psi0 = np.cos(x)
    psi1 = np.sin(x)
    chi0 = -np.sin(x)
    chi1 = np.cos(x)
    xi1 = psi1 - chi1 *1j
    qsca = 0
    gsca = 0
    p = -1

    nstop = int(xstop)
    for n in range(0, nstop):
        en = n + 1
        fn = (2 * en + 1)/(en * (en + 1))

    # for given N, PSI  = psi_n        CHI  = chi_n
    #              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
    #              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
    # Calculate psi_n and chi_n
        psi = (2 * en - 1) * psi1/x - psi0
        chi = (2 * en - 1) * chi1/x - chi0
        xi = psi - chi *1j

    #*** Store previous values of AN and BN for use
    #    in computation of g=<np.cos(theta)>
        if (n > 0):
            an1 = an
            bn1 = bn

    #*** Compute AN and BN:
        an = (d[n]/refrel + en/x) * psi - psi1
        an = an / ((d[n]/refrel + en/x) * xi - xi1)
        bn = (refrel * d[n] + en/x) * psi - psi1
        bn = bn / ((refrel * d[n] + en/x) * xi - xi1)

    #*** Augment sums for Qsca and g=<np.cos(theta)>
        qsca += (2 * en + 1) * (abs(an)**2 + abs(bn)**2)
        gsca += ((2 * en + 1) / (en * (en + 1))) * \
            (np.real(an) * np.real(bn) + np.imag(an) * np.imag(bn))

        if (n > 0):
            gsca += ((en - 1) * (en + 1)/en) * \
            (np.real(an1) * np.real(an) + np.imag(an1) * np.imag(an) + \
             np.real(bn1) * np.real(bn) + np.imag(bn1) * np.imag(bn))

    #*** Now calculate scattering intensity pattern
    #    First do angles from 0 to 90
        pi = np.copy(pi1)
        tau = en * mu * pi - (en + 1) * pi0
        s1_1 += fn * (an * pi + bn * tau)
        s2_1 += fn * (an * tau + bn * pi)

    #*** Now do angles greater than 90 using PI and TAU from
    #    angles less than 90.
    #    P=1 for N=1,3,...% P=-1 for N=2,4,...
    #   remember that we have to reverse the order of the elements
    #   of the second part of s1 and s2 after the calculation
        p = -p
        s1_2 += fn * p * (an * pi - bn * tau)
        s2_2 += fn * p * (bn * pi - an * tau)

        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = psi1 - chi1 *1j

    #*** Compute pi_n for next value of n
    #    For each angle J, compute pi_n+1
    #    from PI = pi_n , PI0 = pi_n-1
        pi1 = ((2 * en + 1) * mu * pi - (en + 1) * pi0) / en
        pi0 = np.copy(pi)

    #*** Have summed sufficient terms.
    #    Now compute QSCA, QEXT, QBACK and GSCA

    # We have to reverse the order of the elements of the second part of s1 and s2
    s1 = np.concatenate((s1_1, s1_2[-2::-1]))
    s2 = np.concatenate((s2_1, s2_2[-2::-1]))
    gsca = 2 * gsca/qsca
    qsca = 2/x**2 * qsca
    qext = 4/x**2 * np.real(s1[0])

    # More common definition of the backscattering efficiency,
    # so that the backscattering cross section really
    # has dimension of length squared
    qback = 4 * ( abs(s1[2*(nang-1)])/x )**2
    #qback = ((abs( s1[2 * nang - 2])/x )**2 )/np.pi  # Old form

    return s1, s2, qext, qsca, qback, gsca


def mie(wvln, nr, ni, r, nang=NANG):
    """
    Compute Mie cross-sections and phase function based
        on Bohren and Huffman theory

    Input
    ------
        wvln (float)    Wavelength (m)
        nr   (float)    Particule real optical index ()
        ni   (float)    Particule real imaginary index ()
        r    (float)    Particule radius (m)
        [nang] (int)    Number of angles for the phase function 
                         (range from 0 to pi/2)

    Outputs
    --------
        qsct  (float)    Scattering cross section (m^-2)
        qext  (float)    Extinction cross section (m^-2)
        qabs  (float)    Absorption cross section (m^-2)
        gg    (float)    Asymmetry parameter
        theta (float[])  Phase function angles (radians)
        P     (float[])  Phase function ()
    """
    Xm = 2. * np.pi * r / wvln
    s1, s2, Qe, Qs, _, gg = mie_bohren_huffman(Xm, complex(nr, ni), nang)
    qsct = Qs * np.pi * r**2
    qext = Qe * np.pi * r**2
    qabs = qext - qsct

    S11 = .5 * (abs(s2)**2 + abs(s1)**2)
    theta = np.linspace(0, np.pi, len(s1))
    norm = .5 * np.trapz(S11*np.sin(theta), x=theta)
    P = S11/norm
    return qsct, qext, qabs, gg, theta, P

