# -*- coding: utf-8 -*-
import numpy as np

from .mie import NANG, mie_bohren_huffman

FORCE = False

def fractals_tomasko_2008(Df, N, Xm, nr, ni, nang=NANG, force=FORCE):
    """
    Compute fractal aerosols scattering based on Tomasko et al. 2008
    empirical model (DOI: 10.1016/j.pss.2007.11.019).

    Input
    ------
        Df      (float) Fractal dimension
        N       (int)   Number of monomers
        Xm      (float) Monomer size parameter (m)
        nr      (float) Particule real optical index
        ni      (float) Particule real imaginary index
        [nang]  (int)   Number of angles for the phase function
                         (range from 0 to pi/2)
        [force] (bool)  Bypass validity checks

    Output
    -------
        Qs_out   Scattering efficiency
        Qs_out   Backscatter efficiency
        Qe_out   Extinction efficiency
        Pij      Funtions which correspond to the (complex) phase functions
    """
    #----------------------------------------
    # Table A1: Single-scattering parameters
    #----------------------------------------
    if not force:
        # Fractal dimension
        if Df != 2:
            raise ValueError("Model tested only for Df = 2 (received Df={:.2f})".format(Df))
        # Number of monomers per aggregate
        if N < 2 or N > 1024:
            raise ValueError("Model tested only for N = 2 - 1024 (received N={})".format(N))
        # Monomer size parameter
        if Xm < 1.e-4 or Xm > 1.5:
            raise ValueError("Model tested only for Xm = 1.e-4 - 1.5 (received Xm={:.2e})".format(Xm))
        # Real part of index of refraction
        if nr < 1.3 or nr > 2:
            raise ValueError("Model tested only for nr = 1.3 - 2 (received nr={:.2f})".format(nr))
        # Imaginary part of index of refraction
        if ni < 0.0 or ni > .7:
            raise ValueError("Model tested only for ni = 0.0 - 0.7 (received ni={:.2f})".format(ni))

    #-----------------------------------------
    # Table A2: Empirical parameters required
    #-----------------------------------------
    # Geometric [Geometric parameters are in units of monomer radius]
    D_cut = 3.194
    # Geometric [Geometric parameters are in units of monomer radius]
    R1 = 1.598
    # Geometric [Geometric parameters are in units of monomer radius]
    R2 = 3.478
    # Geometric [Geometric parameters are in units of monomer radius]
    Rcut = 10000
    # Geometric [Geometric parameters are in units of monomer radius] [Theoretically constrained]
    Rmin = 2
    C_abs_m_1 = 0.606     # Absorption
    E_abs_m_1 = 2.525     # Absorption
    C_abs_m_2 = 1.537     # Absorption
    E_abs_m_2 = 1.273     # Absorption
    C_abs_x_1 = 1.931     # Absorption
    C_abs_x_2 = 1.152     # Absorption
    C_p11_m_1 = 0.054     # P11 vs. P22 offset
    C_p11_m_2 = 1.12      # P11 vs. P22 offset
    E_p11_t_1 = 0.529     # P11 vs. P22 offset
    E_p11_t_2 = 0.531     # P11 vs. P22 offset
    C_p11_m_3 = 0.88      # P11 vs. P22 offset
    E_p11_m_1 = 0.662     # P11 vs. P22 offset
    C_p21_m_1 = 0.06      # P21 for Xm < 1.6 and Xm*M0 < 0.6
    C_p21_m_2 = 0.71      # P21 for Xm < 1.6 and Xm*M0 < 0.6
    E_p21_m_1 = 2.3       # P21 for Xm < 1.6 and Xm*M0 < 0.6
    E_p21_n_1 = 0.25      # P21 for Xm < 1.6 and Xm*M0 < 0.6
    C_p21_ta = 1          # P21 for Xm < 1.6 and Xm*M0 < 0.6
    C_p21_ts = 0          # P21 for Xm < 1.6 and Xm*M0 < 0.6
    C_p21_m_3 = 20        # P21 for Xm < 1.6 and Xm*M0 < 0.6
    C_sca_m_1 = 0.200     # Scattering cross
    C_sca_m_2 = 0.164     # Scattering cross
    C_sca_m_3 = 1.047     # Scattering cross
    C_sca_m_4 = 0.127     # Scattering cross
    C_sca_x_1 = 3.082     # Scattering cross
    C_sca_x_2 = 0.757     # Scattering cross

    # A.2.1. Geometry
    #-----------------
    Rmax = R1 * np.power(N * np.log(Rcut), 1./Df)   # (A.1a)
    R0 = np.arange(Rmin, Rmax, 1/8.)
    if len(R0) < 100:
        R0 = np.linspace(Rmin, Rmax, 100)          # Nb pt > 100

    Nc = (  (1 - np.exp( -np.power(R0/R1, Df)/N ))
          * (1 - np.exp( -np.power(R0/R2, D_cut)))
          + 2./N) / float(1 + 2./N)                # (A.1b)
    F0 = [Nc[0]]
    for ii in range(1, len(R0)-1):
        F0.append(.5 * (Nc[ii+1] - Nc[ii-1]))      # (A.1c)
    F0.append(0)

    # Monomer scattering Mie parameters
    #----------------------------------------------
    s1, s2, Qe, Qs, _, _ = mie_bohren_huffman(Xm, complex(nr, ni), nang)
    Qa = Qe - Qs
    theta = np.linspace(0, np.pi, len(s1))

    S11 = .5 * (np.abs(s2)**2 + np.abs(s1)**2)
    S12 = .5 * (np.abs(s2)**2 - np.abs(s1)**2)
    S33 = .5 * (np.conj(s2)*s1 + s2*np.conj(s1))
    S34 = .5j * (np.conj(s2)*s1 - s2*np.conj(s1))

    norm = .5 * np.trapz(S11*np.sin(theta), x=theta)

    P11_mie = S11 / norm
    P21_mie = S12 / norm             # S12 = S21
    P33_mie = np.real(S33) / norm
    P43_mie = - np.real(S34) / norm  # S34 = S43

    # A.2.2. Monomer scattering Mie
    #----------------------------------------------
    Csca_mon = Qs * np.pi * Xm**2              # (A.2a)
    Cext_mon = Qe * np.pi * Xm**2              # (A.2b)
    Cabs_mon = Qa * np.pi * Xm**2              # (A.2c)

    Ymon = P11_mie * Csca_mon                  # (A.2d)
    Polar_mon = -P21_mie / P11_mie             # (A.2e)
    R33_mon = P33_mie / P11_mie                # (A.2f)
    R43_mon = P43_mie / P11_mie                # (A.2g)
    Ray_11 = .75*(1 + np.cos(theta)**2)        # (A.3a)
    Ray_21 = -.75 * np.sin(theta)**2           # (A.3b)
    Polar_Ray = -Ray_21 / Ray_11               # (A.3c)
    m = complex(nr, ni)                        # (A.3d)
    M0 = abs((m**2 - 1)/(m**2 + 2))            # (A.3e)

    # A.2.3. Coherent scattering and optical depth
    #----------------------------------------------
    DIST = R0 * Xm                                               # (A.4)
    # Total coherent scattering
    Fc = []
    for tt in theta:
        # (A.5) # WARNING: sinc(x) = sin(pi.x)/(pi.x)
        Fi = np.sinc(2 * DIST * np.sin(tt/2.) / np.pi)
        Fc.append(np.sum(np.multiply(Fi, F0)) * (N**2 - N) + N)  # (A.6)

    tau_coef = np.sum(np.divide(F0, DIST**2)) * (N-1)/(4*np.pi)  # (A.7a)

    taue_out = tau_coef * Cext_mon                               # (A.7b)
    taus_out = tau_coef * Csca_mon                               # (A.7c)
    taua_out = tau_coef * Cabs_mon                               # (A.7d)

    # A.2.4: Absorption (step by step correction)
    #---------------------------------------------
    # (A.8 + A.9)
    Cabs = Cabs_mon * N * np.exp(-taua_out)
    corr_abs = 1 + (  C_abs_m_1 * M0**E_abs_m_1
                    + C_abs_m_2 * M0**E_abs_m_2 * np.sin(C_abs_x_1 * Xm)) \
        * np.exp(-C_abs_x_2 * Xm)                                # (A.10a)

    # NOTE: corr_abs is apply later cf A.2.9

    # A.2.5: Single-scattering approximation
    #----------------------------------------
    P22 = Fc * Ymon * np.exp(-taue_out)   # (A.11a)
    P33 = P22 * R33_mon                   # (A.11b)
    P43 = P22 * R43_mon                   # (A.11c)

    # A.2.6. Empirical correction for multiple scattering within aggregate
    #----------------------------------------------------------------------
    depol = C_p11_m_1 * M0**2 / np.power(N-1, 2/3.) * (1 + Polar_Ray) # (A.12a)
    depol_ll = C_p11_m_2 * M0 * np.power(taus_out, E_p11_t_1)         # (A.12b)

    if Xm <= 1.6:
        for ii in range(len(depol)):
            if depol[ii] < depol_ll:
                depol[ii] = depol_ll

    depol = depol * P22[nang-1] * (1-depol[nang-1])  # (A.12c)
    P11 = P22 + depol                                # (A.12d)
    P44 = P33 + depol * (2/np.pi * theta - 1)        # (A.12e)

    # A.2.7. Linear Polarizartion
    #-----------------------------
    Mpol = 1 - C_p21_m_1 * M0**2 / np.sqrt(N-1) \
             - C_p21_m_2 * M0*Xm * E_p21_m_1    \
        * np.exp(-C_p21_ta * taua_out)          \
        * (N-1) * E_p21_n_1                     \
        * np.exp(C_p21_ts * taus_out)  # (A.13a)

    polar_agg = Polar_mon * Mpol       # (A.13b)
    P21 = - P11 * polar_agg            # (A.13c)

    # A.2.8. Scattering cross section
    #---------------------------------
    Csca = .5 * np.trapz(P11*np.sin(theta), x=theta)  # (A.14a)

    P11_out = P11 / Csca
    P22_out = P22 / Csca
    P33_out = P33 / Csca
    P44_out = P44 / Csca
    P21_out = P21 / Csca
    P43_out = P43 / Csca

    corr_sca = 1. + C_sca_m_3                      \
        * (M0-C_sca_m_4) * np.sin(C_sca_x_1 * Xm)  \
        * np.exp(-Xm * C_sca_x_2)                     # (A.14c)

    # A.2.9. Efficiencies
    #---------------------
    Qs_out = Csca / (np.pi * Xm**2 * np.power(N, 2/3.)) * corr_sca  # (A.15a)
    Qa_out = Cabs / (np.pi * Xm**2 * np.power(N, 2/3.)) * corr_abs  # (A.10b + A.15b)
    Qe_out = Qs_out + Qa_out                                       # (A.15c)

    return Qs_out, Qa_out, Qe_out, P11_out, P22_out, P33_out, P44_out, P21_out, P43_out


def fractals(wvln, nr, ni, rm, Df, N, nang=NANG, force=FORCE):
    """
    Compute fractals cross-sections and phase function based
        on Tomasko 2008

    Input
    ------
        wvln    (float) Wavelength (m)
        nr      (float) Particule real optical index
        ni      (float) Particule real imaginary index
        rm      (float) Monomer radius (m)
        Df      (float) Fractal dimension
        N       (int)   Number of monomers
        [nang]  (int)   Number of angles for the phase function
                         (range from 0 to pi/2)
        [force] (bool)  Bypass validity checks

    Outputs
    --------
        qsct  (float)    Scattering cross section (m^-2)
        qext  (float)    Extinction cross section (m^-2)
        qabs  (float)    Absorption cross section (m^-2)
        gg    (float)    Asymmetry parameter
        theta (float[])  Phase function angles (radians)
        P     (float[])  Phase function
    """
    Xm = 2. * np.pi * rm / wvln

    Qs, Qa, Qe, P, _, _, _, _, _ = fractals_tomasko_2008(Df, N, Xm, nr, ni, nang, force)
    qsct = Qs * np.pi * rm**2 * np.power(N, 2/3.)
    qext = Qe * np.pi * rm**2 * np.power(N, 2/3.)
    qabs = Qa * np.pi * rm**2 * np.power(N, 2/3.)
    theta = np.linspace(0, np.pi, len(P))
    gg = None  # <- Not calculated

    return qsct, qext, qabs, gg, theta, P
