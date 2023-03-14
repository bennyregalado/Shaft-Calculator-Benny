# -*- coding: utf-8 -*-
"""Final Shaft Calculation Code

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1vgRQxfS8QrLWF8j1OWqEihTDwx4B_2yB
"""

import numpy as np
def shaft(M,T,sut):
    # All parameters, inputs and outputs are in English Units
    # Input:
    # M: Maximum Moment on the shaft [lbf*in]
    # T: Maximum Torque on the shaft [lbf*in]
    # sut: Ultimate Tensile Strength of the material chosen  [psi]
    # Output:
    # D: Bigger diameter 
    # d: smaller diameter
    # r: filet radius 
    # kt: stress concentration factor, using Table 7-1 an assumption is made
    # kts: stress concentration factor for torsion, using Table 7-1 an assumption is made
    # Kf: stress concentration due to fatigue from Moment
    # Kfs: stress concentration due to fatigue from Torsion
    """
    Assumptions:
    a: factor for calculating ka, for machined/cold-drawn material 
    b: factor for calculating ka, for machined/cold-drawn material 
    n: factor of safety is 2
    D/d = 1.1
    r/d = .05
    kd: temperature modification factor, assuming that the maximum temperature
    is experienced for shaft (~482 F or 250 C)
    ke: reliability factor, assuming 99% reliability
    kf: miscellaneous-effects modification factor, assumption of 1, is only added 
    in special cases. It comes from experience of design engineer.
    """
    # Note: In the absence of information regarding 𝑘𝑑, 𝑘𝑒 and 𝑘𝑓, we approximate 
    # them to be equal to 1.
    # ka: surface condition modification factor 
    # kb: size modification factor, assuming the diameter will be 0.11 <= d <= 2
    # kc: load modification factor, assuming bending
    # kd: temperature modification factor, assuming that the maximum temperature 
    # is experienced for shaft (~482 F or 250 C)
    # ke: reliability factor, assuming 99% reliability 
    # kf: miscellaneous-effects modification factor, assumption of 1, is only added 
    # in special cases. It comes from experience of design engineer.
    # se_p: rotary-beam test specimen endurance limit, assuming sut <= 200 kpsi
    # Fatigue stress concentration calculation factors

    n = 2
    a = 2.7
    b = -.265
    ka = a*(sut*10**-3)**b
    do = 1 # This will be the initial diameter that will be used as a guess
    kb = 0.897*(do)**(-0.107)
    kc = 1
    TF = 482
    kd = 0.975+.432*(10**-3)*TF-0.115*(10**-5)*TF**2+0.104*(10**-8)*TF**3-0.595*(10**-12)*TF**4
    ke = 0.814
    kf = 1
    se_p = .5*sut*10**-3
    # Calculation for the Endurance Limit
    se = ka*kb*kc*kd*ke*kf*se_p 
    kt = 1.7 # Used as an assumption 
    kts = 1.5 # Used as an assumption
    Kf = kt
    Kfs = kts
    # Below is from the Distortion Energy Goodman Equation 
    sigma_ap = Kf*(32*M*10**-3/np.pi)
    sigma_pp = np.sqrt(3)*Kfs*(32*T*10**-3/np.pi)
    do = ((sigma_ap/(se)+sigma_pp/(sut))*n)**(1/3)
    ro=.05*do
    # While loop for final calculations
    z = 1-(1/1.1)
    sqrt_a = 0.246-(3.08*10**-3)*(sut*10**-3)+(1.51*10**-5)*(sut*10**-3)**2-(2.67*10**-8)*(sut*10**-3)**3
    sqrt_as = 0.19-(2.51*10**-3)*(sut*10**-3)+(1.35*10**-5)*(sut*10**-3)**2-(2.67*10**-8)*(sut*10**-3)**3
    error = 1
    error_1 = .01
    while error > error_1:
      y = ro/do
      x = .05/y
      rn = y*do
      # Bending Stress Concentration Factor
      c1 = 0.947+1.206*np.sqrt(x)-0.131*x
      c2 = 0.022-3.405*np.sqrt(x)+0.915*x
      c3 = 0.869-1.777*np.sqrt(x)-0.555*x
      c4 = 0.810-0.422*np.sqrt(x)-0.260*x
      kt = c1+c2*(z)+c3*(z)**2+c4*(z)**3
      # Torsion Stress Concentration Factor
      c11 = 0.905+0.783*np.sqrt(x)-0.075*x
      c22 = -0.437-1.969*np.sqrt(x)+0.553*x
      c33 = 1.557+1.073*np.sqrt(x)-0.578*x
      c44 = -1.061+0.171*np.sqrt(x)+.086*x 
      kts = c11+c22*(z)+c33*(z)**2+c44*(z)**3
      # Calculate fatigue stress concentration factors
      q = 1/(1+(sqrt_a/np.sqrt(rn)))
      qs = 1/(1+(sqrt_as/np.sqrt(rn)))
      # Calculate fatigue stress concentration
      Kf = 1+q*(kt-1)
      Kfs = 1+qs*(kts-1)
      kb = 0.897*(do)**(-0.107)
      # Calculation for the Endurance Limit
      se = ka*kb*kc*kd*ke*kf*se_p
      # Calculate new diameter using the Distortion Energy Goodman Theory
      sigma_ap = Kf*(32*M*10**-3/np.pi)
      sigma_mp = np.sqrt(3)*Kfs*(32*T*10**-3/np.pi)
      dn = ((sigma_ap/(se)+sigma_mp/(sut))*n)**(1/3)
      error = np.abs((do-dn)/dn)
      if error > error_1:
        break 

M = 1260 # Example 7-1 if n = 2
T =  1100 # Example 7-1 if n = 2
sut = 108*10**3 # [psi] # Example 7-1
shaft(M,T,sut)
