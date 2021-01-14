# This module contains functions to estimate a few meteorological parameters.
#
# List of functions: - Theta      -  Potential temperature
#                    - Thetav     -  Virtual potential temperature
#                    - Thetae     -  Equivalent potential temperature
#                    - TD         -  Dew point
#                    - RW         -  Water mass miwing ratio
#                    - Q          -  Specific humidity
#                    - RH         -  Relative humidity
#                    - RHI        -  Relative humidity with respect to ice
#                    - ES         -  Saturation water vapour pressure
#                    - ESI        -  Saturation water vapour pressure w.r.t. ice
#                    - LV         -  Latent heat of vaporization
#                    - NU         -  Air kinematic viscosity 
#
# Author : Virginie Guemas - 2020
################################################################################
import numpy as np
import sys

cp = 1004   # J.K-1.kg-1
cpv= 1860.0 # J.K-1.kg-1
Ra = 287    # J.K-1.kg-1
Rv = 461.4  # J.K-1.kg-1
################################################################################
def Theta(z=None,T=None,P=None,q=None) :
   """ 
   This function computes the potential temperature in Kelvin as a function of:
   - the temperature T in Kelvin,
   - the height z in m,
   or as a function of :
   - the temperature T in Kelvin, 
   - the pressure P in hPa, 
   - optionally, the specific humidity in kg/kg.
   """
   
   if T is not None:
     check_T(T)
   else:
     sys.exit('Whichever method you choose, you need to provide T')

   if z is not None:
     gamma = 0.00975 # K.m-1 (as in Peterson and Renfrew, QJRMS, 2009)
                     # Fairall et al (1996) use 0.0098
                     # Smith (1988) uses 0.01
     theta = T + gamma*z  # Approximation with substantial errors for flights
   elif P is not None:
     if q is None:
       kappa = Ra/cp
     else:
       check_q(q)
       kappa = (Ra*(1-q)+Rv*q)/(cp*(1-q)+cpv*q) 
     theta = T * (1000/P)**kappa
   else:
     sys.exit('theta can be computed from (T,z), from (T,P) or from (T,P,q)')

   return theta
################################################################################
def Thetav(theta,q) :
   """ 
   This function computes the virtual potential temperature in Kelvin as a function of the potential temperature in Kelvin and specific humidity in kg/kg.
   It can also be used to compute the virtual temperature in Kelvin as a function of the temperature in Kelvin and specific humidity in kg/kg.
   """
  
   check_T(theta)
   check_q(q)

   delta = (Rv-Ra)/Ra

   thetav = theta * (1 + delta*q)   # From PV = (ma Ra + mv Rv) T 

   return thetav
################################################################################
def Thetae(T,q,P) :
   """
   This function computes the equivalent potential temperature in Kelvin as a function of the temperature in Kelvin, the specific humidity in kg/kg and the pressure in hPa.
   """

   check_T(T)
   check_q(q)

   thetae = Theta(T = T + LV(T)/cp*q, P = P)        

   return thetae
################################################################################
def TD(T,case=2,rh=None,q=None,P=None) :
   """ 
   This function computes the dew point temperature in Kelvin from the temperature in Kelvin and relative humidity in % or specific humidity in kg/kg and pressure in hPa.
   """

   check_T(T)

   if q is None and rh is None:
     sys.exit('At least one of q or rh should be provided')
   
   if rh is None:
     if P is None:
       sys.exit('If rh is not provided, both q and P should be')
     else:
       check_q(q)
       rh = RH(q,T,P)  
   else:
     check_rh(rh)
     print('The relative humidity provided is used')

   TC = T - 273.15

   if case == 1:
     a = 6.112 # Bolton, 1980, MWR
     b = 17.67 # Bolton, 1980, MWR
     c = 243.5 # Bolton, 1980, MWR
   elif case == 2:
     a = 6.1121 # Buck, 1981, JAMC
     b = np.where(TC>0,17.368,17.966) # Buck, 1981, JAMC
     c = np.where(TC>0,238.88,247.15) # Buck, 1981, JAMC
   elif case == 3:
     a = 6.1078 # meteoblue glossary
     b = np.where(TC<0,17.84362,17.08085) # meteoblue glossary
     c = np.where(TC<0,245.425,234.175) # meteoblue glossary
   else:
     sys.exit('Unknow case')

   gamma = np.log(rh/100) + b*TC/(c+TC)

   Tdew = (c * gamma)/(b - gamma) + 273.15

   return Tdew
################################################################################
def RW(q):
   """
   This function computes the water mass mixing ratio (in kg/kg) as a function of the specific humidity (in kg/kg)
   """

   r = q/(1-q)

   return r
################################################################################
def Q(P,rh=None,Td=None,T=None,case=2):
   """
   This function computes the specific humidity (in kg/kg) as a function of :
   - the dew point temperature Td (in Kelvin), 
   - the pressure P (in hPa)
   or as a function of :
   - the relative humidity rh (in %),
   - the temperature (in Kelvin),
   - the pressure (in hPa). 
   """
 
   if Td is not None :
     check_T(Td)
     e = ES(Td,case=2)
   elif rh is not None and T is not None:
     check_T(T)
     check_rh(rh)
     e = rh * ES(T) / 100
   else:
     sys.exit('q can be computed from (Td,P) or from (rh,T,P)')

   if case == 1:
     q = (e*Ra)/(P*Rv) / (1 - e*(Rv-Ra)/(Rv*P))
   elif case == 2:
     q = (e*Ra)/(P*Rv) # Assumption that R humid air = R dry air
   else:
     sys.exit('Unknown case')
 
   return q
################################################################################
def RH(T,Td=None,q=None,P=None,case=2) :
   """ 
   This function computes the relative humidity in % as a function of:
   - the specific humidity in kg/kg, 
   - the temperature in Kelvin,
   - the pressure in hPa
   or as a function of:
   - the dew point temperature in Kelvin,
   - the temperature in Kelvin.
   """

   check_T(T)
   
   if Td is not None:
     check_T(Td)
     e = ES(Td)
   elif q is not None and P is not None:
     check_q(q)
     TC = T - 273.15 

     if case == 1:
       e = q * P / ((Rv-Ra)/Rv*q + (Ra/Rv))
     elif case == 2:
       e = q * P / (Ra/Rv)  # Assumption that R humid air = R dry air
     else:
       sys.exit('Unknow case')
   
   else:
     sys.exit('RH can be computed from (T,Td) or from (T,Q,P)')

   rh = e / ES(T) * 100

   return rh
################################################################################
def RHI(q,T,P):
   """
   This function computes the relative humidity with respect to ice in % from the specific humidity in kg/kg, the temperature in Kelvin, and the pressure in hPa.
   """

   check_q(q)
   check_T(T)

   e = q * P / (Ra/Rv)

   rhi = e / ESI(T) * 100

   return rhi
################################################################################
def ES(T,case=2):
   """ This function computes the saturation vapour pressure in hPa as a function of the temperature in Kelvin.
   """

   check_T(T)

   TC = T - 273.15

   # August-Roche-Magnus (or Magnus-Tetens or Magnus)
   if case == 1:
     es = 6.112 * np.exp((17.67 * TC)/(TC + 243.5))
   elif case == 2:
     es = 6.1094 * np.exp((17.625 * TC)/(TC + 243.04))
   # Tetens - lower than the others
   elif case == 3:
     es = 6.1078 * np.exp((17.27 * TC)/(TC + 237.3))
   elif case == 4:
   # Buck
     es = 6.1121 * np.exp ((18.678 - TC/234.5)*(TC/(257.14+TC)))
   elif case == 5: 
   # NASA GISS
     es = 6.108 * np.exp(2500000 * (0.00000793252 - 0.002166847/T)) 
   else:
     sys.exit('Unknow case')

   return es
################################################################################
def ESI(T) :
   """
   This function computes the saturation vapour pressure with respect to ice in hPa as a function of the temperature in Kelvin.
   """

   check_T(T)

   TC = T - 273.15

   # NASA GISS
   es = 6.108 * np.exp(2834000 * (0.00000793252 - 0.002166847/T))

   # Tetens
   es = 6.108 * np.exp((21.875 * TC)/(TC+265.5))

   return es
################################################################################
def LV(T,case=1) :
   """
   This function computes the latent heat of vaporisation/condensation in J.kg-1 as a function of temperature in Kelvin.
   """

   if case == 1:
     check_T(T)
     TC = T - 273.15
     Lv = 2500.8 - 2.36*TC + 0.0016*TC**2 - 0.00006*TC**3

   elif case == 2:
     Lv = 2500.0

   elif case == 3: # Fleagle & Businger (1980) p113 used in Fairall et al (1996)
     check_T(T)
     TC = T - 273.15
     Lv = (25 - 0.022274*TC)*10**5

   else:
     sys.exit('Unknown case')

   return Lv*1000
################################################################################
def NU(T,case=2):
   """
   This function computes the air kinematic viscosity as a function of temperature (in Kelvin). 
   """
  
   check_T(T)

   if case == 1: # Fit from Pierre Bouteloup
     nu = - 1.1555*10**(-14)*T**3 + 9.5728*10**(-11)*T**2 + 3.7604*10**(-8)*T -3.4484*10**(-6)

   elif case == 2: # Fit from Andreas (1989) used in Fairal et al (1996)
     TC = T - 273.15
     nu = 1.326*10**(-5)*(1 + 6.542*10**(-3)*TC + 8.301*10**(-6)*TC**2 - 4.840*10**(-9)*T**3)

   else:
     sys.exit('Unknown case')

   return nu
################################################################################
def check_q(q):
   """
   This function checks whether the specific humidity is most probably specified in g/kg or kg/kg. It exits in case the unit is g/kg with an error message.
   """

   if np.max(q) > 0.1 :
     sys.exit('Specific humidity q should be provided in kg/kg')

################################################################################
def check_T(T):
    """
    This function checks whether the temperature is most probably provided in Kelvin or Celsius degrees. It exits in case the unit is Celsius degrees with an error message.
    """

    if (np.min(T)) < 100 :
      sys.exit('Temperature should be provided in Kelvin degrees')

################################################################################
def check_rh(rh):
    """
    This function checks whether the relative humidity is most probably expressed in % or in fraction. It exits in case the unit is fractional with an error message.
    """

    if (np.max(rh)) < 0.1 :
      sys.exit('The relative humidity rh is expressed in %. 0 < rh < 100')

################################################################################
