"""
Functions to compute density, brunt-vaisala frequency...
"""

# from IPython import embed
# import time
import numpy        as np

####################################################################################################
def eos_insitu(temp, sali, alpha=2.0e-4, beta=7.7e-4, **kwargs) :
    """
    Compute the in situ density from potential temperature and
    salinity using a bilinear equation of state : prd = beta * sali -
    alpha * temp
    
    Usage : eos_insitu(temp, sali, alpha=2.0e-4, beta=7.7e-4)
    
    Parameters
    ----------
    temp : 3D or 4D numpy array
       potential temperature in [Celsius] 
    sali  : 3D or 4D numpy array
       salinity in [psu]
    alpha : float
       temperature coefficient in [1/Celsius]
    beta  : float
       salinity coefficient in [1/psu]
    **kwargs : not passed

    Returns
    -------
    prd : 3D or 4D numpy array containin insitu density [no unit]
    """

    print("> eos_instu")
    prd   = beta  * sali - alpha * temp
    return prd
# END eos_insitu
####################################################################################################

####################################################################################################
def eos_insitu_pot(temp, sali, rau0=1035.0, **kwargs) :
    """
    Compute the in situ potential density from potential temperature
    and salinity using a bilinear equation of state :
    prho   = rau0 * ( 1.0 +  prd)
    prd is th in situ density : prd = beta  * sali - alpha * temp

    Usage : eos_insitu_pot(temp, sali, rau0 = 1035.0)

    Parameters
    ----------
    temp  : 3D or 4D numpy array
       potential temperature in [Celsius] 
    sali  : 3D or 4D numpy array
       salinity in [psu]
    rau0  : float
       ref density in [kg/m3]
    **kwargs : passed to eos_insitu
    can be alpha=float or beta=float

    Returns
    -------
    prho : 3D or 4D numpy array containing insitu potential density in
    [kg/m3]
    """

    print("> eos_instu_pot")
    prd  = eos_insitu(temp, sali, **kwargs)
    prho = (1+ prd ) * rau0
    return prho
# END eos_insit_pot
####################################################################################################

####################################################################################################
def eos_bn2(temp, sali, e3w, dim='tzyx', alpha=2.0e-4, beta=7.7e-4) :
    """
    Compute the local Brunt-Vaisala frequency using a linear equation
    of state : 
    grav  = 9.80665     # gravity [m/s2]
    pn2 = grav * (alpha * dk[ temp ] - beta * dk[ sali ] ) / e3w 
    NB : N^2 is set to zero at the first level
    and is never used at this level.

    Usage: eos_bn2(temp, sali, e3w, dim='tzyx', alpha=2.0e-4,
    beta=7.7e-4)

    Parameters
    ----------
    temp : 3D or 4D numpy array
       potential temperature in [Celsius] 
    sal  : 3D or 4D numpy array
       salinity in [psu]
    e3w  : 3D numpy array
       vertical aspect ratio between T point
    dim  : string
       specify the dimension of the fields
    alpha : float
       temperature coefficient in [1/Celsius]
    beta  : float
       salinity coefficient in [1/psu]
    
    Returns
    -------
    pn2 : Brun-Vaisala frequency in [s-2]
    """
    print("> eos_bn2")
    grav = 9.80665     # gravity [m/s2]
    if dim != 'tzyx':
        temp = temp[np.newaxis] 
        sali = sali[np.newaxis] 
    #
    e3w = e3w[np.newaxis]
    pn2 = np.zeros_like(temp) + np.float('nan')
    pn2[:, 1:] = grav * ( alpha * (temp[:, :-1] - temp[:, 1:]) \
                              - beta  * (sali[:, :-1] - sali[:, 1:]) ) / e3w[:, 1:]
    return pn2
# END eos_bn2
####################################################################################################
