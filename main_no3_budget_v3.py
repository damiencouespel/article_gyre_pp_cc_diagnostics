"""
Compute no3 budget in a box
Save outputs as pckl files
"""
from IPython import embed
import time

# import functions_GYRE as fGYRE
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('./my_matplotlibrc.mplstyle')
import pickle

from FGYRE import reading, interpolating, averaging

##################################################
#                PARAMETETERS                    #
##################################################


fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
start_date  = '0166-01-01'
end_date    = '0170-12-31'
time_span=(start_date, end_date)
res = {'CTL1':'1', 'CTL1_NOGM':'1', 'CTL9':'9',  'CTL27':'27', 
       'CC1':'1', 'CC1_NOGM':'1', 'CC9':'9', 'CC27':'27',
       'CTL1_KL2000':'1', 'CTL1_KL500':'1', 'CTL1_KGM500_KL500':'1', 'CTL1_KGM2000_KL2000':'1',
       'CC1_KL2000':'1', 'CC1_KL500':'1', 'CC1_KGM500_KL500':'1', 'CC1_KGM2000_KL2000':'1'}
# res = {'CTL1':'05', 'CTL1_NOGM':'05', 'CTL9':'05', 'CTL27':'05', 
#        'CC1' :'05', 'CC1_NOGM' :'05', 'CC9' :'05', 'CC27' :'05'}
# res = {'CTL1':'1', 'CTL1_NOGM':'1', 'CTL9':'1', 'CTL27':'1', 
#        'CC1' :'1', 'CC1_NOGM' :'1', 'CC9' :'1', 'CC27' :'1'}

print(time_span)

#___________________
# define the box where to compute the budget
lon1, lon2 = None, None
lat1, lat2 = 35, 45

#___________________
# which simulation
sim = ['CC1_KL500']

#___________________
# set outputs files
outd = '/gpfswork/rech/dyk/rdyk004/MY_PYTHON3/PCKL/'
outf = 'main_no3_budget_v3_cc1_kl500_35N_1ymean.pckl'

#___________________
# input files
fsufU   = '_1y_01010101_01701230_grid_U.xml'
fsufV   = '_1y_01010101_01701230_grid_V.xml'
fsufW   = '_1y_01010101_01701230_grid_W.xml'
fsufTr  = '_1y_01010101_01701230_ptrc_T.xml'
fsuftrd = '_1y_01010101_01701230_no3trd.xml'
fsufbio = '_1y_01010101_01701230_diad_T.xml'

output = {'fsufU'   : fsufU  ,
          'fsufV'   : fsufV  ,
          'fsufW'   : fsufW  ,
          'fsufTr'  : fsufTr ,
          'fsuftrd' : fsuftrd,
          'fsufbio' : fsufbio, 
          'start_date' : start_date, 
          'end_date' : end_date, 
          'lon1' : lon1, 'lon2' : lon2, 'lat1':lat1, 'lat2':lat2} 

##################################################
#             LOOP ON SIMULATIONS                #
##################################################

for vsim in sim : 
    
    output[vsim]={}
    
    zwmeshfile = '/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+res[vsim]+'.nc' 
    # zwmeshfile = '/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' # TEST
    zwmesh = reading.read_mesh(zwmeshfile)
    
    output[vsim]['mesh'] = zwmesh

    def cumint(zw, ze3t) :
        # note :  zw is vertical only, idem ze3t
        print("cumulative integral")
        zw = np.cumsum(zw * ze3t) # positive values = flx towards surface
        zw = np.roll(zw, 1) # we shift because first value is 2nd W-point
        zw[0] = 0. # flux at the surface (W-point = 0) is zero
        return zw
    #

    def diag_trend(inp) : 
        zw = averaging.tmean(inp, zwmesh, dim='tzyx', grid='T')
        zw = averaging.ymean(zw, zwmesh, dim='zyx'  , grid='T', ymin=lat1, ymax=lat2, integral = True)
        zw = averaging.xmean(zw, zwmesh, dim='zx'   , grid='T', xmin=lon1, xmax=lon2, integral = True)
        zw = cumint(zw, zwmesh['e3t_0'])
        return zw
    #

    #---------------------------------------
    # COMPUTE OFFLINE ADVECTIVE FLUXES
    #---------------------------------------

    #___________________
    # read Tr
    zwfTr = fdir + vsim + fsufTr
    zw   = reading.read_ncdf('NO3', zwfTr, time=(start_date, end_date))
    zwTr = zw['data']

    #___________________
    # compute uNO3
    zwfU = fdir + vsim + fsufU
    zw  = reading.read_ncdf('vozocrtx', zwfU, time=(start_date, end_date))
    zwU = zw['data']
    uTr = zwU * .5 * (zwTr + np.roll(zwTr, -1, axis=-1)) * zwmesh['umask'][np.newaxis]
    del zwU
    uTr = averaging.tmean(uTr, zwmesh, dim='tzyx', grid='U')
    if lon1 != None : uTr1 = interpolating.xinterpol(uTr, zwmesh, lon1, dim='zyx', grid='U')
    else : uTr1 = uTr[:, :, 0]
    if lon2 != None : uTr2 = interpolating.xinterpol(uTr, zwmesh, lon2, dim='zyx', grid='U')
    else : uTr2 = uTr[:, :, -1]
    del uTr
    uTr1 = averaging.ymean(uTr1, zwmesh, dim='zy' , grid='U', ymin=lat1, ymax=lat2, integral = True)
    uTr2 = averaging.ymean(uTr2, zwmesh, dim='zy' , grid='U', ymin=lat1, ymax=lat2, integral = True)
    uTr1 = cumint(uTr1, zwmesh['e3t_0'])
    uTr2 = cumint(uTr2, zwmesh['e3t_0'])
    output[vsim]['uN1'] = uTr1
    output[vsim]['uN2'] = -uTr2

    #___________________
    # compute vNO3
    zwfV = fdir + vsim + fsufV
    zw  = reading.read_ncdf('vomecrty', zwfV, time=(start_date, end_date))
    zwV = zw['data']
    vTr = zwV * .5 * (zwTr + np.roll(zwTr, -1, axis=-2)) * zwmesh['vmask'][np.newaxis]
    del zwV
    vTr = averaging.tmean(vTr, zwmesh, dim='tzyx', grid='V')
    if lat1 != None : vTr1 = interpolating.yinterpol(vTr, zwmesh, lat1, dim='zyx', grid='V')
    else : vTr1 = vTr[:, 0, :]
    if lat2 != None : vTr2 = interpolating.yinterpol(vTr, zwmesh, lat2, dim='zyx', grid='V')
    else : vTr2 = vTr[:, -1, :]
    del vTr
    vTr1 = averaging.xmean(vTr1, zwmesh, dim='zx' , grid='V', xmin=lon1, xmax=lon2, integral = True)
    vTr2 = averaging.xmean(vTr2, zwmesh, dim='zx' , grid='V', xmin=lon1, xmax=lon2, integral = True)
    vTr1 = cumint(vTr1, zwmesh['e3t_0'])
    vTr2 = cumint(vTr2, zwmesh['e3t_0'])
    output[vsim]['vN1'] = vTr1
    output[vsim]['vN2'] = -vTr2

    #___________________
    # compute wNO3
    zwfW = fdir + vsim + fsufW
    zw  = reading.read_ncdf('vovecrtz', zwfW, time=(start_date, end_date))
    zwW = zw['data']
    wTr = np.roll(zwW, -1, axis=-3) * .5 * (zwTr + np.roll(zwTr, -1, axis=-3)) * zwmesh['wmask'][np.newaxis]
    wTr = np.roll(wTr, 1, axis=-3)
    wTr[:, 0, :, :] = zwW[:, 0, :, :] * .5 * zwTr[:, 0, :, :]
    del zwW
    wTr = averaging.tmean(wTr, zwmesh, dim='tzyx', grid='W')
    wTr = averaging.ymean(wTr, zwmesh, dim='zyx', grid='W', ymin=lat1, ymax=lat2, integral = True)
    wTr = averaging.xmean(wTr, zwmesh, dim='zx' , grid='W', xmin=lon1, xmax=lon2, integral = True)
    output[vsim]['wN'] = wTr

    sfx = uTr1 - uTr2 + vTr1 - vTr2 + wTr
    output[vsim]['sfx'] = sfx

    #---------------------------------------
    # COMPUTE OFFLINE  ADVECTIVE FLUXES GM
    #---------------------------------------

    if vsim in ['CTL1', 'CTL1_KL2000', 'CTL1_KL500', 'CTL1_KGM500_KL500', 'CTL1_KGM2000_KL2000',
                'CC1', 'CC1_KL2000', 'CC1_KL500', 'CC1_KGM500_KL500', 'CC1_KGM2000_KL2000'] :

        #___________________
        # compute ugmNO3
        zwfU = fdir + vsim + fsufU
        zw  = reading.read_ncdf('vozoeivx', zwfU, time=(start_date, end_date))
        zwU = zw['data']
        ugmTr = zwU * .5 * (zwTr + np.roll(zwTr, -1, axis=-1)) * zwmesh['umask'][np.newaxis]
        del zwU
        ugmTr = averaging.tmean(ugmTr, zwmesh, dim='tzyx', grid='U')
        if lon1 != None : ugmTr1 = interpolating.xinterpol(ugmTr, zwmesh, lon1, dim='zyx', grid='U')
        else : ugmTr1 = ugmTr[:, :, 0]
        if lon2 != None : ugmTr2 = interpolating.xinterpol(ugmTr, zwmesh, lon2, dim='zyx', grid='U')
        else : ugmTr2 = ugmTr[:, :, -1]
        del ugmTr
        ugmTr1 = averaging.ymean(ugmTr1, zwmesh, dim='zy' , grid='U', ymin=lat1, ymax=lat2, integral = True)
        ugmTr2 = averaging.ymean(ugmTr2, zwmesh, dim='zy' , grid='U', ymin=lat1, ymax=lat2, integral = True)
        ugmTr1 = cumint(ugmTr1, zwmesh['e3t_0'])
        ugmTr2 = cumint(ugmTr2, zwmesh['e3t_0'])
        output[vsim]['ugmN1'] = ugmTr1
        output[vsim]['ugmN2'] = -ugmTr2

        #___________________
        # compute vgmNO3
        zwfV = fdir + vsim + fsufV
        zw  = reading.read_ncdf('vomeeivy', zwfV, time=(start_date, end_date))
        zwV = zw['data']
        vgmTr = zwV * .5 * (zwTr + np.roll(zwTr, -1, axis=-2)) * zwmesh['vmask'][np.newaxis]
        del zwV
        vgmTr = averaging.tmean(vgmTr, zwmesh, dim='tzyx', grid='V')
        if lat1 != None : vgmTr1 = interpolating.yinterpol(vgmTr, zwmesh, lat1, dim='zyx', grid='V')
        else : vgmTr1 = vgmTr[:, 0, :]
        if lat2 != None : vgmTr2 = interpolating.yinterpol(vgmTr, zwmesh, lat2, dim='zyx', grid='V')
        else : vgmTr2 = vgmTr[:, -1, :]
        del vgmTr
        vgmTr1 = averaging.xmean(vgmTr1, zwmesh, dim='zx' , grid='V', xmin=lon1, xmax=lon2, integral = True)
        vgmTr2 = averaging.xmean(vgmTr2, zwmesh, dim='zx' , grid='V', xmin=lon1, xmax=lon2, integral = True)
        vgmTr1 = cumint(vgmTr1, zwmesh['e3t_0'])
        vgmTr2 = cumint(vgmTr2, zwmesh['e3t_0'])
        output[vsim]['vgmN1'] = vgmTr1
        output[vsim]['vgmN2'] = -vgmTr2
        
        #___________________
        # compute wgmNO3
        zwfW = fdir + vsim + fsufW
        zw  = reading.read_ncdf('voveeivz', zwfW, time=(start_date, end_date))
        zwW = zw['data']
        wgmTr = np.roll(zwW, -1, axis=-3) * .5 * (zwTr + np.roll(zwTr, -1, axis=-3)) * zwmesh['wmask'][np.newaxis]
        wgmTr = np.roll(wgmTr, 1, axis=-3)
        wgmTr[:, 0, :, :] = zwW[:, 0, :, :] * .5 * zwTr[:, 0, :, :]
        del zwW
        wgmTr = averaging.tmean(wgmTr, zwmesh, dim='tzyx', grid='W')
        wgmTr = averaging.ymean(wgmTr, zwmesh, dim='zyx', grid='W', ymin=lat1, ymax=lat2, integral = True)
        wgmTr = averaging.xmean(wgmTr, zwmesh, dim='zx' , grid='W', xmin=lon1, xmax=lon2, integral = True)
        output[vsim]['wgmN'] = wgmTr

        sfxgm = ugmTr1 - ugmTr2 + vgmTr1 - vgmTr2 + wgmTr
        output[vsim]['sfxgm'] = sfxgm
    #

    #---------------------------------------
    # COMPUTE OFFLINE ADVECTIVE TREND
    #---------------------------------------

    #____________________
    # compute udNO3
    zwfU = fdir + vsim + fsufU
    zw  = reading.read_ncdf('vozocrtx', zwfU, time=(start_date, end_date))
    zwU = zw['data'] * zwmesh['umask'][np.newaxis, :, :, :]
    zw = zwmesh['e1u'][np.newaxis, np.newaxis, :, :]
    # Trd(i) = 0.5 * ( u(i) * [ T(i+1) - T(i) ] + u(i-1) * [ T(i) -T(i-1) ] )
    udTr = 0.5 * zwU * ( np.roll(zwTr, -1, axis=-1) - zwTr ) / zw + \
        0.5 * np.roll(zwU, 1, axis=-1) * ( zwTr - np.roll(zwTr, 1, axis=-1) ) / np.roll(zw, 1, axis = -1)
    udTr = udTr * zwmesh['tmask'][np.newaxis, :, :, :]
    del zwU
    udTr = diag_trend(udTr)
    output[vsim]['xadoff'] = - udTr # minus sign so that influx is positive
    del udTr
    
    #____________________
    # compute vdNO3
    zwfV = fdir + vsim + fsufV
    zw  = reading.read_ncdf('vomecrty', zwfV, time=(start_date, end_date))
    zwV = zw['data'] * zwmesh['vmask'][np.newaxis, :, :, :]
    zw = zwmesh['e2v'][np.newaxis, np.newaxis, :, :]
    # Trd(j) = 0.5 * ( v(j) * [ T(j+1) - T(j) ] + v(j-1) * [ T(j) -T(j-1) ] )
    vdTr = 0.5 * zwV * ( np.roll(zwTr, -1, axis=-2) - zwTr ) / zw + \
        0.5 * np.roll(zwV, 1, axis=-2) * ( zwTr - np.roll(zwTr, 1, axis=-2) ) / np.roll(zw, 1, axis = -2)
    vdTr = vdTr * zwmesh['tmask'][np.newaxis, :, :, :]
    del zwV
    vdTr = diag_trend(vdTr)
    output[vsim]['yadoff'] = - vdTr # minus sign so that influx is positive
    del vdTr
    
    #____________________
    # compute wdNO3
    zwfW = fdir + vsim + fsufW
    zw  = reading.read_ncdf('vovecrtz', zwfW, time=(start_date, end_date))
    zwW = zw['data'] * zwmesh['wmask'][np.newaxis, :, :, :]
    zw = zwmesh['e3w'][np.newaxis, :, :, :]
    # Trd(k) = 0.5 * ( w(k) * [ T(k-1) - T(k) ] + w(k+1) * [ T(k) -T(k+1) ] )
    wdTr = 0.5 * zwW * ( np.roll(zwTr, 1, axis=-3) - zwTr ) / zw + \
        0.5 * np.roll(zwW, -1, axis=-3) * ( zwTr - np.roll(zwTr, -1, axis=-3) ) / np.roll(zw, -1, axis = -3)
    wdTr[:, -1, :, :] = 0
    wdTr = wdTr * zwmesh['tmask'][np.newaxis, :, :, :]
    del zwW
    wdTr = diag_trend(wdTr)
    output[vsim]['zadoff'] = - wdTr # minus sign so that influx is positive
    del wdTr
    
    # sum of off adv trd
    output[vsim]['hadoff'] = output[vsim]['xadoff'] + output[vsim]['yadoff']
    output[vsim]['advoff'] = output[vsim]['hadoff'] + output[vsim]['zadoff']

    #---------------------------------------
    # COMPUTE OFFLINE ADVECTIVE TREND GM
    #---------------------------------------

    if vsim in ['CTL1', 'CTL1_KL2000', 'CTL1_KL500', 'CTL1_KGM500_KL500', 'CTL1_KGM2000_KL2000',
                'CC1', 'CC1_KL2000', 'CC1_KL500', 'CC1_KGM500_KL500', 'CC1_KGM2000_KL2000'] : 

        #____________________
        # compute ugmdNO3
        zwfU = fdir + vsim + fsufU
        zw  = reading.read_ncdf('vozoeivx', zwfU, time=(start_date, end_date))
        zwU = zw['data'] * zwmesh['umask'][np.newaxis, :, :, :]
        zw = zwmesh['e1u'][np.newaxis, np.newaxis, :, :]
        # Trd(i) = 0.5 * ( u(i) * [ T(i+1) - T(i) ] + u(i-1) * [ T(i) -T(i-1) ] )
        ugmdTr = 0.5 * zwU * ( np.roll(zwTr, -1, axis=-1) - zwTr ) / zw + \
            0.5 * np.roll(zwU, 1, axis=-1) * ( zwTr - np.roll(zwTr, 1, axis=-1) ) / np.roll(zw, 1, axis = -1)
        ugmdTr = ugmdTr * zwmesh['tmask'][np.newaxis, :, :, :]
        del zwU
        ugmdTr = diag_trend(ugmdTr)
        output[vsim]['xadgmoff'] = - ugmdTr # minus sign so that influx is positive
        del ugmdTr
        
        #____________________
        # compute vgmdNO3
        zwfV = fdir + vsim + fsufV
        zw  = reading.read_ncdf('vomeeivy', zwfV, time=(start_date, end_date))
        zwV = zw['data'] * zwmesh['vmask'][np.newaxis, :, :, :]
        zw = zwmesh['e2v'][np.newaxis, np.newaxis, :, :]
        # Trd(j) = 0.5 * ( v(j) * [ T(j+1) - T(j) ] + v(j-1) * [ T(j) -T(j-1) ] )
        vgmdTr = 0.5 * zwV * ( np.roll(zwTr, -1, axis=-2) - zwTr ) / zw + \
            0.5 * np.roll(zwV, 1, axis=-2) * ( zwTr - np.roll(zwTr, 1, axis=-2) ) / np.roll(zw, 1, axis = -2)
        vgmdTr = vgmdTr * zwmesh['tmask'][np.newaxis, :, :, :]
        del zwV
        vgmdTr = diag_trend(vgmdTr)
        output[vsim]['yadgmoff'] = - vgmdTr  # minus sign so that influx is positive
        del vgmdTr
       
        #____________________
        # compute wgmdNO3
        zwfW = fdir + vsim + fsufW
        zw  = reading.read_ncdf('voveeivz', zwfW, time=(start_date, end_date))
        zwW = zw['data'] * zwmesh['wmask'][np.newaxis, :, :, :]
        zw = zwmesh['e3w'][np.newaxis, :, :, :]
        # Trd(k) = 0.5 * ( w(k) * [ T(k-1) - T(k) ] + w(k+1) * [ T(k) -T(k+1) ] )
        wgmdTr = 0.5 * zwW * ( np.roll(zwTr, 1, axis=-3) - zwTr ) / zw + \
            0.5 * np.roll(zwW, -1, axis=-3) * ( zwTr - np.roll(zwTr, -1, axis=-3) ) / np.roll(zw, -1, axis = -3)
        wgmdTr[:, -1, :, :] = 0
        wgmdTr = wgmdTr * zwmesh['tmask'][np.newaxis, :, :, :]
        del zwW
        wgmdTr = diag_trend(wgmdTr)
        output[vsim]['zadgmoff'] = - wgmdTr # minus sign so that influx is positive
        del wgmdTr

        # sum of trend
        output[vsim]['hadgmoff'] = output[vsim]['xadgmoff'] + output[vsim]['yadgmoff']
        output[vsim]['advgmoff'] = output[vsim]['hadgmoff'] + output[vsim]['zadgmoff']
    #

    #---------------------------------------
    # DIAG ON ONLINE FLUXES/TRENDS
    #---------------------------------------
    
    #____________________
    # Vert diffusive flux

    zwftrd = fdir + vsim + fsuftrd
    zw=reading.read_ncdf('ZDF_NO3', zwftrd, time = (start_date, end_date))
    zdf = diag_trend(zw['data'])
    output[vsim]['zdf'] = zdf

    #____________________
    # Lat diffusive flux

    zwftrd = fdir + vsim + fsuftrd
    zw=reading.read_ncdf('LDF_NO3', zwftrd, time = (start_date, end_date))
    ldf = diag_trend(zw['data'])
    output[vsim]['ldf'] = ldf
    
    #____________________
    # Online advective trends

    zwftrd = fdir + vsim + fsuftrd

    zw=reading.read_ncdf('XAD_NO3', zwftrd, time = (start_date, end_date))
    xad = zw['data']
    zw=reading.read_ncdf('YAD_NO3', zwftrd, time = (start_date, end_date))
    yad = zw['data']
    zw=reading.read_ncdf('ZAD_NO3', zwftrd, time = (start_date, end_date))
    zad = zw['data']
    
    xad = diag_trend(xad)
    yad = diag_trend(yad)
    zad = diag_trend(zad)
    output[vsim]['xadonl'] = xad
    output[vsim]['yadonl'] = yad
    output[vsim]['zadonl'] = zad
    output[vsim]['hadonl'] = xad + yad
    output[vsim]['advonl'] = xad + yad + zad
    
    output[vsim]['dynonl'] = xad + yad + zad + zdf + ldf
    
    #____________________
    # Bio

    zwftrd = fdir + vsim + fsuftrd
    zw = reading.read_ncdf('SMS_NO3', zwftrd, time = (start_date, end_date))
    sms = diag_trend(zw['data'])
    output[vsim]['sms'] = sms
    
    zwfbio = fdir + vsim + fsufbio
    zw=reading.read_ncdf('FNO3PHY', zwfbio, time = (start_date, end_date))
    nwp=zw['data']/3600./24.
    nwp = diag_trend(nwp)
    output[vsim]['nwp'] = nwp
    zw=reading.read_ncdf('FNH4NO3', zwfbio, time = (start_date, end_date))
    nit=zw['data']/3600./24.
    nit = diag_trend(nit)
    output[vsim]['nit'] = nit

    output[vsim]['exp'] = sms - (nit - nwp)
    output[vsim]['src'] = sms + nwp
    
##################################################
#             END LOOP ON SIMULATIONS            #
##################################################

f = open(outd+outf, 'wb')
pickle.dump(output, f)
f.close()
