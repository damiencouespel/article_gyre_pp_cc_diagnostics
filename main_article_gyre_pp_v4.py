from IPython import embed
import time
import pickle

# import functions_GYRE as fGYRE
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('./matplotlibrc_nature_cc.mplstyle')

from FGYRE import reading, interpolating, averaging, \
    plotting, stream_functions


plot_fig1 = False
plot_fig2 = False
plot_fig3 = False
plot_fig4 = True
plot_fig5 = False
plot_fig6 = False
save_values_no3budget = False

unc = 1. # number of std for uncertainty

##################################################
# FIG1
##################################################
if plot_fig1 :

    savefig='fig1.png'

    ####################
    # BSF CTL1
    ####################

    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsufV  = '_1y_01010101_01701230_grid_V.xml'
    start_date  = '0166-01-01'
    end_date    = '0170-12-31'
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc')
    mapBSF={'mesh':zwmesh}
    for vkR in kR : 
        fff = fdir + 'CTL' + vkR + fsufV
        zw  = reading.read_ncdf('vomecrty', fff, time=(start_date, end_date))
        zw = zw['data']
        bsf = stream_functions.bsfv(zw, zwmesh)
        mapBSF['CTL'+vkR] = np.mean(bsf, axis = 0)
    #
    # AVG1
    zw = []
    for vkR in kR : zw.append(mapBSF['CTL'+vkR])
    zw = np.array(zw)
    mapBSF['AVG1'] = np.ma.array(data=np.nanmean(zw, axis=0), mask=mapBSF['CTL1'].mask)
    
    ####################
    # TEMPERATURE SPINUP, CTL, CC
    ####################

    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsufT = '_1y_00010101_01701230_grid_T.xml'
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1']
    kS = ['CTL', 'CC']

    timeT = {}
    for vkR in kR :
        if vkR in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        timeT['R'+vkR]={}
        for vkS in kS :
            inty = 10
            yyy1 = np.arange(1, 170, inty)
            timeT['R'+vkR][vkS] = []
            for zwy in yyy1 : 
                start_date  = '0'+str(zwy)+'-01-01'
                end_date    = '0'+str(zwy+inty-1)+'-12-31'
                fff = fdir + vkS + vkR + fsufT
                # fff = fdir + vkS + '1' + fsufT # TEST
                zw  = reading.read_ncdf('votemper', fff, time=(start_date, end_date))
                zw = zw['data']
                zw = averaging.zmean(zw, zwmesh, dim='tzyx', zmax = 700.)
                zw = averaging.ymean(zw, zwmesh, dim='tyx')
                zw = averaging.xmean(zw, zwmesh, dim='tx')
                timeT['R'+vkR][vkS].extend(list(zw))
            #
        #
    #

    # AVG1, STD1
    timeT['AVG1'], timeT['STD1'], timeT['+STD'], timeT['-STD'] = {}, {}, {}, {}
    for vkS in ['CTL', 'CC'] : 
        zw = []
        kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
        for vkR in kR : zw.append(timeT['R'+vkR][vkS])
        zw = np.array(zw)
        timeT['AVG1'][vkS] = np.nanmean(zw, axis=0)
        timeT['STD1'][vkS] = unc*np.nanstd(zw, axis=0)
        timeT['+STD'][vkS] = timeT['AVG1'][vkS] + timeT['STD1'][vkS]
        timeT['-STD'][vkS] = timeT['AVG1'][vkS] - timeT['STD1'][vkS]
    #
    
    ####################
    # FORCINGS 
    ####################

    forcings = {}
    zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )

    # ATM TEMPERATURE
    zlat = zwmesh['latT'][:, 2]
    zemp_S    = 0.7       # intensity of COS in the South
    zemp_N    = 0.8       # intensity of COS in the North
    zemp_sais = 0.1       #
    zTstar    = 25.6      # intemsity from 25.6 degC at 20 degN !!DC2
    tstarmax = (zTstar+25.6/15.2) * np.cos( np.pi*(zlat-20.) / (33.8*(1+8.3/33.8)*2.) ) 
    tstarmin = (zTstar-25.6/15.2) * np.cos( np.pi*(zlat-20.) / (33.8*(1-8.3/33.8)*2.) ) 
    warming = 0.04*70
    forcings['atmT'] = {'lat':zlat, 'summer':tstarmax, 'winter': tstarmin, 'warming':warming}

    # WIND STRESS
    zlat = zwmesh['latU'][:, 2]
    ztau = 0.105 # mean intensity
    ztau_sais = 0.015 # seasonal oscillation intensity
    ztaunmax = ztau - ztau_sais
    ztaunmin = ztau + ztau_sais
    zjm1 = 1  - 1.
    zphiv1 = 20.
    rad = np.pi / 180 # conversion from degre into radian
    ra  = 6371229.   # earth radius [m]
    zphiv2 = 20 + 30  * 106000/(ra*rad)
    # sinusoidal zonal windstress : - tau0 * COS( 2 * pi * lat / L )
    utaumax = - ztaunmax * np.cos( 2.*np.pi*(zlat-zphiv1) / (zphiv2-zphiv1) )
    utaumin = - ztaunmin * np.cos( 2.*np.pi*(zlat-zphiv1) / (zphiv2-zphiv1) )
    forcings['wind'] = {'lat':zlat, 'summer':utaumax, 'winter':utaumin}
    
    ####################
    # PLOT
    ####################

    infact  = 1/2.54
    ncol, nrow = 3, 2
    fsize = (ncol*5*infact, nrow*5*infact) #(width, height)
    fig, ax   = plt.subplots(nrow, ncol, figsize=fsize, sharey='row') # nrow, ncol

    # BSF
    zax = ax[0, 0]
    X, Y, Z = mapBSF['mesh']['lonV'], mapBSF['mesh']['latV'], mapBSF['AVG1']
    aa = np.linspace(-30, 0, 6)
    bb = np.linspace(0, 30, 6)
    levBSF = np.concatenate((aa[:-1], bb))
    cmapBSF = 'RdBu'
    cfbsf = zax.contourf(X, Y, Z, levels=levBSF, cmap=cmapBSF, extend='both')
    plotting.make_XY_axis(zax, title='(a) CTL1 bar. circ.')

    # WIND
    zax = ax[0, 1]
    wsum = zax.plot(forcings['wind']['summer'], forcings['wind']['lat'], c='orange'     , ls='-', lw=1.5)
    wwin = zax.plot(forcings['wind']['winter'], forcings['wind']['lat'], c='darkturquoise', ls='-', lw=1.5)
    zax.set_title('(b) Wind stress [N/m$^2$]', loc='left')
    zax.tick_params(left=False)
    zax.set_frame_on(False)
    (xmin, xmax) = zax.xaxis.get_view_interval()
    (ymin, ymax) = zax.yaxis.get_view_interval()
    zax.add_artist(plt.Line2D((xmin, xmax), (ymin, ymin), color = 'dimgrey', linewidth=1))
    zax.add_artist(plt.Line2D((xmin, xmax), (ymax, ymax), color = 'dimgrey', linewidth=1))
    x0, y0 = forcings['wind']['summer'][27], forcings['wind']['lat'][27]
    zax.annotate('summer', xy=(x0, y0), xytext=(x0+.01, y0+.2), color='orange')
    x0, y0 = forcings['wind']['winter'][24], forcings['wind']['lat'][24]
    zax.annotate('winter', xy=(x0, y0), xytext=(x0-.08, y0-.4), color='darkturquoise')
    
    # ATM TEMP
    zax = ax[0, 2]
    tsum   = zax.plot(forcings['atmT']['summer'], forcings['atmT']['lat'], c='orange'       , ls='--', dashes=[1.5, .7], lw=1.5)
    twin   = zax.plot(forcings['atmT']['winter'], forcings['atmT']['lat'], c='darkturquoise', ls='--', dashes=[1.5, .7], lw=1.5)
    tsumcc = zax.plot(forcings['atmT']['summer']+forcings['atmT']['warming'], \
                      forcings['atmT']['lat'], c='orange'     , ls='-', lw=1.5)
    twincc = zax.plot(forcings['atmT']['winter']+forcings['atmT']['warming'], \
                      forcings['atmT']['lat'], c='darkturquoise', ls='-', lw=1.5)
    zax.set_title('(c) Pseudo atm. temp. [°C]', loc='left')
    zax.yaxis.set_ticks_position('right')
    zax.yaxis.set_label_position('right')
    zax.set_ylabel('Northward km')
    zax.tick_params(axis='y', length=0)
    zax.set_frame_on(False)
    (xmin, xmax) = zax.xaxis.get_view_interval()
    (ymin, ymax) = zax.yaxis.get_view_interval()
    zax.add_artist(plt.Line2D((xmin, xmax), (ymin, ymin), color = 'dimgrey'))
    zax.add_artist(plt.Line2D((xmin, xmax), (ymax, ymax), color = 'dimgrey'))
    x0, y0 = forcings['atmT']['summer'][21], forcings['atmT']['lat'][21]
    zax.arrow(x0, y0, 0.04*70, 0, length_includes_head=True, width=.1, head_width=1., zorder=30, head_length=1.2)
    zax.annotate('Warming\n in CC simul.', xy=(x0+1, y0+1))
    
    # SPINUP CTL CC
    zax = plt.subplot2grid((nrow, ncol), (1, 0), colspan=3, fig=fig)
    X = np.arange(170)+1
    lines, labels = {}, {}
    ooo=5
    # UNCERTAINTY
    Yu, Yd = timeT['+STD']['CTL'], timeT['-STD']['CTL']
    zax.fill_between(X, Yd, Yu, color='k', alpha=0.1, zorder=ooo)
    ooo+=5
    Yu, Yd = timeT['+STD']['CC'], timeT['-STD']['CC']
    zax.fill_between(X, Yd, Yu, color='k', alpha=0.1, zorder=ooo)
    ooo+=5
    www={'R9':2.5, 'R27':1.5, \
         'AVG1':3.5, '+STD':1, '-STD':1}
    names={'R9':'1/9°', 'R27':'1/27°', \
           'AVG1':'1° avg $\pm$ st. dev.', '+STD':'1°, avg+std', '-STD':'1°, avg-std'}

    for vk in ['AVG1', 'R9', 'R27'] :
        ztim=timeT[vk]
        ll1, = zax.plot(X, ztim['CTL'], c='k', lw=www[vk], zorder=ooo, ls='--', dashes=[1.5,.7])
        ll2, = zax.plot(X, ztim['CC'] , c='k', lw=www[vk], zorder=ooo, ls='-')
        lines[vk]=[ll1, ll2]
        ooo+=5
    #
    zax.set_title('(d) Ocean temperature 0-700 m [°C]', loc='left')
    zax.set_xlim(1, 171)
    zax.xaxis.set_ticks([1, 100, 171])
    zax.xaxis.set_ticklabels(['-100 years', '0', '70 years'])
    zax.locator_params(axis='y', nbins=6)
    zw = zax.get_position()
    zax.set_position([zw.x0, zw.y0-.13, zw.width, 1.2*zw.height])
    hdl=(lines['AVG1'][1], lines['R9'][1], lines['R27'][1])
    lab=(names['AVG1']   , names['R9']   , names['R27']   )
    # hdl = (lines['AVG1'][1], lines['AVG1'][1], lines['AVG1'][1])
    # lab = (names['AVG1']   , names['AVG1']   , names['AVG1']   )
    leg = zax.legend(hdl, lab, loc='best')
    leg._legend_box.align = "left"
    x0, y0 = 75, timeT['AVG1']['CTL'][75]
    zax.annotate('Spin-up', xy=(x0, y0), xytext=(x0+5, y0+.4), \
                 bbox=dict(boxstyle="round", fc="1"), \
                 arrowprops=dict(arrowstyle="-", shrinkA=0))
    x0, y0 = 145, timeT['AVG1']['CC'][145]
    zax.annotate('CC', xy=(x0, y0), xytext=(x0+6, y0-.4), \
                 bbox=dict(boxstyle="round", fc="1"), \
                 arrowprops=dict(arrowstyle="-", shrinkA=0))
    x0, y0 = 157, timeT['AVG1']['CTL'][157]
    zax.annotate('CTL', xy=(x0, y0), xytext=(x0+5, y0+.35), \
                 bbox=dict(boxstyle="round", fc="1"), \
                 arrowprops=dict(arrowstyle="-", shrinkA=0))
 
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    print("Figure saved: ", dirfig+suffig+savefig)
    fig.savefig(dirfig+suffig+savefig)
    plt.close()
#
##################################################
# END FIG1
##################################################

##################################################
# FIG2
##################################################
if plot_fig2 : 

    savefig='fig2.png'
    
    ####################
    # VORTICITY
    ####################

    plotVRT={}
    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsufU  = '_2d_01660101_01701230_grid_U.xml'
    fsufV  = '_2d_01660101_01701230_grid_V.xml'
    start_date  = '0170-01-01'
    end_date    = '0170-01-02'
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    for vkR in kR : 
        if vkR in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        plotVRT['CTL'+vkR]={'mesh':zwmesh}
        fff = fdir + 'CTL' + vkR + fsufU
        # fff = fdir + 'CTL1' + fsufU  # TEST
        zw  = reading.read_ncdf('vozocrtx', fff, time=(start_date, end_date))
        zw = zw['data'].squeeze()
        zwU = zw[0].squeeze()
        fff = fdir + 'CTL' + vkR + fsufV
        # fff = fdir + 'CTL1' + fsufV  # TEST
        zw  = reading.read_ncdf('vomecrty', fff, time=(start_date, end_date))
        zw = zw['data'].squeeze()
        zwV = zw[0].squeeze()
        f = 2* 7.292115083046062e-005 # earth vrticity
        zw = ( np.roll(zwV, 1, axis = 1) - zwV ) / zwmesh['e1f'] - ( np.roll(zwU, 1, axis=0) - zwU ) / zwmesh['e2f']
        plotVRT['CTL'+vkR]['data'] = zw / f * 100
    #
    # AVG1
    plotVRT['AVG1']={'mesh':plotVRT['CTL1']['mesh']}
    zw = []
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    for vkR in kR : zw.append(plotVRT['CTL'+vkR]['data'])
    zw = np.array(zw)
    plotVRT['AVG1']['data'] = np.ma.array(data=np.nanmean(zw, axis=0), mask=plotVRT['CTL1']['data'].mask)
    
    ####################
    # MAP PP
    ####################

    plotNPP={}
    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_pp.xml'
    start_date  = '0166-01-01'
    end_date    = '0170-12-31'
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    for vkR in kR : 
        if vkR in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        plotNPP['CTL'+vkR]={'mesh':zwmesh}
        fff = fdir + 'CTL' + vkR + fsuf
        # fff = fdir + 'CTL1' + fsuf  # TEST
        zw  = reading.read_ncdf('PP', fff, time=(start_date, end_date))
        zw = zw['data']
        zw = averaging.tmean(zw, zwmesh) 
        zw = averaging.zmean(zw, zwmesh, dim='zyx', integral=True) 
        plotNPP['CTL'+vkR]['data'] = zw
    #
    # AVG1
    plotNPP['AVG1']={'mesh':plotNPP['CTL1']['mesh']}
    zw = []
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    for vkR in kR : zw.append(plotNPP['CTL'+vkR]['data'])
    zw = np.array(zw)
    plotNPP['AVG1']['data'] = np.ma.array(data=np.nanmean(zw, axis=0), mask=plotNPP['CTL1']['data'].mask)

    ####################
    # SECTION NO3
    ####################

    plotNO3={}
    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_ptrc_T.xml'
    start_date  = '0166-01-01'
    end_date    = '0170-12-31'
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    for vkR in kR : 
        if vkR in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        plotNO3['CTL'+vkR]={'mesh':zwmesh}
        fff = fdir + 'CTL' + vkR + fsuf
        # fff = fdir + 'CTL1' + fsuf  # TEST
        zw  = reading.read_ncdf('NO3', fff, time=(start_date, end_date))
        zw = zw['data']
        zw = averaging.tmean(zw, zwmesh) 
        zw = averaging.xmean(zw, zwmesh, dim='zyx') 
        plotNO3['CTL'+vkR]['data'] = zw
        #
    #
    # AVG1
    plotNO3['AVG1']={'mesh':plotNO3['CTL1']['mesh']}
    zw = []
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    for vkR in kR : zw.append(plotNO3['CTL'+vkR]['data'])
    zw = np.array(zw)
    plotNO3['AVG1']['data'] = np.ma.array(data=np.nanmean(zw, axis=0), mask=plotNO3['CTL1']['data'].mask)
    
    ####################
    # PLOT
    ####################

    infact  = 1/2.54
    ncol, nrow = 3, 3
    fsize = (ncol*5*infact, nrow*5*infact) #(width, height)
    fig, ax   = plt.subplots(nrow, ncol, figsize=fsize, sharey='row') # nrow, ncol

    #___________________
    # VORTICITY 
    aa   = np.linspace(-20, 0, 6)
    bb   = np.linspace(0, 20, 6)
    levVRT  = np.concatenate((aa[:-1], bb))
    cmapVRT = 'RdBu'
    # CTL1
    zax = ax[0, 0]
    X, Y, Z = plotVRT['AVG1']['mesh']['lonF'], plotVRT['AVG1']['mesh']['latF'], plotVRT['AVG1']['data']
    cfvrt1  = zax.contourf(X, Y, Z, levels=levVRT, cmap=cmapVRT, extend='both')
    plotting.make_XY_axis(zax, title='(a) Surf. vort. CTL1')
    # CTL9
    zax = ax[0, 1]
    X, Y, Z = plotVRT['CTL9']['mesh']['lonF'], plotVRT['CTL9']['mesh']['latF'], plotVRT['CTL9']['data']
    cfvrt9  = zax.contourf(X, Y, Z, levels=levVRT, cmap=cmapVRT, extend='both')
    plotting.make_XY_axis(zax, title='(b) Surf. vort. CTL9')
    # CTL27
    zax = ax[0, 2]
    X, Y, Z = plotVRT['CTL27']['mesh']['lonF'], plotVRT['CTL27']['mesh']['latF'], plotVRT['CTL27']['data']
    cfvrt27 = zax.contourf(X, Y, Z, levels=levVRT, cmap=cmapVRT, extend='both')
    plotting.make_XY_axis(zax, title='(c) Surf. vort. CTL27')

    for zax in ax[0, :].flatten() : zax.label_outer()
    for zax in ax[0, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[0, 1:].flatten() : zax.tick_params(left=False)    
    
    cbtitle = '[% of the planetary vort.]'
    zw = ax[0, -1].get_position()
    cbar_ax = fig.add_axes([zw.x0+zw.width+.02, zw.y0, .0158, zw.height])
    fig.colorbar(cfvrt1, cax=cbar_ax, orientation='vertical', ticklocation='right', label=cbtitle)

    #___________________
    # NPP
    levNPP = np.linspace(0, 4., 11)
    cmapNPP = 'viridis'
    # CTL1
    zax = ax[1, 0]
    X, Y, Z = plotNPP['AVG1']['mesh']['lonT'], plotNPP['AVG1']['mesh']['latT'], plotNPP['AVG1']['data']
    cfnpp1  = zax.contourf(X, Y, Z, levels=levNPP, cmap=cmapNPP, extend='max')
    plotting.make_XY_axis(zax, title='(d) Vert. int. NPP CTL1')
    # CTL9
    zax = ax[1, 1]
    X, Y, Z = plotNPP['CTL9']['mesh']['lonT'], plotNPP['CTL9']['mesh']['latT'], plotNPP['CTL9']['data']
    cfnpp9  = zax.contourf(X, Y, Z, levels=levNPP, cmap=cmapNPP, extend='max')
    plotting.make_XY_axis(zax, title='(e) Vert. int. NPP CTL9')
    # CTL27
    zax = ax[1, 2]
    X, Y, Z = plotNPP['CTL27']['mesh']['lonT'], plotNPP['CTL27']['mesh']['latT'], plotNPP['CTL27']['data']
    cfnpp27 = zax.contourf(X, Y, Z, levels=levNPP, cmap=cmapNPP, extend='max')
    plotting.make_XY_axis(zax, title='(f) Vert. int. NPP CTL27')
    
    for zax in ax[1, 1:].flatten() :
        zax.tick_params(left=False)
        zax.set_ylabel('')
    #
    
    cbtitle = '[mmol-N/m$^2$/d]'
    zw = ax[1, -1].get_position()
    cbar_ax = fig.add_axes([zw.x0+zw.width+.02, zw.y0, .0158, zw.height])
    fig.colorbar(cfnpp1, cax=cbar_ax, orientation='vertical', ticklocation='right', label=cbtitle)

    #___________________
    # NO3
    levNO3 = np.linspace(0, 20., 11)
    cmapNO3 = 'viridis'
    # CTL1
    zax = ax[2, 0]
    X, Y, Z = plotNO3['AVG1']['mesh']['latT'][:, 5], plotNO3['AVG1']['mesh']['depT'], plotNO3['AVG1']['data']
    cfno31  = zax.contourf(X, Y, Z, levels=levNO3, cmap=cmapNO3, extend='max')
    plotting.make_YZ_axis(zax, depmax=800, title='(g) NO3 CTL1')
    # CTL9
    zax = ax[2, 1]
    X, Y, Z = plotNO3['CTL9']['mesh']['latT'][:, 5], plotNO3['CTL9']['mesh']['depT'], plotNO3['CTL9']['data']
    cfno39  = zax.contourf(X, Y, Z, levels=levNO3, cmap=cmapNO3, extend='max')
    plotting.make_YZ_axis(zax, depmax=800, title='(h) NO3 CTL9')
    # CTL27
    zax = ax[2, 2]
    X, Y, Z = plotNO3['CTL27']['mesh']['latT'][:, 5], plotNO3['CTL27']['mesh']['depT'], plotNO3['CTL27']['data']
    cfno327 = zax.contourf(X, Y, Z, levels=levNO3, cmap=cmapNO3, extend='max')
    plotting.make_YZ_axis(zax, depmax=800, title='(i) NO3 CTL27')
    
    for zax in ax[2, :].flatten() :
        zw = zax.get_position()
        zax.set_position([zw.x0, zw.y0, zw.width, .8*zw.height])
        zax.locator_params(axis='y', nbins=5)
    #
    for zax in ax[2, 1:].flatten() :
        zax.tick_params(left=False)
        zax.set_ylabel('')
    #
    
    cbtitle = '[mmol-N/m$^3$]'
    zw = ax[2, -1].get_position()
    cbar_ax = fig.add_axes([zw.x0+zw.width+.02, zw.y0, .0158, zw.height])
    fig.colorbar(cfno31, cax=cbar_ax, orientation='vertical', ticklocation='right', label=cbtitle)

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    print("Figure saved: ", dirfig+suffig+savefig)
    fig.savefig(dirfig+suffig+savefig)
    plt.close()
##################################################
# END FIG2
##################################################


##################################################
# FIG3
##################################################
if plot_fig3 : 
    
    ####################
    # PARAM
    ####################

    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsufGYRE  = '_1y_01010101_01701230_pp.xml'
    var = 'PP'
    savefig='fig3.png'
    zmin, zmax, zint = None, None, True
    ymin, ymax = 35, 45
    sdatemap, edatemap  = '0166-01-01', '0170-12-31'
    sdatetim, edatetim  = '0101-01-01', '0170-12-31'
    unc = 1.

    ####################
    # PREPARE DATA
    ####################
    
    #___________________
    # MAP
    
    mapGYRE={}

    # R1S
    zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_'+'R1'+'.nc')
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    kS = ['CTL', 'CC']
    for vkR in kR : 
        mapGYRE['R'+vkR]={'mesh':zwmesh}
        for vkS in kS :  
            fff = fdir + vkS + vkR + fsufGYRE
            zw  = reading.read_ncdf(var, fff, time=(sdatemap, edatemap))
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx')
            zw = averaging.zmean(zw, zwmesh, grid='T', dim='zyx', integral=zint, zmin=zmin, zmax=zmax)
            mapGYRE['R'+vkR][vkS] = zw
        #
        mapGYRE['R'+vkR]['DELTA'] = mapGYRE['R'+vkR]['CC'] - mapGYRE['R'+vkR]['CTL']
    #
    del zw

    # R9, R27
    kR = ['9', '27']
    kS = ['CTL', 'CC']
    for vkR in kR : 
        zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_'+'R1'+'.nc') # R1 TEST
        mapGYRE['R'+vkR]={'mesh':zwmesh}
        for vkS in kS :  
            fff = fdir + vkS + vkR + fsufGYRE
            # fff = fdir + vkS + '1' + fsufGYRE  # R1 TEST
            zw  = reading.read_ncdf(var, fff, time=(sdatemap, edatemap))
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx')
            zw = averaging.zmean(zw, zwmesh, grid='T', dim='zyx', integral=zint, zmin=zmin, zmax=zmax)
            mapGYRE['R'+vkR][vkS] = zw
        #
        mapGYRE['R'+vkR]['DELTA'] = mapGYRE['R'+vkR]['CC'] - mapGYRE['R'+vkR]['CTL']
    #
    del zw

    #___________________
    # TIMESERIE

    timGYRE = {}

    # R1S
    zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_'+'R1'+'.nc')
    kR = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    kS = ['CTL', 'CC']
    for vkR in kR : 
        timGYRE['R'+vkR]={}
        for vkS in kS :  
            fff = fdir + vkS + vkR+ fsufGYRE
            zw  = reading.read_ncdf(var, fff, time=(sdatetim, edatetim))
            zw = zw['data']
            zw = averaging.zmean(zw, zwmesh, grid='T', dim='tzyx', integral=zint, zmin=zmin, zmax=zmax)
            zw = averaging.ymean(zw, zwmesh, grid='T', dim='tyx', ymin=ymin, ymax=ymax) 
            zw = averaging.xmean(zw, zwmesh, grid='T', dim='tx')
            timGYRE['R'+vkR][vkS] = zw
        #
        timGYRE['R'+vkR]['DELTA'] = timGYRE['R'+vkR]['CC'] - timGYRE['R'+vkR]['CTL']
    #
    del zw

    # AVG1, STD1
    zw = []
    for vkR in kR : zw.append(timGYRE['R'+vkR]['DELTA'])
    zw = np.array(zw)
    timGYRE['AVG1'] = {'DELTA':np.nanmean(zw, axis=0)}
    timGYRE['STD1'] = {'DELTA':unc*np.nanstd(zw, axis=0)}
    timGYRE['+STD'] = {'DELTA':timGYRE['AVG1']['DELTA'] + timGYRE['STD1']['DELTA']}
    timGYRE['-STD'] = {'DELTA':timGYRE['AVG1']['DELTA'] - timGYRE['STD1']['DELTA']}

    # R9, R27
    kR = ['9', '27']
    kS = ['CTL', 'CC']
    for vkR in kR :
        zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+vkR+'.nc')
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_'+'R1'+'.nc') # R1 TEST
        timGYRE['R'+vkR]={}
        for vkS in kS :
            fff = fdir + vkS + vkR + fsufGYRE
            # fff = fdir + vkS + '1' + fsufGYRE # R1 TEST
            zw  = reading.read_ncdf(var, fff, time=(sdatetim, edatetim))
            zw = zw['data']
            zw = averaging.zmean(zw, zwmesh, grid='T', dim='tzyx', integral=zint, zmin=zmin, zmax=zmax)
            zw = averaging.ymean(zw, zwmesh, grid='T', dim='tyx', ymin=ymin, ymax=ymax)
            zw = averaging.xmean(zw, zwmesh, grid='T', dim='tx')
            timGYRE['R'+vkR][vkS] = zw
        #
        timGYRE['R'+vkR]['DELTA'] = timGYRE['R'+vkR]['CC'] - timGYRE['R'+vkR]['CTL']
    #
    del zw
    
    ####################
    # PLOT
    ####################

    infact  = 1/2.54
    ncol, nrow = 3, 3
    fsize = (ncol*5*infact, nrow*5*infact) #(width, height)
    fig, ax   = plt.subplots(nrow, ncol, figsize=fsize, sharex='col', sharey='row') # nrow, ncol


    unit = '[mmol-N/m$^2$/d]'

    # map colorbar
    cbtitle = 'Vert. int. NPP change'+unit
    levmax=1.5
    aa = np.linspace(-levmax, 0, 6)
    bb = np.linspace(0, levmax, 6)
    lev = np.concatenate((aa[:-1], bb))
    cmap = 'PiYG'
    ext = 'both'
    
    #___________________
    # GYRE MAP
    cf={}
    ttls={'R1'               :'k$_{gm}$=1e$^3$, k$_{redi}$=1e$^3$', \
          'R1_KL500'         :'k$_{gm}$=1e$^3$, k$_{redi}$=500'      , \
          'R1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{redi}$=2e$^3$', \
          'R1_KGM500_KL500'  :'k$_{gm}$=500, k$_{redi}$=500'      , \
          'R1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{redi}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmno')
    subnum.reverse()
    def plot_map(vk, zax) : 
        zmap=mapGYRE[vk]
        zwmesh=zmap['mesh']
        X, Y = zwmesh['lonT'], zwmesh['latT']
        Z = zmap['DELTA']
        cf[vk] = zax.contourf(X, Y, Z, cmap=cmap, levels=lev, extend=ext)
        (xmin, xmax) = zax.xaxis.get_view_interval()
        zax.add_artist(plt.Line2D((xmin, xmax), (35, 35), color = 'black', linestyle='--'))
        zax.add_artist(plt.Line2D((xmin, xmax), (45, 45), color = 'black', linestyle='--'))
        plotting.make_XY_axis(zax, title='('+subnum.pop()+') '+ttls[vk])
        zax.annotate(note[vk], xy=(-83, 23), \
                     bbox=dict(boxstyle="round", fc="1"))
        zax.label_outer()
    #
    plot_map('R1_KGM500_KL500', ax[0, 0])
    plot_map('R1', ax[1, 0])
    plot_map('R1_KGM2000_KL2000', ax[2, 0])
    plot_map('R1_KL500' , ax[0, 1])
    plot_map('R1_KL2000', ax[1, 1])
    plot_map('R9', ax[0, 2])
    plot_map('R27', ax[1, 2])
    ax[2, 1].axis('off')
    ax[2, 2].axis('off')
    
    for zax in ax[0, :].flatten()  : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    ax[1, 0].tick_params(bottom=False)
    for zax in ax[1, 1:].flatten() :
        zax.xaxis.set_ticks_position('bottom')
        zax.xaxis.set_label_position('bottom')
        zax.set_xlabel('Eastward km')
    #
    for zax in ax[:, -1].flatten() :
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
        zax.set_ylabel('Northward km')
    #

    cbar_ax = fig.add_axes([0.22, 0.93, 0.6, 0.023])
    fig.colorbar(cf['R9'], cax=cbar_ax, orientation='horizontal', ticklocation='top',\
                 label=cbtitle)

    #___________________
    # TIMESERIES

    # timeseries param
    names={'R9':'1/9°', 'R27':'1/27°', \
           'AVG1':'1° avg. $\pm$ st. dev.'}
    timtitle='('+subnum.pop()+') GSP mean NPP change'
    ylim=(None, None)
    ccc={'R1':'k', 'R9':'k', 'R27':'k', \
         'R1_KL500':'peachpuff', 'R1_KL2000':'paleturquoise', \
         'R1_KGM500_KL500':'palegreen', 'R1_KGM2000_KL2000':'lavender', \
         'AVG1':'k', '+STD':'lightgray', '-STD':'lightgray'}
    www={'R1':3.5, 'R9':2.5, 'R27':1.5, \
         'R1_KL500':3.5, 'R1_KL2000':3.5, \
         'R1_KGM500_KL500':3.5, 'R1_KGM2000_KL2000':3.5, \
         'AVG1':3.5, '+STD':1, '-STD':1}

    axtime = plt.subplot2grid((3, 3), (2, 1), colspan=2, fig=fig)
    X = np.arange(70)+1
    lines, labels = [], []
    ooo=5
    # UNCERTAINTY R1
    Yu, Yd = timGYRE['+STD']['DELTA'], timGYRE['-STD']['DELTA']
    axtime.fill_between(X, Yd, Yu, color='k', alpha=0.1, zorder=ooo)
    ooo+=5
    for vk in ['AVG1', 'R9', 'R27'] :
        ztim=timGYRE[vk]
        Y = ztim['DELTA']
        ll, = axtime.plot(X, Y, c=ccc[vk], lw=www[vk], zorder=ooo)
        lines.append(ll)
        labels.append(names[vk])
        ooo+=5
    #
    
    leg = axtime.legend(lines, labels, loc='best')
    leg._legend_box.align = "left"
    axtime.set_title(timtitle, loc='left')
    axtime.set_xlim(0, 70)
    axtime.xaxis.set_ticks([0, 35, 70])
    axtime.xaxis.set_ticklabels(['0', '35 years', '70 years'])
    axtime.set_ylim(ylim) 
    axtime.set_ylabel(unit)
    axtime.yaxis.set_ticks_position('right')
    axtime.yaxis.set_label_position('right')
    zw = axtime.get_position()
    axtime.set_position([zw.x0, zw.y0, zw.width, 0.8*zw.height])

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END FIG3
##################################################

##################################################
# FIG4
##################################################
if plot_fig4 : 

    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/dyk/rdyk004/MY_PYTHON3/PCKL/'

    data2process = {}
    
    #___________________
    # R1s TOTAL NO3 BUDGET

    for zr1 in r1s : 
        data2process['R'+zr1]={}
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CTL'] = open_pckl(zwin, 'CTL'+zr1)
        zwin = 'main_no3_budget_v3_cc' +zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CC'] = open_pckl(zwin, 'CC'+zr1)
    #
    
    #___________________
    # R1s MEAN NO3 BUDGET

    for zr1 in r1s : 
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N_1ymean.pckl'
        data2process['R'+zr1]['CTL_M'] = open_pckl(zwin, 'CTL'+zr1)
        zwin = 'main_no3_budget_v3_cc' +zr1.lower()+'_35N_1ymean.pckl'
        data2process['R'+zr1]['CC_M'] = open_pckl(zwin, 'CC'+zr1)
    #

    #___________________
    # R9 TOTAL NO3 BUDGET

    data2process['R9']={}
    zwin = 'main_no3_budget_v3_ctl9_35N.pckl'
    data2process['R9']['CTL'] = open_pckl(zwin, 'CTL9')
    zwin = 'main_no3_budget_v3_cc9_35N.pckl'
    data2process['R9']['CC'] = open_pckl(zwin, 'CC9')

    #___________________
    # R9 MEAN NO3 BUDGET

    zwin = 'main_no3_budget_v3_ctl9_35N_onR1_1ymean.pckl'
    data2process['R9']['CTL_M'] = open_pckl(zwin, 'CTL9')
    zwin = 'main_no3_budget_v3_cc9_35N_onR1_1ymean.pckl'
    data2process['R9']['CC_M'] = open_pckl(zwin, 'CC9')

    #___________________
    # R27 TOTAL NO3 BUDGET

    data2process['R27']={}
    
    for zsim in ['CTL', 'CC'] :
        data2process['R27'][zsim] = {}
        zwvars={}
        for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                    'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff', \
                    'xadonl', 'yadonl', 'zadonl', 'hadonl', 'advonl', \
                    'zdf', 'ldf', 'dynonl', \
                    'sms', 'nwp', 'nit', 'exp', 'src'] :
            zwvars[vvv]=[]
        #
        zwsuf = 'main_no3_budget_v3_'+zsim.lower()+'27_35N_y1'
        for yyy in np.arange(66, 71) :
            zwin = zwsuf+str(yyy)+'.pckl'
            zw = open_pckl(zwin, zsim+'27')
            if yyy==66 : data2process['R27'][zsim]['mesh']=zw['mesh']
            for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
        #
        for vvv in zwvars.keys() : data2process['R27'][zsim][vvv]=np.mean(zwvars[vvv], axis=0)
    #

    #___________________
    # R27 MEAN NO3 BUDGET

    zwin = 'main_no3_budget_v3_ctl27_35N_onR1_1ymean.pckl'
    data2process['R27']['CTL_M'] = open_pckl(zwin, 'CTL27')
    zwin = 'main_no3_budget_v3_cc27_35N_onR1_1ymean.pckl'
    data2process['R27']['CC_M'] = open_pckl(zwin, 'CC27')

    ####################
    # PREPARE NO3 BUDGET DATA
    ####################

    data2plot={}

    for zkres in data2process.keys() :
        data2plot[zkres]={}
        for zksim in ['CTL', 'CC'] :
            data2plot[zkres][zksim]={'depW':data2process[zkres][zksim]['mesh']['depW']}
            data2plot[zkres][zksim]['zmix']  = data2process[zkres][zksim]['zdf']
            data2plot[zkres][zksim]['advt']  = data2process[zkres][zksim]['vN1'] + \
                data2process[zkres][zksim]['wN']
            # add ldf to total advection
            data2plot[zkres][zksim]['advt'] += data2process[zkres][zksim]['ldf']
            # add advection by GM velocities to total advection
            if 'sfxgm' in data2process[zkres][zksim].keys() :
                data2plot[zkres][zksim]['advt'] += data2process[zkres][zksim]['vgmN1']
                data2plot[zkres][zksim]['advt'] += data2process[zkres][zksim]['wgmN']
            #
            #if zkres in ['R9','R27'] : # this condition is temporary because no mean R1 yet
            data2plot[zkres][zksim]['advm']  = data2process[zkres][zksim+'_M']['vN1'] + \
                data2process[zkres][zksim+'_M']['wN']
            data2plot[zkres][zksim]['adve'] = data2plot[zkres][zksim]['advt'] - \
                data2plot[zkres][zksim]['advm']
            #
            # if zkres[1:] in r1s : # temporary because no mean R1 yet
            #     data2plot[zkres][zksim]['adve'] = data2process[zkres][zksim]['ldf'] + \
            #         data2process[zkres][zksim]['vgmN1'] + data2process[zkres][zksim]['wgmN']
            #     data2plot[zkres][zksim]['advm'] = data2plot[zkres][zksim]['advt'] - data2plot[zkres][zksim]['adve']
            #
        #
        data2plot[zkres]['DELTA']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['DELTA']['zmix'] = data2plot[zkres]['CC']['zmix'] - \
            data2plot[zkres]['CTL']['zmix']
        data2plot[zkres]['DELTA']['advt'] = data2plot[zkres]['CC']['advt'] - \
            data2plot[zkres]['CTL']['advt']
        data2plot[zkres]['DELTA']['advm'] = data2plot[zkres]['CC']['advm'] - \
            data2plot[zkres]['CTL']['advm']
        data2plot[zkres]['DELTA']['adve'] = data2plot[zkres]['CC']['adve'] - \
            data2plot[zkres]['CTL']['adve']
    #
        
    ####################
    # READ AND PROCESS NO3
    ####################

    sdate = '0166-01-01'
    edate = '0170-12-31'
    fdirN = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          

    for zkres in ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', \
                '1_KGM2000_KL2000', '9', '27'] :
        for zksim in ['CTL', 'CC'] :
            fff = fdirN + zksim + zkres + '_1y_01010101_01701230_ptrc_T.xml'
            zw = reading.read_ncdf('NO3', fff, time = (sdate, edate))
            zwmesh = data2process['R'+zkres][zksim]['mesh']
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx')
            zw = averaging.ymean(zw, zwmesh, dim='zyx', ymin=35, ymax=45)
            zw = averaging.xmean(zw, zwmesh, dim='zx')
            data2plot['R'+zkres][zksim]['no3'] = zw
            data2plot['R'+zkres][zksim]['depT'] = zwmesh['depT']
        #
        data2plot['R'+zkres]['DELTA']['no3'] = data2plot['R'+zkres]['CC']['no3'] - \
            data2plot['R'+zkres]['CTL']['no3']
        data2plot['R'+zkres]['DELTA']['depT'] = zwmesh['depT']
    #

    ####################
    # AVG1, STD1
    ####################

    data2plot['AVG1'] = {}
    data2plot['STD1'] = {}
    for zksim in ['CTL', 'CC', 'DELTA'] :
        data2plot['AVG1'][zksim] = {'depW': data2plot['R'+r1s[0]][zksim]['depW'], \
                                    'depT': data2plot['R'+r1s[0]][zksim]['depT']}
        data2plot['STD1'][zksim] = {'depW': data2plot['R'+r1s[0]][zksim]['depW'], \
                                    'depT': data2plot['R'+r1s[0]][zksim]['depT']}
        for vvv in ['zmix', 'advt', 'no3', 'advm', 'adve'] : 
            zw = []
            for vkr1 in r1s : zw.append(data2plot['R'+vkr1][zksim][vvv])
            zw = np.array(zw)
            data2plot['AVG1'][zksim][vvv] = np.nanmean(zw, axis=0)
            data2plot['STD1'][zksim][vvv] = unc*np.nanstd(zw, axis=0)
        #
    #

    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    fsize   = (15*infact, 10*infact) #(width, height)
    fig, ax = plt.subplots(2, 4, sharey='row', figsize=fsize)
    
    def plotdata(zax, zres, zsim, ztitle) :
        zdat=data2plot[zres][zsim]
        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
        zlzmix, = zax.plot(zfact * zdat['zmix'], zdat['depW'], lw=1.5, c='coral')
        zladvt, = zax.plot(zfact * zdat['advt'], zdat['depW'], lw=1.5, c='dodgerblue')
        zladvm, = zax.plot(zfact * zdat['advm'], zdat['depW'], lw=1.5, c='dodgerblue', ls='--', dashes=[2., 1.])
        zladve, = zax.plot(zfact * zdat['adve'], zdat['depW'], lw=1.5, c='dodgerblue', ls='--', dashes=[1., 2.])
        if zres == 'AVG1' :
            zdatstd=data2plot['STD1'][zsim]
            Xm, Xp = zdat['zmix']-zdatstd['zmix'], zdat['zmix']+zdatstd['zmix']
            zax.fill_betweenx(zdat['depW'], zfact*Xm, zfact*Xp, alpha=0.2, color='coral')
            Xm, Xp = zdat['advt']-zdatstd['advt'], zdat['advt']+zdatstd['advt']
            zax.fill_betweenx(zdat['depW'], zfact*Xm, zfact*Xp, alpha=0.2, color='dodgerblue')
            Xm, Xp = zdat['advm']-zdatstd['advm'], zdat['advm']+zdatstd['advm']
            zax.fill_betweenx(zdat['depW'], zfact*Xm, zfact*Xp, alpha=0.2, color='dodgerblue')
            Xm, Xp = zdat['adve']-zdatstd['adve'], zdat['adve']+zdatstd['adve']
            zax.fill_betweenx(zdat['depW'], zfact*Xm, zfact*Xp, alpha=0.2, color='dodgerblue')
        #
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        zax.set_title(ztitle, loc='left')
        zax.vlines(0, 0, 400, color='grey')
        # zax.label_outer()
        zlines = [zlzmix, zladvt, zladvm, zladve]
        znames = ['Vert. mix.', 'Tot. trp.', 'Mean', 'Eddy']
        return (zlines, znames)
    #

    lw = {'AVG1':3.5, 'R9':2.5, 'R27':1.5}

    lln1 , = ax[0, 0].plot(data2plot['AVG1']['CTL']['no3'] , data2plot['AVG1']['CTL']['depT'], 'k', lw = lw['AVG1'])
    Xm = data2plot['AVG1']['CTL']['no3']-data2plot['STD1']['CTL']['no3']
    Xp = data2plot['AVG1']['CTL']['no3']+data2plot['STD1']['CTL']['no3']
    ax[0, 0].fill_betweenx(data2plot['AVG1']['CTL']['depT'], Xm, Xp, alpha=0.1, color='k')
    # for zr1 in r1s : ax[0, 0].plot(data2plot['R'+zr1]['CTL']['no3'],\
    #                                data2plot['R'+zr1]['CTL']['depT'], 'k', alpha=0.1, lw=1.5)
    lln9 , = ax[0, 0].plot(data2plot['R9']  ['CTL']['no3'] , data2plot['R9']  ['CTL']['depT'], 'k', lw = lw['R9']  )
    lln27, = ax[0, 0].plot(data2plot['R27'] ['CTL']['no3'] , data2plot['R27'] ['CTL']['depT'], 'k', lw = lw['R27'] )
    ax[0, 0].set_title('(a) NO3', loc='left')
    ax[0, 0].set_xlim((0, 25))
    ax[0, 0].xaxis.set_ticks([0, 10, 20])
    ax[0, 0].set_ylabel('Depth [m]')
    
    ll1  = plotdata(ax[0, 1], 'AVG1', 'CTL', '(b) CTL1')
    ll9  = plotdata(ax[0, 2], 'R9'  , 'CTL', '(c) CTL9')
    ll27 = plotdata(ax[0, 3], 'R27' , 'CTL', '(d) CTL27')
    for zax in ax[0, 1:].flatten() : 
        zax.set_xlim((-1.5, 3.5))
        zax.xaxis.set_ticks([-1.5, 0, 1.5, 3])
        zax.tick_params(left=False)
        if zax.is_last_col() :
            zax.set_ylabel('Depth [m]')
            zax.yaxis.set_ticks_position('right')
            zax.yaxis.set_label_position('right')
        #
    #
    
    lldn1 , = ax[1, 0].plot(data2plot['AVG1']['DELTA']['no3'] , data2plot['AVG1']['DELTA']['depT'], 'k', lw = lw['AVG1'])
    Xm = data2plot['AVG1']['DELTA']['no3']-data2plot['STD1']['DELTA']['no3']
    Xp = data2plot['AVG1']['DELTA']['no3']+data2plot['STD1']['DELTA']['no3']
    ax[1, 0].fill_betweenx(data2plot['AVG1']['DELTA']['depT'], Xm, Xp, alpha=0.1, color='k')
    # for zr1 in r1s : ax[1, 0].plot(data2plot['R'+zr1]['DELTA']['no3'],\
    #                                data2plot['R'+zr1]['DELTA']['depT'], 'k', alpha=0.1, lw = 1.5)
    lldn9 , = ax[1, 0].plot(data2plot['R9']  ['DELTA']['no3'] , data2plot['R9']  ['DELTA']['depT'], 'k', lw = lw['R9']  )
    lldn27, = ax[1, 0].plot(data2plot['R27'] ['DELTA']['no3'] , data2plot['R27'] ['DELTA']['depT'], 'k', lw = lw['R27'] )
    ax[1, 0].set_title('(e) NO3 change', loc='left')
    ax[1, 0].set_xlim((-7, 0))
    ax[1, 0].xaxis.set_ticks([-6, -3, 0])
    ax[1, 0].set_ylabel('Depth [m]')
    ax[1, 0].set_xlabel('[mmol-N/m$^3$]')

    lld1  = plotdata(ax[1, 1], 'AVG1', 'DELTA', '(f) CC1 - CTL1')
    lld9  = plotdata(ax[1, 2], 'R9'  , 'DELTA', '(g) CC9 - CTL9')
    lld27 = plotdata(ax[1, 3], 'R27' , 'DELTA', '(h) CC27 - CTL27')
    for zax in ax[1, 1:].flatten() : 
        zax.set_xlabel('[mmolN/m$^2$/d]')
        zax.set_xlim((-1.7, .7))
        zax.xaxis.set_ticks([ -1.4, -.7, 0, .7])
        zax.tick_params(left=False)
        if zax.is_last_col() :
            zax.set_ylabel('Depth [m]')
            zax.yaxis.set_ticks_position('right')
            zax.yaxis.set_label_position('right')
        #
    #
    
    fig.subplots_adjust(hspace=.3)

    hdl = [lln1, lln9, lln27]
    nam = ['CTL1', 'CTL9', 'CTL27']
    ax[0, 0].legend(hdl, nam, handlelength=1., loc='best')
    hdl = lld27[0]
    nam = lld27[1]
    zw = ax[0, 2].get_position()
    legax = fig.add_axes([zw.x0, .92, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=4, handlelength=2., loc='lower center')
    #fig.legend( hdl, nam, handlelength=2., ncol=2, bbox_to_anchor=(.67, .92), loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    savefig='fig4.png'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END FIG4
##################################################

##################################################
# FIG5
##################################################
if plot_fig5 : 

    savefig='fig5.png'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/dyk/rdyk004/MY_PYTHON3/PCKL/'

    data2process = {}
    
    #___________________
    # R1s TOTAL NO3 BUDGET

    for zr1 in r1s : 
        data2process['R'+zr1]={}
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CTL'] = open_pckl(zwin, 'CTL'+zr1)
        zwin = 'main_no3_budget_v3_cc' +zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CC'] = open_pckl(zwin, 'CC'+zr1)
        zwin = 'main_no3_budget_dynXtrc_ctlXcc_r'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CTLxCC'] = open_pckl(zwin, 'R'+zr1)
        zwin = 'main_no3_budget_dynXtrc_ccXctl_r'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CCxCTL'] = open_pckl(zwin, 'R'+zr1)

    #
    
    #___________________
    # R9 TOTAL NO3 BUDGET

    data2process['R9']={}
    zwin = 'main_no3_budget_v3_ctl9_35N.pckl'
    data2process['R9']['CTL'] = open_pckl(zwin, 'CTL9')
    zwin = 'main_no3_budget_v3_cc9_35N.pckl'
    data2process['R9']['CC'] = open_pckl(zwin, 'CC9')
    zwin = 'main_no3_budget_dynXtrc_ctlXcc_r9_35N.pckl'
    data2process['R9']['CTLxCC'] = open_pckl(zwin, 'R9')
    zwin = 'main_no3_budget_dynXtrc_ccXctl_r9_35N.pckl'
    data2process['R9']['CCxCTL'] = open_pckl(zwin, 'R9')

    #___________________
    # R27 TOTAL NO3 BUDGET

    data2process['R27']={}
    
    for zsim in ['CTL', 'CC'] :
        data2process['R27'][zsim] = {}
        zwvars={}
        for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                    'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff', \
                    'xadonl', 'yadonl', 'zadonl', 'hadonl', 'advonl', \
                    'zdf', 'ldf', 'dynonl', \
                    'sms', 'nwp', 'nit', 'exp', 'src'] :
            zwvars[vvv]=[]
        #
        zwsuf = 'main_no3_budget_v3_'+zsim.lower()+'27_35N_y1'
        for yyy in np.arange(66, 71) :
            zwin = zwsuf+str(yyy)+'.pckl'
            zw = open_pckl(zwin, zsim+'27')
            if yyy==66 : data2process['R27'][zsim]['mesh']=zw['mesh']
            for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
        #
        for vvv in zwvars.keys() : data2process['R27'][zsim][vvv]=np.mean(zwvars[vvv], axis=0)
    #

    data2process['R27']['CTLxCC'] = {}
    zwvars={}
    for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff'] :
        zwvars[vvv]=[]
    #
    zwsuf = 'main_no3_budget_dynXtrc_ctlXcc_r27_35N_y1'
    for yyy in np.arange(66, 71) :
        zwin = zwsuf+str(yyy)+'.pckl'
        zw = open_pckl(zwin, 'R27')
        if yyy==66 : data2process['R27']['CTLxCC']['mesh']=zw['mesh']
        for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
    #
    for vvv in zwvars.keys() : data2process['R27']['CTLxCC'][vvv]=np.mean(zwvars[vvv], axis=0)

    data2process['R27']['CCxCTL'] = {}
    zwvars={}
    for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff'] :
        zwvars[vvv]=[]
    #
    zwsuf = 'main_no3_budget_dynXtrc_ccXctl_r27_35N_y1'
    for yyy in np.arange(66, 71) :
        zwin = zwsuf+str(yyy)+'.pckl'
        zw = open_pckl(zwin, 'R27')
        if yyy==66 : data2process['R27']['CCxCTL']['mesh']=zw['mesh']
        for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
    #
    for vvv in zwvars.keys() : data2process['R27']['CCxCTL'][vvv]=np.mean(zwvars[vvv], axis=0)

    ####################
    # PREPARE NO3 BUDGET DATA
    ####################

    data2plot={}

    for zkres in data2process.keys() :
        data2plot[zkres]={}
        for zksim in ['CTL', 'CC', 'CTLxCC', 'CCxCTL'] :
            data2plot[zkres][zksim]={'depW':data2process[zkres][zksim]['mesh']['depW']}
            data2plot[zkres][zksim]['advt']  = data2process[zkres][zksim]['vN1'] + \
                data2process[zkres][zksim]['wN']
            # add advection by GM velocities to total advection
            if 'sfxgm' in data2process[zkres][zksim].keys() :
                data2plot[zkres][zksim]['advt'] += data2process[zkres][zksim]['vgmN1']
                data2plot[zkres][zksim]['advt'] += data2process[zkres][zksim]['wgmN']
            #
        #
        data2plot[zkres]['dVN']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['dVN']['advt'] = data2plot[zkres]['CC']['advt'] - \
            data2plot[zkres]['CTL']['advt']
        data2plot[zkres]['VdN']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['VdN']['advt'] = data2plot[zkres]['CTLxCC']['advt'] - \
            data2plot[zkres]['CTL']['advt']
        data2plot[zkres]['NdV']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['NdV']['advt'] = data2plot[zkres]['CCxCTL']['advt'] - \
            data2plot[zkres]['CTL']['advt']
        data2plot[zkres]['dVdN']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['dVdN']['advt'] = data2plot[zkres]['dVN']['advt'] - \
            data2plot[zkres]['VdN']['advt'] - data2plot[zkres]['NdV']['advt']
    #
        
    ####################
    # AVG1, STD1
    ####################

    data2plot['AVG1'] = {}
    data2plot['STD1'] = {}
    for zksim in ['CTL', 'CC', 'CTLxCC', 'dVN', 'VdN', 'NdV', 'dVdN'] :
        data2plot['AVG1'][zksim] = {'depW': data2plot['R'+r1s[0]][zksim]['depW']}
        data2plot['STD1'][zksim] = {'depW': data2plot['R'+r1s[0]][zksim]['depW']}
        zw = []
        for vkr1 in r1s : zw.append(data2plot['R'+vkr1][zksim]['advt'])
        zw = np.array(zw)
        data2plot['AVG1'][zksim]['advt'] = np.nanmean(zw, axis=0)
        data2plot['STD1'][zksim]['advt'] = unc*np.nanstd(zw, axis=0)
    #

    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    fsize   = (11*infact, 5*infact) #(width, height)
    fig, ax = plt.subplots(1, 3, sharey='row', figsize=fsize)
    
    ccc = {'dVN':'k', 'VdN':'firebrick', 'NdV':'mediumblue', 'dVdN':'.7'}
    nnn = {'dVN':'$\Delta$(u$\cdot$N)', 'VdN':'$\Delta$Nitrate', \
           'NdV':'$\Delta$Circulation', 'dVdN':'Non-linear $\Delta$'}
    ooo = {'dVN':10, 'VdN':15, 'NdV':20, 'dVdN':5}

    def plotdata(zax, zres, ztitle) :
        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
        zlines=[]
        znames = []
        for kpart in ['dVN', 'VdN', 'NdV', 'dVdN'] :
            zw=data2plot[zres][kpart]
            zl, =zax.plot(zfact * zw['advt'], zw['depW'], lw=1.5, c=ccc[kpart], zorder=ooo[kpart])
            zlines.append(zl)
            znames.append(nnn[kpart])
            if zres == 'AVG1' :
                zwstd=data2plot['STD1'][kpart]
                Xm, Xp = zw['advt']-zwstd['advt'], zw['advt']+zwstd['advt']
                zax.fill_betweenx(zw['depW'], zfact*Xm, zfact*Xp, alpha=0.2, color=ccc[kpart], zorder=5)
            #
        #
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        zax.set_title(ztitle, loc='left')
        zax.vlines(0, 0, 400, color='grey')
        # zlines = [zldvn, zlvdn]
        # znames = ['$\Delta$(v$\cdot$NO3)', 'v$_{CTL}$$\cdot$$\Delta$NO3']
        # zlines = [zldvn, zlvdn, zlndv, zldvdn]
        # znames = ['$\Delta$(v$\cdot$NO3)', 'v$_{CTL}$$\cdot$$\Delta$NO3', \
        #           'N$_{CTL}$$\cdot$$\Delta$v', '$\Delta$v$_{CTL}$$\cdot$$\Delta$NO3']
        return (zlines, znames)
    #
    
    lld1  = plotdata(ax[0], 'AVG1', '(a) 1° avg. $\pm$ st. dev.')
    lld9  = plotdata(ax[1], 'R9'  , '(b) 1/9°')
    lld27 = plotdata(ax[2], 'R27' , '(c) 1/27°')
    for zax in ax[:].flatten() : 
        zax.set_xlabel('[mmolN/m$^2$/d]')
        zax.set_xlim((-1.7, .7))
        zax.xaxis.set_ticks([ -1.4, -.7, 0, .7])
        if zax.is_last_col() :
            zax.set_ylabel('Depth [m]')
            zax.yaxis.set_ticks_position('right')
            zax.yaxis.set_label_position('right')
        #
    #
    ax[0].set_ylabel('Depth [m]')
    for zax in ax[1:].flatten() : zax.tick_params(left=False)
    
    hdl = lld27[0]
    nam = lld27[1]
    zw = ax[1].get_position()
    legax = fig.add_axes([zw.x0, .95, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=4, loc='lower center')
    #ax[1].legend( hdl, nam, loc='best')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    print("Figure saved: ", dirfig+suffig+savefig)
    fig.savefig(dirfig+suffig+savefig)
    plt.close()

##################################################
# END FIG5
##################################################

##################################################
# FIG6
##################################################
if plot_fig6 : 

    savefig='fig6.png'
    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    sdate = '0166-01-01'
    edate = '0170-12-31'
    res = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', \
           '9', '27']
    lat1, lat2 = 35, 45

    data2plot = {}
    
    ####################
    # READ AND PROCESS W
    ####################

    fsuf = '_1y_01010101_01701230_grid_W.xml'
    var  = 'vovecrtz'
    data2plot['W GSP'] = {}
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['W GSP'][zksim] = {}
        for zkres in res :
            fff = fdir + zksim + zkres + fsuf
            zw = reading.read_ncdf(var, fff, time = (sdate, edate))
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx', grid='W')
            zw = averaging.ymean(zw, zwmesh, dim='zyx', ymin=lat1, ymax=lat2, grid='W')
            zw = averaging.xmean(zw, zwmesh, dim='zx', grid='W')
            zfact =  3180*1e3 * (2833-1721)*1e3*1e-6  # m/s -> 1e6m3/s=Sv
            data2plot['W GSP'][zksim]['R'+zkres]={'data':zfact*zw, 'dep':zwmesh['depW']}
        #
    #
    #___________________
    # DELTA
    data2plot['W GSP']['DELTA'] = {}
    for zkres in res :
        data2plot['W GSP']['DELTA']['R'+zkres]= \
            {'data':data2plot['W GSP']['CC']['R'+zkres]['data'] - \
             data2plot['W GSP']['CTL']['R'+zkres]['data'], \
             'dep':data2plot['W GSP']['CTL']['R'+zkres]['dep']}
    #

    ####################
    # READ AND PROCESS V
    ####################

    fsuf = '_1y_01010101_01701230_grid_V.xml'
    var  = 'vomecrty'
    data2plot['V 35N'] = {}
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['V 35N'][zksim] = {}
        for zkres in res :
            fff = fdir + zksim + zkres + fsuf
            zw = reading.read_ncdf(var, fff, time = (sdate, edate))
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx', grid='V')
            zw = interpolating.yinterpol(zw, zwmesh, lat1, dim='zyx', grid='V')
            zw = averaging.xmean(zw, zwmesh, dim='zx', grid='V')
            zfact =  3180*1e3*1e-6 # m/s -> 1e6m2/s=Sv/m
            data2plot['V 35N'][zksim]['R'+zkres]={'data':zfact*zw*zwmesh['e3t_0'], 'dep':zwmesh['depV']}
        #
    #
    #___________________
    # DELTA
    data2plot['V 35N']['DELTA'] = {}
    for zkres in res :
        data2plot['V 35N']['DELTA']['R'+zkres]= \
            {'data':data2plot['V 35N']['CC']['R'+zkres]['data'] - \
             data2plot['V 35N']['CTL']['R'+zkres]['data'], \
             'dep':data2plot['V 35N']['CTL']['R'+zkres]['dep']}
    #
    
    ####################
    # AVG1, STD1
    ####################

    r1s = ['R1', 'R1_KL2000', 'R1_KL500', 'R1_KGM500_KL500', 'R1_KGM2000_KL2000']
    for zkvar in data2plot.keys() :
        for zksim in data2plot[zkvar].keys() :
            data2plot[zkvar][zksim]['AVG1']={'dep':data2plot[zkvar][zksim]['R1']['dep']}
            data2plot[zkvar][zksim]['STD1']={'dep':data2plot[zkvar][zksim]['R1']['dep']}
            zw = []
            for zkr1 in r1s : zw.append(data2plot[zkvar][zksim][zkr1]['data'])
            zw = np.array(zw)
            data2plot[zkvar][zksim]['AVG1']['data'] = \
                np.ma.array( data=np.nanmean(zw, axis=0),\
                             mask=data2plot[zkvar][zksim]['R1']['data'].mask )
            data2plot[zkvar][zksim]['STD1']['data'] = \
                np.ma.array( data=unc*np.nanstd(zw, axis=0),\
                             mask=data2plot[zkvar][zksim]['R1']['data'].mask )
        #
    #

    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    fsize   = (6*infact, 8*infact) #(width, height)
    fig, ax = plt.subplots(2, 2, sharey='row', figsize=fsize) # nrow, ncol

    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()
    ttls={'W GSP'   : 'Upw. vel.', \
          'V 35N'   : 'Northw. vel.'}
    units={'W GSP'   : '[Sv]', \
           'V 35N'   : '[Sv]'}
    nnn = {'R1'               :'1°, k$_{gm}$=1e$^3$\nk$_{redi}$=1e$^3$' , \
           'R1_KL500'         :'1°, k$_{gm}$=1e$^3$\nk$_{redi}$=500'       , \
           'R1_KL2000'        :'1°, k$_{gm}$=1e$^3$\nk$_{redi}$=2e$^3$' , \
           'R1_KGM500_KL500'  :'1°, k$_{gm}$=500\nk$_{redi}$=500'       , \
           'R1_KGM2000_KL2000':'1°, k$_{gm}$=2e$^3$\nk$_{redi}$=2e$^3$' , \
           'AVG1':'1° avg $\pm$ st. dev.', 'R9':'1/9°', 'R27':'1/27°'}
    www = {'R1'               :1., \
           'R1_KL500'         :1., \
           'R1_KL2000'        :1., \
           'R1_KGM500_KL500'  :1., \
           'R1_KGM2000_KL2000':1., \
           'AVG1':3.5, 'R9':2.5, 'R27':1.5}
    ccc = {'R1'               :'0', \
           'R1_KL500'         :'0', \
           'R1_KL2000'        :'0', \
           'R1_KGM500_KL500'  :'0', \
           'R1_KGM2000_KL2000':'0', \
           'AVG1':'0', 'R9':'.4', 'R27':'.7'}

    def plotdata(zax, zvar, zsim) :
        zdat=data2plot[zvar][zsim]
        zlines = []
        znames = []
        for zkres in ['AVG1', 'R9', 'R27'] :
            X, Y = zdat[zkres]['data'], zdat[zkres]['dep']
            zl, = zax.plot(X, Y, lw = www[zkres], c=ccc[zkres], ls='-')
            zlines.append(zl)
            znames.append(nnn[zkres])
            if zkres == 'AVG1' :
                zdatstd=zdat['STD1']['data']
                Xm, Xp = X - zdatstd, X + zdatstd
                zax.fill_betweenx(Y, Xm, Xp, alpha=0.1, color='k')
            #
        #
        zax.set_xlabel(units[zvar])
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        if   zsim=='CTL'   :
            zax.set_xlim((-1.5, 5.5))
            zax.set_title('('+subnum.pop()+') '+ttls[zvar], loc='left')
        elif zsim=='DELTA' :
            zax.set_xlim((-1., .5))
            zax.set_title('('+subnum.pop()+') ∆ '+ttls[zvar], loc='left')
        zax.vlines(0, 0, 400, color='grey')
        return (zlines, znames)
    #

    ll = np.zeros_like(ax)
    icol = 0
    for zkvar in data2plot.keys() :
        ll[0, icol] = plotdata( ax[0, icol], zkvar, 'CTL' )
        ll[1, icol] = plotdata( ax[1, icol], zkvar, 'DELTA' )
        icol+=1
    #
    for zax in ax[:, 0].flatten() : zax.set_ylabel('Depth [m]')
    for zax in ax[0, :].flatten() : zax.set_xlabel('')
    for zax in ax[:, 1:-1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #
    fig.subplots_adjust(hspace=.3)

    hdl = ll[0, 0][0]
    nam = ll[0, 0][1]
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 1].get_position()
    nwidth = zw2.x0+zw2.width - zw1.x0
    ny0 = zw1.y0+zw1.height*1.2
    legax = fig.add_axes([zw1.x0, ny0,  nwidth, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=3, loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/dyk/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_v4_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END FIG6
##################################################

##################################################
# SAVE VALUES NO3BUDGET
################################################## 
if save_values_no3budget: 

    savefile='save_values_no3budget.txt'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/dyk/rdyk004/MY_PYTHON3/PCKL/'

    data2process = {}
    
    #___________________
    # R1s TOTAL NO3 BUDGET

    for zr1 in r1s : 
        data2process['R'+zr1]={}
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CTL'] = open_pckl(zwin, 'CTL'+zr1)
        zwin = 'main_no3_budget_v3_cc' +zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CC'] = open_pckl(zwin, 'CC'+zr1)
    #
    
    #___________________
    # R9 TOTAL NO3 BUDGET

    data2process['R9']={}
    zwin = 'main_no3_budget_v3_ctl9_35N.pckl'
    data2process['R9']['CTL'] = open_pckl(zwin, 'CTL9')
    zwin = 'main_no3_budget_v3_cc9_35N.pckl'
    data2process['R9']['CC'] = open_pckl(zwin, 'CC9')

    #___________________
    # R27 TOTAL NO3 BUDGET

    data2process['R27']={}
    
    for zsim in ['CTL', 'CC'] :
        data2process['R27'][zsim] = {}
        zwvars={}
        for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                    'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff', \
                    'xadonl', 'yadonl', 'zadonl', 'hadonl', 'advonl', \
                    'zdf', 'ldf', 'dynonl', \
                    'sms', 'nwp', 'nit', 'exp', 'src'] :
            zwvars[vvv]=[]
        #
        zwsuf = 'main_no3_budget_v3_'+zsim.lower()+'27_35N_y1'
        for yyy in np.arange(66, 71) :
            zwin = zwsuf+str(yyy)+'.pckl'
            zw = open_pckl(zwin, zsim+'27')
            if yyy==66 : data2process['R27'][zsim]['mesh']=zw['mesh']
            for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
        #
        for vvv in zwvars.keys() :
            zw = np.mean(zwvars[vvv], axis=0)
            zmsk = zw!=zw
            data2process['R27'][zsim][vvv]=np.ma.array(data=zw, mask=zmsk)
    #

    ####################
    # PREPARE NO3 BUDGET DATA
    ####################

    data2process2={}

    for zkres in data2process.keys() :
        data2process2[zkres]={}
        for zksim in ['CTL', 'CC'] :
            data2process2[zkres][zksim]={'mesh':data2process[zkres][zksim]['mesh']}
            data2process2[zkres][zksim]['zmix']  = data2process[zkres][zksim]['zdf']
            data2process2[zkres][zksim]['advt']  = data2process[zkres][zksim]['vN1'] + \
                data2process[zkres][zksim]['wN']
            # add ldf to total advection
            data2process2[zkres][zksim]['advt'] += data2process[zkres][zksim]['ldf']
            # add advection by GM velocities to total advection
            if 'sfxgm' in data2process[zkres][zksim].keys() :
                data2process2[zkres][zksim]['advt'] += data2process[zkres][zksim]['vgmN1']
                data2process2[zkres][zksim]['advt'] += data2process[zkres][zksim]['wgmN']
            #
        #
        data2process2[zkres]['DELTA']={'mesh':data2process[zkres]['CTL']['mesh']}
        data2process2[zkres]['DELTA']['zmix'] = data2process2[zkres]['CC']['zmix'] - \
            data2process2[zkres]['CTL']['zmix']
        data2process2[zkres]['DELTA']['advt'] = data2process2[zkres]['CC']['advt'] - \
            data2process2[zkres]['CTL']['advt']
    #
        
    ####################
    # PREPARE PP DATA
    ####################

    fdir  = '/gpfswork/rech/dyk/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_pp.xml'
    var = 'PP'
    ymin, ymax = 35, 45
    sdate, edate  = '0166-01-01', '0170-12-31'

    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    sim  = ['CTL', 'CC']
    for zkres in res :
        if zkres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/dyk/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        for zksim in sim :
            fff = fdir + zksim + zkres + fsuf
            # fff = fdir + zksim + '1' + fsuf # TEST
            zw  = reading.read_ncdf(var, fff, time=(sdate, edate))
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx')/24/3600. # convert to s-1
            zw = averaging.ymean(zw, zwmesh, dim='zyx', ymin=ymin, ymax=ymax, integral=True)
            zw = averaging.xmean(zw, zwmesh, dim='zx', integral=True)
            zw = np.cumsum(zw * zwmesh['e3t_0'])
            zw = np.roll(zw, 1) # we shift because first value is 2nd W-point
            zw[0] = 0. # flux at the surface (W-point = 0) is zero
            data2process2['R'+zkres][zksim]['npp'] = zw
        #
        data2process2['R'+zkres]['DELTA']['npp'] = data2process2['R'+zkres]['CC']['npp'] - \
            data2process2['R'+zkres]['CTL']['npp']
    #
    
    ####################
    # AVG1, STD1
    ####################

    data2process2['AVG1'] = {}
    data2process2['STD1'] = {}
    for zksim in ['CTL', 'CC', 'DELTA'] :
        data2process2['AVG1'][zksim] = {'mesh': data2process2['R'+r1s[0]][zksim]['mesh']}
        data2process2['STD1'][zksim] = {'mesh': data2process2['R'+r1s[0]][zksim]['mesh']}
        for vvv in ['zmix', 'advt', 'npp'] : 
            zw = []
            for vkr1 in r1s : zw.append(data2process2['R'+vkr1][zksim][vvv])
            zw1 = np.nanmean(zw, axis=0)
            zmsk = zw1!=zw1
            data2process2['AVG1'][zksim][vvv] = np.ma.array(data=zw1, mask=zmsk)
            zw1 = unc*np.nanstd(zw, axis=0)
            zmsk = zw1!=zw1
            data2process2['STD1'][zksim][vvv] = np.ma.array(data=zw1, mask=zmsk)
        #
    #
    
    ####################
    # EXTRACT VALUES AT FIXED DEPTH
    ####################

    data2save = {}
    valdep = [25, 50, 100, 200, 400, 4000]

    for idep in valdep : # loop on depth
        data2save[str(idep)]={}
        for zkres in data2process2.keys() : # loop on resolution
            data2save[str(idep)][zkres]={}
            for zksim in ['CTL', 'CC', 'DELTA'] : # loop on simulation
                data2save[str(idep)][zkres][zksim]={}
                zwmesh = data2process2[zkres][zksim]['mesh']
                for zkvar in data2process2[zkres][zksim].keys() : # loop on variables
                    if zkvar != 'mesh' :
                        print(str(idep)+', '+zkres+', '+zksim+', '+zkvar)
                        zw = data2process2[zkres][zksim][zkvar]
                        zw = interpolating.zinterpol(zw, zwmesh, idep, grid='W', dim='z')
                        zw = zw.data
                        zw = zw.flatten()
                        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
                        zw    = np.around(zw[0]*zfact, decimals=2)
                        data2save[str(idep)][zkres][zksim][zkvar] = zw
                    #
                # end loop on variables
            # end loop on simulations
        # end loop on resolution
    # end loop on depth

    ####################
    # SAVE IN TEXT FILE
    ####################

    dirfile = '/gpfswork/rech/dyk/rdyk004/MY_PYTHON3/'
    suffile = 'main_article_gyre_pp_v4_'
    fff = open(dirfile+suffile+savefile, "w+")

    for vdep in data2save.keys() :
        zdat=data2save[vdep]
        fff.write("\n##################################################\n")
        fff.write("               DEPTH: "+str(vdep)+" meters")
        fff.write("\n##################################################\n")
        fff.write("     ")
        fff.write("                   CTL              //")
        fff.write("                   CC               //")
        fff.write("                  DELTA               \n")
        fff.write("             AVG1        R9      R27     //")
        fff.write("        AVG1        R9      R27     //")
        fff.write("        AVG1        R9      R27     ")
        fff.write("\n-------------------------------------------------------------------------------------------------------------------\n")
        fff.write("NPP: ")
        fff.write('    ' + str(zdat['AVG1']['CTL']['npp']) + ' pm ' + str(zdat['STD1']['CTL']['npp']) + '    ')
        fff.write(         str(zdat['R9']  ['CTL']['npp']) + '    ')
        fff.write(         str(zdat['R27'] ['CTL']['npp']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['CC']['npp']) + ' pm ' + str(zdat['STD1']['CC']['npp']) + '    ')
        fff.write(         str(zdat['R9']  ['CC']['npp']) + '    ')
        fff.write(         str(zdat['R27'] ['CC']['npp']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['DELTA']['npp']) + ' pm ' + str(zdat['STD1']['DELTA']['npp']) + '    ')
        fff.write(         str(zdat['R9']  ['DELTA']['npp']) + '    ')
        fff.write(         str(zdat['R27'] ['DELTA']['npp']) + '    ')
        fff.write("\n-------------------------------------------------------------------------------------------------------------------\n")
        fff.write("ZMIX:")
        fff.write('    ' + str(zdat['AVG1']['CTL']['zmix']) + ' pm ' + str(zdat['STD1']['CTL']['zmix']) + '    ')
        fff.write(         str(zdat['R9']  ['CTL']['zmix']) + '    ')
        fff.write(         str(zdat['R27'] ['CTL']['zmix']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['CC']['zmix']) + ' pm ' + str(zdat['STD1']['CC']['zmix']) + '    ')
        fff.write(         str(zdat['R9']  ['CC']['zmix']) + '    ')
        fff.write(         str(zdat['R27'] ['CC']['zmix']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['DELTA']['zmix']) + ' pm ' + str(zdat['STD1']['DELTA']['zmix']) + '    ')
        fff.write(         str(zdat['R9']  ['DELTA']['zmix']) + '    ')
        fff.write(         str(zdat['R27'] ['DELTA']['zmix']) + '    ')
        fff.write("\n-------------------------------------------------------------------------------------------------------------------\n")
        fff.write("ADVT:")
        fff.write('    ' + str(zdat['AVG1']['CTL']['advt']) + ' pm ' + str(zdat['STD1']['CTL']['advt']) + '    ')
        fff.write(         str(zdat['R9']  ['CTL']['advt']) + '    ')
        fff.write(         str(zdat['R27'] ['CTL']['advt']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['CC']['advt']) + ' pm ' + str(zdat['STD1']['CC']['advt']) + '    ')
        fff.write(         str(zdat['R9']  ['CC']['advt']) + '    ')
        fff.write(         str(zdat['R27'] ['CC']['advt']) + '    //')
        fff.write('    ' + str(zdat['AVG1']['DELTA']['advt']) + ' pm ' + str(zdat['STD1']['DELTA']['advt']) + '    ')
        fff.write(         str(zdat['R9']  ['DELTA']['advt']) + '    ')
        fff.write(         str(zdat['R27'] ['DELTA']['advt']) + '    ')
        fff.write("\n-------------------------------------------------------------------------------------------------------------------\n")        
        fff.write("\n")
        fff.write("\n")
    #

    fff.close()
    print("File saved: ", dirfile+suffile+savefile)
    
##################################################
# END SAVE VALUES NO3BUDGET
################################################## 

embed()
