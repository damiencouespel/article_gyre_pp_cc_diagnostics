from IPython import embed
import time
import pickle

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('./matplotlibrc_nature_cc.mplstyle')

from FGYRE import reading, averaging, interpolating, \
    density, stream_functions, plotting

plot_sf1  = False
plot_sf2  = False
plot_sf3  = False
plot_sf3_old  = False
plot_sf4  = False
plot_sf5  = False
plot_sf6  = False
plot_sf7  = False
plot_sf8  = False
plot_sf9  = False
plot_sf10 = False
plot_sf11 = True
plot_sf12 = False

unc = 1.

##################################################
# SF1
##################################################
if plot_sf1 :
    print(">>> plot_sf1 <<<")
    
    savefig = 'sf1.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    iyears = 10
    years  = np.arange(1, 170, iyears)
    data2plot = {}
    
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST

    ####################
    # READ AND PREPARE TEMPERATURE 
    ####################

    data2plot['tem']={}
    fsuf = '_1y_00010101_01701230_grid_T.xml'
    sim  = ['CTL', 'CC']
    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['tem']['R'+vres]={}
        for vsim in sim :
            data2plot['tem']['R'+vres][vsim] = []
            for vyears in years : 
                sdate  = '0'+str(vyears)+'-01-01'
                edate    = '0'+str(vyears+iyears-1)+'-12-31'
                fff = fdir + vsim + vres + fsuf
                zw  = reading.read_ncdf('votemper', fff, time=(sdate, edate))
                zw = zw['data']
                zw = averaging.zmean(zw, zwmesh, dim='tzyx', zmax = 700.)
                zw = averaging.ymean(zw, zwmesh, dim='tyx')
                zw = averaging.xmean(zw, zwmesh, dim='tx')
                data2plot['tem']['R'+vres][vsim].extend(list(zw))
            #
        #
    #

    ####################
    # READ AND PREPARE NO3
    ####################

    data2plot['no3']={}
    fsuf = '_1y_00010101_01701230_ptrc_T.xml'
    sim  = ['CTL', 'CC']
    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['no3']['R'+vres]={}
        for vsim in sim :
            data2plot['no3']['R'+vres][vsim] = []
            for vyears in years : 
                sdate  = '0'+str(vyears)+'-01-01'
                edate    = '0'+str(vyears+iyears-1)+'-12-31'
                fff = fdir + vsim + vres + fsuf
                zw  = reading.read_ncdf('NO3', fff, time=(sdate, edate))
                zw = zw['data']
                zw = averaging.zmean(zw, zwmesh, dim='tzyx', zmax = 700.)
                zw = averaging.ymean(zw, zwmesh, dim='tyx')
                zw = averaging.xmean(zw, zwmesh, dim='tx')
                data2plot['no3']['R'+vres][vsim].extend(list(zw))
            #
        #
    #

    ####################
    # READ AND PREPARE NPP
    ####################

    data2plot['npp']={}
    fsuf = '_1y_00010101_01701230_diad_T.xml'
    sim  = ['CTL', 'CC']
    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['npp']['R'+vres]={}
        for vsim in sim :
            data2plot['npp']['R'+vres][vsim] = []
            for vyears in years : 
                sdate  = '0'+str(vyears)+'-01-01'
                edate    = '0'+str(vyears+iyears-1)+'-12-31'
                fff = fdir + vsim + vres + fsuf
                zw1 = reading.read_ncdf('TNO3PHY', fff, time=(sdate, edate))
                zw = zw1['data']
                zw1 = reading.read_ncdf('TNH4PHY', fff, time=(sdate, edate))
                zw += zw1['data']
                zw = averaging.ymean(zw, zwmesh, dim='tyx')
                zw = averaging.xmean(zw, zwmesh, dim='tx')
                data2plot['npp']['R'+vres][vsim].extend(list(zw))
            #
        #
    #

    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    fsize   = (10*infact, 15*infact) #(width, height)
    fig, ax = plt.subplots(3, 1, sharex='col', figsize=fsize) # nrow, ncol

    ccc={'R1'               :'.5', \
         'R1_KL500'         :'.25', \
         'R1_KL2000'        :'.7', \
         'R1_KGM500_KL500'  :'0.', \
         'R1_KGM2000_KL2000':'.9', \
         'R9':'0', 'R27':'0'}
    www={'R1'               :3.5, \
         'R1_KL500'         :3.5, \
         'R1_KL2000'        :3.5, \
         'R1_KGM500_KL500'  :3.5, \
         'R1_KGM2000_KL2000':3.5, \
         'R9':2.5, 'R27':1.5}
    nnn={'R1'               :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=1e$^3$', \
         'R1_KL500'         :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=500'      , \
         'R1_KL2000'        :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=2e$^3$', \
         'R1_KGM500_KL500'  :'1°, k$_{gm}$=500\nk$_{iso}$=500'      , \
         'R1_KGM2000_KL2000':'1°, k$_{gm}$=2e$^3$\nk$_{iso}$=2e$^3$', \
         'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()
    ttls = {'tem':'Ocean temperature 0-700 m [°C]', \
            'no3':'NO$_3$ concentration 0-700 m [mmol$\,$N$\,$m$^{-3}$]', \
            'npp':'Net primary production [mmol$\,$N$\,$m$^{-2}$d$^{-1}$]'}
    
    def plotdata(zax, zvar) :
        zdat = data2plot[zvar]
        zlines1 = []
        zlines2 = []
        znames = []
        X = np.arange(170)+1
        rrr = ['R1_KGM500_KL500', 'R1_KL500', 'R1' ,'R1_KL2000','R1_KGM2000_KL2000', 'R9', 'R27']
        for zres in rrr :
            zl1, = zax.plot(X, zdat[zres]['CTL'], lw=www[zres], c=ccc[zres], ls='--', dashes=[1.5, .7])
            zl2, = zax.plot(X, zdat[zres]['CC' ], lw=www[zres], c=ccc[zres], ls='-' )
            zlines1.append(zl1)
            zlines2.append(zl2)
            znames.append(nnn[zres])
        #
        zax.set_title('('+subnum.pop()+') '+ttls[zvar], loc='left')
        zax.locator_params(axis='y', nbins=6)
        return (zlines1, zlines2, znames)
    #

    ll = np.zeros_like(ax)
    irow = 0
    for var in ['tem', 'no3', 'npp'] : 
        ll[irow] = plotdata(ax[irow], var)
        irow+=1
    #
    ax[-1].set_xlim(1, 171)
    ax[-1].xaxis.set_ticks([1, 101, 171])
    ax[-1].xaxis.set_ticklabels(['-100 years', '0', '70 years'])

    hdl = ll[0][1]
    nam = ll[0][2]
    zw = ax[1].get_position()
    legax = fig.add_axes([.9, zw.y0, 0.01, zw.height])
    legax.axis('off')
    legax.legend( hdl, nam, loc='center left')

    ax[-1].annotate('Spin-up', xy=(.3, .15), xycoords='axes fraction', \
                 bbox=dict(boxstyle="round", fc="1"))
    hdl = [ ll[0][0][-1], ll[0][1][-1] ]
    nam = [ 'CTL'       , 'CC'         ]
    ax[-1].legend( hdl, nam, edgecolor='k', bbox_to_anchor=(.6, .15), \
                   loc='center left')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig) 
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()


##################################################
# END SF1
##################################################

##################################################
# SF2
##################################################
if plot_sf2 :
    print(">>> plot_sf2 <<<")

    savefig = 'sf2.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}

    ####################
    # BSF 
    ####################

    data2plot['bsf']={}
    fsuf  = '_1y_01010101_01701230_grid_V.xml'
    for vres in res : 
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['bsf']['CTL'+vres]={'lon':zwmesh['lonV'], 'lat':zwmesh['latV']}
        fff = fdir + 'CTL' + vres + fsuf
        zw  = reading.read_ncdf('vomecrty', fff, time=(sdate, edate))
        zw = zw['data']
        zw = stream_functions.bsfv(zw, zwmesh)
        data2plot['bsf']['CTL'+vres]['data'] = np.mean(zw, axis = 0)
    #
    
    ####################
    # NPP
    ####################

    data2plot['npp']={}
    fsuf = '_1y_00010101_01701230_diad_T.xml'
    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['npp']['CTL'+vres]={'lon':zwmesh['lonT'], 'lat':zwmesh['latT']}
        fff = fdir + 'CTL' + vres + fsuf
        zw1 = reading.read_ncdf('TNO3PHY', fff, time=(sdate, edate))
        zw = zw1['data']
        zw1 = reading.read_ncdf('TNH4PHY', fff, time=(sdate, edate))
        zw += zw1['data']
        zw = averaging.tmean(zw, zwmesh, dim='tyx') 
        data2plot['npp']['CTL'+vres]['data'] = zw
    #
    
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 2
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    aa = np.linspace(-30, 0, 6)
    bb = np.linspace(0, 30, 6)
    levbsf = np.concatenate((aa[:-1], bb))
    levnpp = np.linspace(0, 4., 11)
    lev  = {'bsf':levbsf, 'npp':levnpp}
    cmap = {'bsf':'RdBu', 'npp':'viridis'}
    ext  = {'bsf':'both', 'npp':'max'}
    ttls={'CTL1'               :'k$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'CTL1_KL500'         :'k$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'CTL1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'CTL1_KGM500_KL500'  :'k$_{gm}$=500, k$_{iso}$=500'      , \
          'CTL1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'CTL9':'', 'CTL27':''}
    note={'CTL1'               :'1°', \
          'CTL1_KL500'         :'1°', \
          'CTL1_KL2000'        :'1°', \
          'CTL1_KGM500_KL500'  :'1°', \
          'CTL1_KGM2000_KL2000':'1°', \
          'CTL9':'1/9°', 'CTL27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zvar, zres) :
        zdat=data2plot[zvar][zres]
        X, Y, Z = zdat['lon'], zdat['lat'], zdat['data']
        zcf = zax.contourf(X, Y, Z, levels=lev[zvar], cmap=cmap[zvar], extend=ext[zvar])
        ttt='('+subnum.pop()+') '+ttls[zres]
        plotting.make_XY_axis(zax, title=ttt)
        zax.label_outer()
        zax.annotate(note[zres], xy=(-83, 23), \
                     bbox=dict(boxstyle="round", fc="1"))
        return zcf
    #

    cf = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zres in rrr : 
        icol = 0
        for zvar in ['bsf', 'npp'] :
            cf[irow, icol] = plotdata(ax[irow, icol], zvar, 'CTL'+zres)
            icol+=1
        #
        irow+=1
    #

    cbtitle = 'Bar. circ. [Sv]'
    zw = ax[0, 0].get_position()
    cbar_ax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, 0.01])
    fig.colorbar(cf[0, 0], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[-30, 0, 30])

    cbtitle = 'NPP [mmol$\,$N$\,$m$^{-2}$d$^{-1}$]'
    zw = ax[0, 1].get_position()
    cbar_ax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, 0.01])
    fig.colorbar(cf[0, 1], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[0, 2, 4])

    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Northward km')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #
    #fig.subplots_adjust(hspace=.3)

        
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig) 
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()


##################################################
# END SF2
##################################################

##################################################
# SF3 OLD
##################################################
if plot_sf3_old : 
    print(">>> plot_sf3_old <<<")

    savefig = 'sf3_old.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf='_1y_01010101_01701230_diad_T.xml'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}

    ####################
    # MAPS ∆NP ∆RP
    ####################

    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    var={'NP':'TNO3PHY',
         'RP':'TNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['R'+vres] = {'lon':zwmesh['lonT'], 'lat':zwmesh['latT']}
        for vvv in ['NP', 'RP'] :
            zw = {}
            for vsim in ['CTL', 'CC'] : 
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                zw[vsim]  = averaging.tmean(zw2, zwmesh, dim='tyx')
            #
            data2plot['R'+vres]['m'+vvv] = zw['CC'] - zw['CTL']
        #
        zw = data2plot['R'+vres]['mNP'] + data2plot['R'+vres]['mRP']
        for vvv in ['NP', 'RP'] : data2plot['R'+vres]['m'+vvv] =  data2plot['R'+vres]['m'+vvv] / np.abs(zw) * 100
    #

    ####################
    # TIMSERIES ∆NP ∆RP
    ####################

    sdate  = '0101-01-01'
    edate  = '0170-12-31'
    var={'NP':'TNO3PHY',
         'RP':'TNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        for vvv in ['NP', 'RP'] :
            zw = {}
            for vsim in ['CTL', 'CC'] : 
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                zw2 = averaging.ymean(zw2, zwmesh, grid='T', dim='tyx', ymin=35, ymax=45) 
                zw2 = averaging.xmean(zw2, zwmesh, grid='T', dim='tx')
                zw[vsim] = zw2
            #
            data2plot['R'+vres]['t'+vvv] = zw['CC'] - zw['CTL']
        #
    #

    ####################
    # TIMSERIES F-RATIO
    ####################

    sdate  = '0101-01-01'
    edate  = '0170-12-31'
    var={'NP':'TNO3PHY',
         'RP':'TNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        for vsim in ['CTL', 'CC'] :
            zw = {}
            for vvv in ['NP', 'RP'] :
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                zw2 = averaging.ymean(zw2, zwmesh, grid='T', dim='tyx', ymin=35, ymax=45) 
                zw2 = averaging.xmean(zw2, zwmesh, grid='T', dim='tx')
                zw[vvv] = zw2
            #
            data2plot['R'+vres]['f'+vsim] = zw['NP'] / (zw['NP'] + zw['RP'])
        #
    #


    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 4
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharex='col',figsize=fsize)

    ttls={'R1'               :'\nk$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'\nk$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'\nk$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'\nk$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'\nk$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum=subnum+list(map('a'.__add__,subnum))
    subnum.reverse()

    def plotdata(zaxrow, zres) :
        zdat=data2plot[zres]
        #___________________
        # MAPS
        aa  = np.linspace(-100, 0, 6)
        bb  = np.linspace(0, 100, 6)
        lev = np.concatenate((aa[:-1], bb))
        X, Y = zdat['lon'], zdat['lat']
        # NP
        zcfNP = zaxrow[0].contourf(X, Y, zdat['mNP'], levels=lev, cmap='PiYG', extend='both')
        ttt='('+subnum.pop()+') NPP from NO$_3$'+ttls[zres]
        plotting.make_XY_axis(zaxrow[0], title=ttt)
        zaxrow[0].label_outer()
        zaxrow[0].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        # RP
        zcfRP = zaxrow[1].contourf(X, Y, zdat['mRP'], levels=lev, cmap='PiYG', extend='both')
        ttt='('+subnum.pop()+') NPP from NH$_4$'+ttls[zres]
        plotting.make_XY_axis(zaxrow[1], title=ttt)
        zaxrow[1].label_outer()
        zaxrow[1].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        # TIMESERIES NP RP
        X = np.arange(70)+1
        Y1 = zdat['tRP']
        Y2 = Y1 + zdat['tNP']
        zllRP   =  zaxrow[2].fill_between(X, Y1    , color='0')
        zllNP   =  zaxrow[2].fill_between(X, Y1, Y2, color='.5')
        ttt='('+subnum.pop()+') ' + ttls[zres]
        zaxrow[2].set_title(ttt, loc='left')
        zaxrow[2].set_ylim((-1., 0.1))
        zaxrow[2].locator_params(axis='y', nbins=4)
        zaxrow[2].annotate(note[zres], xy=(7, -.85), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        # TIMESERIES F-RATIO
        X = np.arange(70)+1
        Y1 = zdat['fCTL']
        Y2 = zdat['fCC']
        zllfCT,  =  zaxrow[3].plot(X, Y1, color='0', ls='--', dashes=[1.5, .7])
        zllfCC,  =  zaxrow[3].plot(X, Y2, color='0', ls='-' )
        ttt='('+subnum.pop()+') ' + ttls[zres]
        zaxrow[3].set_title(ttt, loc='left')
        zaxrow[3].set_ylim((0, 1))
        zaxrow[3].locator_params(axis='y', nbins=4)
        zaxrow[3].annotate(note[zres], xy=(7, .136), \
                        bbox=dict(boxstyle="round", fc="1"))
        # zaxrow[3].annotate(note[zres], xy=(7, .527), \
        #                 bbox=dict(boxstyle="round", fc="1"))
        return (zcfNP, zcfRP, zllNP, zllRP, zllfCT, zllfCC)
    #

    cf = np.zeros_like(ax[:, 0])
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zres in rrr :
        cf[irow] = plotdata(ax[irow], 'R'+zres)
        irow+=1
    #
    fig.subplots_adjust(hspace=.4)

    cbtitle = 'Percent of NPP change '
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 1].get_position()
    nwidth = zw1.width*1.5
    nx0 = (zw1.x0+zw1.width+zw2.x0)/2 - nwidth/2
    ny0 = zw1.y0+zw1.height*1.4
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.01])
    fig.colorbar(cf[0][0], cax=cbar_ax, orientation='horizontal', ticklocation='top', label=cbtitle)

    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -2].flatten() :
        zax.set_ylabel('[mmolN/m$^2$/d]')
        zw = zax.get_position()
        zax.set_position([zw.x0+.09, zw.y0, zw.width, zw.height])
    #
    for zax in ax[:, -1].flatten() :
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
        zw = zax.get_position()
        zax.set_position([zw.x0+.11, zw.y0, zw.width, zw.height])
    #
    ax[-1, -2].set_xlim(1, 71)
    ax[-1, -2].xaxis.set_ticks([1, 71])
    ax[-1, -2].xaxis.set_ticklabels(['0', '70 years'])
    ax[-1, -1].set_xlim(1, 71)
    ax[-1, -1].xaxis.set_ticks([1, 71])
    ax[-1, -1].xaxis.set_ticklabels(['0', '70 years'])

    hdl=(cf[0][3], cf[0][2])
    nam=('from NH$_4$', 'from NO$_3$')  
    zw = ax[0, -2].get_position()
    legax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, .01])
    legax.axis('off')
    leg=legax.legend( hdl, nam, loc='lower center', title='Mean NPP change:')
    leg._legend_box.align = "left"

    hdl=(cf[0][4], cf[0][5])
    nam=('CTL', 'CC')  
    zw = ax[0, -1].get_position()
    legax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, .01])
    legax.axis('off')
    leg=legax.legend( hdl, nam, loc='lower center', title='Mean f-ratio:')
    leg._legend_box.align = "left"

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF3 OLD
##################################################

##################################################
# SF3
##################################################
if plot_sf3 : 
    print(">>> plot_sf3 <<<")

    savefig = 'sf3.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf='_1y_01010101_01701230_diad_T.xml'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}

    ####################
    # MAPS ∆NP ∆RP
    ####################

    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    var={'NP':'FNO3PHY',
         'RP':'FNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        data2plot['R'+vres] = {'lon':zwmesh['lonT'], 'lat':zwmesh['latT']}
        for vvv in ['NP', 'RP'] :
            zw = {}
            for vsim in ['CTL', 'CC'] : 
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                # remove nitrification from NP and add to RP
                zwnit = reading.read_ncdf('FNH4NO3', fff, time=(sdate, edate))
                zwnit = zwnit['data']
                if   vvv == 'NP' : zw2 = zw2-zwnit
                elif vvv == 'RP' : zw2 = zw2+zwnit
                else : exit('vvv must be NP or RP')
                # DC
                zw2  = averaging.tmean(zw2, zwmesh, dim='tzyx')
                zw[vsim]  = averaging.zmean(zw2, zwmesh, dim='zyx', zmax=100, integral=True)
                zw2 = averaging.zmean(zw2, zwmesh, grid='T', dim='zyx', zmax=100, integral=True) 
                zw2 = averaging.ymean(zw2, zwmesh, grid='T', dim='yx', ymin=35, ymax=45) 
                zw2 = averaging.xmean(zw2, zwmesh, grid='T', dim='x')                
                zwnit = averaging.tmean(zwnit, zwmesh, dim='tzyx')
                zwnit = averaging.zmean(zwnit, zwmesh, grid='T', dim='zyx', zmax=100, integral=True) 
                zwnit = averaging.ymean(zwnit, zwmesh, grid='T', dim='yx', ymin=35, ymax=45) 
                zwnit = averaging.xmean(zwnit, zwmesh, grid='T', dim='x')
                print(">>> "+vsim+vres)
                print(vvv+": "+str(zw2))
                print("NIT: "+str(zwnit))
            #
            data2plot['R'+vres]['m'+vvv] = zw['CC'] - zw['CTL']
        #
        zw = data2plot['R'+vres]['mNP'] + data2plot['R'+vres]['mRP']
        for vvv in ['NP', 'RP'] : data2plot['R'+vres]['m'+vvv] =  data2plot['R'+vres]['m'+vvv] / np.abs(zw) * 100
    #

    ####################
    # TIMSERIES ∆NP ∆RP
    ####################

    sdate  = '0101-01-01'
    edate  = '0170-12-31'
    var={'NP':'FNO3PHY',
         'RP':'FNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        for vvv in ['NP', 'RP'] :
            zw = {}
            for vsim in ['CTL', 'CC'] : 
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                # DC
                zwnit = reading.read_ncdf('FNH4NO3', fff, time=(sdate, edate))
                if   vvv == 'NP' : zw2 = zw2-zwnit['data']
                elif vvv == 'RP' : zw2 = zw2+zwnit['data']
                else : exit('vvv must be NP or RP')
                # DC
                zw2 = averaging.zmean(zw2, zwmesh, grid='T', dim='tzyx', zmax=100, integral=True) 
                zw2 = averaging.ymean(zw2, zwmesh, grid='T', dim='tyx', ymin=35, ymax=45) 
                zw2 = averaging.xmean(zw2, zwmesh, grid='T', dim='tx')
                zw[vsim] = zw2
            #
            data2plot['R'+vres]['t'+vvv] = zw['CC'] - zw['CTL']
        #
    #

    ####################
    # TIMSERIES F-RATIO
    ####################

    sdate  = '0101-01-01'
    edate  = '0170-12-31'
    var={'NP':'FNO3PHY',
         'RP':'FNH4PHY'}

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
        for vsim in ['CTL', 'CC'] :
            zw = {}
            for vvv in ['NP', 'RP'] :
                fff = fdir + vsim + vres + fsuf
                zw2  = reading.read_ncdf(var[vvv], fff, time=(sdate, edate))
                zw2  = zw2['data']
                # DC
                zwnit = reading.read_ncdf('FNH4NO3', fff, time=(sdate, edate))
                if   vvv == 'NP' : zw2 = zw2-zwnit['data']
                elif vvv == 'RP' : zw2 = zw2+zwnit['data']
                else : exit('vvv must be NP or RP')
                # DC
                zw2 = averaging.zmean(zw2, zwmesh, grid='T', dim='tzyx', zmax=100, integral=True) 
                zw2 = averaging.ymean(zw2, zwmesh, grid='T', dim='tyx', ymin=35, ymax=45) 
                zw2 = averaging.xmean(zw2, zwmesh, grid='T', dim='tx')
                zw[vvv] = zw2
            #
            data2plot['R'+vres]['f'+vsim] = zw['NP'] / (zw['NP'] + zw['RP'])
            print(">>> "+vsim+vres)
            print("f-ratio: "+str(data2plot['R'+vres]['f'+vsim][-1]))
        #
    #


    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 4
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharex='col',figsize=fsize)

    ttls={'R1'               :'\nk$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'\nk$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'\nk$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'\nk$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'\nk$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum=subnum+list(map('a'.__add__,subnum))
    subnum.reverse()

    def plotdata(zaxrow, zres) :
        zdat=data2plot[zres]
        #___________________
        # MAPS
        aa  = np.linspace(-100, 0, 6)
        bb  = np.linspace(0, 100, 6)
        lev = np.concatenate((aa[:-1], bb))
        X, Y = zdat['lon'], zdat['lat']
        # NP
        zcfNP = zaxrow[0].contourf(X, Y, zdat['mNP'], levels=lev, cmap='PiYG', extend='both')
        ttt='('+subnum.pop()+') New NPP'+ttls[zres]
        plotting.make_XY_axis(zaxrow[0], title=ttt)
        zaxrow[0].label_outer()
        zaxrow[0].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        # RP
        zcfRP = zaxrow[1].contourf(X, Y, zdat['mRP'], levels=lev, cmap='PiYG', extend='both')
        ttt='('+subnum.pop()+') Reg NPP'+ttls[zres]
        plotting.make_XY_axis(zaxrow[1], title=ttt)
        zaxrow[1].label_outer()
        zaxrow[1].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        # TIMESERIES NP RP
        X = np.arange(70)+1
        Y1 = zdat['tRP']
        Y2 = Y1 + zdat['tNP']
        zllRP   =  zaxrow[2].fill_between(X, Y1    , color='0')
        zllNP   =  zaxrow[2].fill_between(X, Y1, Y2, color='.5')
        ttt='('+subnum.pop()+') ' + ttls[zres]
        zaxrow[2].set_title(ttt, loc='left')
        zaxrow[2].set_ylim((-1., 0.1))
        zaxrow[2].locator_params(axis='y', nbins=4)
        zaxrow[2].annotate(note[zres], xy=(7, -.85), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        # TIMESERIES F-RATIO
        X = np.arange(70)+1
        Y1 = zdat['fCTL']
        Y2 = zdat['fCC']
        zllfCT,  =  zaxrow[3].plot(X, Y1, color='0', ls='--', dashes=[1.5, .7])
        zllfCC,  =  zaxrow[3].plot(X, Y2, color='0', ls='-' )
        ttt='('+subnum.pop()+') ' + ttls[zres]
        zaxrow[3].set_title(ttt, loc='left')
        zaxrow[3].set_ylim((0, 1))
        zaxrow[3].locator_params(axis='y', nbins=4)
        zaxrow[3].annotate(note[zres], xy=(7, .136), \
                        bbox=dict(boxstyle="round", fc="1"))
        # zaxrow[3].annotate(note[zres], xy=(7, .527), \
        #                 bbox=dict(boxstyle="round", fc="1"))
        return (zcfNP, zcfRP, zllNP, zllRP, zllfCT, zllfCC)
    #

    cf = np.zeros_like(ax[:, 0])
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zres in rrr :
        cf[irow] = plotdata(ax[irow], 'R'+zres)
        irow+=1
    #
    fig.subplots_adjust(hspace=.4)

    cbtitle = 'Percent of NPP change '
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 1].get_position()
    nwidth = zw1.width*1.5
    nx0 = (zw1.x0+zw1.width+zw2.x0)/2 - nwidth/2
    ny0 = zw1.y0+zw1.height*1.4
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.01])
    fig.colorbar(cf[0][0], cax=cbar_ax, orientation='horizontal', ticklocation='top', label=cbtitle)

    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -2].flatten() :
        zax.set_ylabel('[mmolN/m$^2$/d]')
        zw = zax.get_position()
        zax.set_position([zw.x0+.09, zw.y0, zw.width, zw.height])
    #
    for zax in ax[:, -1].flatten() :
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
        zw = zax.get_position()
        zax.set_position([zw.x0+.11, zw.y0, zw.width, zw.height])
    #
    ax[-1, -2].set_xlim(1, 71)
    ax[-1, -2].xaxis.set_ticks([1, 71])
    ax[-1, -2].xaxis.set_ticklabels(['0', '70 years'])
    ax[-1, -1].set_xlim(1, 71)
    ax[-1, -1].xaxis.set_ticks([1, 71])
    ax[-1, -1].xaxis.set_ticklabels(['0', '70 years'])

    hdl=(cf[0][3], cf[0][2])
    nam=('New NPP', 'Reg NPP')  
    zw = ax[0, -2].get_position()
    legax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, .01])
    legax.axis('off')
    leg=legax.legend( hdl, nam, loc='lower center', title='Mean NPP change:')
    leg._legend_box.align = "left"

    hdl=(cf[0][4], cf[0][5])
    nam=('CTL', 'CC')  
    zw = ax[0, -1].get_position()
    legax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, .01])
    legax.axis('off')
    leg=legax.legend( hdl, nam, loc='lower center', title='Mean f-ratio:')
    leg._legend_box.align = "left"

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF3
##################################################

##################################################
# SF4
##################################################
if plot_sf4 :
    print(">>> plot_sf4 <<<")

    savefig = 'sf4.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_ptrc_T.xml'
    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}
    
    ####################
    # PREPARE NO3
    ####################

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        data2plot['R'+vres]={'lat':zwmesh['latT'][:, 5], 'dep':zwmesh['depT']}
        for vsim in ['CTL', 'CC'] : 
            fff = fdir + vsim + vres + fsuf
            zw  = reading.read_ncdf('NO3', fff, time=(sdate, edate))
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh) 
            zw = averaging.xmean(zw, zwmesh, dim='zyx') 
            data2plot['R'+vres][vsim] = zw
        #
        data2plot['R'+vres]['DELTA'] = data2plot['R'+vres]['CC'] - data2plot['R'+vres]['CTL']
    #
    
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*4*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    aa  = np.linspace(-4, 0, 6)
    bb  = np.linspace(0, 4, 6)
    levD = np.concatenate((aa[:-1], bb))
    levC = np.linspace(0, 20., 11)
    lev  = {'CTL':levC, 'CC':levC, 'DELTA':levD}
    cmap = {'CTL':'viridis', 'CC':'viridis', 'DELTA':'PiYG'}
    ext  = {'CTL':'max', 'CC':'max', 'DELTA':'both'}
    ttls={'R1'               :'\nk$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'\nk$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'\nk$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'\nk$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'\nk$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zres, zsim) :
        zdat=data2plot[zres]
        X, Y, Z = zdat['lat'], zdat['dep'], zdat[zsim]
        zcf = zax.contourf(X, Y, Z, levels=lev[zsim], cmap=cmap[zsim], extend=ext[zsim])
        ttt='('+subnum.pop()+') '+zsim+' '+ttls[zres]
        plotting.make_YZ_axis(zax, depmax=800, title=ttt)
        zax.label_outer()
        zax.annotate(note[zres], xy=(23, 700), \
                     bbox=dict(boxstyle="round", fc="1"))
        return zcf
    #

    cf = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zres in rrr :
        icol = 0
        for zsim in ['CTL', 'CC', 'DELTA'] :
            cf[irow, icol] = plotdata(ax[irow, icol], 'R'+zres, zsim)
            icol+=1
        #
        irow+=1
    #


    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1:].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #
    fig.subplots_adjust(hspace=.4)

    cbtitle = 'NO$_3$ [mmol$\,$N$\,$m$^{-3}$]'
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 1].get_position()
    nwidth = zw1.width*1.5
    nx0 = (zw1.x0+zw1.width+zw2.x0)/2 - nwidth/2
    ny0 = zw1.y0+zw1.height*1.4
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.01])
    fig.colorbar(cf[0,0], cax=cbar_ax, orientation='horizontal', ticklocation='top', label=cbtitle)

    cbtitle = '∆NO$_3$ [mmol$\,$N$\,$m$^{-3}$]'
    zw = ax[0, -1].get_position()
    ny0 = zw.y0+zw.height*1.4
    cbar_ax = fig.add_axes([zw.x0, ny0, zw.width, 0.01])
    fig.colorbar(cf[0, -1], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[-4, 0, 4])

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig) 
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()
#
##################################################
# END SF4
##################################################

##################################################
# SF5
##################################################
if plot_sf5 : 
    print(">>> plot_sf5 <<<")

    savefig='sf5.png'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/eee/rdyk004/MY_PYTHON3/PCKL/'

    data2process = {}
    
    #___________________
    # R1s TOTAL NO3 BUDGET

    for zr1 in r1s : 
        data2process['R'+zr1]={}
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N.pckl'
        data2process['R'+zr1]['CTL'] = open_pckl(zwin, 'CTL'+zr1)
    #
    
    #___________________
    # R1s MEAN NO3 BUDGET

    for zr1 in r1s : 
        zwin = 'main_no3_budget_v3_ctl'+zr1.lower()+'_35N_1ymean.pckl'
        data2process['R'+zr1]['CTL_M'] = open_pckl(zwin, 'CTL'+zr1)
    #

    #___________________
    # R9 TOTAL NO3 BUDGET

    data2process['R9']={}
    zwin = 'main_no3_budget_v3_ctl9_35N.pckl'
    data2process['R9']['CTL'] = open_pckl(zwin, 'CTL9')

    #___________________
    # R9 MEAN NO3 BUDGET

    zwin = 'main_no3_budget_v3_ctl9_35N_onR1_1ymean.pckl'
    data2process['R9']['CTL_M'] = open_pckl(zwin, 'CTL9')

    #___________________
    # R27 TOTAL NO3 BUDGET

    data2process['R27']={}
    
    
    data2process['R27']['CTL'] = {}
    zwvars={}
    for vvv in ['uN1', 'uN2', 'vN1', 'vN2', 'wN', 'sfx', \
                'xadoff', 'yadoff', 'zadoff', 'hadoff', 'advoff', \
                'xadonl', 'yadonl', 'zadonl', 'hadonl', 'advonl', \
                'zdf', 'ldf', 'dynonl', \
                'sms', 'nwp', 'nit', 'exp', 'src'] :
        zwvars[vvv]=[]
    #
    zwsuf = 'main_no3_budget_v3_ctl27_35N_y1'
    for yyy in np.arange(66, 71) :
        zwin = zwsuf+str(yyy)+'.pckl'
        zw = open_pckl(zwin, 'CTL27')
        if yyy==66 : data2process['R27']['CTL']['mesh']=zw['mesh']
        for vvv in zwvars.keys() : zwvars[vvv].append(zw[vvv])
    #
    for vvv in zwvars.keys() : data2process['R27']['CTL'][vvv]=np.mean(zwvars[vvv], axis=0)

    #___________________
    # R27 MEAN NO3 BUDGET

    zwin = 'main_no3_budget_v3_ctl27_35N_onR1_1ymean.pckl'
    data2process['R27']['CTL_M'] = open_pckl(zwin, 'CTL27')

    ####################
    # PREPARE NO3 BUDGET DATA
    ####################

    data2plot={}

    for zkres in data2process.keys() :
        data2plot[zkres]={}
        #___________________
        # TOTAL
        data2plot[zkres]['TOTAL']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['TOTAL']['vN'] = data2process[zkres]['CTL']['vN1']
        data2plot[zkres]['TOTAL']['wN'] = data2process[zkres]['CTL']['wN']
        # add advection by GM velocities to total advection
        if 'sfxgm' in data2process[zkres]['CTL'].keys() :
            data2plot[zkres]['TOTAL']['vN'] += data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['TOTAL']['wN'] += data2process[zkres]['CTL']['wgmN']
        #
        data2plot[zkres]['TOTAL']['ldf']  = data2process[zkres]['CTL']['ldf']
        data2plot[zkres]['TOTAL']['adv'] = data2plot[zkres]['TOTAL']['vN'] + \
            data2plot[zkres]['TOTAL']['wN'] + data2plot[zkres]['TOTAL']['ldf']
        #___________________
        # MEAN
        data2plot[zkres]['MEAN']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['MEAN']['vN'] = data2process[zkres]['CTL_M']['vN1']
        data2plot[zkres]['MEAN']['wN'] = data2process[zkres]['CTL_M']['wN']
        data2plot[zkres]['MEAN']['adv'] = data2plot[zkres]['MEAN']['vN'] + \
            data2plot[zkres]['MEAN']['wN']
        #___________________
        # EDDY
        data2plot[zkres]['EDDY']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['EDDY']['vN']  = data2plot[zkres]['TOTAL']['vN']  - data2plot[zkres]['MEAN']['vN']
        data2plot[zkres]['EDDY']['wN']  = data2plot[zkres]['TOTAL']['wN']  - data2plot[zkres]['MEAN']['wN']
        data2plot[zkres]['EDDY']['ldf'] = data2process[zkres]['CTL']['ldf']
        data2plot[zkres]['EDDY']['adv'] = data2plot[zkres]['TOTAL']['adv'] - data2plot[zkres]['MEAN']['adv']
        if 'sfxgm' in data2process[zkres]['CTL'].keys() :
            data2plot[zkres]['EDDY']['vNgm'] = data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['EDDY']['wNgm'] = data2process[zkres]['CTL']['wgmN']
        #
    #
        
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    ccc = {'adv':'k', 'vN':'darkorange', 'wN':'royalblue', 'ldf':'limegreen',\
           'vNgm':'k', 'wNgm':'k'}
    www = {'adv':1.5, 'vN':1.5, 'wN':1.5, 'ldf':1.5, 'vNgm':.5, 'wNgm':.5}
    sss = {'adv':'-', 'vN':'-', 'wN':'-', 'ldf':'-', 'vNgm':':', 'wNgm':':'}
    aaa = {'adv':.3, 'vN':1, 'wN':1, 'ldf':1, 'vNgm':1, 'wNgm':1}
    nnn = {'adv':'Total', 'vN':'Meri. adv.', 'wN':'Vert. adv.', 'ldf':'Iso. mix.', \
           'vNgm':'GM adv.', 'wNgm':'GM adv.'}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    ttls={'R1'               :'k$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'k$_{gm}$=1e$^3$, k$_{iso}$=500', \
          'R1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'k$_{gm}$=500, k$_{iso}$=500', \
          'R1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    def plotdata(zax, zres, zpart) :
        zdat=data2plot[zres][zpart]
        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
        zlines = []
        znames = []
        Y = zdat.pop('depW')
        for zkproc in zdat.keys() :
            zl, = zax.plot(zfact * zdat[zkproc], Y, lw=1.5, c=ccc[zkproc], ls=sss[zkproc], alpha=aaa[zkproc])
            zax.annotate(note[zres], xy=(-1.3, 60), \
                         bbox=dict(boxstyle="round", fc="1"))
            zlines.append(zl)
            znames.append(nnn[zkproc])
        #
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        zax.set_title('('+subnum.pop()+') '+zpart.upper()+'\n'+ttls[zres], loc='left')
        zax.vlines(0, 0, 400, color='grey')
        return (zlines, znames)
    #

    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    ll = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zkr in rrr :
        icol = 0
        for zkpart in ['TOTAL', 'MEAN', 'EDDY'] :
            ll[irow, icol] = plotdata(ax[irow, icol], 'R'+zkr, zkpart)
            icol+=1
        #
        irow+=1
    #
    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, 0].flatten() : zax.set_ylabel('Depth [m]')
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    for zax in ax[-1, :].flatten() : 
        zax.set_xlim((-1.5, 3.5))
        zax.set_xlabel('mmol$\,$N$\,$m$^{-2}$d$^{-1}$')
        zax.xaxis.set_ticks([-1.5, 0, 1.5, 3])
    #    
    fig.subplots_adjust(hspace=.5)

    hdl = ll[0, -1][0][:-1]
    nam = ll[0, -1][1][:-1]
    zw = ax[0, 1].get_position()
    legax = fig.add_axes([zw.x0, .91, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=5, handlelength=1., loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig) 
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF5
##################################################

##################################################
# SF6
##################################################
if plot_sf6 : 
    print(">>> plot_sf6 <<<")

    savefig='sf6.png'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/eee/rdyk004/MY_PYTHON3/PCKL/'

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
        #___________________
        # TOTAL
        data2plot[zkres]['TOTAL']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['TOTAL']['vN'] = data2process[zkres]['CC']['vN1'] - data2process[zkres]['CTL']['vN1']
        data2plot[zkres]['TOTAL']['wN'] = data2process[zkres]['CC']['wN']  - data2process[zkres]['CTL']['wN']
        # add advection by GM velocities to total advection
        if 'sfxgm' in data2process[zkres]['CTL'].keys() :
            data2plot[zkres]['TOTAL']['vN'] += data2process[zkres]['CC']['vgmN1'] - data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['TOTAL']['wN'] += data2process[zkres]['CC']['wgmN']  - data2process[zkres]['CTL']['wgmN'] 
        #
        data2plot[zkres]['TOTAL']['ldf'] = data2process[zkres]['CC']['ldf'] - data2process[zkres]['CTL']['ldf']
        data2plot[zkres]['TOTAL']['adv'] = data2plot[zkres]['TOTAL']['vN'] + \
            data2plot[zkres]['TOTAL']['wN'] + data2plot[zkres]['TOTAL']['ldf'] 
        #___________________
        # MEAN
        data2plot[zkres]['MEAN']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['MEAN']['vN'] = data2process[zkres]['CC_M']['vN1'] - data2process[zkres]['CTL_M']['vN1']
        data2plot[zkres]['MEAN']['wN'] = data2process[zkres]['CC_M']['wN']  - data2process[zkres]['CTL_M']['wN'] 
        data2plot[zkres]['MEAN']['adv'] = data2plot[zkres]['MEAN']['vN'] + \
            data2plot[zkres]['MEAN']['wN']
        #___________________
        # EDDY
        data2plot[zkres]['EDDY']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['EDDY']['vN']  = data2plot[zkres]['TOTAL']['vN']  - data2plot[zkres]['MEAN']['vN']
        data2plot[zkres]['EDDY']['wN']  = data2plot[zkres]['TOTAL']['wN']  - data2plot[zkres]['MEAN']['wN']
        data2plot[zkres]['EDDY']['ldf'] = data2process[zkres]['CC']['ldf'] - data2process[zkres]['CTL']['ldf']
        data2plot[zkres]['EDDY']['adv'] = data2plot[zkres]['TOTAL']['adv'] - data2plot[zkres]['MEAN']['adv']
        if 'sfxgm' in data2process[zkres]['CTL'].keys() :
            data2plot[zkres]['EDDY']['vNgm'] = data2process[zkres]['CC']['vgmN1'] - data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['EDDY']['wNgm'] = data2process[zkres]['CC']['wgmN']  - data2process[zkres]['CTL']['wgmN'] 
        #
    #
        
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    ccc = {'adv':'k', 'vN':'darkorange', 'wN':'royalblue', 'ldf':'limegreen',\
           'vNgm':'k', 'wNgm':'k'}
    www = {'adv':1.5, 'vN':1.5, 'wN':1.5, 'ldf':1.5, 'vNgm':.5, 'wNgm':.5}
    sss = {'adv':'-', 'vN':'-', 'wN':'-', 'ldf':'-', 'vNgm':':', 'wNgm':':'}
    aaa = {'adv':.3, 'vN':1, 'wN':1, 'ldf':1, 'vNgm':1, 'wNgm':1}
    nnn = {'adv':'Total', 'vN':'Meri. adv.', 'wN':'Vert. adv.', 'ldf':'Iso. mix.', \
           'vNgm':'GM adv.', 'wNgm':'GM adv.'}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    ttls={'R1'               :'k$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'k$_{gm}$=1e$^3$, k$_{iso}$=500', \
          'R1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'k$_{gm}$=500, k$_{iso}$=500', \
          'R1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    def plotdata(zax, zres, zpart) :
        zdat=data2plot[zres][zpart]
        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
        zlines = []
        znames = []
        Y = zdat.pop('depW')
        for zkproc in zdat.keys() :
            zl, = zax.plot(zfact * zdat[zkproc], Y, lw=1.5, c=ccc[zkproc], ls=sss[zkproc], alpha=aaa[zkproc])
            zlines.append(zl)
            znames.append(nnn[zkproc])
            zax.annotate(note[zres], xy=(-1.4, 80), \
                         bbox=dict(boxstyle="round", fc="1"))
        #
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        zax.set_title('('+subnum.pop()+') '+zpart.upper()+'\n'+ttls[zres], loc='left')
        zax.vlines(0, 0, 400, color='grey')
        return (zlines, znames)
    #

    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()
    ll = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zkr in rrr : 
        icol = 0
        for zkpart in ['TOTAL', 'MEAN', 'EDDY'] :
            ll[irow, icol] = plotdata(ax[irow, icol], 'R'+zkr, zkpart)
            icol+=1
        #
        irow+=1
    #
    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, 0].flatten() : zax.set_ylabel('Depth [m]')
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    for zax in ax[-1, :].flatten() : 
        zax.set_xlabel('mmol$\,$N$\,$m$^{-2}$d$^{-1}$')
        zax.set_xlim((-1.7, .7))
        zax.xaxis.set_ticks([ -1.4, -.7, 0, .7])
    #    
    fig.subplots_adjust(hspace=.5)

    hdl = ll[0, -1][0][:-1]
    nam = ll[0, -1][1][:-1]
    zw = ax[0, 1].get_position()
    legax = fig.add_axes([zw.x0, .91, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=5, handlelength=1., loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF6
##################################################

##################################################
# SF7
##################################################
if plot_sf7 : 
    print(">>> plot_sf7 <<<")

    savefig='sf7.png'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']

    ####################
    # READ NO3 BUDGET DATA
    ####################
    
    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #    

    fdir = '/gpfswork/rech/eee/rdyk004/MY_PYTHON3/PCKL/'

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
        # vert. adv.
        data2plot[zkres]['vert. adv.']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['vert. adv.']['dvN'] = data2process[zkres]['CC']['wN'] - \
            data2process[zkres]['CTL']['wN']
        data2plot[zkres]['vert. adv.']['vdN'] = data2process[zkres]['CTLxCC']['wN'] - \
            data2process[zkres]['CTL']['wN']
        data2plot[zkres]['vert. adv.']['Ndv'] = data2process[zkres]['CCxCTL']['wN'] - \
            data2process[zkres]['CTL']['wN']
        # meri. adv.
        data2plot[zkres]['meri. adv.']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['meri. adv.']['dvN'] = data2process[zkres]['CC']['vN1'] - \
            data2process[zkres]['CTL']['vN1']
        data2plot[zkres]['meri. adv.']['vdN'] = data2process[zkres]['CTLxCC']['vN1'] - \
            data2process[zkres]['CTL']['vN1']
        data2plot[zkres]['meri. adv.']['Ndv'] = data2process[zkres]['CCxCTL']['vN1'] - \
            data2process[zkres]['CTL']['vN1']
        # GM
        if 'sfxgm' in data2process[zkres]['CTL'].keys() :
            # vert
            data2plot[zkres]['vert. adv.']['dvN'] += data2process[zkres]['CC']['wgmN'] - \
                data2process[zkres]['CTL']['wgmN']
            data2plot[zkres]['vert. adv.']['vdN'] += data2process[zkres]['CTLxCC']['wgmN'] - \
                data2process[zkres]['CTL']['wgmN']
            data2plot[zkres]['vert. adv.']['Ndv'] += data2process[zkres]['CCxCTL']['wgmN'] - \
                data2process[zkres]['CTL']['wgmN']
            # meri. adv.
            data2plot[zkres]['meri. adv.']['dvN'] += data2process[zkres]['CC']['vgmN1'] - \
                data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['meri. adv.']['vdN'] += data2process[zkres]['CTLxCC']['vgmN1'] - \
                data2process[zkres]['CTL']['vgmN1']
            data2plot[zkres]['meri. adv.']['Ndv'] += data2process[zkres]['CCxCTL']['vgmN1'] - \
                data2process[zkres]['CTL']['vgmN1']
        #
        data2plot[zkres]['vert. adv.']['dvdN'] = data2plot[zkres]['vert. adv.']['dvN'] - \
            data2plot[zkres]['vert. adv.']['vdN'] - data2plot[zkres]['vert. adv.']['Ndv']
        data2plot[zkres]['meri. adv.']['dvdN'] = data2plot[zkres]['meri. adv.']['dvN'] - \
            data2plot[zkres]['meri. adv.']['vdN'] - data2plot[zkres]['meri. adv.']['Ndv']
        # total
        data2plot[zkres]['total']={'depW':data2process[zkres]['CTL']['mesh']['depW']}
        data2plot[zkres]['total']['dvN'] = data2plot[zkres]['vert. adv.']['dvN'] + \
            data2plot[zkres]['meri. adv.']['dvN']
        data2plot[zkres]['total']['vdN'] = data2plot[zkres]['vert. adv.']['vdN'] + \
            data2plot[zkres]['meri. adv.']['vdN']
        data2plot[zkres]['total']['Ndv'] = data2plot[zkres]['vert. adv.']['Ndv'] + \
            data2plot[zkres]['meri. adv.']['Ndv']
        data2plot[zkres]['total']['dvdN'] = data2plot[zkres]['vert. adv.']['dvdN'] + \
            data2plot[zkres]['meri. adv.']['dvdN']
    #
    
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    ccc = {'dvN':'k', 'vdN':'firebrick', 'Ndv':'mediumblue', 'dvdN':'.7'}
    nnn = {'dvN':'$\Delta$($\mathbf{u}\cdot$N)', 'vdN':'$\Delta$Nitrate', \
           'Ndv':'$\Delta$Circulation', 'dvdN':'Non-linear $\Delta$'}
    ooo = {'dvN':10, 'vdN':15, 'Ndv':20, 'dvdN':5}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    ttls={'R1'               :'k$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'k$_{gm}$=1e$^3$, k$_{iso}$=500', \
          'R1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'k$_{gm}$=500, k$_{iso}$=500', \
          'R1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zres, zpart) :
        zdat=data2plot[zres][zpart]
        zfact = 3600*24. / ( 3180*1e3 * (2833-1721)*1e3 ) # mmol/s -> mmol/m2/d
        zlines = []
        znames = []
        Y = zdat.pop('depW')
        for zkcontrib in zdat.keys() :
            zl, = zax.plot(zfact * zdat[zkcontrib], Y, lw=1.5, c=ccc[zkcontrib], zorder=ooo[zkcontrib])
            zlines.append(zl)
            znames.append(nnn[zkcontrib])
        #
        zax.set_ylim((400, 0))
        zax.yaxis.set_ticks([0, 100, 200, 300, 400])
        zax.set_title('('+subnum.pop()+') '+zpart.upper()+'\n'+ttls[zres], loc='left')
        zax.vlines(0, 0, 400, color='grey')
        zax.annotate(note[zres], xy=(-1.4, 80), \
                     bbox=dict(boxstyle="round", fc="1"))
        return (zlines, znames)
    #

    ll = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zkr in rrr :
        icol = 0
        for zkpart in ['total', 'vert. adv.', 'meri. adv.'] :
            ll[irow, icol] = plotdata(ax[irow, icol], 'R'+zkr, zkpart)
            icol+=1
        #
        irow+=1
    #
    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1].flatten() : zax.tick_params(left=False)
    for zax in ax[:, 0].flatten() : zax.set_ylabel('Depth [m]')
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    for zax in ax[-1, :].flatten() : 
        zax.set_xlabel('mmol$\,$N$\,$m$^{-2}$d$^{-1}$')
        zax.set_xlim((-1.7, .7))
        zax.xaxis.set_ticks([ -1.4, -.7, 0, .7])
    #    
    fig.subplots_adjust(hspace=.5)

    hdl = ll[0, -1][0]
    nam = ll[0, -1][1]
    zw = ax[0, 1].get_position()
    legax = fig.add_axes([zw.x0, .91, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=5, handlelength=2., loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF7
##################################################

##################################################
# SF8
##################################################
if plot_sf8 : 
    print(">>> plot_sf8 <<<")

    savefig='sf8.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    sdate = '0166-01-01'
    edate = '0170-12-31'
    res = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', \
           '9', '27']
    lat1, lat2 = 35, 45

    data2plot = {}
    
    ####################
    # READ AND PROCESS NO3
    ####################

    fsuf = '_1y_01010101_01701230_ptrc_T.xml'
    var  = 'NO3'
    data2plot['NO3 GSP'] = {}
    data2plot['NO3 35N'] = {}
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['NO3 GSP'][zksim] = {}
        data2plot['NO3 35N'     ][zksim] = {}
        for zkres in res :
            fff = fdir + zksim + zkres + fsuf
            zw = reading.read_ncdf(var, fff, time = (sdate, edate))
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx')
            zw1 = averaging.ymean(zw, zwmesh, dim='zyx', ymin=lat1, ymax=lat2)
            zw1 = averaging.xmean(zw1, zwmesh, dim='zx')
            data2plot['NO3 GSP'][zksim]['R'+zkres]={'data':zw1, 'dep':zwmesh['depT']}
            zw1 = interpolating.yinterpol(zw, zwmesh, lat1, dim='zyx')
            zw1 = averaging.xmean(zw1, zwmesh, dim='zx')
            data2plot['NO3 35N'][zksim]['R'+zkres]={'data':zw1, 'dep':zwmesh['depT']}
        #
    #
    #___________________
    # DELTA
    data2plot['NO3 GSP']['DELTA'] = {}
    data2plot['NO3 35N']['DELTA'] = {}
    for zkres in res :
        data2plot['NO3 GSP']['DELTA']['R'+zkres]= \
            {'data':data2plot['NO3 GSP']['CC']['R'+zkres]['data'] - \
             data2plot['NO3 GSP']['CTL']['R'+zkres]['data'], \
             'dep':data2plot['NO3 GSP']['CTL']['R'+zkres]['dep']}
        data2plot['NO3 35N']['DELTA']['R'+zkres]= \
            {'data':data2plot['NO3 35N']['CC']['R'+zkres]['data'] - \
             data2plot['NO3 35N']['CTL']['R'+zkres]['data'], \
             'dep':data2plot['NO3 35N']['CTL']['R'+zkres]['dep']}
    #

    ####################
    # READ AND PROCESS W
    ####################

    fsuf = '_1y_01010101_01701230_grid_W.xml'
    var  = 'vovecrtz'
    zfact = 3600*24.
    data2plot['W GSP'] = {}
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['W GSP'][zksim] = {}
        for zkres in res :
            fff = fdir + zksim + zkres + fsuf
            zw = reading.read_ncdf(var, fff, time = (sdate, edate))
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx', grid='W')
            zw = averaging.ymean(zw, zwmesh, dim='zyx', ymin=lat1, ymax=lat2, grid='W')
            zw = averaging.xmean(zw, zwmesh, dim='zx', grid='W')
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
    zfact = 3600*24.
    data2plot['V 35N'] = {}
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['V 35N'][zksim] = {}
        for zkres in res :
            fff = fdir + zksim + zkres + fsuf
            zw = reading.read_ncdf(var, fff, time = (sdate, edate))
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            zw = zw['data']
            zw = averaging.tmean(zw, zwmesh, dim='tzyx', grid='V')
            zw = interpolating.yinterpol(zw, zwmesh, lat1, dim='zyx', grid='V')
            zw = averaging.xmean(zw, zwmesh, dim='zx', grid='V')
            data2plot['V 35N'][zksim]['R'+zkres]={'data':zfact*zw, 'dep':zwmesh['depV']}
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
    # READ AND PROCESS N2
    ####################

    fsuf = '_1y_01010101_01701230_grid_T.xml'
    varT = 'votemper'
    varS = 'vosaline'
    data2plot['N2 GSP'] = {}
    zfact = 1e5
    #___________________
    # CTL CC
    for zksim in ['CTL', 'CC'] :
        data2plot['N2 GSP'][zksim] = {}
        for zkres in res :
            if zkres in ['9', '27'] :
                zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+zkres+'.nc')
            else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
            # zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' ) # TEST
            fff = fdir + zksim + zkres + fsuf
            zw  = reading.read_ncdf(varT, fff, time = (sdate, edate))
            zwT = zw['data']
            zw  = reading.read_ncdf(varS, fff, time = (sdate, edate))
            zwS = zw['data']
            zw  = density.eos_bn2(zwT, zwS, zwmesh['e3w'])
            zw  = averaging.tmean(zw, zwmesh, dim='tzyx')
            zw  = averaging.ymean(zw, zwmesh, dim='zyx', ymin=lat1, ymax=lat2)
            zw  = averaging.xmean(zw, zwmesh, dim='zx')
            data2plot['N2 GSP'][zksim]['R'+zkres]={'data':zfact*zw, 'dep':zwmesh['depT']}
        #
    #
    #___________________
    # DELTA
    data2plot['N2 GSP']['DELTA'] = {}
    for zkres in res :
        data2plot['N2 GSP']['DELTA']['R'+zkres]= \
            {'data':data2plot['N2 GSP']['CC']['R'+zkres]['data'] - \
             data2plot['N2 GSP']['CTL']['R'+zkres]['data'], \
             'dep':data2plot['N2 GSP']['CTL']['R'+zkres]['dep']}
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
    fsize   = (14*infact, 8*infact) #(width, height)
    fig, ax = plt.subplots(2, 5, sharey='row', figsize=fsize) # nrow, ncol

    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()
    ttls={'NO3 GSP' : 'NO$_3$ GSP', \
          'NO3 35N' : 'NO$_3$ 35°N', \
          'W GSP'   : 'Upw. vel.', \
          'V 35N'   : 'Northw. vel.', \
          'N2 GSP'  : 'Stratif.'}
    units={'NO3 GSP' : '[mmol$\,$N$\,$m$^{-3}$]', \
           'NO3 35N' : '[mmol$\,$N$\,$m$^{-3}$]', \
           'W GSP'   : '[m$\,$d$^{-1}$]', \
           'V 35N'   : '[m$\,$d$^{-1}$]', \
           'N2 GSP'  : '[s$^{-2}$]'}
    nnn = {'R1'               :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=1e$^3$' , \
           'R1_KL500'         :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=500'       , \
           'R1_KL2000'        :'1°, k$_{gm}$=1e$^3$\nk$_{iso}$=2e$^3$' , \
           'R1_KGM500_KL500'  :'1°, k$_{gm}$=500\nk$_{iso}$=500'       , \
           'R1_KGM2000_KL2000':'1°, k$_{gm}$=2e$^3$\nk$_{iso}$=2e$^3$' , \
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
        if   zsim=='CTL'   : zax.set_title('('+subnum.pop()+') '+ttls[zvar], loc='left')
        elif zsim=='DELTA' : zax.set_title('('+subnum.pop()+') ∆ '+ttls[zvar], loc='left')
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
    zw = ax[0, 2].get_position()
    legax = fig.add_axes([zw.x0, .92, zw.width, .01])
    legax.axis('off')
    legax.legend( hdl, nam, ncol=3, handlelength=2., loc='lower center')
    
    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF8
##################################################

##################################################
# SF9
##################################################

if plot_sf9 : 
    print(">>> plot_sf9 <<<")

    savefig = 'sf9.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf='_1m_01010101_01701230_grid_T.xml'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    data2plot={}

    ####################
    # PREPARE MAPS AND DISTRIBUTION OF MLD
    ####################

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        lon, lat = zwmesh['lonT'], zwmesh['latT']
        data2plot['R'+vres] = { 'MAP':{'lon':lon, 'lat':lat}, 'DIS':{} }
        zw = {}
        for vsim in ['CTL', 'CC'] : 
            fff = fdir + vsim + vres + fsuf
            zw2  = reading.read_ncdf('somxl010', fff, time=(sdate, edate))
            zw2  = zw2['data']
            #___________________
            # MAPS OF MAXIMUM MLD
            zw3 = []
            for yy in range(5) : zw3.append(np.nanmax(zw2[yy:yy+12], axis = 0))
            zw2       = np.nanmean(zw3, axis=0)
            zw[vsim]  = zw2
            #___________________
            # SELECT MLD VALUES FOR DISTRIBUTION
            zw2 = np.where( lat>=45 , zw2, float('nan') )
            data2plot['R'+vres]['DIS'][vsim] = zw2
        #
        data2plot['R'+vres]['MAP']['CTL']   = zw['CTL']
        data2plot['R'+vres]['MAP']['DELTA'] = zw['CC'] - zw['CTL']
    #

    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*3*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharex='col',figsize=fsize)

    ttls={'R1'               :'k$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'k$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'k$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'k$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'k$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zaxrow, zres) :
        zdat=data2plot[zres]
        #___________________
        # MAPS
        X, Y = zdat['MAP']['lon'], zdat['MAP']['lat']
        # CTL
        lev = np.linspace(0, 800, 11)
        zcfctl = zaxrow[0].contourf(X, Y, zdat['MAP']['CTL'], levels=lev, cmap='plasma_r', extend='max')
        ttt='('+subnum.pop()+') '+ttls[zres]
        plotting.make_XY_axis(zaxrow[0], title=ttt)
        zaxrow[0].label_outer()
        zaxrow[0].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        # DELTA
        aa = np.linspace(-400, 0, 6)
        bb = np.linspace(0, 400, 6)
        levD = np.concatenate((aa[:-1], bb))
        zcfdelta = zaxrow[1].contourf(X, Y, zdat['MAP']['DELTA'], levels=levD, cmap='RdBu', extend='both')
        ttt='('+subnum.pop()+') '+ttls[zres]
        plotting.make_XY_axis(zaxrow[1], title=ttt)
        zaxrow[1].label_outer()
        zaxrow[1].tick_params(left=False)
        zaxrow[1].annotate(note[zres], xy=(-83, 23), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        # DISTRIBUTIONS
        kwctl = {'bins':50, 'range':(1, 2501), 'density':True, 'color':'royalblue'}
        kwcc  = {'bins':50, 'range':(1, 2501), 'density':True, 'color':'firebrick', 'rwidth':.5}
        zlhctl = zaxrow[2].hist(zdat['DIS']['CTL'].flatten(), **kwctl)
        zlhcc  = zaxrow[2].hist(zdat['DIS']['CC'].flatten() , **kwcc )
        ttt='('+subnum.pop()+') '+ttls[zres]
        zaxrow[2].set_ylabel('Probability')
        zaxrow[2].yaxis.set_ticks_position('right')
        zaxrow[2].yaxis.set_label_position('right')
        zaxrow[2].set_title(ttt, loc='left')
        zaxrow[2].set_ylim((0, 0.007))
        zaxrow[2].locator_params(axis='y', nbins=4)
        zaxrow[2].annotate(note[zres], xy=(1600, .005), \
                        bbox=dict(boxstyle="round", fc="1"))
        #___________________
        return (zcfctl, zcfdelta, zlhctl, zlhcc)
    #

    cf = np.zeros_like(ax[:, 0])
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    # rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '1', '1']# TEST
    for zres in rrr :
        cf[irow] = plotdata(ax[irow], 'R'+zres)
        irow+=1
    #
    ax[-1, -1].set_xlabel('Metres')

    cbtitle = 'MLD, CTL simulation [m]'
    zw = ax[0, 0].get_position()
    ny0 = zw.y0+zw.height*1.4
    cbar_ax = fig.add_axes([zw.x0, ny0, zw.width, 0.01])
    fig.colorbar(cf[0][0], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[0, 400, 800])
    cbtitle = 'MLD change [m]'
    zw = ax[0, 1].get_position()
    ny0 = zw.y0+zw.height*1.4
    cbar_ax = fig.add_axes([zw.x0, ny0, zw.width, 0.01])
    fig.colorbar(cf[0][1], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[-400, 0, 400])

    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)

    hdl=(cf[0][-1][2][0], cf[0][-2][2][0])
    nam=('CTL', 'CC')  
    zw = ax[0, -1].get_position()
    legax = fig.add_axes([zw.x0, zw.y0+zw.height*1.4, zw.width, .01])
    legax.axis('off')
    leg=legax.legend( hdl, nam, loc='lower center', title='Max MLD distribution', ncol=2)

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF9
##################################################

##################################################
# SF10
##################################################

if plot_sf10 : 
    print(">>> plot_sf10 <<<")

    savefig = 'sf10.png'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_grid_T.xml'
    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}
    
    ####################
    # PREPARE N2
    ####################

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        data2plot['R'+vres]={'lat':zwmesh['latT'][:, 5], 'dep':zwmesh['depT']}
        for vsim in ['CTL', 'CC'] : 
            fff = fdir + vsim + vres + fsuf
            zw  = reading.read_ncdf('votemper', fff, time=(sdate, edate))
            zwT = zw['data']
            zw  = reading.read_ncdf('vosaline', fff, time=(sdate, edate))
            zwS = zw['data']
            zw  = density.eos_bn2(zwT, zwS, zwmesh['e3w'])
            zw = averaging.tmean(zw, zwmesh) 
            zw = averaging.xmean(zw, zwmesh, dim='zyx') 
            zfact = 1e5
            data2plot['R'+vres][vsim] = zw*zfact
        #
        data2plot['R'+vres]['DELTA'] = data2plot['R'+vres]['CC'] - data2plot['R'+vres]['CTL']
    #
    
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*4*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    aa  = np.linspace(-5, 0, 6)
    bb  = np.linspace(0, 5, 6)
    levD = np.concatenate((aa[:-1], bb))
    levC = np.linspace(0, 30., 11)
    lev  = {'CTL':levC, 'CC':levC, 'DELTA':levD}
    cmap = {'CTL':'viridis', 'CC':'viridis', 'DELTA':'RdBu_r'}
    ext  = {'CTL':'max', 'CC':'max', 'DELTA':'both'}
    ttls={'R1'               :'\nk$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'\nk$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'\nk$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'\nk$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'\nk$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zres, zsim) :
        zdat=data2plot[zres]
        X, Y, Z = zdat['lat'], zdat['dep'], zdat[zsim]
        zcf = zax.contourf(X, Y, Z, levels=lev[zsim], cmap=cmap[zsim], extend=ext[zsim])
        ttt='('+subnum.pop()+') '+zsim+' '+ttls[zres]
        plotting.make_YZ_axis(zax, depmax=800, title=ttt)
        zax.label_outer()
        zax.annotate(note[zres], xy=(23, 700), \
                     bbox=dict(boxstyle="round", fc="1"))
        return zcf
    #

    cf = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    # rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '1', '1'] # TEST
    for zres in rrr :
        icol = 0
        for zsim in ['CTL', 'CC', 'DELTA'] :
            cf[irow, icol] = plotdata(ax[irow, icol], 'R'+zres, zsim)
            icol+=1
        #
        irow+=1
    #


    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1:].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #
    fig.subplots_adjust(hspace=.4)

    cbtitle = 'Brunt-Vaiasala frequency [1e$^{-5}$s$^{-1}$]'
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 1].get_position()
    nwidth = zw1.width*1.5
    nx0 = (zw1.x0+zw1.width+zw2.x0)/2 - nwidth/2
    ny0 = zw1.y0+zw1.height*1.4
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.01])
    fig.colorbar(cf[0,0], cax=cbar_ax, orientation='horizontal', ticklocation='top', label=cbtitle)

    cbtitle = 'Change [1e$^{-5}$s$^{-1}$]'
    zw = ax[0, -1].get_position()
    ny0 = zw.y0+zw.height*1.4
    cbar_ax = fig.add_axes([zw.x0, ny0, zw.width, 0.01])
    fig.colorbar(cf[0, -1], cax=cbar_ax, orientation='horizontal', ticklocation='top', \
                 label=cbtitle, ticks=[-5, 0, 5])

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig)
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF10
##################################################

##################################################
# SF11
##################################################
if plot_sf11 : 
    print(">>> plot_sf11 <<<")

    savefig='sf11.png'
    
    ####################
    # PARAM
    ####################

    fdir = '/gpfswork/rech/eee/rdyk004/MY_PYTHON3/PCKL/'

    inctl = 'main_no3_budget_horizontal_ctl27_z400m_y1'
    incc  = 'main_no3_budget_horizontal_cc27_z400m_y1'
    inctl_m = 'main_no3_budget_horizontal_ctl27_z400m_onR1_1ymean.pckl'
    incc_m  = 'main_no3_budget_horizontal_cc27_z400m_onR1_1ymean.pckl'
    avgfactor = 27

    ####################
    # PREPARE DATA
    ####################

    def open_pckl(zwin, zwsim) : 
        fff = open(fdir+zwin, 'rb')
        zw = pickle.load(fff)
        return zw[zwsim]
    #     
    
    def open_mul_pckl(zinpckl, zsimname) :
        zw={}
        var=['strd']
        for vvv in var : zw[vvv] = []
        for yyy in np.arange(66, 71) : 
            zname = zinpckl + str(yyy) + '.pckl'
            zw2 = open_pckl(zname, zsimname)
            for vvv in var : 
                zw[vvv].append(zw2[vvv])
            #
        #
        zw3 = {}
        zw3['mesh'] = zw2['mesh']
        for vvv in var : zw3[vvv] = np.ma.array(np.mean(zw[vvv], axis = 0))
        return zw3
    #
    def avg_onRX(zin, zoutmesh, zfactor) :
        zfactor = int(zfactor)
        zw = zin[1:-1, 1:-1] # remove edge
        # only on T grid
        zout = np.zeros(( zw.shape[0]//zfactor + 2, zw.shape[1]//zfactor + 2 ))
        # reshape to get subarray of shape zfactor*zfactor
        zw = np.reshape(zw, (zout.shape[0]-2, zfactor, zout.shape[1]-2, zfactor))
        zw = np.nanmean(zw, axis = (1, 3)) # mean the subarray
        zout[1:-1, 1:-1] = zw # fill output
        # fill corner
        zout[ 0,  0] = zin[ 0,  0]
        zout[ 0, -1] = zin[ 0, -1]
        zout[-1,  0] = zin[-1,  0]
        zout[-1, -1] = zin[-1, -1]
        # fill edge
        zout[1:-1,    0] = zin[zfactor//2+1::zfactor,                   0]
        zout[1:-1,   -1] = zin[zfactor//2+1::zfactor,                  -1]
        zout[   0, 1:-1] = zin[                  0, zfactor//2+1::zfactor]
        zout[  -1, 1:-1] = zin[                 -1, zfactor//2+1::zfactor]
        zout = np.ma.array(zout)
        zout.mask=[x!=1 for x in zoutmesh['tmask']]
        return zout
    #
    data2plot={}
    zw = open_pckl(inctl_m, 'CTL27')
    data2plot['CTL MEAN'] = zw['strd']
    data2plot['lat'] = zw['mesh']['latT']
    data2plot['lon'] = zw['mesh']['lonT']
    zwmesh = zw['mesh']
    zw = open_pckl(incc_m, 'CC27')
    data2plot['CC MEAN'] = zw['strd']
    zw = open_mul_pckl(inctl, 'CTL27')
    zw = zw['strd']
    data2plot['CTL TOTAL'] = avg_onRX(zw, zwmesh, 27)
    zw = open_mul_pckl(incc, 'CC27')
    zw = zw['strd']
    data2plot['CC TOTAL'] = avg_onRX(zw, zwmesh, 27)
    data2plot['CTL EDDY'] = data2plot['CTL TOTAL'] - data2plot['CTL MEAN']
    data2plot['CC EDDY']  = data2plot['CC TOTAL']  - data2plot['CC MEAN']

    for kkk in ['TOTAL', 'MEAN', 'EDDY'] : data2plot['DELTA '+kkk] = \
        data2plot['CC '+kkk]  - data2plot['CTL '+kkk]

    ####################
    # POT
    ####################

    infact  = 1/2.54
    ncol, nrow = 3, 3
    fsize = (ncol*5*infact, nrow*5*infact) #(width, height)
    fig, ax   = plt.subplots(nrow, ncol, figsize=fsize, sharex='col', sharey='row') # nrow, ncol

    X, Y = data2plot['lon'], data2plot['lat']
    aa = np.linspace(-4, 0, 6)
    bb = np.linspace(0, 4, 6)
    lev = np.concatenate((aa[:-1], bb))
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zplot) :
        zfact = 3600*24. # mmol/m2/s -> mmol/m2/d
        Z = zfact * data2plot[zplot]
        zcf = zax.contourf(X, Y, Z, levels=lev, cmap='PiYG', extend='both')
        ttt='('+subnum.pop()+') '+zplot
        plotting.make_XY_axis(zax, title=ttt)
        zax.label_outer()
        return zcf
    #

    cf = np.zeros_like(ax)
    irow = 0
    for vsim in ['CTL', 'CC', 'DELTA'] :
        icol = 0
        for vpart in ['TOTAL', 'MEAN', 'EDDY'] :
            cf[irow, icol] = plotdata(ax[irow, icol], vsim + ' ' +vpart)
            icol+=1
        #
        irow+=1
    #
        
    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1:].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Northward km')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #

    cbtitle = 'NO$_3$ adv. fluxes divergence [mmol$\,$N$\,$m$^{-2}$d$^{-1}$]'
    zw0 = ax[0, 0].get_position()
    zw2 = ax[0, 2].get_position()
    nx0 = zw0.x0 + 0.5*zw0.width
    nwidth = zw2.x0 + 0.5*zw2.width - nx0
    ny0 = zw0.y0+zw0.height*1.2
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.02])
    fig.colorbar(cf[0, 0], cax=cbar_ax, orientation='horizontal', \
                 label=cbtitle, ticklocation='top')

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    print("Figure saved: ", dirfig+suffig+savefig)
    fig.savefig(dirfig+suffig+savefig)
    plt.close()

##################################################
# END SF11
##################################################

##################################################
# SF12
##################################################
if plot_sf12 :
    print(">>> plot_sf12 <<<")

    savefig = 'sf12.png'
    savefile='sf12.txt'
    fdir  = '/gpfswork/rech/eee/rdyk004/GYRE_XML/'          
    fsuf  = '_1y_01010101_01701230_grid_V.xml'
    sdate  = '0166-01-01'
    edate  = '0170-12-31'
    r1s = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000']
    res = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '9', '27']
    # res  = ['1', '1_KL2000', '1_KL500', '1_KGM500_KL500', '1_KGM2000_KL2000', '1', '1'] # TEST
    data2plot={}
    
    ####################
    # PREPARE MSF
    ####################

    for vres in res :
        if vres in ['9', '27'] : zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R'+vres+'.nc')
        else :  zwmesh = reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )
        data2plot['R'+vres]={'lat':zwmesh['latW'][:, 5], 'dep':zwmesh['depW']}
        for vsim in ['CTL', 'CC'] : 
            fff = fdir + vsim + vres + fsuf
            zw  = reading.read_ncdf('vomecrty', fff, time=(sdate, edate))
            zw = zw['data']
            if vres in r1s : 
                zw1  = reading.read_ncdf('vomeeivy', fff, time=(sdate, edate))
                zw += zw1['data']
            #
            zw = stream_functions.msf(zw, zwmesh, dim='tzyx')
            zw = averaging.tmean(zw, zwmesh, dim='tzy') 
            data2plot['R'+vres][vsim] = zw
        #
        data2plot['R'+vres]['DELTA'] = data2plot['R'+vres]['CC'] - data2plot['R'+vres]['CTL']
    #
    
    ####################
    # AVG1, STD1
    ####################
    
    avg1={'mesh':reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )}
    std1={'mesh':reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R1.nc' )}
    for vsim in ['CTL', 'CC', 'DELTA'] : 
        zw= []
        for vres in r1s : zw.append(data2plot['R'+vres][vsim])
        zw = np.array(zw)
        avg1[vsim] = np.ma.array( data=np.nanmean(zw, axis=0),\
                                  mask=data2plot['R1']['CTL'].mask )
        std1[vsim] = np.ma.array( data=unc*np.nanstd(zw, axis=0),\
                                  mask=data2plot['R1']['CTL'].mask )
    #

    dirfile = '/gpfswork/rech/eee/rdyk004/MY_PYTHON3/'
    suffile = 'main_article_gyre_pp_supfig_v2_'
    fff = open(dirfile+suffile+savefile, "w+")

    fff.write("\n##################################################")
    fff.write("\n##################################################")
    fff.write("\n##                                              ##")
    fff.write("\n##    Meri. flux at 35N, 100-400 metres depth   ##")
    fff.write("\n##                                              ##")
    fff.write("\n##################################################")
    fff.write("\n##################################################")

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                      AVG1 pm STD1                    ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=avg1['mesh']
        zw1    = interpolating.yinterpol(avg1[vsim], zwmesh, 35, dim='zy')
        zw2    = np.around(interpolating.zinterpol(zw1, zwmesh, 100, dim='z'), decimals=2)
        zw3    = np.around(interpolating.zinterpol(zw1, zwmesh, 400, dim='z'), decimals=2)
        zw1    = interpolating.yinterpol(std1[vsim], zwmesh, 35, dim='zy')
        zw2std = np.around(interpolating.zinterpol(zw1, zwmesh, 100, dim='z'), decimals=2)
        zw3std = np.around(interpolating.zinterpol(zw1, zwmesh, 400, dim='z'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4  = np.float(zw2), np.float(zw3), np.float(zw4)
        zw2std, zw3std = np.float(zw2std), np.float(zw3std)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('100m: '+str(zw2)+' pm '+str(zw2std)+'   ---   400m: '+\
                  str(zw3)+' pm '+str(zw3std)+'   ---   diff: '+str(zw4))
    #

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                            R9                        ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R9.nc' )
        zdat=data2plot['R9']
        zw1 = interpolating.yinterpol(zdat[vsim], zwmesh, 35, dim='zy')
        zw2 = np.around(interpolating.zinterpol(zw1, zwmesh, 100, dim='z'), decimals=2)
        zw3 = np.around(interpolating.zinterpol(zw1, zwmesh, 400, dim='z'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4 = np.float(zw2), np.float(zw3), np.float(zw4)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('100m: '+str(zw2)+'   ---   400m: '+str(zw3)+'   ---   diff: '+str(zw4))
    #

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                           R27                        ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R27.nc' )
        zdat=data2plot['R27']
        zw1 = interpolating.yinterpol(zdat[vsim], zwmesh, 35, dim='zy')
        zw2 = np.around(interpolating.zinterpol(zw1, zwmesh, 100, dim='z'), decimals=2)
        zw3 = np.around(interpolating.zinterpol(zw1, zwmesh, 400, dim='z'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4 = np.float(zw2), np.float(zw3), np.float(zw4)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('100m: '+str(zw2)+'   ---   400m: '+str(zw3)+'   ---   diff: '+str(zw4))
    #

    fff.write("\n\n")
    fff.write("\n##################################################")
    fff.write("\n##################################################")
    fff.write("\n##                                              ##")
    fff.write("\n##    Vert. flux at 400 metres depth, 35-45N    ##")
    fff.write("\n##                                              ##")
    fff.write("\n##################################################")
    fff.write("\n##################################################")

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                      AVG1 pm STD1                    ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=avg1['mesh']
        zw1    = interpolating.zinterpol(avg1[vsim], zwmesh, 400, dim='zy')
        zw2    = np.around(interpolating.yinterpol(zw1, zwmesh, 35, dim='y'), decimals=2)
        zw3    = np.around(interpolating.yinterpol(zw1, zwmesh, 45, dim='y'), decimals=2)
        zw1    = interpolating.zinterpol(std1[vsim], zwmesh, 400, dim='zy')
        zw2std = np.around(interpolating.yinterpol(zw1, zwmesh, 35, dim='y'), decimals=2)
        zw3std = np.around(interpolating.yinterpol(zw1, zwmesh, 45, dim='y'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4  = np.float(zw2), np.float(zw3), np.float(zw4)
        zw2std, zw3std = np.float(zw2std), np.float(zw3std)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('35N: '+str(zw2)+' pm '+str(zw2std)+'   ---   45N: '+\
                  str(zw3)+' pm '+str(zw3std)+'   ---   diff: '+str(zw4))
    #

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                            R9                        ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R9.nc' )
        zdat=data2plot['R9']
        zw1 = interpolating.zinterpol(zdat[vsim], zwmesh, 400, dim='zy')
        zw2 = np.around(interpolating.yinterpol(zw1, zwmesh, 35, dim='y'), decimals=2)
        zw3 = np.around(interpolating.yinterpol(zw1, zwmesh, 45, dim='y'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4 = np.float(zw2), np.float(zw3), np.float(zw4)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('35N: '+str(zw2)+'   ---   45N: '+str(zw3)+'   ---   diff: '+str(zw4))
    #

    fff.write("\n")
    fff.write("\n##################################################\n")
    fff.write("                           R27                        ")
    fff.write("\n##################################################\n")
    for vsim in ['CTL', 'CC', 'DELTA'] :
        zwmesh=reading.read_mesh('/gpfswork/rech/eee/rdyk004/MESH/mesh_mask_R27.nc' )
        zdat=data2plot['R27']
        zw1 = interpolating.zinterpol(zdat[vsim], zwmesh, 400, dim='zy')
        zw2 = np.around(interpolating.yinterpol(zw1, zwmesh, 35, dim='y'), decimals=2)
        zw3 = np.around(interpolating.yinterpol(zw1, zwmesh, 45, dim='y'), decimals=2)
        zw4 = np.around(zw3-zw2, decimals=2)
        zw2, zw3, zw4 = np.float(zw2), np.float(zw3), np.float(zw4)
        fff.write('\n---------- '+vsim+' ----------\n')
        fff.write('35N: '+str(zw2)+'   ---   45N: '+str(zw3)+'   ---   diff: '+str(zw4))
    #

    fff.write('\n')
    fff.close()
    print("File saved: ", dirfile+suffile+savefile)
    
    ####################
    # PLOT DATA
    ####################

    infact  = 1/2.54
    nrow, ncol = 7, 3
    fsize   = (ncol*4*infact, nrow*3*infact) #(width, height)
    fig, ax = plt.subplots(nrow, ncol, sharey='row', sharex='col', figsize=fsize)

    aa  = np.linspace(-5, 0, 6)
    bb  = np.linspace(0, 5, 6)
    lev = np.concatenate((aa[:-1], bb))
    cmap='RdBu'
    ext = 'both'
    ttls={'R1'               :'\nk$_{gm}$=1e$^3$, k$_{iso}$=1e$^3$', \
          'R1_KL500'         :'\nk$_{gm}$=1e$^3$, k$_{iso}$=500'      , \
          'R1_KL2000'        :'\nk$_{gm}$=1e$^3$, k$_{iso}$=2e$^3$', \
          'R1_KGM500_KL500'  :'\nk$_{gm}$=500, k$_{iso}$=500'      , \
          'R1_KGM2000_KL2000':'\nk$_{gm}$=2e$^3$, k$_{iso}$=2e$^3$', \
          'R9':'', 'R27':''}
    note={'R1'               :'1°', \
          'R1_KL500'         :'1°', \
          'R1_KL2000'        :'1°', \
          'R1_KGM500_KL500'  :'1°', \
          'R1_KGM2000_KL2000':'1°', \
          'R9':'1/9°', 'R27':'1/27°'}
    subnum=list('abcdefghijklmnopqrstuvwxyz')
    subnum.reverse()

    def plotdata(zax, zres, zsim) :
        zdat=data2plot[zres]
        X, Y, Z = zdat['lat'], zdat['dep'], zdat[zsim]
        zcf = zax.contourf(X, Y, Z, levels=lev, cmap=cmap, extend=ext)
        ttt='('+subnum.pop()+') '+zsim+' '+ttls[zres]
        plotting.make_YZ_axis(zax, depmax=800, title=ttt)
        zax.vlines(35, 100, 400, color='k', linestyle='-', linewidth=1.)
        zax.label_outer()
        zax.annotate(note[zres], xy=(23, 700), \
                     bbox=dict(boxstyle="round", fc="1"))
        return zcf
    #

    cf = np.zeros_like(ax)
    irow = 0
    rrr = ['1_KGM500_KL500', '1_KL500', '1' ,'1_KL2000','1_KGM2000_KL2000', '9', '27']
    for zres in rrr :
        icol = 0
        for zsim in ['CTL', 'CC', 'DELTA'] :
            cf[irow, icol] = plotdata(ax[irow, icol], 'R'+zres, zsim)
            icol+=1
        #
        irow+=1
    #


    for zax in ax[:-1, :].flatten() : zax.tick_params(bottom=False)
    for zax in ax[:, 1:].flatten() : zax.tick_params(left=False)
    for zax in ax[:, -1].flatten() :
        zax.set_ylabel('Depth [m]')
        zax.yaxis.set_ticks_position('right')
        zax.yaxis.set_label_position('right')
    #
    fig.subplots_adjust(hspace=.4)

    cbtitle = 'Meridional stream function [Sv]'
    zw1 = ax[0, 0].get_position()
    zw2 = ax[0, 2].get_position()
    nx0 = zw1.x0+zw1.width
    nwidth = zw2.x0 - nx0
    ny0 = zw1.y0+zw1.height*1.4
    cbar_ax = fig.add_axes([nx0, ny0, nwidth, 0.01])
    fig.colorbar(cf[0,0], cax=cbar_ax, orientation='horizontal', ticklocation='top', label=cbtitle)

    for zax in fig.axes : zax.tick_params(color='dimgrey')
    dirfig='/gpfswork/rech/eee/rdyk004/MY_FIG_PYTHON3/'
    suffig='main_article_gyre_pp_supfig_v2_'
    fig.savefig(dirfig+suffig+savefig) 
    print("Figure saved: ", dirfig+suffig+savefig)
    plt.close()
#
##################################################
# END SF12
##################################################

