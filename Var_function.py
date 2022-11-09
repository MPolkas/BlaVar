#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 18:52:18 2020

@author: markos polkas
"""

import os
import numpy as np
from scipy import interpolate as inpl
from scipy import integrate as intg
from scipy import stats as st
from scipy.optimize import curve_fit
import scipy.optimize as optmz
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from  matplotlib import pyplot as plt, image as mpimg, gridspec, colors 
from astropy.io import ascii
import sys



def load_input_information_on_objects(filename='nofile'):

    global objects,redsh, Dists, Rbb, SSC, BBextra, Gamma2_arr, factbb_arr, Tbb2_arr, Rbb2_arr
    #Bibliograpgical Data for objects
    objects=['3C273','3C279','PKS2155304','Mrk421']
    redsh=[0.158,0.536,0.0, 0.031]
    Rbb=[6.3e17,5e17,10e12, 1e12] #Radii of the BB radiation zone #UNITS;cm
    SSC=[False,False,True, True]
    BBextra=[True,False,False,False]
    #Second BB into the model
    Gamma2_arr=[1,1,1, 1] #Deboosting/Boosting for observation
    factbb_arr=[10,500,'-', '-'] #(L_disk-L_BLR)/L_BLR
    Tbb2_arr=[11000,20000,'-', '-'] #T_disk #UNITS; K
    Rbb2_arr=[10**15.89,10**15.62,'-', '-'] #disk radii #UNITS;cm
    
    #allows for an input of 8 x number of objects .txt file if you have a list of objects
    if filename!='nofile':
        redsh, Rbb, SSC, BBextra, Gamma2_arr, factbb_arr, Tbb2_arr, Rbb2_arr = np.loadtxt(filename).T
    
    from astropy.cosmology import FlatLambdaCDM 
    cosmo=FlatLambdaCDM(H0=70, Om0=0.3)
    Dists=cosmo.luminosity_distance(redsh).value #UNITS;Mpc
    del cosmo
    del FlatLambdaCDM
    


"""Constants"""
def load_constants():
    global G, fortio, c,h,kb,hbar,eV,me,mp,sigmaT,aBB,aBBeV,sigmaBB,Jy,pc,Msun,Rsun,Lsun,Ledd
    G=6.67259*10**(-8)
    fortio=4.80*10**(-10)
    c = 2.997925*10**10
    h= 6.625*10**(-27)
    kb=1.38*10**(-16)
    hbar=h/2/3.1415
    eV = 1.602192*10**(-12)
    me=9.109*10**(-28)
    mp = 1.672614*10**(-24)
    sigmaT=0.665245*10**(-24)
    aBB = 7.564*10**(-15)
    aBBeV = 7.564*10**(-15)/eV
    sigmaBB=aBB*c/4
    Jy=10**-23
    pc = 3.0856*10**18  
    Msun = 1.989*10**33
    Rsun = 6.96*10**10
    Lsun=3.9*10**33
    Ledd=3.2*10**4*Lsun*10**8 #M/10^8 Msun

"""Functions"""
def logparabola_g(x, logA, alpha, E_eV= 500*1e9 , beta=0.7):
    x0 = E_eV  * (eV/h)
    return 10**logA*(x/x0)**(alpha-beta*np.log(x/x0)+2.0)

def logparabola_x(x, logA, alpha, E_eV= 1e3 , beta=0.38):
    x0 = E_eV  * (eV/h)
    return 10**logA* (x/x0)**(alpha-beta*np.log(x/x0)+2.0)


def surface_density(m1, m2,lims):
        if lims:
            xmin, xmax, ymin, ymax= lims
        else:
            xmin, xmax, ymin, ymax=min(m1), max(m1), min(m2), max(m2)
        X, Y = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]                                                     
        positions = np.vstack([X.ravel(), Y.ravel()])                                                       
        values = np.vstack([m1, m2])                                                                        
        kernel = st.gaussian_kde(values)                                                                 
        Z = np.reshape(kernel(positions).T, X.shape)
        return X, Y, Z
    
def BB(x,T,lext,Rbb):
        Ivbb=3*x+np.log10((2*me**3*c**4/h**2)/(np.exp(10**x*me*c**2/kb/T)-1)+10**-100)
        opacity=(lext*me*c**2/sigmaT/R)/aBB/T**4
        #make it dimensionelss / scaling
        Ibb=x+np.log10(me*c**2/h)+Ivbb+np.log10(opacity*4*np.pi*Rbb**2)\
            +np.log10(sigmaT/(4*np.pi*R*me*c**3))
        return Ibb
    
def surface_density_per_energy(m1, m2,lims):
            if lims:
                #xmin, xmax=lims[0:2]
                xmin, xmax, =min(m1), max(m1),
                ymin, ymax= lims[2::]
            else:
                xmin, xmax, ymin, ymax=min(m1), max(m1), min(m2), max(m2)
            X, Y = np.mgrid[xmin:xmax:2j, ymin:ymax:50j]
            hist=np.histogram(m2,bins=list(Y[0])+[ymax],density=True)[0]
            Z =[hist,hist] #uniform in x-axis
            return X, Y, Z        
        
#FIT POWER LAW
def line(x,a):
    return a*x

def corrplot_withspots(xarr, yarr, tobs, axmain, axcbar=None, labs=['',''],
                       title ='',  clr = 'b', XLAB=True , YLAB=True, lvls=[],
                       SPOTS = False, LW=1.5, mscorr=2.0, color='k', corrfont = 14.0, 
                       color_continuous = True, color_dic = {}, axcorr=[]):
    #CORRELATION X-rays VHEband
    CONTOUR=False
    if len(lvls)!=0: CONTOUR=True
    xlab , ylab= labs[0] , labs[1]
    axmain.plot(np.log10(xarr),np.log10(yarr),'-',color=clr, linewidth=LW, alpha=0.3) 
    for i,t in enumerate(tobs):
        color=(0.0+max(0.0,abs(tobs[i]-tobs[0])/(tobs[-1]-tobs[0])),  
                      0.75,
                       0.75-min(0.75,abs(tobs[i]-tobs[0])/(tobs[-1]-tobs[0])))
        if not color_continuous: color = color_dic[str(int(max(tobs[i], 56393)))]
        axmain.plot(np.log10(xarr[i]),np.log10(yarr[i]),'o', color = color, 
                  ms=mscorr,zorder=0) #oldci
    vals = np.ones((24,4))*0.75
    vals[:,0] = np.linspace(0.0, 0.75,24)
    vals[:,2] = np.linspace(0.75, 0.0, 24)
    if axcbar!=None and color_continuous:
        cbar=plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=tobs[0],vmax=tobs[-1]), cmap=colors.ListedColormap(vals)),\
            ax=[axcbar],location='left', extend='both')
        axcbar.annotate(r'$t_{obs}-t_0 [days]$', (0.25, 0.3), rotation=90,size=15)
        axcbar.get_xaxis().set_visible(False)
        axcbar.get_yaxis().set_visible(False)
        axcbar.spines['top'].set_visible(False)
        axcbar.spines['right'].set_visible(False)
        axcbar.spines['bottom'].set_visible(False)
        axcbar.spines['left'].set_visible(False)
        axcbar.tick_params(labelsize=10)
    incl=optmz.curve_fit(line,np.log10(xarr/np.mean(xarr)),np.log10(yarr/np.mean(yarr)))[0][0]
    aa,bb,cc,dd=[axmain.get_xlim()[0],axmain.get_xlim()[1],axmain.get_ylim()[0],axmain.get_ylim()[1]]
    ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
    axlim1, axlim2 = min(np.mean(np.log10(xarr)),np.mean(np.log10(yarr)))-3.0 , min(np.mean(np.log10(xarr)),np.mean(np.log10(yarr))) +3.0
    axmain.plot(np.linspace(axlim1,axlim2),np.linspace(axlim1,axlim2),'k-', label= r'lin' )
    axmain.plot(np.linspace(axlim1,axlim2),incl*(np.linspace(axlim1,axlim2)+9.5)-9.5,'k-.',label=r'p='+str(round(incl,2)))
    if CONTOUR:      
        X,Y,Z=surface_density(np.log10(xarr),np.log10(yarr),[])
        axmain.contour(X,Y,Z,colors=['k']*len(lvls),linestyles='dashed',\
                              linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
    axmain.legend(fontsize='small')
    axmain.set_xlabel(xlab,fontsize=corrfont)
    axmain.set_ylabel(ylab,fontsize=corrfont)
    if len(title):
        axmain.set_title(title,fontsize=corrfont+2)
    #axmain.axis([-11.0, -8.0 , -11.0 , -8.0])
    axmain.xaxis.set_minor_locator(AutoMinorLocator())
    axmain.yaxis.set_minor_locator(AutoMinorLocator())
    axmain.tick_params(axis='y',labelsize=corrfont-3)
    axmain.tick_params(axis='x',labelsize=corrfont-3)
    if len(axcorr)==0: axcorr= [ ae, be ,ce ,de]
    axmain.axis(axcorr)
    #    #plot scale errors
    #    xerr, yerr =  errx, errg
    #    plt.errorbar(axcorr[0]+0.3*(axcorr[1]-axcorr[0]),axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #    plt.annotate('obs error',(axcorr[0]+0.1*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr))

all_mjd_days = list(np.linspace(5e4,7e4,int(2e4+1)).astype(int))
#import scipy.stats as stats
"""%%%%%%%%%%%%%%%%  SETTINGS  %%%%%%%%%%%%%%%%%%% """
def Var_analysisVar_analysis(obj, route_save, ALLTOOLS=False , SAVE=True,   imtyp='.png', ELECTRONS=False,
             SED=  'on', TIMELAPSE=False ,  #SED 'on'(lines)/'con'(contour)/'no'(skip)
             TIMECURVES=True,  CORR=    False ,TIMELAPSE_CORR=False, #plot timecurves,correlation plot(multiple frames)
             TIMELAGS=False, DCFcalc= False, SFcalc= False , STINGRAY=False, #time-lag analysis
             CALC_THEO_INDICES=True,   #calculate the rate of change on the parameter provided, 
                                       # to study the effectiveness of your POWER selection
             STATIONARY = [False, all_mjd_days, 'lowstate' ,  50.0, 4.9e15, './','None'],
	     EXTRA =  [False, [] , 'extrastate',150.0, 4.9e15/2**5.0,'./','None'],   feedback_VHE_on_extrastate =   False,
	     TIMEDEPENDENT_DATA_ON_PLOT = [True, './data/seds/']):
    load_constants()
    #Please, walk through settings down to END SETTINGS for selecting the parameters for each function of the script
    #The user is requeste to have multiple one fort.85 or multiple fort_*index*.85 
    # under the directory with name /'NameObject'/'Namevar'_... for the current script to operate
    #Function Switches , for parameters of functions see below
    #ALLTOOLS=False #except SFcalc
    #SED=     'on' #'no', 'on'(overplot), 'con' (contour)  vFv Diagram [CAUTION!]only up to *limSEDs*, set below
    #TIMECURVES=True #plot TC diagrams and Color(t)
    #CORR=    True #SWITCH for flux flux and Color diagrams
    #QUICK_TIMELAGS=False  #Quick method to calculate timelag from flares (#PEAKS)
    #DCFcalc= False #Calculate or not DCF // For 10-yr long sims DCFcalc could take >1hour
    #SFcalc=  False #Structure Function
    #STINGRAY=False # Time-dependent timelags quick method
    #SAVE=    True #saves diagrammes (for timelapse below 'saveframe') 
    #TIMELAPSE=False # Run a timelapse and print hundreds of frames to create a GIF/video later on or observe the evolution live.
    #TIMELAPSE_CORR=False
    STATIONARY_STATE ,days_with_lowstate ,name_stat_state , delta_stat , R_stat , route_S, nameflc_S= STATIONARY #saved background NT state to add on variable spectrum// second zone
    EXTRA_STATE , days_with_extrastate,  name_extra_state , delta_extra , R_extra, route_E, nameflc_E = EXTRA#saved background NT state to add on variable spectrum// second zone
    TIME_DEPENDENT_DATA ,route_timelapse = TIMEDEPENDENT_DATA_ON_PLOT
    
    #define your own colors depending on your segmentation
    colors_daysplit = {'56392':'black', '56393': 'darkblue', '56394': 'blue' , '56395': 'cyan'  , '56396': 'limegreen' ,
                       '56397': 'tab:olive'  , '56398': 'tab:orange' , '56399' : 'orangered',
                       '56400':'red','56401':'red'}

    #directories
    if route[-1]!='/':route+='/'
    route_main=  route.split('Results')[0]
    routesave = route  #where to print images
    routedata=route_main+'data/'#route.split('/BlaVar/')[0]+'/BlaVar/Var/' #directory of data points for SED
    routeSS=route_main+'Steady_States/'
    routelc=route_main+'Lightcurve_Simulation/' # routelc/Nameobj/ is the directory of real LC
    
    #you can choose your own comparison data light-curves here depending on what you used
    nameflcs=['NUSTAR_LC_7-30keV_3cols.txt' , 'MAGIC_LC_3cols.txt','VERITAS_3cols.txt'] # name(s) of the file of the real LC
    ir_colors = ['darkred','darkblue','tab:purple']

    dt_ELAPSE= 0.05
    
    """Selection of Object and Features"""

    
    #VARYING VARIABLE
    units={'le':'','B':'G','lext':'','gmax':'','theta':'^o','Gamma':'', 'p':'','delta':''}
    namevar=route.split('/')[-2].split('_')[0]
    if '+' in namevar: 
        namevar_sec =  namevar.split('+')[1]
        namevar= namevar.split('+')[0]
        DOUBLE_VAR, nunits_sec = True ,units[namevar_sec]
    else:
        DOUBLE_VAR, namevar_sec = False, 'None' 
    nunits=units[namevar]
    POWER = route.split('_')[-1].replace('/','') #route.split('^')[1].split('_')[-1].replace('/','')
    if namevar not in ['le','B','gmax','lext','p']:
        #namevar=input("Give Name of Var:\t")
        al = os.listdir('./')
        att = [os.path.getmtime('./'+a) for a in al]
        def att(a, al=al,att=att):  return att[al.index(a)]
        namevar = [i[9] for i in sorted(al, key=att) if 'code_var_' in i][-1]
    if namevar=='delta':
        print('Please use delta_var.py')
    

    
    
    """Variation"""
    #timesteps to plot as stored in fort.81 (warning! double times during crashing not removed yet)
    n1=1 #min 1
    nmax=-1  #maximum number of iterations , -1 for length of array
    initt=0 #MJD  : converts given value of t_0 in MJD
    #Timing analysis on segments for long light curves
    SEGnum=1 #Default is 1 for full LC analysis // Use greater number to bootstrap results,
                                                #splitting the time-interval into segments
    """Timelaspe"""
    pausetime=0.01 #seconds
    pastlength, lnwdth1, lnwdth2 =20 , 0.2 , 0.1 #previous states plotted in diagram half with lnwdth1 and half lnwdth2
    saveframes= True #only saves total SED frames into .jpgs
    monochromVIR, radio, gammarays, Xrays= True , True , True, True #only if TIMELAPSE on
    
    """Plot frames in SED with specific time interval"""
    SAMPLE_FRAMES=False #to sample and plot less states in the SED='on' option
    dtsmpl=100 # t_cross (check your time - interval for selection) 
    
    """Timecurves and Bands"""
    ADD_BB=   True  # add the black body spectrums of disk and BLR
    addnoise= True  # add Poisson noiser at optical bands
    tbin= -1 # bin output (days) -1 to use code output binning (<0.25 tcross / delta)
    tbin_S = 1 
    
    """Correlations"""
    SPOTS=False #spots appear / Connect dots if false
    mscorr=3.0 #markersize of scatter plot points
    CONTOUR=False #plot contour lines of 1 sigma included points
    corrfont=16 #size of x,ylabel+2 tile -3 ticks
    legsize='x-small' #legend size  medium, large x-large etc
    colcont='k'#'tab:orange' #color of contour in all plots
    lvls=[1/2] # list of values,  fraction of max density in CPs, to contour on
    
    """LightCurve Settings"""
    #observational
    realLC=     True  #Load the real LC UNEDITED for making histograms 
    #indobs  <-1.0 index of F(epsilon) for modelling observations, not used if band data is loaded
    indobs, indobs_S, indobs_E =    None, -1.5, -1.0 
    REALHIST=   False #add histogram near to the SED of all states
    
    #model light curves
    bands=['VHEband','xband','optical'] #more to less energetic
    #We cannot analysis 9 bands x 9 bands, you can choose here a combination that works for you
    optband=1 #which band of optical to keep from array of created optical channels based on code bins
             #each optical channel created has width 0.4 (log nu) and usually matches,
             #print 'opt + index', e.g. opt1 to find the central frequency of the optical band generated
    ixband=1 # 0:soft 3-7 keV  / 1:medium band 7-30 2:hard 30-80 keV / -1 Xtreme UV 10-124 eV "bands" list
    iVHEband=2 # 0soft, 1medium , 2hard bands  for vFv single data point at selected band for counts
    addtcross=0.0 #to compare escaped photons with input parameter variation (default 1)
    timespan=50 #days to zoom in in timecurve  (create zoome in window to study short-term variability)
    
    #average errorbars from LCs #average magnitude error from optical LCs #useful for fractional variability calc
    erro , errx, errg = 0.01, 0.01, 0.01
    
    """PLOTTING"""
    """SED settings, SED='on'"""
    limSEDs=500 #select nmax-n1 less than this limit to plot multiple SEDs
    MULTIcolor=True #Plot data points with different color for every interval of 2 years
    diconec={'le':'m','B':'c','lext':'y','gmax':'m', 'p':'g'}
    onecolor=diconec[namevar] # for mixed color and 1.0 transparancy
    #Indexing legend up-right
    lgd2opts={'le':'','B':'','lext':'','theta':'','p':'','delta':''}
    
    """Contour settings, if SED='con' """
    cmps={'le':plt.cm.RdPu,'B':plt.cm.GnBu,'delta':plt.cm.YlGn,\
          'lext':plt.cm.YlOrRd,'gmax':plt.cm.PuBu, 'p':plt.cm.PuBu,'delta':plt.cm.RdPu }
    nbx=500 #x-histogram
    nby=250 #y-histogram / check data-reduction parameter THRESHOLD for setting value
    Nlvs=10
    toplvl=0.3
    
    #x,y labels/units
    xlab='v'
    xoptions=['mec2','Hz','GeV']
    yoptions=['le','vLv','vFv','radio','Fv','Xrays','gamma'] #radio -> Jy / Xrays -> counts /sec / Gamma -> counts/hour
    yunits={'le':' ','vLv':'erg/s','vFv':'erg/s/cm^2','radio':'Jy',\
           'Fv':'erg/s/cm^2/Hz','Xrays':'counts/s/cm^2','gamma':'counts/hour/cm^2'}
    #Select here from above
    xaxis=xoptions[1]# CAUTION! indexing of the list begins from zero
    yaxis=yoptions[2]# CAUTION! indexing of the list begins from zero    
    #MULTIcolor Data Points 
    data_name=routedata+'/'+obj+".ascii"
    dtpts=[]
    colors =['tomato','gold','limegreen','purple']
    if MULTIcolor:
       dtpts.append(routedata+obj+'/'+obj+'_08-10.ascii')
       dtpts.append(routedata+obj+'/'+obj+'_10-12.ascii')
       dtpts.append(routedata+obj+'/'+obj+'_12-14.ascii')
       dtpts.append(routedata+obj+'/'+obj+'_14-16.ascii')
       legend_pts=['\'08-\'10','\'10-\'12','\'12-\'14','\'14-\'18']
       form=['o','.','s','+']
       fillstyle=['full','full','full','none']
       #dtpts.append(obj+'_16-18.ascii')
    
    """Color"""
    CI=     True      #calculates color if 'on', else omit all color diagrams and info
    cibnda , cibndb = 'B', 'J'  #bluer, redder color / choose from below
    optdic={'B':10**14.837, 'V':10**14.74, 'R':10**14.67, 'I':10**14.575,\
           'J':10**14.383, 'H':10**14.267, 'K':10**14.18} #optical filters
    FWHMdic={'B': 0.2, 'V': 0.16, 'R': 0.25, 'I': 0.19,\
             'J': 0.12, 'H': 0.18, 'K': 0.16} #FWHM foroptical filters
    ZPdic={'B':4000, 'V':3580, 'R':2971, 'I':2405,\
           'J':1565, 'H':1039, 'K':647.6} #zero point for optical filters (wavelength Angstrom)
    
    #all arrays must be delivered in frequencies
    range_optical=[10**14.17,10**14.837]
    range_EUV=[10*(eV/h) ,124*(eV/h)] #eV/h = 2.41e14
    range_VHE=[0.25*10**12*(eV/h),2.5*10**12*(eV/h)] #eV/h = 2.41e14
    range_VHE1=[0.2*10**12*(eV/h),0.4*10**12*(eV/h)] #eV/h = 2.41e14
    range_VHE2=[0.4*10**12*(eV/h),0.8*10**12*(eV/h)] #eV/h = 2.41e14 
    range_VHE3=[0.8*10**12*(eV/h),2.5*10**12*(eV/h)] #eV/h = 2.41e14 
    counts_labels= ['#/cm^2/s 0.2-0.4TeV', '#/cm^2/s 0.4-0.8TeV', '#/cm^2/s >0.8TeV']
    range_xsoft=[3e3*(eV/h),7e3*(eV/h)] #eV/h = 2.41e14
    range_xmed=[7e3*(eV/h),30e3*(eV/h)] #eV/h = 2.41e14
    range_xhard=[30e3*(eV/h),80e3*(eV/h)] #eV/h = 2.41e14
    MULTIoptical=   True # True-> generate multiple optical functions and study one of them with width 0.4,
                         # False-> integrate over the optical band 
    Nint=10000 #intergration bins (defaults 10^4)
    #for integration of logparabola
    x1, x2 = range_xmed[0] , 30e3*2.41e14
    beta_x , E_ev_x = 0.4 ,1e3
    x1, x2 = range_VHE[0] , 2.41e14*1e12
    beta_g , E_ev_g= 0.79 , 500e9
    
    
    
    """Timelags"""
    ####QUICK METHOD is commended out
    ####PEAKS #deactivated in this version // yields timelags using peaks or minima
    QUICK_TIMELAGS=False #SWITCH/ does not control Stingray
    num_peaks=25 #top peaks to compare for timelags
    peaktype='maxima' #'maxima' or 'minima'
    lendel=50 #delete -lendel, lendel pointsaround the peaks, to remove it
    
    #Discrete Correlation Function and Structure Function, same properties for structure function
    taubin=0.01
    lentau=50*taubin # days|| integer: -before/ + after zero timelag to calculate
    
    #Timelag per Frequency
    segment_size=20 #days #depends on Run Time Interval
    
    """Data Reduction"""
    nu_cutoff=13 #energy in mec^2 to cutoff plotting 
    HEIGHT=1 #height of y-axis (threshold) in logscale measured below the minimum value of flux
    THRESHOLD=25 #orders of magnituded for data reduction from maximum value
    """%%%%%%%%%%%%%%%%%%    END OF SETTINGS   %%%%%%%%%%%%%%"""
    
    #
    #
    #
    # Activate all tools besides SED
    if ALLTOOLS:   CORR, TIMECURVES, CI , DCFcalc, STINGRAY =[True]*5
    
    """The fort files to use"""
    nmlist=os.listdir(route)
    names=[]
    cnt=0
    for nm in nmlist:
        if 'fort_' in nm and '.85' in nm and 'steady' not in nm:
            cnt=cnt+1 #.85, .89 files        
            names.append(nm.replace('.85',''))
    names.sort()
    if len(names)==0:
        names.append('fort')
    #names=['fort']
    
    
    """INPUT Observable Parameters"""
        
    # The object + the name of the variable parameter is detected from the named directory
    load_input_information_on_objects()
    #obj=os.getcwd().split('/')[-3]
    if obj not in objects:
        print(objects)
        oo=int(input("\n Give Number of object [0 , 1 ..] from following list:\t"))
        #oo = len(objects)-1
        #obj=objects[oo]
    oo=objects.index(obj)
    
    D=Dists[oo]*10**6*pc #"Distance of Source in Mpc"
    zsh=redsh[oo] #redshift
    dilate=1-((((1+zsh)**2-1)/((1+zsh)**2+1))**2)**2 #dilation factor, https://www.scirp.org/journal/paperinformation.aspx?paperid=65598
    SSConly=SSC[oo] #if SSC or EC modelling is used
    SECONDBB=BBextra[oo] #add an additional BB on spectrum (usually DISK) but,
                             # that is not involved in external photon fields of code
    Gamma2 , factbb ,Tbb2, Rbb2 = Gamma2_arr[oo], factbb_arr[oo], Tbb2_arr[oo], Rbb2_arr[oo]
    
    cn4=open(routeSS+'/code_new4.inp','r')
    lines=cn4.readlines()
    cn4.close()
    
    tend=float(lines[0].split()[-1].replace('\n',''))
    nsteps= float(lines[0].split()[-2])
    R=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
    B=float(lines[2].split()[1]) #dimension of source in cm
    p=float(lines[4].split()[3])
    loggmin=float(lines[4].split()[1])
    loggmax=float(lines[4].split()[2])
    logle=float(lines[4].split()[4])
    delta=float(lines[7].split()[0]) #set zero if used the above for computation
    Gamma=float(lines[8].split()[0])
    T=float(lines[5].split()[1])
    lext=float(lines[5].split()[2])
    tcross=R/c/3600/24 #days jet frame
    tcr=tcross*dilate/delta
    
    #second/slow component
    if STATIONARY_STATE:
        cn4=open(routeSS+'code_new4_slow.inp','r')
        lines=cn4.readlines()
        cn4.close()
        tend_S=float(lines[0].split()[-1].replace('\n',''))
        nsteps_S= float(lines[0].split()[-2])
        R_S=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
        B_S=float(lines[2].split()[1]) #dimension of source in cm
        p_S=float(lines[4].split()[3])
        loggmin_S=float(lines[4].split()[1])
        loggmax_S=float(lines[4].split()[2])
        logle_S=float(lines[4].split()[4])
        delta_S=float(lines[7].split()[0]) #set zero if used the above for computation
        Gamma_S=float(lines[8].split()[0])
        T_S=float(lines[5].split()[1])
        lext_S=float(lines[5].split()[2])
        tcross_S=R_S/c/3600/24 #days jet frame
        tcr_S=tcross_S*dilate/delta_S
        
    #second/slow component
    if EXTRA_STATE:
        cn4=open(routeSS+'/code_new4_fast.inp','r')
        lines=cn4.readlines()
        cn4.close()
        tend_E=float(lines[0].split()[-1].replace('\n',''))
        nsteps_E= float(lines[0].split()[-2])
        R_E=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
        B_E=float(lines[2].split()[1]) #dimension of source in cm
        p_E=float(lines[4].split()[3])
        loggmin_E=float(lines[4].split()[1])
        loggmax_E=float(lines[4].split()[2])
        logle_E=float(lines[4].split()[4])
        delta_E=float(lines[7].split()[0]) #set zero if used the above for computation
        Gamma_E=float(lines[8].split()[0])
        T_E=float(lines[5].split()[1])
        lext_E=float(lines[5].split()[2])
        tcross_E=R_E/c/3600/24 #days jet frame
        tcr_E=tcross_E*dilate/delta_E
        
    if feedback_VHE_on_extrastate:
        extra_tc  = np.loadtxt(routelc+'residuals_MAGIC.txt')
        extra_tc[:,1][extra_tc[:,1]<0.1] =  0.0
        extra_tc[:,1] = (1.0+np.sign(extra_tc[:,1])*extra_tc[:,1]**0.25)*delta_S*2.0
       
    
    def read_state(filename, NUM_OUT = 3, THRESHOLD=15, nu_cutoff=13, xaxis=xaxis, yaxis=yaxis, yunits= yunits, D=D, zsh=zsh, R=R,
                   delta = delta,  SSConly = SSConly, T = T, lext = lext, Rext = Rbb, Gamma = Gamma,
                   Tbb2 = Tbb2 ,Rbb2 = Rbb2, factbb = factbb , Gamma2=Gamma2):    
       """Reading Files"""
       n5=filename+'.85'
       file85=open(n5,'r')
       f85=file85.readlines()
       file85.close()
       eps85=[]
       I85=[]
    
       LEN=int(len(f85)/NUM_OUT-1)
       for y in f85:
           x=y.split()
           z=x[0].strip()
           w=x[1].replace('\n','').strip()
           if w=='' or w=='NAN':
               z='-100'
               w='-100'
           eps85.append(z)
           I85.append(w)
       eps85=np.array(eps85[-LEN::]).astype(float)
       I85=np.array(I85[-LEN::]).astype(float)
       
       """REMOVE UNWANTED POINTS"""    
       """Remove Extremely High Energy Gamma Rays/ Right Cut-off of Frequencies""" 
       zz=np.where(eps85>nu_cutoff,eps85,0)
       zzz=zz.nonzero()
       for i in zzz:
           eps85[i]=-100
           I85[i]=-100
       
       """Remove Extremely High Outlier Points"""
       #THRESHOLD=15 #orders of magnitude of allowed difference of consecutive points 
       eps85=eps85[I85>max(I85)-THRESHOLD]
       I85=I85[I85>max(I85)-THRESHOLD]
       
       """BLACK BODY ADDITION"""
       eps85bb=eps85
       Ibb=eps85
       if not SSConly:
           if len(eps85bb)<200:
               eps85bb=np.linspace(min(eps85),max(eps85),200)
           Ibb=BB(eps85bb,T/Gamma,lext/Gamma**2,Rext)+np.log10(Gamma) #due to \nu boosting that wasn't used while computing nu *Fnu
           if SECONDBB:
               lext2=1/Gamma**2*factbb*(Rext/Rbb2)**2*lext
               Ibb2=BB(eps85bb,Tbb2,lext2,Rbb2)+np.log10(Gamma)
               xbb2=eps85bb+np.log10((me*c**2/h))
               ybb2=np.log10(10**Ibb2/D**2*R*me*c**3/sigmaT)
       
       """TRANSFORMATION to the frame of reference of the Observer"""
       if xaxis=='mec2':
           x85=eps85+np.log10(delta/(1+zsh))
           xbb=eps85bb
       if xaxis=='Hz':
           x85=eps85+np.log10(delta*(me*c**2/h)/(1+zsh))
           xbb=eps85bb+np.log10((me*c**2/h))
       if xaxis=='GeV':
           x85=np.log10(delta*10**eps85/eV/10**9*me*c**2/(1+zsh))
           xbb=np.log10(10**eps85bb/eV/10**9*me*c**2)
       if yaxis=='le':
           y85=I85+np.log10(delta**4)
           ybb=Ibb
       if yaxis=='vFv':
           y85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3)
           ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT)
           
       if yaxis=='Fv':
           y85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3)-eps85-np.log10(delta*(me*c**2/h))
           ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT)-eps85bb-np.log10((me*c**2/h))
       if yaxis=='vLv':
           y85=np.log10(10**I85*delta**4*(4*np.pi)*R*me*c**3/sigmaT/3)
           ybb=np.log10(10**Ibb**R*me*c**3/sigmaT)
       if yaxis=='radio':
           y85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3/Jy)-eps85-np.log10(delta*(me*c**2/h))
           ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT/Jy)-eps85bb-np.log10((me*c**2/h)) 
       if yaxis=='Xrays':
           y85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3)-eps85-np.log10(delta)
           ybb=np.log10(10**Ibb/D**2*R*c/sigmaT)-eps85bb
       if yaxis=='gamma':
           y85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3*3600)-eps85-np.log10(delta)
           ybb=np.log10(10**Ibb/D**2*R*c/sigmaT*3600)-eps85bb
           
       """RECREATION Of the Total Thermal + Non-Thermal SED""" 
       if not SSConly:
           x85nn=x85[y85>(max(y85)-THRESHOLD)] #without noise and unwanted points for interpolation to work
           y85nn=y85[y85>(max(y85)-THRESHOLD)]
           funcy85=inpl.interp1d(x85nn,y85nn)
           #x81tbb=np.concatenate([x85,xbb]) #overpose NOT addition
           #y81tbb=np.concatenate([y85,ybb]) #overpose NOT addition
           x81tbb=xbb[xbb<max(x85nn)][xbb[xbb<max(x85nn)]>min(x85nn)]
           y81tbb=np.log10(10**ybb[xbb<max(x85nn)][xbb[xbb<max(x85nn)]>min(x85nn)]+10**funcy85(x81tbb))
           if SECONDBB:
                   funcybb2=inpl.interp1d(xbb2,ybb2)
                   y81tbb=np.log10(10**funcybb2(x81tbb)+10**y81tbb)
           tbb=open('fort.81tbb','w')
           for i,j in zip(x81tbb,y81tbb):
                   tbb.write(str(i)+'\t'+str(j)+'\n')
           tbb.close()
       else:
            x81tbb=x85
            y81tbb=y85
       return x85, y85, x81tbb, y81tbb
    """
    %
    %
    %
    """
    #%%
    """Observational Data"""
    ##Intergrating intervals for optical, X-rays and gamma-ray bands
    lnopt=np.linspace(range_optical[0],range_optical[1],Nint)
    lnx=np.linspace(range_xsoft[0],range_xsoft[1],Nint)
    lnxm=np.linspace(range_xmed[0],range_xmed[1],Nint)
    lnxx=np.linspace(range_xhard[0],range_xhard[1],Nint)
    #counts when integrated
    

    #time-dependent data input
    if TIME_DEPENDENT_DATA:
        llf = [[float(i.split('_')[1]),float(i.split('_')[2].replace('.txt',''))] for i in os.listdir(route_timelapse) if 'sed_' in i]
        llf.sort()
        nlf  = np.array(llf) #-llf[0][0]
        MULTIcolor = False
        mjdday = int(np.sum(llf[0])/2)
    
    LOGdata = True
    points=ascii.read(data_name)
    Cunitx=0 #4 Hz 
    Cunity=0 #for erg/cm^2/s
    v_pts=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
    vFv_pts=np.array([points["col2"]])[0]-np.ones(1)*Cunity
    #errorv=abs(points["col4"])
    errorvFv=(np.array(abs(points["col3"]))+ np.array(abs(points["col4"])))/2
    if LOGdata:
        v_pts = np.log10(v_pts)
        vFv_pts = np.log10(vFv_pts)
        errorvFv = errorvFv/vFv_pts/np.log(10)
    #average measurement errors
    #erro=np.mean(errorvFv[v_pts<np.log10(10**15.2)][v_pts[v_pts<np.log10(10**15.2)]>np.log10(10**14.0)])
    #errx=np.mean(errorvFv[v_pts<np.log10(range_xhard[1])][v_pts[v_pts<np.log10(range_xhard[1])]>np.log10(range_xsoft[0])])
    #errg=np.mean(errorvFv[v_pts<np.log10(range_VHE[1])][v_pts[v_pts<np.log10(range_VHE[1])]>np.log10(range_VHE[0])])
    
    if MULTIcolor:
      for i in range(len(dtpts)):
        points=ascii.read(dtpts[i])
        locals()['v_pts'+str(i)]=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
        locals()['vFv_pts'+str(i)]=np.array([points["col3"]])[0]-np.ones(1)*Cunity
        locals()['errorv'+str(i)]=0.009
        locals()['errorvFv'+str(i)]=np.array([points["col4"]])[0]
    
    func_obs=inpl.interp1d(v_pts,vFv_pts)
    
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    """
    %
    %
    %
    %
    """
    """READING FILES"""
    #read fort files .85 (.81 ) , steady state .85 and .89
    epsfile, Ifile, gammafile, Nfile, filetime, filevar = [],[],[],[],[],[] #blob
    epsfile_S, Ifile_S, gammafile_S, Nfile_S, filetime_S, filevar_S = [],[],[],[],[],[] #slow-varying stationary blob (variations are external, or stationary per given segment)
    epsfile_E, Ifile_E, gammafile_E, Nfile_E, filetime_E, filevar_E = [],[],[],[],[],[] #fast-moving blob (variations induced on-fly of the others)
    R_all , B_all, p_all, logle_all, loggmin_all, loggmax_all, delta_all  = [],[],[],[],[],[],[]
    R_all_S , B_all_S, p_all_S, logle_all_S, loggmin_all_S, loggmax_all_S, delta_all_S  = [],[],[],[],[],[],[]
    R_all_E, B_all_E, p_all_E, logle_all_E, loggmin_all_E, loggmax_all_E, delta_all_E  = [],[],[],[],[],[],[]
    lens, lens_S , lens_E, lens89,  lens89_S, lens89_E, filevar_sec = [],[],[],[],[],[],[]
    n0 , n0_S  , n0_E= len(filetime) , len(filetime_S) , len(filetime_E)
    tfirst ,  tlast =  [] , []   
    tfirst_dic = np.load('tfirst.npy', allow_pickle=True).all()
    tlast_dic = np.load('tlast.npy', allow_pickle=True).all()
    diffcolgamma = np.loadtxt(routelc+'diff_colgamma.txt')
    diffcolgamma[:,1] = smooth(diffcolgamma[:,1],20)
    m = [i  not in [56396,56397] for  i in list(diffcolgamma[:,0].astype(int))]
    diffcolgamma[:,1][m] = 0.0
    # tfirst_dic_new = {}
    # tlast_dic_new = {}
    # for t in tfirst_dic.keys():
    #     v = tfirst_dic[t]
    #     tn = t
    #     if '93' not in t:
    #         tn = str(int(t)+1)
    #     tfirst_dic_new[tn] = v
    #     tlast_dic_new[tn] = v
    # tfirst_dic , tlast_dic = tfirst_dic_new , tlast_dic_new
    #read multi.log to correct for times of different runs
    #IF INCOMPATIBILITY WITH INPUT RUNTIME ADN Tmax
    #if len(names)>1:
    #    mlog=open(route+'multi.log','r')
    #    tms=mlog.readlines()
    #    mlog.close()
    #    tms=np.array(tms).astype(float)
    for name in names:
        #Read fort.81 or .85
        fphotons=open(route+name+'.85','r')
        f85=fphotons.readlines()
        fphotons.close()    
        LEN,   LEN_S , LEN_E = 0 , 0  , 0
        for ln in f85:
            if float(ln.split("  ")[1])==0. and LEN!=0:
               break
            LEN+=1
                
        n0+=int(len(f85)/LEN)
        
        if len(f85)/LEN-int(len(f85)/LEN)!=0: print('File '+name+' has some rows missing. Please check output')
        
        cn4=open(routeSS+name.split('_')[1]+'_new4.inp','r')
        lines=cn4.readlines()
        cn4.close()
    
        tend=float(lines[0].split()[-1].replace('\n',''))
        nsteps= float(lines[0].split()[-2])
        R=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
        B=float(lines[2].split()[1]) #dimension of source in cm
        p=float(lines[4].split()[3])
        loggmin=float(lines[4].split()[1])
        loggmax=float(lines[4].split()[2])
        logle=float(lines[4].split()[4])
        delta=float(lines[7].split()[0]) #set zero if used the above for computation
        Gamma=float(lines[8].split()[0])
        T=float(lines[5].split()[1])
        lext=float(lines[5].split()[2])
        tcross=R/c/3600/24 #days jet frame
        tcr=tcross*dilate/delta
    
        for y in f85:
                x=y.split()
                z=x[0].strip()
                w=x[1].replace('\n','').strip()
                if w in ['','-INF']:
                    w='-100'        
                if w=='NAN':
                    w='-100'
                epsfile.append(z)
                Ifile.append(w)   
                if w=='-':
                    print(name)
        
        if ELECTRONS:
            for ln in f89:
               if float(ln.split("  ")[1])==0. and LEN89!=0:
                  break
               LEN89+=1

                   
            felectrons=open(route+name+'.89','r')
            f89=felectrons.readlines()
            felectrons.close()
            LEN89,  LEN89_S, LEN89_E = 0 , 0 , 0 
                   
            for y in f89:
                x=y.split()
                z=x[0].strip()
                w=x[1].replace('\n','').strip()
                if w in ['','-INF']:
                    w='-100'        
                if w=='NAN':
                    w='-100'
                gammafile.append(z)
                Nfile.append(w)
                if w=='-': print(name)
                    
            
        timebin=tend/nsteps    
        #read fort.55
        name55=name+'.55' #name55='fakeTC.txt' #read analytically all the input fake variables
        varfile=open(route+name55,'r')
        for y in varfile:
            x=y.split()
            z=x[0]
            w=x[1].replace('\n','')
            if DOUBLE_VAR:
                ww=x[2].replace('\n','')
                filevar_sec.append(ww)
            filetime.append(z)
            filevar.append(w)
            R_all.append(R) 
            B_all.append(B) 
            p_all.append(p) 
            logle_all.append(logle) 
            loggmin_all.append(loggmin) 
            loggmax_all.append(loggmax) 
            delta_all.append(delta) 
            lens.append(LEN)
            if ELECTRONS:lens89.append(LEN89)
            tfirst.append(tfirst_dic[name.split('_')[1]])
            tlast.append(tlast_dic[name.split('_')[1]])
        varfile.close()
                    
        #slow state
        if STATIONARY_STATE:
            #Read fort.81 or .85
            fphotons=open(route_S+'fort_0.85','r')
            f85_S=fphotons.readlines()
            fphotons.close()
           
            for ln in f85_S:
                if float(ln.split("  ")[1])==0. and LEN_S!=0:
                   break
                LEN_S+=1
            
            n0_S+=int(len(f85_S)/LEN_S)
            #second/slow component
            cn4=open(routeSS+'lowstate_'+name.split('_')[1]+'_new4.inp','r')
            lines=cn4.readlines()
            cn4.close()
            tend_S=float(lines[0].split()[-1].replace('\n',''))
            nsteps_S= float(lines[0].split()[-2])
            R_S=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
            B_S=float(lines[2].split()[1]) #dimension of source in cm
            p_S=float(lines[4].split()[3])
            loggmin_S=float(lines[4].split()[1])
            loggmax_S=float(lines[4].split()[2])
            logle_S=float(lines[4].split()[4])
            delta_S=float(lines[7].split()[0]) #set zero if used the above for computation
            Gamma_S=float(lines[8].split()[0])
            T_S=float(lines[5].split()[1])
            lext_S=float(lines[5].split()[2])
            tcross_S=R_S/c/3600/24 #days jet frame
            tcr_S=tcross_S*dilate/delta_S
            
            for y in f85_S:
                    x=y.split()
                    z=x[0].strip()
                    w=x[1].replace('\n','').strip()
                    if w in ['','-INF']:
                        w='-100'        
                    if w=='NAN':
                        w='-100'
                    epsfile_S.append(z)
                    Ifile_S.append(w)  
            if ELECTRONS:
                felectrons=open(route_S+'fort_0.89','r')
                f89_S=felectrons.readlines()
                felectrons.close()
                
                for ln in f89_S:
                   if float(ln.split("  ")[1])==0. and LEN89_S!=0:
                      break
                   LEN89_S+=1
                   
                for y in f89_S:
                        x=y.split()
                        z=x[0].strip()
                        w=x[1].replace('\n','').strip()
                        if w in ['','-INF']:
                            w='-100'        
                        if w=='NAN':
                            w='-100'
                        gammafile_S.append(z)
                        Nfile_S.append(w)
                    
            timebin_S=tend_S/nsteps_S    
            #read fort.55
            name55_S='fort_0.55' #name55='fakeTC.txt' #read analytically all the input fake variables
            #varfile=open(route_S+name55_S,'r') use_the_main run varfile
            varfile=open(route+name55,'r') 
            for y in varfile:
                x=y.split()
                z=x[0]
                w=x[1].replace('\n','')
                filetime_S.append(z)
                filevar_S.append(w)
                R_all_S.append(R_S) 
                B_all_S.append(B_S) 
                p_all_S.append(p_S) 
                logle_all_S.append(logle_S) 
                loggmin_all_S.append(loggmin_S) 
                loggmax_all_S.append(loggmax_S) 
                delta_all_S.append(delta_S) 
                lens_S.append(LEN_S)
                if ELECTRONS:lens89_S.append(LEN89_S)
            varfile.close()
            
            
        #Read fort.81 or .85
        if EXTRA_STATE:
            fphotons=open(route_E+'fort_0.85','r')
            f85_E=fphotons.readlines()
            fphotons.close()            
            for ln in f85_E:
                if float(ln.split("  ")[1])==0. and LEN_E!=0:
                   break
                LEN_E+=1
            if int(name.split('_')[1].replace('half','')) in days_with_extrastate:
                n0_E+=int(len(f85_E)/LEN_E)
                #second/slow component
                cn4=open(routeSS+'extrastate_'+name.split('_')[1]+'_new4.inp','r')
                lines=cn4.readlines()
                cn4.close()
                tend_E=float(lines[0].split()[-1].replace('\n',''))
                nsteps_E= float(lines[0].split()[-2])
                R_E=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
                B_E=float(lines[2].split()[1]) #dimension of source in cm
                p_E=float(lines[4].split()[3])
                loggmin_E=float(lines[4].split()[1])
                loggmax_E=float(lines[4].split()[2])
                logle_E=float(lines[4].split()[4])
                delta_E=float(lines[7].split()[0]) #set zero if used the above for computation
                Gamma_E=float(lines[8].split()[0])
                T_E=float(lines[5].split()[1])
                lext_E=float(lines[5].split()[2])
            else:
                R_E , B_E, p_E, logggmin_E , loggmax_E, logle_E, delta_E ,Gamma_E, T_E, lext_E = [np.nan,]*10
            tcross_E=R_E/c/3600/24 #days jet frame
            tcr_E=tcross_E*dilate/delta_E
            for y in f85_E:
                    x=y.split()
                    z=x[0].strip()
                    w=x[1].replace('\n','').strip()
                    if w in ['','-INF']:
                        w='-100'        
                    if w=='NAN':
                        w='-100'
                    epsfile_E.append(z)
                    Ifile_E.append(w)  
            if ELECTRONS:
                felectrons=open(route_E+'fort_0.89','r')
                f89_E=felectrons.readlines()
                felectrons.close()
                
                for ln in f89_E:
                   if float(ln.split("  ")[1])==0. and LEN89_E!=0:
                      break
                   LEN89_E+=1
                   
                for y in f89_E:
                        x=y.split()
                        z=x[0].strip()
                        w=x[1].replace('\n','').strip()
                        if w in ['','-INF']:
                            w='-100'        
                        if w=='NAN':
                            w='-100'
                        gammafile_E.append(z)
                        Nfile_E.append(w)
            timebin_E=tend_E/nsteps_E    
            #read fort.55
            name55_E='fort_0.55' #name55='fakeTC.txt' #read analytically all the input fake variables
            #varfile=open(route_E+name55_E,'r') #use the main-run varfile
            varfile=open(route+name55,'r') 
            for y in varfile:
                x=y.split()
                z=x[0]
                w=x[1].replace('\n','')
                filetime_E.append(z)
                filevar_E.append(w)
                R_all_E.append(R_E) 
                B_all_E.append(B_E) 
                p_all_E.append(p_E) 
                logle_all_E.append(logle_E) 
                loggmin_all_E.append(loggmin_E) 
                loggmax_all_E.append(loggmax_E) 
                delta_all_E.append(delta_E) 
                lens_E.append(LEN_E)
                if ELECTRONS:lens89_E.append(LEN89_E)
            varfile.close()
                    
    epsfile=np.array(epsfile).astype(float)
    Ifile=np.array(Ifile).astype(float)
    filetime=np.array(filetime).astype(float) 
    filevar=np.array(filevar).astype(float)
    filevar_sec=np.array(filevar_sec).astype(float)
    deps=round(epsfile[-1]-epsfile[-2],5) 
    lens  = np.array(lens)
    R_all , B_all, p_all, logle_all, loggmin_all, loggmax_all, delta_all  = np.array(R_all) , np.array(B_all) , np.array(p_all) , np.array(logle_all) , np.array(loggmin_all) , np.array( loggmax_all) , np.array(delta_all)
    tcross_all=R_all/c/3600/24 #days jet frame
    tcr_all=tcross_all*dilate/delta_all   
    tfirst, tlast =  np.array(tfirst) , np.array(tlast)
    if ELECTRONS:
        gammafile=np.array(gammafile).astype(float)
        Nfile=np.array(Nfile).astype(float)
        lens89 = np.array(lens89)
    if STATIONARY_STATE:
        epsfile_S =np.array(epsfile_S).astype(float)
        Ifile_S =np.array(Ifile_S).astype(float)
        filetime_S =np.array(filetime_S).astype(float) 
        filevar_S =np.array(filevar_S).astype(float)
        deps_S =round(epsfile_S[-1]-epsfile_S[-2],5) 
        if ELECTRONS:
            gammafile_S =np.array(gammafile_S).astype(float)
            Nfile_S =np.array(Nfile_S).astype(float)
            lens89_S = np.array(lens89_S)
        lens_S  = np.array(lens_S)
        R_all_S , B_all_S, p_all_S, logle_all_S, loggmin_all_S, loggmax_all_S, delta_all_S  = np.array(R_all_S) , np.array(B_all_S) , np.array(p_all_S) , np.array(logle_all_S) , np.array(loggmin_all_S) , np.array( loggmax_all_S) , np.array(delta_all_S)
        tcross_all_S=R_all_S/c/3600/24 #days jet frame
        tcr_all_S=tcross_all_S*dilate/delta_all_S 
    if EXTRA_STATE:
        epsfile_E =np.array(epsfile_E).astype(float)
        Ifile_E =np.array(Ifile_E).astype(float)
        filetime_E =np.array(filetime_E).astype(float) 
        filevar_E =np.array(filevar_E).astype(float)
        deps_E =round(epsfile_E[-1]-epsfile_E[-2],5)
        if ELECTRONS:
            gammafile_E =np.array(gammafile_E).astype(float)
            Nfile_E =np.array(Nfile_E).astype(float) 
            lens89_E= np.array(lens89_E)
        lens_E  = np.array(lens_E)
        R_all_E , B_all_E, p_all_E, logle_all_E, loggmin_all_E, loggmax_all_E, delta_all_E  = np.array(R_all_E) , np.array(B_all_E) , np.array(p_all_E) , np.array(logle_all_E) , np.array(loggmin_all_E) , np.array( loggmax_all_E) , np.array(delta_all_E)
    n0 , n0_S  , n0_E= len(filetime) , len(filetime_S) , len(filetime_E)
    """ Real LC used for that object/ check info file in Lightcurve Simulation"""
    if realLC:
            ir_names=  []
            #for VHE
            band_rf_x=np.logspace(np.log10(range_xsoft[0]),np.log10(range_xhard[1]))
            band_rf_vhe= np.logspace(np.log10(range_VHE[0]),np.log10(range_VHE[1]),5)
            band_rf_fermi= np.logspace(0.1e9*eV/h,np.log10(300e9*eV/h),5)
            obsbands = [band_rf_x, band_rf_vhe, band_rf_vhe]
            for ir,nameflc in enumerate(nameflcs):
                obst,obsf = [], []
                ir_names.append(nameflc.replace('_3cols.txt',''))
                realfile=routelc+nameflc    
                rlines=open(realfile)
                for r in rlines:
                    obst.append(float(r.split()[0]))
                    obsf.append(float(r.split()[1]))
                rlines.close()
                if ir==0:
                    t0  = obst[0]
                rt=np.array(obst)
                rt = rt #-t0
                rf=np.log10(np.array(obsf))
                locals()['rt_'+str(ir)],   locals()['rf_'+str(ir)] = rt, rf
                
                """ The analysis below is if importing an integrated flux with unknown spectrum (indobs !=None)"""
                if indobs!=None:
                    locals()['obs_flux'+str(ir)]=[]
                    bands_rf = obsbands[ir]
                    #band_rf=np.logspace(np.log10(range_VHE[0]),np.log10(range_VHE[1]),5) #5 X monochromatic
                    if indobs!=-1.:
                        Af=(1+indobs)*10**rf/(range_VHE[1]**(1+indobs)-range_VHE[0]**(1+indobs)) #monochromatic
                    if indobs==-1.:
                        Af=10**rf/np.log(range_VHE[1]/range_VHE[0]) #monochromatic 
                    Af=10**rf  #integral
                        
                    normrf=np.mean(Af)#3.75e-10
                    if oo==2:
                        normrf=3.017823941751141e-10
                    renorm=normrf/np.mean(Af)
                    Af=renorm*Af #renormalize
                    for feps in band_rf:
                            locals()['obs_flux'+str(ir)].append(np.log10(Af*feps**(1+indobs)))
                    bands_rf = np.log10(bands_rf)
    
            if STATIONARY_STATE and nameflc_S!='None':   
                realfile=routelc+nameflc_S
                rlines=open(realfile)
                rt_S, rf_S =[], []
                title = rlines.readline()
                for r in rlines: 
                    rt_S.append(float(r.split()[0]))
                    rf_S.append(float(r.split()[1]))
                rlines.close()
                t0_S = rt_S[0]
                rt_S=np.array(rt_S)
                rt_S = rt_S - t0_S
                rf_S=np.log10(np.array(rf_S))
                if indobs_S!=None:
                    indobs=indobs_S
                    band_rf= band_rf_fermi #5 X monochromatic
                    if indobs!=-1.:
                        Af=(1+indobs)*10**rf/(range_VHE[1]**(1+indobs)-range_VHE[0]**(1+indobs)) #monochromatic
                    if indobs==-1.:
                        Af=10**rf/np.log(range_VHE[1]/range_VHE[0]) #monochromatic 
                    Af_S=10**rf_S  #integral
                    normrf_S=np.mean(Af_S)#3.75e-10
                    renorm=normrf_S/np.mean(Af_S)
                    Af_S=renorm*Af_S #renormalize
                    obs_flux_S=[]
                    for feps in band_rf:
                            obs_flux_S.append(np.log10(Af*feps**(1+indobs)))
                    band_rf_S = np.log10(band_rf)
            if EXTRA_STATE and nameflc_E!='None':   
                realfile=routelc+nameflc_E
                rlines=open(realfile)
                rt_E, rf_E =[], []
                title = rlines.readline()
                for r in rlines: 
                    rt_E.append(float(r.split()[0]))
                    rf_E.append(float(r.split()[1]))
                rlines.close()
                t0_E = rt_E[0]
                rt_E=np.array(rt_E)
                if 'fakeTC' in name55:
                    rt_E=rt_E+filetime_E[list(filevar_E).index(max(filevar_E))]*tcr_E\
                       -rt_E[rf.index(max(rf))]+addtcross*tcr_E
                else:
                    rt_E=rt_E+filetime_E[list(filevar_E).index(max(filevar_E))]\
                      -rt[rf_E.index(max(rf_E))]+addtcross*tcr_E
                rf_E=np.log10(np.array(rf_E))
        
                if indobs_E!=None:
                    indobs=indobs_E
                    obs_flux_E=[]
                    band_rf=np.logspace(np.log10(range_VHE3[0]),np.log10(range_VHE3[1]),3) #5 X monochromatic
                    if indobs!=-1.:
                        Af=(1+indobs)*10**rf/(range_VHE[1]**(1+indobs)-range_VHE[0]**(1+indobs)) #monochromatic
                    if indobs==-1.:
                        Af=10**rf/np.log(range_VHE[1]/range_VHE[0]) #monochromatic 
                    Af_E=10**rf_E  #integral
                    normrf_E=np.mean(Af_E)#3.75e-10
                    renorm=normrf_E/np.mean(Af_E)
                    Af_E=renorm*Af_E #renormalize
                    for feps in band_rf:
                            obs_flux_E.append(np.log10(Af*feps**(1+indobs)))
                    band_rf_E = np.log10(band_rf)
    
    
    addt , addt_S, addt_E = 0.0 , 0.0 ,  0.0
    for it in range(len(filetime)):
        if filetime[it]<1e-5 and it>1:
            addt=filetime[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
            # the time of steady state (conjoint part of two TC)
        filetime[it]=filetime[it]+addt
    print('Runtime calculated from adding up all fakeTC_*.txt files:\t',\
          str(max(filetime*tcr)))
    if STATIONARY_STATE:
        for it in range(len(filetime_S)):
            if filetime_S[it]<1e-5 and it>1:
                addt_S=filetime_S[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
                # the time of steady state (conjoint part of two TC)
            filetime_S[it] += addt_S
        print('Also Stationary file option is on: :\t',\
                  str(max(filetime_S*tcr_S)))
        X_S=[]
        for i in range(len(filetime_S)-1):
            X_S.append(np.sum((filetime_S[1:i+1]-filetime_S[0:i])*tcr_all_S[0:i])-t0+t0_S)
            
    if EXTRA_STATE:
        for it in range(len(filetime_E)):
            if filetime_E[it]<1e-5 and it>1:
                addt_E=filetime_E[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
                # the time of steady state (conjoint part of two TC)
            filetime_E[it] += addt_E
    print('Check if this(these) number(s) agrees with the initialized runtime.\n Edit manually if gaps are introduced during the run')
    #%%
    if namevar=='B':
            filevar=10**filevar
            if STATIONARY_STATE: filevar_S=10**filevar_S      
    if namevar_sec=='le':
        filevar_sec = np.log10(filevar_sec)
        
    nmax_S , nmax_E = nmax , nmax
    if nmax==-1:
            nmax=n0-1#len(names) #-3 if added SS at the end/ times 0. 0. 5.
            if STATIONARY_STATE:nmax_S=n0_S-1#len(names) #-3 if added SS at the end/ times 0. 0. 5.
            if EXTRA_STATE:     nmax_E=n0_E-1#len(names) #-3 if added SS at the end/ times 0. 0. 5.
    
    """ STEADY STATE"""
    #STEADY STATE
    epsst,Ist= [], []
    filest=open(routeSS+name.replace('fort_','')+'.85','r')
    fsteady=filest.readlines()
    filest.close()
    LENst = int(len(fsteady)/3)
    for y in fsteady:
                x=y.split()
                z=x[0].strip()
                w=x[1].replace('\n','').strip()
                if w=='':
                    w='-100'
                if w=='NAN':
                    w='-100'
                epsst.append(z)
                Ist.append(w)
    epsst=np.array(epsst).astype(float)[-LENst+1::]
    Ist=np.array(Ist).astype(float)[-LENst+1::]
    epsst=epsst[Ist>max(Ist)-THRESHOLD]
    Ist=Ist[Ist>max(Ist)-THRESHOLD]
    xsteady=np.log10(delta*10**epsst*(me*c**2/h)/(1+zsh))  #observing (SED) units
    ysteady=np.log10(10**Ist*delta**4/D**2*R*me*c**3/sigmaT/3)  #observing (SED) units
    
    #useful for time-dependent secondary state only
    # if STATIONARY_STATE:
    #    filest=open(routeSS+name+'.85','r')
    #    fsteady=filest.readlines()
    #    filest.close()
    #   LENst_S = int(len(fsteady)/3)
    #     epsst,Ist= [], []
    #     for y in fsteady_S:
    #                 x=y.split()
    #                 z=x[0].strip()
    #                 w=x[1].replace('\n','').strip()
    #                 if w=='':
    #                     w='-100'
    #                 if w=='NAN':
    #                     w='-100'
    #                 epsst.append(z)
    #                 Ist.append(w)
    #     epsst_S=np.array(epsst).astype(float)[-LENst_S+1::]
    #     Ist_S=np.array(Ist).astype(float)[-LENst_S+1::]
    #     epsst_S=epsst_S[Ist_S>max(Ist_S)-THRESHOLD]
    #     Ist_S=Ist_S[Ist_S>max(Ist_S)-THRESHOLD]
    #     xsteady_S=np.log10(delta_S*10**epsst_S*(me*c**2/h)/(1+zsh))  #observing (SED) units
    #     ysteady_S=np.log10(10**Ist_S*delta_S**4/D**2*R_S*me*c**3/sigmaT/3)  #observing (SED) units
        
    if nbx==-1:
        nbx=LEN-1
    xsed=np.linspace(min(epsfile[1:LEN]),max(epsfile[1:LEN]),nbx)
    
    """BLACK BODIES (BB)"""
    #Correct for 3c273 (lessen the BB from bibliographical values)
    if oo==0: Gamma=2.36*Gamma
    
    eps85bb, Ibb =epsst ,epsst
    if not SSConly:
        if len(eps85bb)<200:
            eps85bb=np.linspace(min(epsst),max(epsst),200)
            Rext=Rbb[oo]
            Ibb=BB(eps85bb,T/Gamma,lext/Gamma**2,Rext)+np.log10(Gamma)
            xbb=eps85bb+np.log10((me*c**2/h))
            ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT)  #observing (SED) units
            bbmodel=inpl.interp1d(10**xbb,10**ybb)  #model BB
            #save BB in file
            BBfile=open('BB_p'+str(POWER)+'.txt','w')
            for i,j in zip(xbb,ybb):
                BBfile.write(str(i)+'\t'+str(j)+'\n')
            BBfile.close()
            if SECONDBB:
                    Ibb2=BB(eps85bb,Tbb2*Gamma2,(Gamma2)**2/Gamma**2*\
                            factbb*(Rext/Rbb2)**2*lext,Rbb2)+np.log10(Gamma)
                    xbb2=eps85bb+np.log10((me*c**2/h)/Gamma2)
                    ybb2=np.log10(10**Ibb2*Gamma2**-2/D**2*R*me*c**3/sigmaT)
                    bb2model=inpl.interp1d(10**xbb2,10**ybb2)
            #save BB2 in file
            BB2file=open('BB2_p'+str(POWER)+'.txt','w')
            for i,j in zip(xbb2,ybb2):
                BB2file.write(str(i)+'\t'+str(j)+'\n')
            BB2file.close()
                    
    
    """Generate optical Bands"""
    #optical bands: one bin and one after out of the range of optical
    iopt=0
    for ix in range(len(xsteady)-1):
        if xsteady[ix]>np.log10(range_optical[0]):
            locals()['opt'+str(iopt)]=10**xsteady[ix-1]
            iopt=iopt+1
        if xsteady[ix]>np.log10(range_optical[1]):
            locals()['opt'+str(iopt)]=10**xsteady[ix]
            break
        
    
    """ITERATIONS"""
    #PLOTTING SETUP
    if (SED!='no' or SAMPLE_FRAMES) and TIMELAPSE:
       lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),\
                                  round(min(vFv_pts)-HEIGHT),int(max(vFv_pts)+HEIGHT)]
       if REALHIST:
            size1=[10,5]
            if SED=='con':
                size1=[12,5]
            fig1,ax1=plt.subplots(1,2,figsize=size1,num=1, \
                                  sharey=True,gridspec_kw={'width_ratios': [3.5, 1]})
            ax10=ax1[0]
            plt.subplots_adjust(wspace=0.05)
            ax10.axis(lims)
            ax10.set_xlabel(r'$log\,('+xlab+'_{obs})\;\;'+' ['+xaxis+']$',fontsize=15,labelpad=-1)
            ax10.set_ylabel(r'$log\,('+yaxis+'_{obs})\;\;['+yunits[yaxis]+']$',fontsize=15)
                #ax10.set_title(obj+'  SED',fontsize=18)
       else:
            size2=[8,5]
            if SED=='con':
                size2=[10,5]
            fig1,ax1=plt.subplots(num=1,figsize=size2)
            ax10=ax1
            plt.axis(lims)
            plt.xlabel(r'$log\,('+xlab+'_{obs})\;\;['+xaxis+']$',fontsize=15,labelpad=-1)
            plt.ylabel(r'$log\,('+yaxis+'_{obs})\;\;['+yunits[yaxis]+']$',fontsize=15)
            #plt.title(obj+'  SED',fontsize=18)
                
                
    if SAMPLE_FRAMES:
        SED , TIMELAPSE , sflgd ='no' , False, []
    
    if SSConly:
        lgd2opts['theta']='(c)'
    lgd2=lgd2opts[namevar]
    ##set arrays for iteration
    legs=[]
    time , tobs, optical , xsoft , xmed,  xhard , VHE, VHE_counts, gbrnum ,alphax , VHEband_counts3,  VHEband_counts2,  VHEband_counts1, alphag,VHE_normfit, VHEband1 ,VHEband2 , VHEband3, VHEband_vFv3,\
    time_S , tobs_S, optical_S , xsoft_S ,xmed_S, xhard_S , VHE_S,VHE_counts_S, gbrnum_S ,alphax_S , VHEband_counts3_S,  VHEband_counts2_S,  VHEband_counts1_S,VHE_normfit_S, alphag_S, VHEband1_S, VHEband2_S, VHEband3_S, VHEband_vFv3_S,\
    time_E , tobs_E, optical_E , xsoft_E ,xmed_E, xhard_E , VHE_E, VHE_counts_E, gbrnum_E ,alphax_E ,VHEband_counts3_E,  VHEband_counts2_E, VHEband_counts1_E, VHE_normfit_E, alphag_E, VHEband1_E, VHEband2_E, VHEband3_E, VHEband_vFv3_E,\
    =  [np.zeros(nmax+1) for i in range(3*19)] #18 common tracing properties + for the main blob
    if DOUBLE_VAR:fkvar_sec = np.zeros(nmax+1)
    time_prev, time_prev_S = 0.0 , 0.0
    
    if MULTIoptical:
        for sopt in range(iopt+1):
                    locals()['optical'+str(sopt)]=np.zeros(nmax+1)
                    locals()['optical_S'+str(sopt)]=np.zeros(nmax_S+1)
                    locals()['optical_E'+str(sopt)]=np.zeros(nmax_E+1)
    if 'fakeTC' in name55:
        fktime, fkvar , tcr_arr =np.zeros(nmax+1) , np.zeros(nmax+1), np.zeros(nmax+1)
        if STATIONARY_STATE:fktime_S, fkvar_S , tcr_arr_S=np.zeros(nmax_S+1) , np.zeros(nmax_S+1) , np.zeros(nmax_S+1)
        if EXTRA_STATE:fktime_E, fkvar_E , tcr_arr_E=np.zeros(nmax_E+1) , np.zeros(nmax_E+1) , np.zeros(nmax_E+1)
    else:
        fktime, fkvar, tcr_arr = filetime, filevar, tcr_all
        if STATIONARY_STATE:fktime_S, fkvar_S, tcr_arr_S = filetime_S, filevar_S , tcr_all_S
        if EXTRA_STATE:fktime_E, fkvar_E, tcr_arr_E = filetime_E, filevar_E , np.ones(len(filevar_E))*tcr_E
    for pst in range(pastlength):
        locals()['xpst'+str(pst)] , locals()['ypst'+str(pst)] = np.array([]) , np.array([])
        locals()['xpst_S'+str(pst)] , locals()['ypst_S'+str(pst)] = np.array([]) , np.array([])
        locals()['xpst_E'+str(pst)] , locals()['ypst_E'+str(pst)] = np.array([]) , np.array([])
    #Write the SED output
    seds=open(route+'SEDs_p'+str(POWER)+'.txt','a')
    #tbb=open('fort.81tbb','w')
    if SAMPLE_FRAMES:
        frms=open(route+'FRAMES_p'+str(POWER)+'.txt','w')
    
    #THE MAIN LOOP
    ttemp=0.0
    prev_len=0
    tobs_prev= tfirst_dic[names[0].split('_')[1]]
    for nj in range(n1,nmax):
        LEN = lens[nj-n1+1] 
        eps85=epsfile[np.sum(lens[0:nj-n1+1]):np.sum(lens[0:nj-n1+1])+LEN]
        time[nj-n1+1]=eps85[0] #single fort.85 file
        I85=Ifile[np.sum(lens[0:nj-n1+1]):np.sum(lens[0:nj-n1+1])+LEN]
        if ELECTRONS:
            LEN89 = lens89[nj-n1+1]
            g89=gammafile[np.sum(lens89[0:nj-n1+1]):np.sum(lens89[0:nj-n1+1])+LEN89] #electrons later section
            N89=Nfile[np.sum(lens89[0:nj-n1+1]):np.sum(lens89[0:nj-n1+1])+LEN89]     #electrons later section
        R , B, p , logle, loggmin, loggmax, delta = R_all[nj-n1+1] , B_all[nj-n1+1] , p_all[nj-n1+1], logle_all[nj-n1+1], loggmin_all[nj-n1+1] , loggmax_all[nj-n1+1], delta_all[nj-n1+1]
        tinit , tfin = tfirst[nj-n1+1] , tlast[nj-n1+1]
        tcross=R/c/3600/24 #days jet frame
        tcr=tcross*dilate/delta
        tcr_arr[nj-n1+1]  = tcr 
        if time[nj-n1+1] not in time[0:nj-n1] and time[nj-n1+1]>0:
            ttemp += (time[nj-n1+1]  - time_prev)*tcr
        if STATIONARY_STATE:
            if nj>n1:
                nj_S =list(X_S).index(min(X_S, key=lambda t:abs(ttemp-t)))
            else:
                nj_S=0
            nj_S = nj  #non-time-dependent slow-fast component (the time-dependent must be commented out)
            LEN_S = lens_S[nj_S-n1+1]
            R_S , B_S , p_S  , logle_S , loggmin_S , loggmax_S , delta_S  = R_all_S[nj_S-n1+1] , B_all_S[nj_S-n1+1] , p_all_S[nj_S-n1+1], logle_all_S[nj_S-n1+1], loggmin_all_S[nj_S-n1+1] , loggmax_all_S[nj_S-n1+1], delta_all_S[nj_S-n1+1]
            tcross_S=R_S/c/3600/24 #days jet frame
            tcr_S=tcross_S*dilate/delta_S
            tcr_arr_S[nj_S-n1+1] = tcr_S
            eps85_S=epsfile_S[int((nj_S)*LEN_S):int((nj_S+1)*LEN_S)]
            I85_S=Ifile_S[int((nj_S)*LEN_S):int((nj_S+1)*LEN_S)]
            if ELECTRONS:
                LEN89_S = lens89_S[nj_S-n1+1]
                g89_S=gammafile_S[int((nj_S)*LEN89_S):int((nj_S+1)*LEN89_S)] #electrons later section
                N89_S=Nfile_S[int((nj_S)*LEN89_S):int((nj_S+1)*LEN89_S)] #electrons later section
            time_S[nj-n1+1]=eps85_S[0]#single fort.85 file
            time_prev_S = time_S[nj-n1+1]
        
        if EXTRA_STATE:
            nj_E = nj #it is not time-dependent itself (it depends on the conditions of time-dependent main blob)
            R_E, B_E , p_E  , logle_E , loggmin_E , loggmax_E , delta_E  = R_all_E[nj_E-n1+1] , B_all_E[nj_E-n1+1] , p_all_E[nj_E-n1+1], logle_all_E[nj_E-n1+1], loggmin_all_E[nj_E-n1+1] , loggmax_all_E[nj_E-n1+1], delta_all_E[nj_E-n1+1]
            tcross_E=R_E/c/3600/24 #days jet frame
            tcr_E=tcross_E*dilate/delta_E
        
        if len(names)>1: time[nj-n1+1]=filetime[nj-n1+1] 
        if time[nj-n1+1]<time[nj-n1]: tobs_prev =  tinit #it's because the segment has changed and t=0.0
        if time[nj-n1+1] not in time[0:nj-n1] and time[nj-n1+1]>0:        
            eps85 , I85 =eps85[1::] , I85[1::]                   
            """DATA REDUCTION"""    
            """Remove Extremely High Energy Gamma Rays/ Right Cut-off of Frequencies"""
            zz=np.where(eps85>nu_cutoff,eps85,0)
            zzz=zz.nonzero()
            for i in zzz: I85[i]=-100 #will be removed by THRESHOLD      
            """Cut off noise """
            bm=max(I85)-THRESHOLD
            eps85=eps85[I85>bm]
            I85=I85[I85>bm]            
            if ELECTRONS:
                g89, N89 = g89[1::], N89[1::]  #electrons later sectionn
                zz=np.where(g89>loggmax+3,g89,0)
                zzz=zz.nonzero()
                for i in zzz:
                    N89[i]=-100  #will be removed by THRESHOLD
                if len(g89)==0:
                    g89=g89n
                    N89=g89n**-2
                Nm=max(N89)-THRESHOLD/2
                g89=g89[N89>Nm]
                N89=N89[N89>Nm]
            if STATIONARY_STATE:
                eps85_S , I85_S =eps85_S[1::] , I85_S[1::]
                zz=np.where(eps85_S>nu_cutoff,eps85_S,0)
                zzz=zz.nonzero()
                for i in zzz: I85[i]=-100 #will be removed by THRESHOLD    
                eps85_S=eps85_S[I85_S>bm]
                I85_S=I85_S[I85_S>bm]
                if ELECTRONS:
                    g89_S, N89_S = g89_S[1::], N89_S[1::]  #electrons later sectionn
                    g89_S=g89_S[N89_S>Nm] #same cut as for main electrons
                    N89_S=N89_S[N89_S>Nm] #same
            
            if len(eps85)>1 and len(I85)>1:
                """Add extrapolated values at high energies """#for plotting reasons// avoid truncation
                eps85=np.append(eps85,(2*eps85[-1]-eps85[-2]))
                I85=np.append(I85,max(I85)-THRESHOLD)
                if STATIONARY_STATE:
                    I85_S=np.append(I85_S,max(I85_S)-THRESHOLD)
                    eps85_S=np.append(eps85_S,(2*eps85_S[-1]-eps85_S[-2]))
                    
                """Observed SED transformation"""
                x85=eps85+np.log10(delta*(me*c**2/h)/(1+zsh)) #observing (SED) units
                y85=I85+np.log10(delta**4/D**2*R*me*c**3/sigmaT/3)  #observing (SED) units
                ymodel=inpl.interp1d(x85, y85, kind='quadratic') #model as a function for the above arrays
                xtot, ytot, ymodel_tot = x85, y85, ymodel
                funcytot=inpl.interp1d(xtot,ytot,kind='slinear')                   
    
                ##Variable parameter
                if 'fakeTC' not in name55 or len(names)==1:
                    fktime[nj-n1+1]=filetime[nj-1]
                    fkvar[nj-n1+1]=filevar[nj-1]
                    if DOUBLE_VAR:fkvar_sec[nj-n1+1]=filevar_sec[nj-1]
                    if STATIONARY_STATE:
                        fktime_S[nj_S-n1+1]=filetime_S[nj_S-1]
                        fkvar_S[nj_S-n1+1]=filevar_S[nj_S-1]
                tobs[nj-n1+1] = round((time[nj-n1+1]- time_prev)*tcr,3)+tobs_prev
                time_prev = time[nj-n1+1]
                tobs_prev = tobs[nj-n1+1] 
                
                if STATIONARY_STATE and int(tinit) in days_with_lowstate:
                    #if time-dependent stationary state (and not stationary within a segment) uncomment below
                    #x85_S=eps85_S+np.log10(delta_S*(me*c**2/h)/(1+zsh)) #observing (SED) units
                    #y85_S=I85_S+np.log10(delta_S**4/D**2*R_S*me*c**3/sigmaT/3)  #observing (SED) units
                    #ymodel_S=inpl.interp1d(10**x85_S,10**y85_S) #model as a function for the above arrays
                    #xtot=x85_S[x85_S<max(x85)][x85_S[x85_S<max(x85)]>min(x85)]
                    #ytot = np.log10(ymodel(10**xtot))#+ymodel_S(10**xtot)) #in stationary_state claw now
                    #ymodel_tot  = inpl.interp1d(10**xtot,10**ytot) 
                    #tobs_S[nj-n1+1]=round(time_S[nj-n1+1]*tcr_S,3) #for channging stationary state
                    for i in range(len(nlf)):
                        if  tobs[nj-n1+1]>nlf[i][0]  and tobs[nj-n1+1]<=nlf[i][1]:
                            mjdday = int(np.sum(llf[i])/2)
                    xstat, ystat  , xstattbb  , ystattbb  = read_state(routeSS+name_stat_state+'_'+str(mjdday), lext = lext, T=T, delta=delta_S, R=R_S)                   
                    funcystat = inpl.interp1d(xstattbb, ystattbb)
                    xtota = xtot[xtot<max(xstattbb)][xtot[xtot<max(xstattbb)]>min(xstattbb)]
                    ytota = np.log10(10**ytot[xtot<max(xstattbb)][xtot[xtot<max(xstattbb)]>min(xstattbb)] + 10**funcystat(xtota))
                    xtot, ytot = xtota , ytota
                if EXTRA_STATE and not np.isnan(delta_E) and int(tinit) in days_with_extrastate:
                    if feedback_VHE_on_extrastate:
                        njex =list(extra_tc[:,0]).index(min(extra_tc[:,0], key=lambda t:abs(tobs[nj-n1+1]-t)))
                        boost_extrastate = extra_tc[njex,1]#(1.0+diffcolgamma[i][1]) #+1.0
                        delta_E_temp =  boost_extrastate
                    else:
                        for i in range(len(nlf)):
                            if  tobs[nj-n1+1]>nlf[i][0]  and tobs[nj-n1+1]<=nlf[i][1]:
                                boost_extrastate = (1.0+diffcolgamma[i][1]) #+1.0
                                delta_E_temp =  delta_E * boost_extrastate
                                #boost_extrastate =  10**(max(diffcolgamma[i][1]*0.25,0.0))
                    xex, yex  , xexbb  , yexbb  = read_state(routeSS+name_extra_state+'_'+str(int(tinit)), lext = lext, T=T, delta=delta_E_temp, R=R_E)                   
                    funcyex = inpl.interp1d(xexbb, yexbb)
                    xtotb = xtot[xtot<max(xexbb)][xtot[xtot<max(xexbb)]>min(xexbb)]
                    ytotb = np.log10(10**ytot[xtot<max(xexbb)][xtot[xtot<max(xexbb)]>min(xexbb)] +10**funcyex(xtotb))
                    xtot , ytot =xtotb, ytotb
                ymodel_tot  = inpl.interp1d(xtot, ytot, kind='quadratic')
                
                """BLACK BODY ADDITION"""    
                if not SSConly:
                    x81tbb=xbb[xbb<max(xtot)][xbb[xbb<max(xtot)]>min(xtot)]
                    y81tbb=np.log10(10**ybb[xbb<max(xtot)][xbb[xbb<max(xtot)]>min(xtot)]\
                                            +10**funcy85(x81tbb))
                    if SECONDBB:
                            funcybb2=inpl.interp1d(xbb2,ybb2)
                            y81tbb=np.log10(10**funcybb2(x81tbb)+10**y81tbb)
                    #for i,j in zip(x81tbb,y81tbb):
                           # tbb.write(str(i)+'\t'+str(j)+'\n')
                else:
                    x81tbb = xtot
                    y81tbb = ytot
                funcy81tbb=inpl.interp1d(x81tbb,y81tbb)
    
                ##Intergrate in bands (need adjustment)
                lng=np.linspace(range_VHE[0],min(range_VHE[1],10**max(xtot)),Nint)
                lng1=np.linspace(range_VHE1[0],min(range_VHE1[1],10**max(xtot)),Nint)
                lng2=np.linspace(range_VHE2[0],min(range_VHE2[1],10**max(xtot)),Nint)
                lng3=np.linspace(range_VHE3[0],min(range_VHE3[1],10**max(xtot)),Nint)
                if not MULTIoptical:
                    if min(lnopt)>10**min(xtot):
                        optical[nj-n1+1]=intg.simps(10**ymodel_tot(lnopt)/lnopt,lnopt)
                    else:
                        optical[nj-n1+1]=optical[nj-n1+1]
                else:
                    for sopt in range(iopt+1):
                        if locals()['opt'+str(sopt)]>10**min(xtot):
                            locals()['optical'+str(sopt)][nj-n1+1]=10**ymodel_tot(np.log10(locals()['opt'+str(sopt)]))*np.log10(opt2/opt1)
                        else:
                            locals()['optical'+str(sopt)][nj-n1+1]=locals()['optical'+str(sopt)][nj-n1]
                            print('Error with soft O/IR at:\t'+str(nj))
                if min(lnx)>10**min(xtot):
                    xsoft[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lnx))/lnx,lnx)
                else:
                    xsoft[nj-n1+1]=xsoft[nj-n1]
                    print('Error with soft X-rays at:\t'+str(nj))
                if min(lnxm)>10**min(xtot):
                    xmed[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lnxm))/lnxm,lnxm)
                else:
                    xmed[nj-n1+1]=xmed[nj-n1]
                    print('Error with medium X-rays at:\t'+str(nj))
                if min(lnxx)>10**min(xtot):
                    xhard[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lnxx))/lnxx,lnxx)
                else:
                    xhard[nj-n1+1]=xhard[nj-n1]
                    print('Error with hard X-rays at step:\t'+str(nj))
                #integrate vFv/v or vFv/v^2 in the energy range selected
                VHE[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lng))/lng,lng)
                VHEband1[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lng1))/lng1,lng1)
                VHEband2[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lng2))/lng2,lng2)
                VHEband3[nj-n1+1]=intg.simps(10**ymodel_tot(np.log10(lng3))/lng3,lng3)
                VHE_counts[nj-n1+1] = intg.simps(10**ymodel_tot(np.log10(lng))/lng**2,lng)/h
                VHEband_counts1[nj-n1+1] = intg.simps(10**ymodel_tot(np.log10(lng1))/lng1**2,lng1)/h
                VHEband_counts2[nj-n1+1] = intg.simps(10**ymodel_tot(np.log10(lng2))/lng2**2,lng2)/h
                VHEband_counts3[nj-n1+1] = intg.simps(10**ymodel_tot(np.log10(lng3))/lng3**2,lng3)/h
                VHEband_vFv3[nj-n1+1] = 10**ymodel_tot(np.log10(range_VHE3[0]))
                resultx = curve_fit(logparabola_x, lnx, 10**ymodel_tot(np.log10(lnx)), p0 = [-8.0,-1.5]) 
                xsolx,dx = resultx
                alphax[nj-n1+1] = xsolx[1] #np.log10(ymodel_tot(x1)/x1**2/10**ymodel_tot(x2)*x2**2)/np.log10(x1/x2)
                try:
                    resultg = curve_fit(logparabola_g, lng, 10**ymodel_tot(np.log10(lng)), p0 = [-8.0,-0.5])
                    xsolg,dx = resultg
                except:
                    print('\n\n Error at curve-fitting step step'+str(nj))
                    pass
                VHE_normfit[nj-n1+1] = 10**xsolg[0]
                alphag[nj-n1+1] = xsolg[1]#np.log10(10**ymodel_tot(x1)/x1**2/10**ymodel_tot(x2)*x2**2)/np.log10(x1/x2)]
                print('Step of integration:\t '+str(nj)+'/'+str(nmax))
                
                """SAVE in .txt files the model"""
                if SED=='con':
                    xa=xsed+np.log10(delta*(me*c**2/h)/(1+zsh)) #observing (SED) units
                    xb=xa[xa>min(xtot)+0.01]
                    xc=xb[xb<max(xtot)-0.01]
                    sedzip = zip(xc,ymodel_tot(xc))
                else:
                    sedzip = zip(xtot,ytot)
                seds.write(str(tobs[nj-n1+1])+'\n')
                dt_temp = tobs[nj-n1+1]- tobs[nj-n1] 
                snapshot=open(route+'SEDmodel_{}_^{}_{}-{}.txt'.format(namevar+DOUBLE_VAR*('-'+namevar_sec),POWER,round(tobs[nj-n1]+dt_temp/2,3),round(tobs[nj-n1+1]+dt_temp/2,3)),'a')
                for i,j in sedzip:
                    seds.write(str(i)+'\t'+str(j)+'\n') #logscale
                    snapshot.write(str(10**i)+'\t'+str(10**j)+'\n') #linear scale, no time title
                snapshot.close()
    
                """PLOTTING"""
                if TIMELAPSE:
                        fig1,ax1  = plt.subplots(num=1)
                        plt.xlabel(r'$'+xlab+'_{obs}\;\;'+' ['+xaxis+']$',fontsize=15)
                        plt.ylabel(r'$'+yaxis+'_{obs}\;\;['+yunits[yaxis]+']$',fontsize=15)
                        plt.plot(x85,y85,'b-',linewidth=1.0,label='fast')
                        #steady secondary components (not time-dependent input)
                        #plt.plot(x85_S,y85_S,'g-',linewidth=1.0,label='slow')
                        #plt.plot(x85_E,y85_E,'g-',linewidth=1.0,label='fast')
                        for i in range(len(nlf)):
                            if  tobs[nj-n1+1]>nlf[i][0] -dt_ELAPSE and tobs[nj-n1+1]<=nlf[i][1] -dt_ELAPSE:
                                mjdday = int(np.sum(llf[i])/2)
                        if STATIONARY_STATE:
                            plt.plot(xstattbb,ystattbb,'c--',linewidth=0.5, label='slow NT')
                        if EXTRA_STATE and int(tinit) in days_with_extrastate:
                            plt.plot(xexbb,yexbb,'r--',linewidth=1.0,label='fast NT')
                        if not SSConly:
                            plt.plot(xbb,ybb,'b--',linewidth=0.5,label='BLR')
                        if SECONDBB:
                            plt.plot(xbb2,ybb2,'b-.',linewidth=0.5,label='Disk')
                        plt.plot(xtot,ytot,'k-',label='Total t='+str(round(tobs[nj-n1+1],3))+' d',lw=2)
                        plt.legend(framealpha=0.5,loc='upper left')
                        plt.title(obj+'  SED',fontsize=18)
                        lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),\
                                  round(min(vFv_pts)-HEIGHT),int(max(vFv_pts)+HEIGHT)]
                        plt.axis(lims)
                        if MULTIcolor:
                            for i in range(len(dtpts)):
                                plt.errorbar(locals()['v_pts'+str(i)],locals()['vFv_pts'+str(i)],\
                                         locals()['errorv'+str(i)],locals()['errorvFv'+str(i)],\
                                         elinewidth=1.5,capsize=2.5,capthick=0.5,\
                                         fmt=form[i],color=colors[i],ms=3.5)
                        else:
                            dcol = colors[0]
                            if TIME_DEPENDENT_DATA:
                                for i in range(len(nlf)):
                                    if  tobs[nj-n1+1]>nlf[i][0] -dt_ELAPSE and tobs[nj-n1+1]<=nlf[i][1] -dt_ELAPSE:
                                        mjdday = int(np.sum(llf[i])/2)
                                        namedata  = route_timelapse+'sed_'+str(llf[i][0])+'_'+str(llf[i][1])+'.txt'
                                        tdata = np.loadtxt(namedata)
                                        v_pts  = tdata[:,0]
                                        vFv_pts  = tdata[:,1]
                                        errorvFv  = tdata[:,2]
                                        if LOGdata:
                                            v_pts = np.log10(v_pts)
                                            vFv_pts = np.log10(vFv_pts)
                                            errorvFv = errorvFv/10**vFv_pts/np.log(10)
                                dcol = colors_daysplit[str(mjdday)]
                            plt.errorbar(v_pts[errorvFv!=0.0],vFv_pts[errorvFv!=0.0],yerr=errorvFv[errorvFv!=0.0],elinewidth=1.5,\
                                capsize=2.5,capthick=0.5,fmt='.',color = dcol , ecolor=dcol,ms=3.5)
                            plt.errorbar(v_pts[errorvFv==0.0], vFv_pts[errorvFv==0.0], yerr=0.0,elinewidth=1.5,\
                                capsize=2.5,capthick=0.5,fmt='v',color = dcol , ecolor=dcol,ms=5.5)
                        # plt.legend(legend_pts)
                        #TRACEBACK CURVES
                        xpst0=x85
                        ypst0=y85                    
                        for pst in range(pastlength-1,int(pastlength/2),-1):
                              if len(locals()['xpst'+str(pst)])!=0:
                               plt.figure(1)
                               plt.plot(locals()['xpst'+str(pst)],locals()['ypst'+str(pst)],'b-',linewidth=lnwdth2)
        
                        for pst in range(int(pastlength/2),0,-1):
                              if len(locals()['xpst'+str(pst)])!=0:
                               plt.figure(1)
                               plt.plot(locals()['xpst'+str(pst)],locals()['ypst'+str(pst)],'b-',linewidth=lnwdth1)
                        xtemp=[xpst0]
                        ytemp=[ypst0]
                        for pst in range(1,pastlength+1):
                            xtemp.append(locals()['xpst'+str(pst-1)])
                            ytemp.append(locals()['ypst'+str(pst-1)])
                        for pst in range(pastlength):
                            locals()['xpst'+str(pst)]=xtemp[pst]
                            locals()['ypst'+str(pst)]=ytemp[pst]
                            
                        
                        #TRACEBACK CURVES
                        # xpst_S0=x85_S
                        # ypst_S0=y85_S                   
                        # for pst in range(pastlength-1,int(pastlength/2),-1):
                        #        if len(locals()['xpst_S'+str(pst)])!=0:
                        #         plt.figure(1)
                        #         plt.plot(locals()['xpst_S'+str(pst)],locals()['ypst_S'+str(pst)],'g-',linewidth=lnwdth2)
         
                        # for pst in range(int(pastlength/2),0,-1):
                        #        if len(locals()['xpst_S'+str(pst)])!=0:
                        #         plt.figure(1)
                        #         plt.plot(locals()['xpst_S'+str(pst)],locals()['ypst_S'+str(pst)],'g-',linewidth=lnwdth1)
                        # xtemp=[xpst_S0]
                        # ytemp=[ypst_S0]
                        # for pst in range(1,pastlength+1):
                        #      xtemp.append(locals()['xpst_S'+str(pst-1)])
                        #      ytemp.append(locals()['ypst_S'+str(pst-1)])
                        # for pst in range(pastlength):
                        #      locals()['xpst_S'+str(pst)]=xtemp[pst]
                        #      locals()['ypst_S'+str(pst)]=ytemp[pst]
                        #saveframe
                        plt.xlim(9.0,28.0)
                        if saveframes:
                           plt.figure(1).savefig(routesave+'vFv_'+str(nj)+'.jpg',bbox_inches="tight")
                        if pausetime>0.0:
                            plt.pause(pausetime)
                        plt.clf()
                elif SED=='on':
                    if REALHIST:
                        ax10.plot(x85,y85,'-',color=onecolor,alpha=0.5,linewidth=0.6-0.1*(nmax//1000))
                    else:
                        plt.figure(1)
                        if onecolor:
                            plt.plot(x85,y85,'-',color=onecolor,alpha=0.5,linewidth=0.6-0.1*(nmax//1000))
                            #plt.plot(np.log10(lng),np.log10(logparabola_g(lng,xsolg[0], xsolg[1])),lw=2)
                            #plt.plot(np.log10(lnx),np.log10(logparabola_x(lnx,xsolx[0], xsolx[1])),lw=2)
                        else:
                                plt.plot(x85,y85,'-',linewidth=0.5)
                                
                    #plt.ylim(-10.0,-4.0)
        
                ##MULTI WINDOW DIAGRAM      
                if TIMELAPSE:
                  dimh = int(1*monochromVIR +1*radio +1*Xrays +1*gammarays)
                  figall=plt.figure(2, figsize=[12,5])
                  figall.suptitle(obj+'  t='+str(round(tobs[nj-n1+1],2))+' days',fontsize=16)
                  addim=0
                  if monochromVIR:
                    
                    nu85=np.log10(delta*10**eps85*(me*c**2/h)/(1+zsh))
                    Fv85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3)-nu85
                    
                    #plt.figure(6)
                    if dimh==1: 
                        plt.subplot()
                    else:
                        plt.subplot(int(dimh//1.5),2,addim)
                        addim+=1
                    plt.xlabel(r'$\nu_{obs}\;[Hz]$')
                    plt.ylabel(r'$ F\nu_{obs}\; [erg/cm^2/Hz/s]$')
                    plt.plot(nu85,Fv85,'y-')
                    plt.title('VIR')
                    if Multicolor:
                      for i in range(len(dtpts)):
                        xp=locals()['v_pts'+str(i)]
                        yp=locals()['vFv_pts'+str(i)]-locals()['v_pts'+str(i)]
                        #erxp=locals()['errorv'+str(i)]
                        eryp=locals()['errorvFv'+str(i)]
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
                    else:
                        xp=v_pts
                        yp=vFv_pts-v_pts
                        #eryp=errorv
                        eryp=errorvFv
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt='r.',ecolor='red',ms=3.5)
                    limlow=12.0
                    limhigh=16.5
                    rangey=   yp[xp<limhigh][xp[xp<limhigh]>limlow]
                    if len(rangey)>1:
                        minyx , maxyx = min(rangey)-HEIGHT, max(rangey)+HEIGHT
                        lims=[limlow,limhigh,  minyx , maxyx]
                        plt.axis(lims)
                    else:
                        plt.xlim(limlow,limhigh)
                        plt.ylim(minyx , maxyx)
                    plt.axis(lims)  
                    
                    
                  if radio:        
                    nu5=np.log10(delta*10**eps85*(me*c**2/h)/(1+zsh))
                    R85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3/Jy)-nu5
                    
                    #plt.figure(10)
                    if dimh==1: 
                        plt.subplot()
                    else:
                        plt.subplot(int(dimh//1.5),2,addim)
                        addim+=1
                    plt.xlabel(r'$\nu_{obs}\;\;[Hz]$')
                    plt.ylabel(r'$F\nu\;\; [Jy]$')
                    plt.plot(nu5,R85,'r-')
                    plt.title('Radio')
                    if Multicolor:
                      for i in range(len(dtpts)):
                        xp=locals()['v_pts'+str(i)]
                        yp=locals()['vFv_pts'+str(i)]+np.log10(1/Jy)-locals()['v_pts'+str(i)]
                        #erxp=locals()['errorv'+str(i)]
                        eryp=locals()['errorvFv'+str(i)]
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
                    else:
                        xp=v_pts
                        yp=vFv_pts+np.log10(1/Jy)-v_pts
                        #eryp=errorv
                        eryp=errorvFv
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt='r.',ecolor='red',ms=3.5)   
                    limlow=8.0
                    limhigh=12.5
                    rangey=   yp[xp<limhigh][xp[xp<limhigh]>limlow]
                    if len(rangey)>1:
                        minyx , maxyx = min(rangey)-HEIGHT, max(rangey)+HEIGHT
                        lims=[limlow,limhigh,  minyx , maxyx]
                        plt.axis(lims)
                    else:
                        plt.xlim(limlow,limhigh)
                        plt.ylim(minyx , maxyx)
                    plt.axis(lims)            
                              
                  if Xrays:
                      
                    eps5=np.log10(delta*10**eps85/(1+zsh))
                    xx5 = eps5+np.log10(1/eV/10**3*me*c**2)
                    X85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3)-eps5
                    #plt.figure(8)
                    if dimh==1: 
                        plt.subplot()
                    else:
                        plt.subplot(int(dimh//1.5),2,addim)
                        addim+=1
                    plt.xlabel(r'$E_{obs}$ [keV]',size=15)
                    plt.ylabel(r'$F_X$  counts [#/cm$^2$/s]', size=15)
                    plt.plot(xx5,X85,'b-',lw=2.0)
                    plt.title('X-rays')
                    if MULTIcolor:
                      for i in range(len(dtpts)):
                        xp=locals()['v_pts'+str(i)]+np.log10(h/eV/10**3*me*c**2)
                        yp=locals()['vFv_pts'+str(i)]+np.log10(1/h)-locals()['v_pts'+str(i)]
                        #erxp=locals()['errorv'+str(i)]
                        eryp=locals()['errorvFv'+str(i)]
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt=form[i],color=colors[i],ms=6.5)
                    else:
                        xp=v_pts+np.log10(h/eV/10**3)
                        yp=vFv_pts+np.log10(1/h)-v_pts
                        #eryp=errorv
                        eryp=errorvFv
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt='r.',ecolor='red',ms=6.5)
                    limlow= -0.5 #un
                    limhigh=1.75
                    rangey=   yp[xp<limhigh][xp[xp<limhigh]>limlow]
                    if len(rangey)>1:
                        minyx , maxyx = min(rangey)-HEIGHT, max(rangey)+HEIGHT
                        lims=[limlow,limhigh,  minyx , maxyx]
                        plt.axis(lims)
                    else:
                        plt.xlim(limlow,limhigh)
                        plt.ylim(minyx , maxyx)
                    
                  if gammarays:
                    eps5=np.log10(delta*10**eps85/(1+zsh))
                    tev5=eps5+np.log10(1/eV/10**12*me*c**2)
                    T85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3)-eps5
                    
                    #plt.figure(12)
                    if dimh==1: 
                        plt.subplot()
                    else:
                        plt.subplot(int(dimh//1.5),2,addim)
                        addim+=1
                    plt.xlabel(r'$E_{obs}$ [TeV]', size=15)
                    plt.ylabel(r'$F_\gamma$  [#/cm$^2$s]',size=15)
                    plt.plot(tev5,T85,'k-',lw=2.0)
                    plt.title('$\gamma$ - rays')
                    if MULTIcolor:
                      for i in range(len(dtpts)):
                        xp=locals()['v_pts'+str(i)]+np.log10(h/eV/10**12)
                        yp=locals()['vFv_pts'+str(i)]+np.log10(1/h)-locals()['v_pts'+str(i)]
                        #erxp=locals()['errorv'+str(i)]
                        eryp=locals()['errorvFv'+str(i)]
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt=form[i],color=colors[i],ms=6.5)
                    else:
                        xp=v_pts+np.log10(h/eV/10**12)
                        yp=vFv_pts+np.log10(1/h)-v_pts
                        #erxp=errorv
                        eryp=errorvFv
                        plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
                                      capthick=0.5,fmt='r.',ecolor='red',ms=6.5)
                    xp=v_pts+np.log10(h/eV/10**12)
                    yp=vFv_pts+np.log10(1/h)-v_pts
                    limlow=-1.0
                    limhigh=0.5
                    rangey=   yp[xp<limhigh][xp[xp<limhigh]>limlow]
                    if len(rangey)>1:
                        minyg , maxyg = min(rangey)-HEIGHT, max(rangey)+HEIGHT
                        lims=[limlow,limhigh,  minyg , maxyg]
                        plt.axis(lims)
                    else:
                        plt.xlim(limlow,limhigh)  
                        plt.ylim(minyg , maxyg)
        
                if TIMELAPSE:
                    plt.pause(pausetime)
                    figall.subplots_adjust(hspace=0.6, wspace=0.5)
                    if saveframes:
                        figall.savefig(routesave+'bands_'+str(nj)+'.jpg',bbox_inches="tight")
                    plt.clf()
                    
                if SAMPLE_FRAMES:
                    nnw=nj+0 #offset to catch specific point
                    if nnw%int(dtsmpl/tcr)==0:
                #    if nnw in [238,470,661,856,1119,1320]:
                        plt.figure(1)
                        xn85=np.linspace(min(x85),max(x85),300)
                        plt.plot(xn85,funcy85(xn85),'-',linewidth=1.3)
                        sflgd.append(str(round(nnw*tcr+initt,1))+' days')
                        print('Frame index to check in lightcurves/ statitstics:\t '+str(nj))
                        if not SSConly:
                            frms.write(str(time[nj-n1+1])+'\n')
                            for i,j in zip(x81tbb,y81tbb):
                                frms.write(str(i)+'\t'+str(j)+'\n')
        
                    
                """ ELECTRONS"""     
                if ELECTRONS:               
                    """Cooling Breaking Calc"""
                    if namevar in ['B','lext'] and len(g89)>1:
                        dif1=N89[1::]-N89[0:-1]
                        dif1mod=inpl.interp1d(g89[1::],dif1)
                        g89n=np.linspace(g89[1],g89[-1],50)
                        dif1=dif1mod(g89n)
                       ###Using 2nd derivative for finding breaking point (alter loop below though dif1->dif2)
                       #dif2=np.array(dif1[0:-1]-dif1[1::])
                       #dif2mod=inpl.interp1d(g89[1:-1],dif2)
                       #dif2=dif2mod(g89n)
                       # dec89=0.2 #declination of dN89/dg (change of slope) for which to assume the start of the breaking
                       # c89=0.81 #offset to higher gamma for actual gbr value (value found from calibration)
                       #c89 is found by the difference of the theoretical and analytical value:
                        for idif in dif1:
                            if idif<(3-2*p)*0.1:
                                gbrnum[nj-n1+1]=g89n[list(dif1).index(idif)]#+c89
                                if gbrnum[nj-n1+1]>g89[-1]:
                                    gbrnum[nj-n1+1]=g89[-1]
                                break
                            if fkvar[nj-n1+1]>locals()[namevar]:
                                gbrnum[nj-n1+1]=g89n[0]
                            else:
                                gbrnum[nj-n1+1]=g89n[-1]
                    
                    if TIMELAPSE:
                        plt.figure(89)
                        plt.plot(g89,N89+(p-2)*g89)
                        if namevar in ['B','lext']:
                            plt.plot(gbrnum[nj-n1+1]*np.ones(50),np.linspace(-10.0,0.0,50),'k--')
                        mingobs=loggmin-1.0
                        if mingobs<0.: mingobs=0.
                        plt.axis([mingobs,loggmax+0.5,-10.,0.0])
                        plt.xlabel(r' $log\,\gamma$')
                        plt.ylabel(r'$\gamma^{p-2}N(\gamma)$')
                        plt.legend([r' $t= '+str(round(tobs[nj-n1+1],1))+'\; days$',r'$\gamma^{num}_{br}$'])
                        plt.pause(pausetime)
                        if saveframes:
                            plt.figure(89).savefig(routesave+'electrons_'+str(nj)+'.jpg')
                        plt.clf()
        else:
            time[nj-n1+1]=0.0
        
    #tbb.close()
    seds.close()
    
    #remove zeros
    if MULTIoptical:
        for sopt in range(iopt+1):
                    locals()['optical'+str(sopt)]=locals()['optical'+str(sopt)][tobs>0]
    else:
        optical=optical[tobs>0]
    VHE, VHEband1,  VHEband2, VHEband3,VHE_normfit , alphag =VHE[tobs>0], VHEband1[tobs>0],  VHEband2[tobs>0], VHEband3[tobs>0] ,VHE_normfit[tobs>0] , alphag[tobs>0]
    VHE_counts, VHEband_counts1, VHEband_counts2, VHEband_counts3 = VHE_counts[tobs>0], VHEband_counts1[tobs>0], VHEband_counts1[tobs>0], VHEband_counts3[tobs>0]
    xsoft , xmed , xhard , alphax =xsoft[tobs>0] , xmed[tobs>0] ,xhard[tobs>0], alphax[tobs>0]
    if 'fakeTC' not in name55:
        tcr_arr = tcr_arr[tobs>0]
        fktime=fktime[tobs>0]+initt/tcr_arr
        fkvar=fkvar[tobs>0]
        if DOUBLE_VAR:fkvar_sec=fkvar_sec[tobs>0]
    if namevar=='B':
        lb=fkvar**2/8/np.pi/me/c**2*sigmaT*R
        gbr=np.log10(3/4/(lb+lext))
        if ELECTRONS:gbrnum=gbrnum[tobs>0]
    elif namevar=='lext':
        lb=B**2/8/np.pi/me/c**2*sigmaT*R
        gbr=np.log10(3/4/(fkvar+lb))
        if ELECTRONS:gbrnum=gbrnum[tobs>0]
    tobs=tobs[tobs>0]+initt
    Time=max(time)+initt/tcr
    
    """X/gamma-rays selection"""
    
    if ixband==2:
        xband=xhard
        range_xband=range_xhard
    elif ixband==1:
        xband=xmed
        range_xband=range_xmed
    else:
        xband=xsoft
        range_xband=range_xsoft
        
    
    if iVHEband==0:
        VHEband  = VHEband1
        VHEband_counts  = VHEband_counts1
        range_gband=range_VHE1
    elif iVHEband==1:
        VHEband  = VHEband2
        VHEband_counts = VHEband_counts2
        range_gband=range_VHE2
    elif iVHEband==2:
        VHEband  = VHEband3
        VHEband_counts = VHEband_counts3
        range_gband=range_VHE3
        
    """optical filter selection"""
    #optical filter selection
    if MULTIoptical:
        #optband=1 #from array of created optical channels based on code bins
        optfq=str(round(np.log10(locals()['opt'+str(optband)]),1))
        optical=locals()['optical'+str(optband)]
    else:
        optfq=''
    
    
    """Add Noise and BB radiation to optical band"""
    if ADD_BB and not SSConly:
        noise=np.ones(len(optical))
        if addnoise:
                noise=np.random.poisson((100)**2,len(optical))/(100)**2
        lnopt=np.linspace(range_optical[0],range_optical[1],Nint)
        optical=optical+noise*(intg.simps(bbmodel(lnopt)/lnopt,lnopt)+intg.simps(bb2model(lnopt)/lnopt,lnopt))
        lnx=np.linspace(range_xband[0],range_xband[1],Nint)
        xband=xband+intg.simps(bbmodel(lnx)/lnx,lnx)+intg.simps(bb2model(lnx)/lnx,lnx)
        if namevar=='lext':
                optical=optical+noise*(intg.simps(bbmodel(lnopt)/lnopt,lnopt)*fkvar/np.mean(fkvar)+intg.simps(bb2model(lnopt)/lnopt,lnopt))
                xband=xband+intg.simps(bbmodel(lnx)/lnx,lnx)*fkvar/np.mean(fkvar)+intg.simps(bb2model(lnx)/lnx,lnx)    
    ci=[]
    filterA, filterB= [], []
    for t in range(len(tobs)):
        Aa,Bb=[],[]
        for sopt in range(iopt+1):
            Aa.append(locals()['opt'+str(sopt)]) #nuFnu
            Bb.append(locals()['optical'+str(sopt)][t]/np.log10(opt2/opt1)) #nuFnu
        optmodel=inpl.interp1d(Aa,Bb)
        af=optmodel(optdic[cibnda])*FWHMdic[cibnda]
        bf=optmodel(optdic[cibndb])*FWHMdic[cibndb]
        filterA.append(-2.5*np.log10(af/(optdic[cibnda]*ZPdic[cibnda]*Jy*FWHMdic[cibnda])))
        filterB.append(-2.5*np.log10(bf/(optdic[cibnda]*ZPdic[cibndb]*Jy*FWHMdic[cibnda])))
        ci.append(filterA[-1]-filterB[-1])
    filterA=np.array(filterA).astype(float)
    filterB=np.array(filterB).astype(float)
    ci=np.array(ci).astype(float)
    #label color
    cilbl='('+cibnda+' - '+cibndb+') '      
    #%%
    """NEW BINNING for output LCs /  Default -1 for code's output binning"""
    #NEW LCs save as %nane%2
    #defaults if tbin = -1
    arn=arn=['xsoft','xmed','xhard','xband','VHE','VHE_counts','VHEband','VHEband_counts',\
         'VHEband_counts1', 'VHEband_counts2','VHEband_counts3','VHEband_vFv3',\
         'VHEband1','VHEband2','VHEband3','VHE_normfit',\
         'VHEband_counts', 'optical','fkvar','fktime'] +DOUBLE_VAR*['fkvar_sec']
    if CI:
        arn.append('ci')
        arn.append('filterA')
        arn.append('filterB')
        arn.append('alphax')
        arn.append('alphag')
    if namevar in ['B','lext']:
        arn.append('gbr')
        if ELECTRONS:arn.append('gbrnum')
    arn.append('tobs')
    for ar in arn:
        locals()[ar+'B']=locals()[ar]
    #change if different binning selected
    binsB = len(tobs)
    if tbin!=-1:
        dt = tobs[1::]-tobs[0:-1]
        dtold = max(dt[dt>0.0])
        #slices = np.linspace(0,max(tobs), int((max(tobs)-1)/tbin)+1, True).astype(np.int)
        binsB=int((max(tobs)-min(tobs))/tbin)+1
        pcntg=[] #for generating randomness in optical obs while binning
        for ar in arn:
            locals()[ar+'B']=st.binned_statistic(tobs,locals()[ar],statistic='mean',bins=binsB)[0]
            #sample randomly from the observations of one day for optical band to create opticalB , ciB
            if ar in ['optical','ci']:
                stmin=st.binned_statistic(tobs,locals()[ar],statistic='min',bins=binsB)[0]
                stmax=st.binned_statistic(tobs,locals()[ar],statistic='max',bins=binsB)[0]
                if len(pcntg)==0:
                    pcntg=np.random.random(binsB)
                locals()[ar+'B']=stmin+(stmax-stmin)*pcntg 
            #periodic condition for averaging, else we end up with huge last value and neglegible first value of timeseries
            #firlas=np.mean(locals()[ar])
            #locals()[ar][0]=firlas
            #locals()[ar][-1]=firlas
        taubin=tbin
        segment_size=int(tbin/dtold*segment_size)
        tbin=round(tobsB[1]-tobsB[0],2)

    
    if SAVE:
        """"Write file""" #does not consider new binning
        np.savetxt('MAGIC_p'+str(POWER)+'.txt', np.c_[tobsB,VHEB])
        np.savetxt('MAGIC_counts_p'+str(POWER)+'.txt', np.c_[tobsB,VHE_countsB])
        np.savetxt('MAGIC_0.2-0.4'+str(POWER)+'.txt', np.c_[tobsB,VHEband1B])
        np.savetxt('MAGIC_0.4-0.8'+str(POWER)+'.txt', np.c_[tobsB,VHEband2B])
        np.savetxt('MAGIC_0.8plus'+str(POWER)+'.txt', np.c_[tobsB,VHEband3B])
        np.savetxt('MAGIC_counts_0.2-0.4'+str(POWER)+'.txt', np.c_[tobsB,VHEband_counts1B])
        np.savetxt('MAGIC_counts_0.4-0.8'+str(POWER)+'.txt', np.c_[tobsB,VHEband_counts2B])
        np.savetxt('MAGIC_counts_0.8plus'+str(POWER)+'.txt', np.c_[tobsB,VHEband_counts3B])
        np.savetxt('MAGIC_normfit_p'+str(POWER)+'.txt', np.c_[tobsB,VHE_normfitB])
        if SAMPLE_FRAMES:
            frms.write('STEADY\n')
            for i,j in zip(xsteady,ysteady):
                frms.write(str(i)+'\t'+str(j)+'\n')
            frms.close()
        np.savetxt('optical_p'+str(POWER)+'.txt',np.c_[tobsB,locals()['opticalB']])
        np.savetxt('X_SOFT_p'+str(POWER)+'.txt', np.c_[tobsB,xsoftB])
        np.savetxt('X_MED_p'+str(POWER)+'.txt', np.c_[tobsB,xmedB])
        np.savetxt('X_HARD_p'+str(POWER)+'.txt', np.c_[tobsB,xhardB])
        if CI:
            np.savetxt('CI_p'+str(POWER)+'.txt', np.c_[tobsB,filterAB,filterBB,ciB])
            np.savetxt('alphax_p'+str(POWER)+'.txt', np.c_[tobsB,alphaxB])
            np.savetxt('alphag_p'+str(POWER)+'.txt', np.c_[tobsB,alphagB])
        np.savetxt('fakeTCB_p'+str(POWER)+'.txt', np.c_[tobsB,fkvarB])
        if DOUBLE_VAR: np.savetxt('fakeTCB_sec_p'+str(POWER)+'.txt', np.c_[tobsB,fkvar_secB])
        
        print('TIME CURVES CREATED AND SAVED')
    #%% continue with plotting
    if not TIMELAPSE and (SED!='no' or SAMPLE_FRAMES):
            """Contour Plot"""
            plt.figure(1)
            if SED=='con':
            # 2-D Histogram --> to be commented out
                seds=open('SEDs_p'+str(POWER)+'.txt','r')
                ss, xseds, yseds=[],[],[]
                for j in seds:
                    ss.append(j)
                for j in ss:
                    if len(j.split())>1:
                          xseds.append(j.split()[0])
                          yseds.append(j.split()[1])
                seds.close()
                xseds=np.array(xseds).astype(float)
                yseds=np.array(yseds).astype(float)
                xbin=xseds[3]-xseds[2]
    
                rangecon=[lims[0:2],lims[2::]]
                heatmap, xedges, yedges =np.histogram2d(xseds,\
                                yseds,bins=[nbx,nby],range=rangecon,density=True)
                extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                #cmap=cmps[namevar]
                cmap=plt.cm.GnBu
                cmap.set_bad(color='white')
                cmap.set_under(color='white')
                #create contour levels
                lvs=np.linspace(0.0,toplvl,Nlvs).round(decimals=2)
                maxlvs=lvs[-1]
                
                c10=ax10.contourf(heatmap.T,lvs,cmap=cmap,extent=extent,extend='both')#,interpolation='bilinear')
                cbar=plt.colorbar(c10,ax=[ax10],location='left',pad=0.15)
                #cbar.ax.set_yticklabels(np.round(100-lvs*100,3))
                cbar.set_label('Number Density of States (norm=1)',size=14,rotation=90)
                plt.axis(lims)
    
            
            ax10.plot(xsteady,ysteady,'k--',linewidth=2.0,alpha=0.7,label='Steady')
            if MULTIcolor:                     
                for i in range(len(dtpts)):
                    ax10.errorbar(locals()['v_pts'+str(i)],locals()['vFv_pts'+str(i)],\
                              yerr=locals()['errorvFv'+str(i)],\
                              elinewidth=1.5,capsize=2.5,capthick=0.5,ecolor=colors[i],\
                              fmt=form[i],color=colors[i],ms=5.0,fillstyle=fillstyle[i],\
                              label=legend_pts[i])
            else:
                ax10.errorbar(v_pts,vFv_pts,yerr=errorvFv,elinewidth=1.5,capsize=2.5,\
                             fmt='r.',ecolor='red',ms=3.5,label='Obs')
            if not SSConly:
                ax10.plot(xbb,ybb,'r--',linewidth=0.9,label='BLR')
            if SECONDBB:
                ax10.plot(xbb2,ybb2,'r-.',linewidth=0.9,label='Disk')
            lg4=ax10.legend( loc='upper left',framealpha=0.5)
            if SAMPLE_FRAMES and nmax/dtsmpl<10:
                    ax10.add_artist(plt.legend(sflgd,framealpha=0.5))
            if lgd2!=-1:
                    ax10.annotate(lgd2, xy=(0.9,0.9), xycoords='axes fraction', size=14,weight='bold')
            ax10.add_artist(lg4)
            ax10.xaxis.set_minor_locator(AutoMinorLocator())
            ax10.yaxis.set_minor_locator(plt.LinearLocator(4*(lims[3]-lims[2])+1))
            if namevar=='lext':
                xTbb=np.log10(kb*T/me/c**2*np.array([0.3,9.0]))+np.log10(delta*(me*c**2/h)/(1+zsh))
                ax10.fill_between(np.linspace(xTbb[0], xTbb[1]),-20,20, facecolor='red', alpha=0.2, hatch='/',edgecolor='red')
                
            if 'on' in SED:
                    #plot scale errors
                wbd1, wbd2, wbdx, wbdf  =  FWHMdic[cibndb], FWHMdic[cibnda],\
                    np.log10(range_xhard[1]/range_xhard[0])/2,\
                    np.log10(range_VHE[1]/range_VHE[0])/2
                ax10.errorbar(np.log10(optdic[cibndb]),lims[3]-.9,xerr=wbd1,\
                              color='red',ms=0.0,capsize=3.5)
                ax10.annotate(cibndb,(np.log10(optdic[cibndb])-0.25,lims[3]-0.8),fontsize=13)
                ax10.errorbar(np.log10(optdic[cibnda]),lims[3]-0.9,xerr=wbd2,\
                              color='c',ms=0.0,capsize=3.5)
                ax10.annotate(cibnda,(np.log10(optdic[cibnda])-0.25,lims[3]-0.8),fontsize=13)
                ax10.errorbar(np.mean(np.log10(range_xsoft)),lims[3]-0.9,\
                              xerr=wbdx,color='grey',ms=0.0,capsize=3.5)
                ax10.annotate('2-10 keV',(np.log10(range_xsoft[0]),lims[3]-0.8),fontsize=13)
                ax10.errorbar(np.mean(np.log10(range_VHE)),lims[3]-0.9,\
                              xerr=wbdf,color='blue',ms=0.0,capsize=3.5)
                ax10.annotate('MAGIC',(np.log10(range_VHE[0]),lims[3]-0.8),fontsize=13)
    
    
            if namevar=='lext':
                #2-D Histogram of BLR / or BB in general
                xbbt=np.log10(10**xbb[ybb>-30])
                ybbt=np.log10(10**ybb[ybb>-30])
                var1=min(fkvarB)
                var2=max(fkvarB)
                ax10.plot(xbbt,ybbt+np.log10(var1/np.mean(fkvarB)),'m:',linewidth=0.8)
                ax10.plot(xbbt,ybbt+np.log10(var2/np.mean(fkvarB)),'m:',linewidth=0.8)
                ax10.fill_between(xbbt,np.log10(var1/np.mean(fkvarB))+ybbt,ybbt+np.log10(var2/np.mean(fkvarB)), facecolor='purple', alpha=0.2,edgecolor='purple')
            ax10.tick_params(axis='x', labelsize=12)
            ax10.tick_params(axis='y', labelsize=12)
            ax10.tick_params(which='major', width=1.2, length=7)
            ax10.tick_params(which='minor', width=1.1, length=4)
            for axis in ['top', 'bottom', 'left', 'right']:
                ax10.spines[axis].set_linewidth(1.1)
            
            if REALHIST:
                #HISTOGRAM
                ax1[1].hist(np.log10(VHE),bins=50,orientation='horizontal',density=True,color=onecolor,alpha=0.5,label='sim')
                ax1[1].hist(obs_flux,bins=50,orientation='horizontal',density=True,color='tab:orange',alpha=0.5,label='obs')
                ax1[1].hist(obs_flux_S,bins=50,orientation='horizontal',density=True,color='tab:orange',alpha=0.5,label='obs Stat')
                #ax1[1].set_title('Simulation/Observations',size=16)
                ax1[1].set_xlabel(r'High Energy  PDF')
                ax1[1].legend(framealpha=0.5)
                ax1[1].axis([0,max(max(np.histogram(np.log10(VHE),density=True)[0]),\
                   max(np.histogram(obs_flux,density=True)[0]))*1.1,lims[2],lims[3]])
    if SAVE:
        titlenums='_'+str(n1)+'-'+str(nmax)
        if tbin!=-1:
            titlenums=titlenums+'_tbin'+str(tbin)
        if (nmax-n1)<limSEDs and SED!='no':
            plt.figure(1).savefig(routesave+SED.strip('on')+'vFv_'+namevar+imtyp,bbox_inches='tight')
        elif (nmax-n1)>limSEDs and SED!='no':
            print('ERROR: Too many SEDs(over '+str(limSEDs)+') to introduce into one diagram')
            plt.figure(1).savefig(routesave+SED.strip('on')+'vFv_'+namevar+imtyp,bbox_inches='tight')
    plt.show()
    #%% SEGMENTATION
    fkts= []
    for tt in ['','B']:
        #set Time Curves to segment initial and 1-day binned
        tx=locals()['tobs'+tt]
        opt, xs,xm, xh , xb, fer , fer1 , fer2, fer3, fkvr, fkt = locals()['optical'+tt], locals()['xsoft'+tt], locals()['xmed'+tt],\
                                locals()['xhard'+tt], locals()['xband'+tt], locals()['VHE'+tt], \
                                locals()['VHEband1'+tt], locals()['VHEband2'+tt], locals()['VHEband3'+tt],\
                                locals()['fkvar'+tt] , locals()['fktime'+tt] 
        
        if DOUBLE_VAR: fkts = locals()['fkvar_sec'+tt]        
        if CI:
            cii, fiA, fiB = locals()['ci'+tt], locals()['filterA'+tt], locals()['filterB'+tt],  
        if namevar in ['B','lext']:
            ggb = locals()['gbr'+tt]
            if ELECTRONS: ggbn = locals()['gbrnum'+tt]
        seg, topr= 0 , 0
        lenseg=int(np.ceil(len(tx)/SEGnum))
        for to in range(1,len(tx)+1):
            if to%lenseg==0 or to==len(tx):
                seg+=1
                #initially 1 tcross resolved TCs
                locals()['TCs'+tt+str(seg)]=[tx[topr:to],opt[topr:to],xs[topr:to],xm[topr:to],xh[topr:to],\
                       xb[topr:to],fer[topr:to], fer1[topr:to], fer2[topr:to], fer3[topr:to],fkvr[topr:to], fkt[topr:to]]+DOUBLE_VAR*[fkts[topr:to]]
                if CI:
                  locals()['TCs'+tt+str(seg)].append(cii[topr:to])
                  locals()['TCs'+tt+str(seg)].append(fiA[topr:to])
                  locals()['TCs'+tt+str(seg)].append(fiB[topr:to])
        
                if namevar in ['B','lext']:
                  locals()['TCs'+tt+str(seg)].append(ggb[topr:to])
                  if ELECTRONS: locals()['TCs'+tt+str(seg)].append(ggbn[topr:to])
                topr=to
    
    
    for seg in range(1,SEGnum+1):
        tobs, optical ,xsoft , xmed, xhard ,xband ,VHE , VHEband1 ,VHEband2 ,VHEband3 ,fkvar ,fktime = locals()['TCs'+str(seg)][0:12]
        tobsB,opticalB,xsoftB,xmedB, xhardB,xbandB,VHEbandB, VHEband1B ,VHEband2B ,VHEband3B ,fkvarB , fktimeB = locals()['TCsB'+str(seg)][0:12]
        iii = 0
        if DOUBLE_VAR:
            fkvar_sec = TCs1[12]
            fkvar_secB = TCsB1[12]
            iii = 1
        if CI:
            ci ,filterA ,filterB  =locals()['TCs'+str(seg)][12+iii:15+iii]
            ciB,filterAB,filterBB =locals()['TCsB'+str(seg)][12+iii:15+iii]
        if namevar in ['B','lext']:
            gbr = locals()['TCs'+str(seg)][15+iii]
            gbrB = locals()['TCsB'+str(seg)][15+iii]
            if ELECTRONS: gbrnum,  gbrnumB  = locals()['TCs'+str(seg)][16+iii] , locals()['TCsB'+str(seg)][16+iii]
        if SEGnum>1:
              titlenums +='_seg'+str(seg)
              print('SEGMENT\t\t'+str(seg))
    #%%"""PLOTTING TIMECURVES"""
        nmv=namevar
        fktm, fkvr , nunit = tobsB, fkvar , units[nmv]
        fkvarrenorm  = (np.mean(VHEband*0.01)/np.mean(fkvr))
        if DOUBLE_VAR:
            nmv_sec='log'+namevar_sec
            fktm_sec, fkvr_sec , nunit_sec = tobsB, fkvar_sec, units[nmv_sec] 
            fkvarrenorm_sec  = (np.mean(VHEband*0.02)/np.mean(fkvr_sec))
        # ftmin, fvarmin =(fktm[list(fkvr).index(min(fkvr))]+addtcross/np.log(2))*tcr_all[list(fkvr).index(max(fkvr))]+t0,\
        #               np.log10(min(fkvr)*np.mean(VHEband*2)/np.mean(fkvr))
        # ftmax, fvarmax =(fktm[list(fkvr).index(max(fkvr))]+addtcross/np.log(2))*tcr_all[list(fkvr).index(max(fkvr))]+t0,\
        #               np.log10(max(fkvr)*np.mean(VHEband*2)/np.mean(fkvr))
                      
        ftmin, fvarmin, ftmax, fvarmax  = fktm[list(fkvr).index(min(fkvr))] , fkvar[list(fkvr).index(min(fkvr))],  fktm[list(fkvr).index(max(fkvr))] , fkvar[list(fkvr).index(max(fkvr))],
        if TIMECURVES:  
            #TIMECURVES
            tplot, tplotB = tobs, tobsB
            figlc , ax= plt.subplots(num=77,figsize=(20,8),sharex=True)
            gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
            axlc = plt.subplot(gs[0])
            axresid  = plt.subplot(gs[1])
            #plt.plot(tplotB,np.log10(opticalB),'r-', label='O/IR')
            for rr in range(ir,-1,-1):
                axlc.plot( locals()['rt_'+str(rr)],   locals()['rf_'+str(rr)],'s',label= ir_names[rr],\
                          color=ir_colors[rr], ms=4.0, markerfacecolor='white')
            axlc.plot(tplotB,np.log10(VHEband1B),'m.',label=counts_labels[0])
            #plt.plot(tplotB,np.log10(VHEband2B),'r.',label=counts_labels[1])
            #plt.plot(tplotB,np.log10(VHEband3B),'m.',label=counts_labels[2])
            axlc.plot([min(days_with_extrastate),max(days_with_extrastate)], [-8.5,-8.5],color='tab:purple')
            axlc.text(np.sum(days_with_extrastate)/len(days_with_extrastate)-1.0, -8.35, 'super-fast component NOT assumed here',color='tab:purple')
            axlc.plot(tplotB,np.log10(VHEband_countsB),'c.',label=counts_labels[iVHEband])
            #plt.plot(tplotB,np.log10(VHE_normfitB*0.25*eV*1e12),'r.',label='normfit')
            axlc.plot(tplotB,np.log10(xbandB),'y.',label='{}-{} keV sim'.format(round(range_xband[0]*h/eV/1e3,1),round(range_xband[1]*h/eV/1e3,1)))
            axlc.plot( fktm, 
                      #(fktm+addtcross/np.log(2))*tcr_arr,
                          np.log10(fkvr*fkvarrenorm),'x', label = nmv+'(t), renorm.', lw=1.5, alpha=0.5)
            if DOUBLE_VAR:axlc.plot( fktm, 
                      #(fktm+addtcross/np.log(2))*tcr_arr,
                          np.log10(fkvr_sec*fkvarrenorm_sec), 'x', label = nmv_sec+'(t), renorm.', lw=1.5, alpha=0.5)
            # plt.plot(fktm, 
            #           #(fktm+addtcross/np.log(2))*tcr_arr,
            #                   np.log10(fkvr*fkvarrenorm),'-.', label = '$B(t)/<B> x <F_{\gamma}>$', lw=1.5, alpha=0.5)
            axlc.plot(ftmin,np.log10(fvarmin*fkvarrenorm),'b*', label=r'min $ '+nmv+' = '+\
              str(round(min(fkvr),int(-np.log10(min(fkvr)))+1))+nunit+'$')
            axlc.plot(ftmax,np.log10(fvarmax*fkvarrenorm),'r*',ms=7.5, label = r'max ${'+nmv+'} = '+str(round(max(fkvr),int(-np.log10(max(fkvr)))+1))+nunit+' $',)
            axlc.set_title('Mrk421 April 2013 Full Flare',fontsize=16)
            axlc.set_ylabel(r'$F_{obs}$[erg/cm$^2$/sec]',fontsize=14) 
            axlc.legend(fontsize=legsize,markerscale=2,framealpha=0.5,loc='upper left')
            limya20=np.floor(min(min(np.log10(VHEbandB[5::])),\
                                  min(np.log10(xband[5::])),np.log10(max(fkvr[5::])*\
                                      np.mean(VHEband_countsB[5::]*2)/np.mean(fkvr[5::]))))
            limyb20=np.ceil(max(max(np.log10(VHEband_countsB[5::])),\
                                max(np.log10(xband[5::])),np.log10(max(fkvr[5::])*\
                                    np.mean(VHEband_countsB[5::]*2)/np.mean(fkvr[5::]))))
            #plt.axis([initt,Time*tcr,limya20,limyb20])
            axlc.set_ylim(-12.0,-8.0)
            #axlc.xaxis.set_minor_locator(plt.LinearLocator(5*int(Time+1)))
            axlc.xaxis.set_minor_locator(AutoMinorLocator())
            axlc.yaxis.set_minor_locator(AutoMinorLocator())
            plt.subplot(2,1,2)
            miti ,mati = max(min(tplotB),min(rt_0),min(rt_1)) , min(max(tplotB),max(rt_0),max(rt_1)) 
            tresid = tplotB[tplotB>miti][tplotB[tplotB>miti]<mati]
            modxmed_res = inpl.interp1d(tresid, np.log10(xbandB[tplotB>miti][tplotB[tplotB>miti]<mati]))
            modmag_res = inpl.interp1d(tresid, np.log10(VHEband_countsB[tplotB>miti][tplotB[tplotB>miti]<mati]))
            miti ,mati = max(min(tresid),min(rt_0),min(rt_1)) , min(max(tresid),max(rt_0),max(rt_1)) 
            tresid_obs0 = rt_0[rt_0>miti][rt_0[rt_0>miti]<mati]
            tresid_obs1 = rt_1[rt_1>miti][rt_1[rt_1>miti]<mati]
            modrf0_res = inpl.interp1d(tresid_obs0, rf_0[rt_0>miti][rt_0[rt_0>miti]<mati])
            modrf1_res = inpl.interp1d(tresid_obs1, rf_1[rt_1>miti][rt_1[rt_1>miti]<mati])
            rel_diff_0  =  (10**modxmed_res(tresid_obs0) - 10**modrf0_res(tresid_obs0))/10**modrf0_res(tresid_obs0)
            rel_diff_0_rev = (10**modrf0_res(tresid_obs0) -10**modxmed_res(tresid_obs0))/10**modxmed_res(tresid_obs1)
            rel_diff_1 = (10**modmag_res(tresid_obs1) - 10**modrf1_res(tresid_obs1))/10**modrf1_res(tresid_obs1)
            rel_diff_2 = (10**modrf1_res(tresid_obs1) -10**modmag_res(tresid_obs1))/10**modmag_res(tresid_obs1)
            if not feedback_VHE_on_extrastate: 
                np.savetxt(routelc+'residuals_'+namevar+'_NUSTAR.txt', np.c_[tresid_obs0, rel_diff_0_rev])
                np.savetxt(routelc+'residuals_'+namevar+'_MAGIC.txt', np.c_[tresid_obs1, rel_diff_2])
                np.savetxt(routelc+'current_residuals.txt', np.c_[tresid_obs1, rel_diff_2]) #save to use last 
            axresid.plot(tresid_obs0,rel_diff_0*100,'ro' ,lw=2,label='res X-rays %')
            axresid.plot(tresid_obs1, rel_diff_1*100,'bo',lw=2,label='res MAGIC %')
            axresid.fill_between(tresid_obs1, -100, 100,color='green', alpha=0.3,label='100%')
            axresid.fill_between(tresid_obs1, -10, 10,color='green', alpha=0.3,label='10%')
            axresid.legend(fontsize='xx-small')
            axresid.set_xlabel(r'$t_{obs}$[ $days$ ]',fontsize=14)
            axresid.set_ylabel(r'(model- obs) / obs [%]')
            axresid.set_xlim(list(axlc.axis()[0:2]))
            axresid.set_yscale('symlog')
            if SAVE:  
                for imtypp in ['.pdf',imtyp]:
                    figlc.savefig(routesave+'LCs_'+namevar\
                                       +titlenums+imtypp,bbox_inches='tight')
                #plt.figure(24).savefig(routesave+'LCs_'+namevar+'_'+str(round(timespan))+'d_'+str(tbin)+imtyp,bbox_inches='tight')
                # if CI:
                #     plt.figure(19).savefig(routesave+'CI_TC_'\
                #               +namevar+titlenums+imtyp,bbox_inches='tight')
            
            """Gamma Break""" #numerical , theoretical (varying process) + color
            if namevar in ['B','lext']:
                fig90,ax90=plt.subplots(num=90,figsize=(15,5))
                if ELECTRONS: ax90.plot(tplot,gbrnum,label=r'numerical $\gamma_{br}$')
                ax90.plot(tplot,gbr,label=r'analytical $\gamma_{br}$')
                ax90.plot(tplot,np.ones(len(tplot))*loggmin,'k--',linewidth=0.5,\
                          label=r'$\gamma_{min}='+str(round(loggmin,2))+'$')
                ax90.plot(tplot,np.ones(len(tplot))*loggmax,'k-.',linewidth=0.5,\
                          label=r'$\gamma_{max}='+str(round(loggmax,2))+'$')
                ax90.plot(tplot,alphax+np.mean(gbr-alphax),label=r'$\alpha_x$')
                ax90.axis([min(tplot),max(tplot),loggmin-0.2,loggmax+0.2])
                plt.xlabel(r'$t_{ obs}$ [ $days$ ]',fontsize=14)
                plt.ylabel(r'$log\,\gamma_{br}\; / C.I. $',fontsize=14)
                plt.title(r'Cooling Break $\gamma_{break}$    '+obj,fontsize=16)
                plt.legend(fontsize=legsize,markerscale=2,framealpha=0.5)
    
                if SAVE:
                    for imtypp in ['.pdf',imtyp]:
                        plt.figure(90).savefig(routesave+'gbr_TC_'\
                              +namevar+titlenums+imtypp,bbox_inches='tight')   
        # #%% HISTOGRAMS / PDFs
        # if TIMECURVES:
        #     plt.figure(6)
        #     plt.hist(np.log10(VHEbandB),bins=50,density=True,color='b',alpha=0.5,label='gamma-rays')
        #     plt.hist(np.log10(xsoftB),bins=50,density=True,color='tab:orange',alpha=0.5,label='2-10 keV')
        #     plt.hist(np.log10(xhardB),bins=50,density=True,color='grey',alpha=0.5,label='10-80 keV')
        #     plt.title('Probability Density Functions',size=16)
        #     plt.ylabel(r'PDFs',size=14)
        #     plt.xlabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]',size=14)
        #     plt.legend()
        #     if SAVE:
        #         plt.figure(6).savefig('XraysPDFs'+imtyp,bbox_inches='tight')
        #%% CORRELATIONS
        if CORR:
            unt = 'erg'
            LW=2.0
            mscorr=5.0
            ADD_SPOTS=         True
            color_continuous = False
            
            #plot correlation of parameters
            if DOUBLE_VAR:                
                plt.figure(figsize=[6,6])
                tpast='56392'
                for i,t in enumerate(tobs):
                    tt = str(int(max(tobs[i], 56393)))
                    color = colors_daysplit[tt]
                    if tt==tpast:
                        plt.plot(fkvar[i], fkvar_sec[i],'s',color=color,ms=4)
                    else:
                        plt.plot(fkvar[i], fkvar_sec[i],'s',ms=4,color=color,label=tt)
                    tpast =tt
                plt.legend(fontsize='small')
                plt.xlabel(namevar,size=16)
                plt.ylabel(namevar_sec,size=16)
                plt.savefig('corr'+namevar+'-'+namevar_sec+imtyp)
                
            #MOSAIC
            ii =  binsB-1
            fig23, ax23 =plt.subplots(3,3, num=23, figsize=[12,12])
            gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 
            for i , zzz in enumerate(zip([xsoftB[0:ii],xmedB[0:ii],xhardB[0:ii]],[r'$\log_{10}\,F_{3-7 keV}$ [erg/cm$^2$/s]', r'$\log_{10}\,F_{7-30 keV}$ [erg/cm$^2$/s]', r'$\log_{10}\,F_{30-80 keV}$ [erg/cm$^2$/s]'])):
                for j, yyy in enumerate(zip([VHEband1B[0:ii], VHEband2B[0:ii], VHEband3B[0:ii]],[r'$\log_{10}\,F_{>250 GeV}$ ['+unt+'/cm$^2$/s]', r'$\log_{10}\,F_{0.4-0.8TeV}$ ['+unt+'/cm$^2$/s]' , r'$\log_{10}\,F_{>0.8TeV}$ ['+unt+'/cm$^2$/s]'])):
                    arrx , lblx  = zzz
                    arry , lbly  = yyy
                    corrplot_withspots(arrx,arry,tobsB[0:ii], ax23[j,i], labs=[lblx,lbly],XLAB= j==2, YLAB= i==0, corrfont=10.0, 
                                       color_continuous=False, color_dic = colors_daysplit)
            # ax23 = plt.subplot(gs[0])
            # axcbar = plt.subplot(gs[1])
            #titlexg = r' X-rays vs $\gamma$-rays'
            #plt.suptitle(titlexg,size=16)
            plt.savefig('X-gamma_mosaic_'+titlenums+imtyp,bbox_inches="tight")
            plt.savefig('X-gamma_mosaic_'+titlenums+'.pdf',bbox_inches="tight")
            
            #plt.close('all')
            listcorrplots=[]
            if TIMELAPSE_CORR:
                listcorrplots = range(1,binsB-1)
            else:
                listcorrplots = [binsB-1-12]
            ars = [xsoftB[0:ii],VHEbandB[0:ii],10**alphagB[0:ii], 10**alphax[0:ii]] # slopes to be linearized in function
            labs = [r'$\log_{10}\,F_{3-7 keV}$ [erg/cm$^2$/s]', r'$\log_{10}\,F_{>250 GeV}$ ['+unt+'/cm$^2$/s]' ,r'$\alpha_{\gamma}$', r'$\alpha_{X}$']
            axs = [ [-10.0 , -8.8  ] , [ -10.0 , -8.8  ],  [ -2.1 ,-1.1 ] ,  [ -1.7 ,-0.9 ]]
            for ii in listcorrplots:
                    for i , zzz,  in enumerate(zip(ars[0:2],labs[0:2], axs[0:2])):
                        for j ,yyy,  in enumerate(zip(ars,labs,axs)):
                            vara = ['x','g','ag','ax'][j]
                            varb = ['x','g'][i]
                            if vara!=varb:
                                arrx , lblx , axx = zzz
                                arry , lbly , axy = yyy
                                axcorr = axx + axy
                                fig,ax = plt.subplots(figsize=[5,5] , num=i+j*10)
                                corrplot_withspots(arrx,arry,tobsB[0:ii], ax, labs=[lblx,lbly],XLAB= True, YLAB=True, corrfont=corrfont, 
                                           color_continuous=False, axcorr = axcorr, color_dic = colors_daysplit)
                                if SAVE:
                                    titlenums='_'+str(n1)+'-'+str(ii)
                                    fig.savefig(routesave+'corr_'+namevar+'_'+vara+'-'+varb+titlenums\
                                         +'.jpg',bbox_inches='tight')       
                                plt.close()
                                
                                                            
                    if TIMELAPSE_CORR:
                         print('step correlation plots:'+str(ii)+'/'+str(listcorrplots[-1]))
                         plt.close('all')     
                    
            #Other correlation plots
            other_labs= ['O/IR '+optfq,
                         r'$\log_{10}\,F_{3-7 keV}$ [erg/cm$^2$/s]',     r'$\log_{10}\,F_{7-30 keV}$ [erg/cm$^2$/s]',            r'$\log_{10}\,F_{30-80 keV}$ [erg/cm$^2$/s]',
                         r'$\log_{10}\,F_{0.2-0.4 TeV}$ erg/cm$^2$/s]',  r'$\log_{10}\,F_{0.4-0.8 TeV}$ erg/cm$^2$/s]',           r'$\log_{10}\,F_{>0.8 TeV} erg/cm$^2$/s]', 
                         r'$\log_{10}\,F_{0.2-0.4 TeV}$ [#/cm$^2$/s]',   r'$\log_{10}\,F_{0.4-0.8 TeV}$ [#/cm$^2$/s]',       r'$\log_{10}\,F_{>0.8 TeV} [#/cm$^2$/s]']
                         
            titlenums='_'+str(n1)+'-'+str(nmax)
            corrleg = 'corr_'+namevar+DOUBLE_VAR*('+'+namevar_sec)
            #optical -gamma-rays
            fig,ax = plt.subplots(figsize=[5,5] , num=101)
            vara, varb='o','g'
            corrplot_withspots(opticalB,VHEbandB,tobsB, ax, labs=[other_labs[0],other_labs[4+iVHEband]],XLAB= True, YLAB=True, lvls=lvls,
                               corrfont=corrfont,  color_continuous=True, axcorr =[])
            if SAVE:   fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 
    
            
            #optical - X-rays
            fig,ax = plt.subplots(figsize=[5,5] , num=102)
            vara, varb='o','x'
            corrplot_withspots(opticalB,xbandB,tobsB, ax, labs=[other_labs[0],other_labs[1+ixband]],XLAB= True, YLAB=True, lvls=lvls,
                               corrfont=corrfont,  color_continuous=True, axcorr =[])
                    
            #gamma-rays g-break
            if ELECTRONS:
                fig,ax = plt.subplots(figsize=[5,5] , num=111)
                vara, varb='g','gbr'
                corrplot_withspots(VHEbandB,10**gbrnumB,tobsB, ax, labs=[other_labs[4+VHEband],other_labs[1]],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr =[])
                if SAVE:   fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 
        
                        
                #X-rays g-break
                fig,ax = plt.subplots(figsize=[5,5] , num=112)
                vara, varb='x','gbr'
                corrplot_withspots(xbandB,10**gbrnumB,tobsB, ax, labs=[other_labs[1+ixband],r'$log_10\,\gamma_{br,cool}$'],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr =[])
                if SAVE:fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 
        
                #X-rays alpha g-break
                fig,ax = plt.subplots(figsize=[5,5] , num=113)
                vara, varb='alphax','gbr'
                corrplot_withspots(10**alphaxB,10**gbrnumB,tobsB, ax, labs=[labs[-1],r'$log_10\,\gamma_{br,cool}$'],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr =[])
                if SAVE:  fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 
                        
                #X-rays alpha g-break
                fig,ax = plt.subplots(figsize=[5,5] , num=114)
                vara, varb='alphag','gbr'
                corrplot_withspots(10**alphagB,10**gbrnumB,tobsB, ax, labs=[labs[-2],r'$log_10\,\gamma_{br,cool}$'],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr =[])
                if SAVE:fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 

            #color variations
            if CI:
                optvci=cibndb
                gci=r'\log_10 (_{>0.2 TeV}/ <F_{>0.2 TeV}>)'
                mag_opt=filterBB #select x axis from optical bands and make it color
                mag_g=np.log10(VHEbandB/np.mean(VHEbandB)) #select x axis from optical bands and make it color
                
                fig,ax = plt.subplots(figsize=[5,5] , num=121)
                vara, varb='o','ci'
                corrplot_withspots(10**mag_opt,10**ciB,tobsB, ax, labs=[r'${}$'.format(optvci), r'${}$'.format(cibndb)],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr = [])
                if SAVE:   fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight') 

                fig,ax = plt.subplots(figsize=[5,5] , num=122)
                vara, varb='g','ci'
                corrplot_withspots(10**mag_g,10**ciB,tobsB, ax, labs=[r'${}$'.format(gci), r'${}$'.format(cibndb)],XLAB= True, YLAB=True, lvls=lvls,
                                   corrfont=corrfont,  color_continuous=True, axcorr = [])
                if SAVE:  fig.savefig(routesave+corrleg+'_'+vara+'-'+varb+titlenums +imtyp,bbox_inches='tight')   
        
        #%%    
        """QUICK_TIMELAGS"""
        """Find timelag from peaks or minima// QUICK METHOD"""
        """QUICK_TIMELAGS"""
        """Find timelag from peaks or minima// QUICK METHOD"""
        if QUICK_TIMELAGS:
           smod=inpl.interp1d(tobs,optical,kind='quadratic')
           fmod=inpl.interp1d(tobs,VHEband,kind='quadratic')
           xmod=inpl.interp1d(tobs,xband,kind='quadratic')
           tobopt=np.linspace(min(tobs),max(tobs),2*len(tobs))
           
           crit = tcr*10*(tobopt[1]-tobopt[0])**-1 #how many tcross can two energies differ in emission given we have fast escape?
           
           if not len(tobs)>num_peaks*lendel*2.0 :
               print('Won\'t calculate quick peaks because of lack of data for the demanded number of peaks')
           elif peaktype in ['minima','maxima']:
               opt=smod(tobopt)
               fm=fmod(tobopt)
               xm=fmod(tobopt)
               inds, indf, indx = [], [] , [] 
               for i in range(num_peaks):
                   if peaktype=='maxima':
                       x=list(opt).index(max(opt))
                       y=list(fm).index(max(fm))
                       z=list(xm).index(max(xm))
                   elif peaktype=='minima':
                       x=list(opt).index(min(opt))
                       y=list(fm).index(min(fm))
                       z=list(xm).index(min(xm))
                   else: break
                   inds.append(x)
                   indf.append(y)
                   indx.append(z)
                   opt=np.concatenate([opt[lendel:max(x-lendel,lendel)],opt[min(x+lendel,len(opt)-1-lendel):len(opt)-1-lendel]])
                   fm=np.concatenate([fm[lendel:max(y-lendel,lendel)],fm[min(y+lendel,len(opt)-1-lendel):len(opt)-1-lendel]])
                   xm=np.concatenate([xm[lendel:max(z-lendel,lendel)],xm[min(z+lendel,len(opt)-1-lendel):len(opt)-1-lendel]])
               
               timelags=[[],[],[]] # select below the delayed band in loop
               for i in range(len(inds)) :
                   if abs(inds[i]-indf[i])<crit :timelags[0].append(round(tobopt[indf[i]]-tobopt[inds[i]],3))
                   if abs(indf[i]-indx[i])<crit :timelags[2].append(round(tobopt[indf[i]]-tobopt[indx[i]],3))
                   if abs(inds[i]-indx[i])<crit:timelags[1].append(round(tobopt[indx[i]]-tobopt[inds[i]],3))
                       
               print('\n Calculated timelags from 25 flares [delayed - firstly observed]:\n',timelags[0],'\n',timelags[1],'\n',timelags[2],\
                     '\n Gamma Rays after optical (Av.+/- std):\t '+str(round(np.mean((timelags[0])),3))+'+/-'+str(round(np.std(timelags[0]),3))+' days',\
                     '\n Xrays after optical (Av. +/- std):\t '+str(round(np.mean((timelags[1])),3))+'+/-'+str(round(np.std(timelags[1]),3))+' days',\
                     '\n Xrays after Gammarays (Av. +/- std):\t '+str(round(np.mean((timelags[2])),3))+'+/-'+str(round(np.std(timelags[2]),3))+' days')
                   
                   
        """ DCF """
        if DCFcalc:
                #CALCULATING
                dtau=tobs[3]-tobs[2]
                if taubin<dtau:
                    taubin=dtau
                if dtau<1.0:
                    taubin=dtau
                taurange=np.linspace(-lentau,lentau,2*int(lentau/taubin)+1)
                tlags, DCFs=[] , []
                lnbd=len(bands)
                for bnd in range(lnbd):
                    if bnd==0:
                        band1=bands[-1]
                        band2=bands[0]
                        bands[-1]=bands[0] #change sequence of comparison for visual reasons  less -> more energetic
                        bands[0]=band1
                    a=np.log10(globals()[bands[bnd-1]+'B']) #+B --> use new binning if altered
                    b=np.log10(globals()[bands[bnd]+'B']) #use new binning
                    amean , bmean = np.mean(a) ,np.mean(b)
                    astd , bstd = np.std(a) , np.std(b)
                    tlag , DCF =[] , []
                    step=0
                    for tau in taurange:
                        udcf=0
                        M=0
                        for i in range(len(a)-1):
                         aa=a[i]
                         for j in range(len(b)-1):
                            bb=b[j]
                            dt=tobsB[i]-tobsB[j]
                            if dt<tau+taubin/2 and dt>tau-taubin/2:
                                udcf=udcf+(aa-amean)*(bb-bmean)/astd/bstd 
                                M=M+1
                        if M!=0:
                            tlag.append(tau)
                            DCF.append(udcf/M)
                            print('Step  '+str(step+1)+' / '+str(len(taurange))+'     bands pair compared:  '+str(bnd+1)+' / '+str(len(bands)))
                            step=step+1 #count steps
                    tlags.append(np.array(tlag).astype(float))
                    DCFs.append(np.array(DCF).astype(float))
                    if bnd==0:
                        bands[-1]=band1
                        bands[0]= band2            
                #PLOTTING
                fig25,ax25=plt.subplots(1,1,num=25)
                for bnd in range(lnbd):
                    plt.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+'-'+bands[bnd-1])
                plt.legend()
                plt.xlabel(r'$\tau\;\;days$',size=14)
                plt.ylabel(r'$DCF$($\tau_{ obs}$)',size=14)
                plt.title(obj+'    Timelags',size=16)
                ll=len(tlags[bnd])
                plt.plot(np.linspace(-lentau,lentau,50),0*np.linspace(-1,1,50),'k--',linewidth=0.5)
                plt.plot(0*np.linspace(-1,1,50),np.linspace(-1.,1,50),'k--',linewidth=0.5)
                plt.axis([-lentau,lentau,-1,1])
                
                pospeak=list(abs(DCFs[0])).index(max(abs(DCFs[0])))
                tpk=tlags[0][pospeak]
                pk=DCFs[0][pospeak] #peak of optical VHEband DCF
                if abs(tpk)<2.0*taubin:
                    # location for the zoomed portion  in fig 20
                    ax34= plt.axes([0.65,0.16, 0.25, 0.25])# units unity for axis scales
                    ax34.set_xlabel(r'$\tau_{ obs}$ [days]',size=10)
                    ax34.tick_params(axis='y',labelright='on',labelleft='off',right='on',direction='in',labelsize=10,pad=3)
                    ax34.tick_params(axis='x',labelbottom='on',bottom='on',direction='in',labelsize=10,pad=-12)
                    #ax34.tick_params(axis='x',labelsize=9,pad=-15)
                    #plt.title('Peak')
                    # insert the zoomed figure
                    plt.setp(ax34)
                    for bnd in range(lnbd):
                        plt.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+'-'+bands[bnd-1])
                    if pk>0:
                        ax34.axis([tpk-lentau/4.0,tpk+lentau/4.0,pk*0.85,min(1.0,pk*1.05)])
                    else:
                        ax34.axis([tpk-lentau/4.0,tpk+lentau/4.0,max(pk*1.05,-1.0),pk*0.85])
                    from mpl_toolkits.axes_grid1.inset_locator import mark_inset
                    mark_inset(ax25, ax34, loc1=1, loc2=2, fc="none", ec="0.5")
                plt.show()
                if SAVE:
                    fig25.savefig(routesave+'tlags_'+namevar+titlenums+imtyp,bbox_inches='tight')
                    np.savetxt('DCFs_p'+str(POWER)+'.txt', np.c_[np.concatenate(tlags, axis=0 ),np.concatenate(DCFs, axis=0 )])
        """ SF """
        if SFcalc:
                logtau=np.logspace(np.log10((tobs[1]-tobs[0])/2.),np.log10(tobsB[-1]/2.))#SF
                fsf=logtau[2]/logtau[1]
                sflags, SFs= [], []
                fig27=plt.figure(27)
                lnbd=len(bands)
                errbands=[erro,errx,errg]
                minSF=1.0
                maxSF=1.0
                for bnd in range(lnbd):
                    a=-2.5*np.log10(locals()[bands[bnd-1]+'B']) #mags
                    b=a
                    sigma2=(errbands[bnd])**2
                    amean,bmean,astd,bstd= 1. , 1. ,1. ,1.
                    sflag=[]
                    SF=[]
                    step=0
                    for tau in logtau:
                        usf=0
                        us12,um1,um2=0,0,0
                        M=0
                        for i in range(len(a)-1):
                         aa=a[i]
                         for j in range(len(b)-1):
                            bb=b[j]
                            dt=tobsB[i]-tobsB[j]
                            if dt>tau/fsf and dt<tau*fsf:
                                usf=usf+abs(aa-bb) #SF
                                M=M+1
                        if M!=0:
                            sflag.append(tau)
                            if sigma2<(usf/M)**2:
                                SF.append(np.sqrt(np.pi/8*((usf/M)**2-sigma2)))
                            else:
                                #SF.append(np.sqrt(np.pi/8)*usf/M)
                                SF.append(1e-5)                     
                        print('Step  '+str(step+1)+' / '+str(len(logtau))+\
                              '    SF for band:  '+str(bnd+1)+' / '+str(len(bands)))
                        step=step+1 #count steps
                    sflags.append(np.array(sflag).astype(float))
                    SFs.append(np.array(SF).astype(float))
                    plt.loglog(sflags[bnd],SFs[bnd],'-',label=bands[bnd])
                    minSF=min(min(SF),minSF)
                    maxSF=max(max(SF),maxSF)
                plt.legend()
                plt.xlabel(r'$\tau_{ obs}\;\;[days]$',size=14)
                plt.ylabel(r'$SF$($\tau_{ obs}$)',size=14)
                plt.title(obj+'  Structure Function',size=16)
                plt.axis([min(logtau),max(logtau),0.95*minSF,1.05*maxSF])
                plt.show()
                
                if SAVE:
                    plt.savefig(routesave+'SF_'+namevar+titlenums+imtyp,bbox_inches='tight')
        """Stingray/ Timelags per frequency range"""
        if STINGRAY:
            from stingray import Lightcurve, AveragedCrossspectrum #timelags
            long_times=tobsB
            nn=0
            figstg,axstg=plt.subplots(3,1,num=26,sharex='col')
            axstg[0].set_title(obj+'    Frequency-dependent lags',size=16)
            for bnd in range(len(bands)):
                long_signal_1=locals()[bands[bnd-1]+'B']  #+B--> if new binning 
                long_signal_2=locals()[bands[bnd]+'B']    
                long_lc1 = Lightcurve(long_times, long_signal_1)
                long_lc2 = Lightcurve(long_times, long_signal_2)
                avg_cs = AveragedCrossspectrum(long_lc1, long_lc2,segment_size)
                freq_lags, freq_lags_err = avg_cs.time_lag()
                freq_lags=-freq_lags #change notation to match DCF
                #PLOT
                ax=axstg[nn]
                ax.hlines(0, avg_cs.freq[0], avg_cs.freq[-1],\
                          color='black', linestyle='dashed', lw=2)
                ax.errorbar(avg_cs.freq, freq_lags, yerr=freq_lags_err,fmt=".", lw=1,\
                            color=colors[nn],label=bands[bnd]+' - '+bands[bnd-1])
                ax.legend()
                axax=ax.axis()
                ax.axis([axax[0],avg_cs.freq[-1]/2,axax[2],axax[3]]) 
                #ax.set_ylabel("Time lag (days)",labelpad=10)
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(axis='x', labelsize=12)
                ax.tick_params(axis='y', labelsize=12)
                ax.tick_params(which='major', width=1.5, length=7)
                ax.tick_params(which='minor', width=1.5, length=4)
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(1.5)
                nn=nn+1
            axstg[1].set_ylabel("Observed Time lag (days)",size=13)
            ax.set_xlabel("Frequency (days^-1)",labelpad=3,size=13)
            plt.show()
            if SAVE:
                figstg.savefig(routesave+'tlags_freq_'+namevar+titlenums+imtyp,bbox_inches='tight')    
        #%%
        """TIMELAGS with Extrapolation (sort of) using neighbour points of DCF"""
        if 'DCFs_p'+str(POWER)+'.txt' in nmlist:
            alldcf=open(route+'DCFs_p'+str(POWER)+'.txt','r')
            dfile = alldcf.readlines()
            alldcf.close()
            lendf = int(len(dfile)/3)
            tlags, DCFs=[ ] , [ ]
            for bnd in range(len(bands)):
                df=dfile[bnd*lendf:lendf*(bnd+1)]
                tlg, dcf = [] , []
                for i in df:
                    tlg.append(i.split()[0])
                    dcf.append(i.split()[1].replace('\n',''))
                tlags.append(np.array(tlg).astype(float))
                DCFs.append(np.array(dcf).astype(float))
    
        if DCFcalc:
            #PLOTTING
            fig25,ax25=plt.subplots(1,1,num=25)
            for bnd in range(lnbd):
                plt.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+'-'+bands[bnd-1])
            plt.legend()
            plt.xlabel(r'$\tau\;\;days$',size=14)
            plt.ylabel(r'$DCF$($\tau_{ obs}$)',size=14)
            plt.title(obj+'    Timelags',size=16)
            ll=len(tlags[bnd])
            plt.plot(np.linspace(-lentau,lentau,50),0*np.linspace(-1,1,50),'k--',linewidth=0.5)
            plt.plot(0*np.linspace(-1,1,50),np.linspace(-1.,1,50),'k--',linewidth=0.5)
            plt.axis([-lentau,lentau,-1,1])
            for bnd in range(len(bands)):
                tlag=tlags[bnd]
                DCF=list(DCFs[bnd])
                if np.mean(DCF)>0:
                    indexm=DCF.index(max(DCF))
                else:
                    indexm=DCF.index(min(DCF))
                tlag_max=tlag[indexm]
                ex1=1/abs(DCF[indexm-1]-DCF[indexm])
                ex2=1/abs(DCF[indexm+1]-DCF[indexm])
                tlag_extra=((ex1*tlag[indexm-1]+ex2*tlag[indexm+1])/(ex1+ex2)+tlag[indexm])/2
                print('\nDCF created for timelag of '+bands[bnd]+' over '+bands[bnd-1])
                print('Most strongly detected timelag (peak) at:\t'+str(tlag_max)+' days')
                print('Maximum with extrapolation method:\t'+str(tlag_extra)+' days')
            #for segment analysis
            locals()['tlag'+str(seg)] = tlag
            locals()['DCFs'+str(seg)] = DCFs
                
        """Fractional Variability"""
        errx = 0.01
        FVin=np.sqrt(np.var(fkvar)/np.mean(fkvar)**2)
        FVo=np.sqrt(np.var(opticalB)/np.mean(opticalB)**2-(10**erro-1)**2)
        FVx1=np.sqrt(np.var(xsoftB)/np.mean(xsoftB)**2-(10**errx-1)**2)
        FVx2=np.sqrt(np.var(xmedB)/np.mean(xmedB)**2-(10**errx-1)**2)
        FVx3=np.sqrt(np.var(xhardB)/np.mean(xhardB)**2-(10**errx-1)**2)
        FVg1=np.sqrt(np.var(VHEband1B)/np.mean(VHEband1B)**2-(10**errg-1)**2)
        FVg2=np.sqrt(np.var(VHEband2B)/np.mean(VHEband2B)**2-(10**errg-1)**2)
        FVg3=np.sqrt(np.var(VHEband3B)/np.mean(VHEband3B)**2-(10**errg-1)**2)
        FVx , FVg = (FVx1+FVx2+FVx3)/3, (FVg1+FVg2+FVg3)/3
    
        
        Nn=len(VHEbandB)
        dFVo=np.sqrt(FVo**2+(((2/Nn)**0.5*(10**erro-1))+(((10**erro-1)/Nn)**0.5*2*FVo)**2)**0.5)-FVo
        dFVx1=np.sqrt(FVx1**2+(((2/Nn)**0.5*(10**errx-1))+(((10**errx-1)/Nn)**0.5*2*FVx1)**2)**0.5)-FVx1
        dFVx2=np.sqrt(FVx2**2+(((2/Nn)**0.5*(10**errx-1))+(((10**errx-1)/Nn)**0.5*2*FVx2)**2)**0.5)-FVx2
        dFVx3=np.sqrt(FVx3**2+(((2/Nn)**0.5*(10**errx-1))+(((10**errx-1)/Nn)**0.5*2*FVx3)**2)**0.5)-FVx3
    
        dFVg1=np.sqrt(FVg1**2+(((2/Nn)**0.5*(10**errg-1))+(((10**errg-1)/Nn)**0.5*2*FVg1)**2)**0.5)-FVg1
        dFVg2=np.sqrt(FVg2**2+(((2/Nn)**0.5*(10**errg-1))+(((10**errg-1)/Nn)**0.5*2*FVg2)**2)**0.5)-FVg2
        dFVg3=np.sqrt(FVg3**2+(((2/Nn)**0.5*(10**errg-1))+(((10**errg-1)/Nn)**0.5*2*FVg3)**2)**0.5)-FVg3
        dFVx , dFVg = (dFVx1+dFVx2+dFVx3)/3, (dFVg1+dFVg2+dFVg3)/3
    
        if np.isnan(FVg1):
            FVg1 = np.sqrt(np.var(VHEband1B))/np.mean(VHEband1B)
            dFVg1 = FVg1
        if np.isnan(FVg2):
            FVg2 = np.sqrt(np.var(VHEband2B))/np.mean(VHEband2B)
            dFVg2 = FVg2
        if np.isnan(FVg3):
            FVg3 = np.sqrt(np.var(VHEband3B))/np.mean(VHEband3B)
            dFVg3 = FVg3
            
        if dFVo>FVo: dFVo = FVo
        if dFVx1>FVx1: dFVx1 = FVx1
        if dFVx2>FVx2: dFVx2 = FVx2
        if dFVx3>FVx3: dFVx3 = FVx3
        if dFVg1>FVg1: dFVg1 = FVg1
        if dFVg2>FVg2: dFVg2 = FVg2
        if dFVg3>FVg3: dFVg3 = FVg3
        
        plt.figure()
        plt.errorbar(14.5, FVo, yerr = dFVo , xerr=0.2)
        plt.errorbar(np.mean(np.log10(range_xsoft)), FVx1, yerr = dFVx1, xerr=0.36)
        plt.errorbar(np.mean(np.log10(range_xmed)), FVx2, yerr = dFVx2, xerr=0.63)
        plt.errorbar(np.mean(np.log10(range_xhard)), FVx3, yerr = dFVx3, xerr=0.42)
        plt.errorbar(np.mean(np.log10(range_VHE1)), FVg1, yerr = dFVg1, xerr=0.3)
        plt.errorbar(np.mean(np.log10(range_VHE2)), FVg2, yerr = dFVg2, xerr=0.3)
        plt.errorbar(np.mean(np.log10(range_VHE3)), FVg3, yerr = dFVg3, xerr=0.2)
        plt.xlabel('log Frequency',size=16)
        plt.ylabel('Fractional Variability')
        plt.savefig('FV.png',bbox_inches="tight")
        print( 'Fractional Variabilities (daily binned LCs) \n optical band:\t'+str(FVo)+'\t+/-'+str(dFVo)\
              +'\n X-rays band 1:\t'+str(FVx1)+'\t+/-'+str(dFVx1)\
              +'\n VHE band 1:\t'+str(FVg1)+'\t+/-'+str(dFVg1))
        
        #SAVE Segment statistics and properties
        locals()['FVin'+str(seg)] =FVin
        locals()['FVo'+str(seg)], locals()['dFVo'+str(seg)]=FVo, dFVo
        locals()['FVx'+str(seg)], locals()['dFVx'+str(seg)]=FVx, dFVx1
        locals()['FVg'+str(seg)], locals()['dFVg'+str(seg)]=FVg, dFVg1
        if SEGnum>1:
            TIMECURVES=False
            plt.pause(1)
            plt.close('all')
    #%%COMPARE SEGMENTS
    if SEGnum>1:
        fig71, ax71 = plt.subplots(num=71)
        for seg in range(1,SEGnum):
            FVin=locals()['FVin'+str(seg)]
            ax71.plot(FVin, locals()['FVo'+str(seg)],'r.')
            ax71.plot(FVin, locals()['FVx'+str(seg)],'gx')
            ax71.plot(FVin, locals()['FVg'+str(seg)],'b+')    
        ax71.set_title('Segment Size: '+str(lenseg)+' observational days',size=corrfont+2)
        ax71.set_xlabel(r'$FV_{in}$',size=corrfont)
        ax71.set_ylabel(r'$FV_{sim}$',size=corrfont)
        ax71.legend(['O/IR',r'$X-$rays', r'$\gamma-$rays'])
        fig71.savefig('FVg_FVin_'+obj+imtyp,bbox_inches='tight')
        
        
        #BOOTSTRAP DCFs
        if DCFcalc:
            sDCF , xDCF, fDCF = np.zeros(len(taurange)) ,np.zeros(len(taurange)) ,np.zeros(len(taurange))
            for seg in range(1,SEGnum):
                for tl in range(0,len(taurange)):
                    ttl=taurange[tl]
                    ind=list(locals()['tlag'+str(seg)]).index(ttl)
                    sDCF[tl]+=locals()['DCFs'+str(seg)][0][ind]
                    xDCF[tl]+=locals()['DCFs'+str(seg)][1][ind]
                    fDCF[tl]+=locals()['DCFs'+str(seg)][2][ind]
            
            DCFa, DCFb, DCFc = sDCF/SEGnum, xDCF/SEGnum, fDCF/SEGnum
            fig72, ax72 = plt.subplots(num=72)
            DCFs= [DCFa, DCFb, DCFc] #the bootstrapped DCFs
            for bnd in range(lnbd):
                            plt.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+'-'+bands[bnd-1])
            plt.legend()
            plt.xlabel(r'$\tau_{ obs}\;\;days$',size=14)
            plt.ylabel(r'$DCF$($\tau_{ obs}$)',size=14)
            plt.title(obj+'    Timelags',size=16)
            ll=len(tlags[bnd])
            plt.plot(np.linspace(-lentau,lentau,50),0*np.linspace(-1,1,50),'k--',linewidth=0.5)
            plt.plot(0*np.linspace(-1,1,50),np.linspace(-1.,1,50),'k--',linewidth=0.5)
            plt.axis([-lentau,lentau,-1,1])  
            if SAVE:fig72.savefig('DCFsboot_'+str(SEGnum)+imtyp,bbox_inches='tight')
        
    #%%
    CALC_THEO_INDICES=False
    """synchro"""
    if namevar=='B':
        fakeB=inpl.interp1d(fktime*tcr,fkvar)
        try:
            B=fakeB(tobs)
        except ValueError:
            B=fkvar
        lb=B**2/8/np.pi/me/c**2*sigmaT*R
        tcsyn=3/4/lb*tcr #days obs frame
        tcmax=tcsyn/10**loggmin
        tcmin=tcsyn/10**loggmax #days obs frame
        vssa=(((p+1)/32/np.pi**2)**2*(3*fortio/2/np.pi/me**3/c)*delta**(-1)*(2.8*10**6*B*10**(2*loggmax))**(p-3)*\
        (4*np.pi*max(vFv_pts[v_pts<18])*D**2)**2*B/R**4)**(1/(p+4))/(1+zsh)
        gssa=0.5*np.log10(vssa/2.8/10**6/B/delta*(1+zsh))
        
        gopt=0.5*np.log10(opt1/2.8/10**6/delta/B)
        gx=0.5*np.log10(range_xsoft[0]/2.8/10**6/delta/B*(1+zsh))
        #
        tcopt=tcsyn/10**gopt
        tcx=tcsyn/10**gx
    #%
    #%
    #%
    #%
    #"""
    if CALC_THEO_INDICES:
    ##actual unrenormalized histogram
        fakevar=inpl.interp1d(fktime*tcr,fkvar)
        try:
            vec=fakevar(tobs)
        except ValueError:
            vec=fkvar
        Nbin=50
        vect=[np.log10(vec),np.log10(optical),np.log10(xband),np.log10(VHEband)]
        vecnames=[namevar,'optical','xband','VHEband']
        
        Nb2=25
        bvec,hvec=np.histogram(vect[0],bins=Nb2)
        bopt,hopt=np.histogram(vect[1],bins=Nb2)
        bx,hx=np.histogram(vect[2],bins=Nb2)
        bf,hf=np.histogram(vect[3],bins=Nb2)
            
        
        import scipy.signal as sg
        pvec=hvec[list(bvec).index(max(bvec))+1]
        popt=hopt[list(bopt).index(max(bopt))+1]
        px=hx[list(bx).index(max(bx))+1]
        pf=hf[list(bf).index(max(bf))+1]
        
        wvec=sg.peak_widths(bvec,[list(bvec).index(max(bvec))])[0][0]*(hvec[1]-hvec[0])
        wopt=sg.peak_widths(bopt,[list(bopt).index(max(bopt))])[0][0]*(hopt[1]-hopt[0])
        wx=sg.peak_widths(bx,[list(bx).index(max(bx))])[0][0]*(hx[1]-hx[0])
        wf=sg.peak_widths(bf,[list(bf).index(max(bf))])[0][0]*(hf[1]-hf[0])
        
        #CALC_THEO_INDICES of variability F ~ (B)**s_i
        sopt=round(np.log10(10**wopt)/np.log10(10**wvec),2)
        sx=round(np.log10(10**wx)/np.log10(10**wvec),2)
        sf=round(np.log10(10**wf)/np.log10(10**wvec),2)
        print('CALC_THEO_INDICES for variability:\n optical:',sopt,'\t X:',sx,' \t VHEband:',sf)
        
        plt.figure(56)
        plt.plot((hvec[1::]-pvec),bvec,label=vecnames[0])
        plt.plot((hopt[1::]-popt),bopt,label=vecnames[1])
        plt.plot((hx[1::]-px),bx,label=vecnames[2])
        plt.plot((hf[1::]-pf),bf,label=vecnames[3])
        plt.xlabel('mag',fontsize=14)
        plt.ylabel('N',fontsize=14)
        plt.title(obj+'   PDFs',size=16)
        plt.legend()
        
    #    #Histogram
    #    plt.figure(57)
    #    plt.hist(vect[0]-pB,alpha=0.8,label=vecnames[0])
    #    plt.hist(vect[1]-popt,alpha=0.65,label=vecnames[1])
    #    plt.hist(vect[2]-px,alpha=0.5,label=vecnames[2])
    #    plt.hist(vect[3]-pf,alpha=0.4,label=vecnames[3])
    #    plt.xlabel('mag',fontsize=14)
    #    plt.ylabel('N',fontsize=14)
    #    plt.title(obj+'   PDFs',size=16)
    #    plt.legend()
        
        #Momentary index
        plt.figure(99)
        plt.xlabel(r'log '+namevar)
        plt.ylabel(r'$log(\,F_{bamd}/<F_{band}>)\,/\,log(\,'+namevar+'/<'+namevar+'>)$')
        plt.title(obj+'Momentary Index of Correlation',size=16)
        indxx=0.001 #fraction of mean value for down limit
        nos1=10**np.mean(np.log10((optical)))
        nog1=10**np.mean(np.log10((VHEband)))
        nox1=10**np.mean(np.log10((xband)))
        noB=10**np.mean(np.log10(B))
        B=np.array(B).astype(float)
        plt.plot(np.log10(B[B>indxx*np.mean(B)]),np.log10(optical[B>indxx*np.mean(B)]/nos1)/np.log10(B[B>indxx*np.mean(B)]/noB),'r.',label='optical')
        plt.plot(np.log10(B[B>indxx*np.mean(B)]),np.log10(VHEband[B>indxx*np.mean(B)]/nog1)/np.log10(B[B>indxx*np.mean(B)]/noB),'b.',label=r'$\gamma$-rays')
        plt.axis([min(vect[0]),max(vect[0]),-0.5,2.0])
        plt.legend()
        
        if namevar=='B':
            vi=0.3*(range_optical[0]+range_optical[1])
            vssa=10**12.5
            vssa=10**12.5
            Bl1=48*(R/10**16)**-0.5*(10**loggmax)**-0.5
            Bl2=48**(4/3)*(R/10**16)**(-2/3)*(vssa/3/10**6/delta*(1+zsh))**(-1/3)
            Bl3=48*(R/10**16)**-0.5*(10**loggmin)**-0.5
            Bi1=48**(4/3)*(R/10**16)**(-2/3)*(vi/3/10**6/delta*(1+zsh))**(-1/3)
            Bi2=vssa/3e6/delta/10**(2*loggmin)*(1+zsh)
            Bi3=vi/3e6/delta/10**(2*loggmin)*(1+zsh)
            Bi4=vi/3e6/delta/10**(2*loggmax)*(1+zsh)
            print('\nCharact B limits:\t',
                  round(Bl1,int(-np.log10(Bl1))+2),\
        
                  round(Bl2,int(-np.log10(Bl2))+2),\
                  round(Bl3,int(-np.log10(Bl3))+2),\
                  round(Bi1,int(-np.log10(Bi1))+2),\
                  round(Bi2,int(-np.log10(Bi2))+2),\
                  round(Bi3,int(-np.log10(Bi3))+2))
            
            p0=-(3*p-1)/6*np.heaviside(B-max(Bi3,Bi2),1) #added after loop
            p1=(p-2)/2
            p2=2*(p-2)
            #p3=(p+1)/2 -3/2*np.heaviside(B-Bi3,1) added in loop
            p4=(p+1)/2
            
            pvar=[]
            for b in B:
                if b>Bl3:
                    pvar.append(round(p1,2))
                if min(Bl3,Bl2)<b<Bl3:
                    pvar.append(round(p2,2))
                if Bl1<b<min(Bl2,Bl3):
                    p3=(p+1)/2
                    if b>Bi3:
                        p3=p3-3/2
                    pvar.append(round(p3,2))
                if b<Bl1:
                    pvar.append(round(p4,2))
            pvar=p0+np.array(pvar)
            fig101=plt.figure(101,figsize=[10.5,4])
            plt.plot(np.log10(B),pvar,'r.',label='analytical')
            plt.xlabel(r'$B\;$[G]',fontsize=14)
            plt.ylabel(r'$s_{optical}=\frac{dlog\,F_{optical}}{dlog\,B}$',fontsize=14)
            plt.title(r''+obj+'\t ($p='+str(p)+'$)',fontsize=16)
            
            
            plt.plot(np.linspace(min(np.log10(B))-0.5,max(np.log10(B))+0.5),np.ones(50)*sopt,'k-',label='numerical')
            plt.axis([min(np.log10(B)-0.1),max(np.log10(B)+0.1),min(pvar)-0.1,max(pvar)+0.1])
            fig101.legend()
            if SAVE:
                fig101.savefig('sopt_'+obj+imtyp,bbox_inches='tight')
