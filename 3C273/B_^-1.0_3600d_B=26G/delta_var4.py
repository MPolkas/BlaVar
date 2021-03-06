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
import scipy.optimize as optmz
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from astropy.io import ascii
from stingray import Lightcurve, AveragedCrossspectrum #timelags
#import scipy.stats as stats
"""%%%%%%%%%%%%%%%%  SETTINGS  %%%%%%%%%%%%%%%%%%% """
#Please, walk through settings down to END SETTINGS for selecting the parameters for each function of the script
#The user is requeste to have multiple one fort.85 or multiple fort_*index*.85 
# under the directory with name /'NameObject'/'Namevar'_... for the current script to operate
#Function Switches , for parameters of functions see below .
ALLTOOLS='no' #except SFcalc
SED='con' #'no', 'on'(overplot), 'con' (contour)  vFv Diagram [CAUTION!]only up to *limSEDs*, set below
TIMELAPSE='no' # Run a timelapse and print hundreds of frames to create a GIF/video later on or observe the evolution live.
CORR='no' #SWITCH for flux flux and Color diagrams
#TIMELAGS='no'  #Quick method to calculate timelag from flares (#PEAKS)
DCFcalc='no' #Calculate or not DCF // For 10-yr long sims DCFcalc could take >1hour
SFcalc='no' #Structure Function
STINGRAY='no' # Time-dependent timelags quick method
SAVE='on' #saves diagrammes (for timelapse below 'saveframe') 

#directories
route=os.getcwd()+'/'

routeimage=route  #where to print images
routef=route.split('/BlaVar/')[0]+'/BlaVar/Var/' #directory of data points for SED
fermifile=route.split('/BlaVar/')[0]+'/BlaVar/Lightcurve_Simulation/' # fermifile/Nameobj/ is the directory of real Fermi LC
nameflc='myfile_new.txt' # name of the file of the real Fermi LC
imtyp='.png' #.eps .png  save figures in that format


"""Selection of Object and Features"""
# The object + the name of the variable parameter is detected from the named directory
objects=['3C273','3C279','PKS2155304','J1544-3']
namevar=route.split('/')[-2].split('_')[0]
POWER=float(route.split('^')[1].split('_')[0])
deltapower=4
if namevar not in ['le','B','gmax','lext']:
    namevar=input("Give Name of Var:\t")
obj=os.getcwd().split('/')[-2]
if obj not in objects:
    print(objects)
    oo=int(input("\n Give Number of object [0 , 1 ..] from following list:\t"))
    obj=objects[oo]
oo=objects.index(obj)

units={'le':'','B':'G','lext':'','gmax':'','theta':'^o','Gamma':''}
nunits=units[namevar]

#Bibliograpgical Data for objects
redsh=[0.158, 0.536, 0.0, 0.171]
Dists=[749,3115,543, 847]
Rbb=[6.3e17,5e17,10e12] #Radii of the BB radiation zone
SSC=['no','no','on','on'] 
BBextra=['on','on','no','on']
#Second BB into the model
secondBB=BBextra[oo] #add an additional BB on spectrum (usually DISK) but that is not involved in external photon fields of code
Gamma2=[1,1,1,1][oo] #Deboosting/Boosting for observation
factbb=[10,500,'-','-'][oo] #(L_disk-L_BLR)/L_BLR
Tbb2=[11000,20000,'-','-'][oo] #T_disk
Rbb2=[10**15.89,10**15.62,'-','-'][oo] #disk radii



"""Variation"""
#timesteps to plot as stored in fort.81 (warning! double times during crashing not removed yet)
n1=1 #min 1
nmax=-1#maximum number of iterations , -1 for length of array
initt=0 #MJD  : converts given value of t_0 in MJD
#Timing analysis on segments for long light curves
SEGnum=1 #Default is 1 for full LC analysis // Use greater number to bootstrap results,
                                            #splitting the time-interval into segments

"""Timelaspe"""
pausetime=0.05 #seconds
pastlength, lnwdth1, lnwdth2 =20 , 0.2 , 0.1 #previous states plotted in diagram half with lnwdth1 and half lnwdth2
saveframes='on' #only saves total SED frames into .jpgs
monochromVIR, radio, gammarays, Xrays='on', 'on', 'on', 'on'  #only if TIMELAPSE on
"""Plot frames in SED with specific time interval"""
sampleframes='no' #to sample and plot less states in the SED='on' option
dtsmpl=100 # t_cross (check your time - interval for selection)


"""Timecurves and Bands"""
TIMECURVES='no' #plot TC diagrams and Color(t)
addBB='on'      # add the black body spectrums of disk and BLR
addnoise='on'   # add Poisson noiser at SMARTS bands
tbin=1 # bin output (days) -1 to use code output binning (<0.25 tcross / delta)


"""Fermi LC"""
realLC='no'  #Load the real Fermi-LAT LC UNEDITED for making histograms 
indobs=-1.0 #<-1 index of F(epsilon) for modelling fermi observations
normrf=-1 #-1 if normalized Real LC is used
if oo==2:
    normrf=3.017823941751141e-10
if oo==3:
    normrf=1e-9
realhist='no' #add histogram near to the SED of all states


bands=['fermi','xband','smarts'] #more to less energetic

smband=1 #which band of SMARTs to keep from array of created SMARTS channels based on code bins
         #each SMARTS channel created has width 0.4 (log nu) and usually matches,
         #print 'sm + index', e.g. sm1 to find the central frequency of the SMARTS band generated
ixband=0 # 0:soft 2-10 keV  / 1:hard 10-80 keV / -1 Xtreme UV 10-124 eV "bands" list
addtcross=1 #to compare escaped photons with input parameter variation (default 1)
timespan=50 #days to zoom in in timecurve  (create zoome in window to study short-term variability)


"""PLOTTING"""
"""SED settings, SED='on'"""
limSEDs=500 #select nmax-n1 less than this limit to plot multiple SEDs
Multicolor='no' #Plot data points with different color for every interval of 2 years
diconec={'le':'m','B':'c','lext':'y','gmax':'m'}
onecolor=diconec[namevar] #'no' for mixed color and 1.0 transparancy
#Indexing legend up-right
lgd2opts={'le':'(a)','B':'(b)','lext':'(c)','theta':'(d)'}

"""Contour settings, if SED='con' """
cmps={'le':plt.cm.RdPu,'B':plt.cm.GnBu,'delta':plt.cm.YlGn,\
      'lext':plt.cm.YlOrRd,'gmax':plt.cm.PuBu}
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
plt.close('all')

#Multicolor Data Points 
data_name=routef+obj+'/'+obj+".ascii"
dtpts=[]
if Multicolor=='on':
   dtpts.append(routef+obj+'/'+obj+'_08-10.ascii')
   dtpts.append(routef+obj+'/'+obj+'_10-12.ascii')
   dtpts.append(routef+obj+'/'+obj+'_12-14.ascii')
   dtpts.append(routef+obj+'/'+obj+'_14-16.ascii')
   legend_pts=['\'08-\'10','\'10-\'12','\'12-\'14','\'14-\'18']
   colors =['tomato','gold','limegreen','purple']
   form=['o','.','s','+']
   fillstyle=['full','full','full','none']
   #dtpts.append(obj+'_16-18.ascii')

"""Color"""
CI='on' #calculates color if 'on', else omit all color diagrams and info
cibnda , cibndb = 'B', 'J'  #bluer, redder color / choose from below
SMdic={'B':10**14.837, 'V':10**14.74, 'R':10**14.67, 'I':10**14.575,\
       'J':10**14.383, 'H':10**14.267, 'K':10**14.18} #SMARTs filters
FWHMdic={'B': 0.2, 'V': 0.16, 'R': 0.25, 'I': 0.19,\
         'J': 0.12, 'H': 0.18, 'K': 0.16} #FWHM forSMARTS filters
ZPdic={'B':4000, 'V':3580, 'R':2971, 'I':2405,\
       'J':1565, 'H':1039, 'K':647.6} #zero point for SMARTS filters
range_smarts=[10**14.17,10**14.837]
range_fermi=[0.1*10**9*(2.4184*10**14),300*10**9*(2.4184*10**14)] #eV/h = 2.41e14 
range_EUV=[10*(2.4184*10**14) ,124*(2.4184*10**14)] #eV/h = 2.41e14
range_xsoft=[2e3*(2.4184*10**14),10e3*(2.4184*10**14)] #eV/h = 2.41e14
range_xhard=[1e4*(2.4184*10**14),8e4*(2.4184*10**14)] #eV/h = 2.41e14
SMARTSmulti='on' # 'on' generate multiple smarts functions and study one of them with width 0.4,
                 # 'no' integrate over the SMARTS band
Nint=10000 #intergration bins (defaults 10^4)

"""Correlations"""
SPOTS='on' #spots 'on' / Connect dots 'no'
mscorr=3.0 #markersize of scatter plot points
CONTOUR='on' #plot contour lines of 1 sigma included points
corrfont=16 #size of x,ylabel +2 tile -3 ticks
legsize='x-large' #legend size  medium, large x-large etc
colcont='k'#'tab:orange' #color of contour in all plots
lvls=[1/2] # list of values,  fraction of max density in CPs, to contour on

"""Timelags"""
#QUICK METHOD is commended out
#PEAKS #deactivated in this version // yields timelags using peaks or minima
#TIMELAGS='no' #SWITCH/ does not control Stingray
#num_peaks=25 #top peaks to compare for timelags
#peaktype='maxima' #'maxima' or 'minima'
#lendel=50 #delete -lendel, lendel pointsaround the peaks, to remove it

#Discrete Correlation Function and Structure Function, same properties for structure function
taubin=1.0
lentau=50 # days|| integer: -before/ + after zero timelag to calculate

#Timelag per Frequency
segment_size=20 #days #depends on Run Time Interval

"""Data Reduction"""
nu_cutoff=8 #energy in mec^2 to cutoff plotting 
HEIGHT=2 #height of y-axis (threshold) in logscale measured below the minimum value of flux
THRESHOLD=15 #orders of magnituded for data reduction from maximum value

"""%%%%%%%%%%%%%%%%%%    END OF SETTINGS   %%%%%%%%%%%%%%"""

#
#
#
# Activate all tools besides SED
if ALLTOOLS=='on': 
    CORR, TIMECURVES, CI , DCFcalc, STINGRAY =['on']*5

"""Constants"""
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

"""The fort files to use"""
nmlist=os.listdir(route)
names=[]
cnt=0
for nm in nmlist:
    if 'fort_' in nm and '.85' in nm and 'steady' not in nm:
        cnt=cnt+1 #.85, .89 files        
for ia in range(cnt):
    names.append('fort_'+str(ia))
if len(names)==0:
    names.append('fort')
#names=['fort']


"""INPUT Observable Parameters"""
cn4=open(routef+obj+'/code_new4.inp','r')
lines=cn4.readlines()
cn4.close()

tend=float(lines[0].split()[-1].replace('\n',''))
nsteps=float(lines[0].split()[-2])
R=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
B=float(lines[2].split()[1]) #dimension of source in cm
p=float(lines[4].split()[3])
loggmin=float(lines[4].split()[1])
loggmax=float(lines[4].split()[2])
logle=float(lines[4].split()[4])
delta0=float(lines[7].split()[0]) #set zero if used the above for computation
Gamma=float(lines[8].split()[0])
T=float(lines[5].split()[1])
lext=float(lines[5].split()[2])
D=Dists[oo]*10**6*pc #"Distance of Source in Mpc"
zsh=redsh[oo] #redshift
dilate=1-((((1+zsh)**2-1)/((1+zsh)**2+1))**2)**2 #dilation factor, https://www.scirp.org/journal/paperinformation.aspx?paperid=65598
SSConly=SSC[oo] #if SSC or EC modelling is used
tcross=R/c/3600/24 #days jet frame
tcr=tcross*dilate/delta0

"""DEFINITIONS"""
if CONTOUR=='on' or SED=='con':
#    def surface_density_per_energy(m1, m2,lims):
#                if lims:
#                    #xmin, xmax=lims[0:2]
#                    xmin, xmax, =min(m1), max(m1),
#                    ymin, ymax= lims[2::]
#                else:
#                    xmin, xmax, ymin, ymax=min(m1), max(m1), min(m2), max(m2)
#                X, Y = np.mgrid[xmin:xmax:2j, ymin:ymax:50j]
#                hist=np.histogram(m2,bins=list(Y[0])+[ymax],density=True)[0]
#                Z =[hist,hist] #uniform in x-axis
#                return X, Y, Z        
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
if SSConly=='no':
        def BB(x,T,lext,Rbb):
            Ivbb=3*x+np.log10((2*me**3*c**4/h**2)/(np.exp(10**x*me*c**2/kb/T)-1)+10**-100)
            opacity=(lext*me*c**2/sigmaT/R)/aBB/T**4
            #make it dimensionelss / scaling
            Ibb=x+np.log10(me*c**2/h)+Ivbb+np.log10(opacity*4*np.pi*Rbb**2)\
                +np.log10(sigmaT/(4*np.pi*R*me*c**3))
            return Ibb
"""
%
%
%
"""
#%%
"""Observational Data"""
points=ascii.read(data_name)
Cunitx=0 #4 Hz 
Cunity=0 #for erg/cm^2/s
v_pts=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
vFv_pts=np.array([points["col2"]])[0]-np.ones(1)*Cunity
#errorv=abs(points["col4"])
errorvFv=np.array(abs(points["col4"]-points["col3"])/2)
#average measurement errors
#erro=np.mean(errorvFv[v_pts<np.log10(10**15.2)][v_pts[v_pts<np.log10(10**15.2)]>np.log10(10**14.0)])
#errx=np.mean(errorvFv[v_pts<np.log10(range_xhard[1])][v_pts[v_pts<np.log10(range_xhard[1])]>np.log10(range_xsoft[0])])
#errg=np.mean(errorvFv[v_pts<np.log10(range_fermi[1])][v_pts[v_pts<np.log10(range_fermi[1])]>np.log10(range_fermi[0])])
#average errorbars from LCs
erro=[0.00307/2.5, 0.00307/2.5, 0.005/2.5, 0.005/2.5][oo] #average magnitude error from SMARTS LCs
errx=np.mean(errorvFv[v_pts<np.log10(range_xhard[1])][v_pts[v_pts<np.log10(range_xhard[1])]>np.log10(range_xsoft[0])])
errg=[0.1456, 0.1456 , 0.26115 , 0.26115][oo] # average relevant error for error < norms cases
if Multicolor=='on':
  for i in range(len(dtpts)):
    points=ascii.read(dtpts[i])
    globals()['v_pts'+str(i)]=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
    globals()['vFv_pts'+str(i)]=np.array([points["col3"]])[0]-np.ones(1)*Cunity
    globals()['errorv'+str(i)]=0.009
    globals()['errorvFv'+str(i)]=np.array([points["col4"]])[0]

func_obs=inpl.interp1d(v_pts,vFv_pts)

"""
%
%
%
%
"""
"""READING FILES"""
#number of bins per time-step (use first run)
name=names[0]
fphotons=open(route+name+'.85','r')
f85=fphotons.readlines()
fphotons.close()
fst=open(route+'steady.85','r')
fsteady=fst.readlines()
fst.close()
felectrons=open(route+name+'.89','r')
f89=felectrons.readlines()
felectrons.close()

LEN, LENst,  LEN89 = 0 , 0 , 0
for ln in f85:
    if float(ln.split("  ")[1])==0. and LEN!=0:
       break
    LEN+=1
       
for ln in fsteady:
    if float(ln.split("  ")[1])==0. and LENst!=0:
       break
    LENst+=1

for ln in f89:
   if float(ln.split("  ")[1])==0. and LEN89!=0:
      break
   LEN89+=1
#read fort files .85 (.81 ) , steady state .85 and .89
n0=0
epsfile, Ifile, gammafile, Nfile, filetime, filevar = [],[],[],[],[],[]
#read multi.log to correct for times of different runs
#IF INCOMPATIBILITY WITH INPUT RUNTIME ADN Tmax
#if len(names)>1:
#    mlog=open(route+'multi.log','r')
#    tms=mlog.readlines()
#    mlog.close()
#    tms=np.array(tms).astype(np.float)

for name in names:
    #Read fort.81 or .85
    fphotons=open(route+name+'.85','r')
    f85=fphotons.readlines()
    fphotons.close()
    felectrons=open(route+name+'.89','r')
    f89=felectrons.readlines()
    felectrons.close()
    
    n0+=int(len(f85)/LEN)
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
    
    timebin=tend/nsteps    
    #read fort.55
    name55=name+'.55' #name55='fakeTC.txt' #read analytically all the input fake variables
    varfile=open(route+name55,'r')
    for y in varfile:
        x=y.split()
        z=x[0]
        w=x[1].replace('\n','')
        filetime.append(z)
        filevar.append(w)
    varfile.close()
epsfile=np.array(epsfile).astype(np.float)
Ifile=np.array(Ifile).astype(np.float)
gammafile=np.array(gammafile).astype(np.float)
Nfile=np.array(Nfile).astype(np.float)
filetime=np.array(filetime).astype(np.float) 
filevar=np.array(filevar).astype(np.float)
deps=round(epsfile[-1]-epsfile[-2],5) 

addt=0
for it in range(len(filetime)):
    if filetime[it]<1e-5 and it>1:
        addt=filetime[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
        # the time of steady state (conjoint part of two TC)
    filetime[it]=filetime[it]+addt
    
print('Runtime calculated from adding up all fakeTC_*.txt files:\t',\
      str(max(filetime*tcross/delta0)))
print('Check if this number agrees with the initialized runtime.\n Edit manually if gaps are introduced during the run')

    
if namevar=='B':
        filevar=10**filevar
        
if nmax==-1:
        nmax=n0-len(names) #-3 if added SS at the end/ times 0. 0. 5.

""" STEADY STATE"""
#STEADY STATE
epsst,Ist= [], []
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
epsst=np.array(epsst).astype(np.float)[-LENst+1::]
Ist=np.array(Ist).astype(np.float)[-LENst+1::]
epsst=epsst[Ist>max(Ist)-THRESHOLD]
Ist=Ist[Ist>max(Ist)-THRESHOLD]
xsteady=np.log10(delta0*10**epsst*(me*c**2/h)/(1+zsh))  #observing (SED) units
ysteady=np.log10(10**Ist*delta0**4/D**2*R*me*c**3/sigmaT/3)  #observing (SED) units
if nbx==-1:
    nbx=LEN-1
xsed=np.linspace(min(epsfile[1:LEN]),max(epsfile[1:LEN]),nbx)

"""BLACK BODIES (BB)"""
#Correct for 3c273 (lessen the BB from bibliographical values)
if oo==0:
    Gamma=2.36*Gamma

eps85bb, Ibb =epsst ,epsst
if SSConly=='no':
    if len(eps85bb)<200:
        eps85bb=np.linspace(min(epsst),max(epsst),200)
        Rext=Rbb[oo]
        Ibb=BB(eps85bb,T/Gamma,lext/Gamma**2,Rext)+np.log10(Gamma)
        xbb=eps85bb+np.log10((me*c**2/h))
        ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT)  #observing (SED) units
        bbmodel=inpl.interp1d(10**xbb,10**ybb)  #model BB
        #save BB in file
        BBfile=open('BB.txt','w')
        for i,j in zip(xbb,ybb):
            BBfile.write(str(i)+'\t'+str(j)+'\n')
        BBfile.close()
        if secondBB=='on':
                Ibb2=BB(eps85bb,Tbb2*Gamma2,(Gamma2)**2/Gamma**2*\
                        factbb*(Rext/Rbb2)**2*lext,Rbb2)+np.log10(Gamma)
                xbb2=eps85bb+np.log10((me*c**2/h)/Gamma2)
                ybb2=np.log10(10**Ibb2*Gamma2**-2/D**2*R*me*c**3/sigmaT)
                bb2model=inpl.interp1d(10**xbb2,10**ybb2)
        #save BB2 in file
        BB2file=open('BB2.txt','w')
        for i,j in zip(xbb2,ybb2):
            BB2file.write(str(i)+'\t'+str(j)+'\n')
        BB2file.close()
                

"""Generate SMARTs Bands"""
#SMART bands: one bin and one after out of the range of SMARTS
ism=0
for ix in range(len(xsteady)-1):
    if xsteady[ix]>np.log10(range_smarts[0]):
        globals()['sm'+str(ism)]=10**xsteady[ix-1]
        ism=ism+1
    if xsteady[ix]>np.log10(range_smarts[1]):
        globals()['sm'+str(ism)]=10**xsteady[ix]
        break
    
""" Real Fermi-LAT LC used for that object/ check info file in Lightcurve Simulation"""
if realLC=='on':
        realfile=fermifile+obj+'/'+nameflc     
        rlines=open(realfile)
        rt, rf =[] ,[]
        for r in rlines:
            rt.append(float(r.split()[0]))
            rf.append(float(r.split()[1]))
        rlines.close()
        rt=np.array(rt)
        if name55=='fakeTC.txt':
            rt=rt+filetime[list(filevar).index(max(filevar))]*tcr\
               -rt[rf.index(max(rf))]+addtcross*tcr
        else:
            rt=rt+filetime[list(filevar).index(max(filevar))]\
              -rt[rf.index(max(rf))]+addtcross*tcr
        rf=np.log10(np.array(rf))
        obsg=[10**np.mean([np.log10(range_fermi[0]),np.log10(range_fermi[1])])] #integral
        #obsg=np.logspace(np.log10(range_fermi[0]),np.log10(range_fermi[1]),5) #5 X monochromatic
        #if indobs!=-1.:
            #Af=(1+indobs)*10**rf/(range_fermi[1]**(1+indobs)-range_fermi[0]**(1+indobs)) #monochromatic
        #if indobs==-1.:
            #Af=10**rf/np.log(range_fermi[1]/range_fermi[0]) #monochromatic 
        Af=10**rf  #integral
            
        if normrf==-1:
            normrf=np.mean(Af)#3.75e-10
        renorm=normrf/np.mean(Af)
        Af=renorm*Af #renormalize
        obsf=[]
        for feps in obsg:
                obsf.append(np.log10(Af*feps**(1+indobs)))
        obsg=np.log10(obsg)

"""ITERATIONS"""
#PLOTTING SETUP
if (SED!='no' or sampleframes=='on') and TIMELAPSE!='on':
   lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),\
                              round(min(vFv_pts)-HEIGHT),int(max(vFv_pts)+HEIGHT)]
   if realhist=='on':
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
            
            
if sampleframes=='on':
    SED='no'
    TIMELAPSE='no'
    sflgd=[]

if SSConly=='on':
    lgd2opts['theta']='(c)'
lgd2=lgd2opts[namevar]
##set arrays for iteration
legs=[]
time , tobs, smarts , xsoft , xhard , fermi, gbrnum = np.zeros(nmax+1), np.zeros(nmax+1), np.zeros(nmax+1), np.zeros(nmax+1), np.zeros(nmax+1), np.zeros(nmax+1), np.zeros(nmax+1)

if SMARTSmulti=='on':
    for ssm in range(ism+1):
                globals()['smarts'+str(ssm)]=np.zeros(nmax+1)
if name55!='fakeTC.txt':
    fktime, fkvar =np.zeros(nmax+1) , np.zeros(nmax+1)
else:
    fktime, fkvar= filetime, filevar

for pst in range(pastlength):
    globals()['xpst'+str(pst)] , globals()['ypst'+str(pst)] = np.array([]) , np.array([])

#Write the SED output
seds=open('SEDs.txt','w')
fermidat=open('FERMI.txt','w')
#tbb=open('fort.81tbb','w')
if sampleframes=='on':
    frms=open('FRAMES.txt','w')


#THE MAIN LOOP
for nj in range(n1+1,nmax):
    eps85=epsfile[int((nj)*LEN):int((nj+1)*LEN)]
    I85=Ifile[int((nj)*LEN):int((nj+1)*LEN)]
    g89=gammafile[int((nj)*LEN89):int((nj+1)*LEN89)] #electrons later section
    N89=Nfile[int((nj)*LEN89):int((nj+1)*LEN89)] #electrons later section
    time[nj-n1+1]=eps85[0] #single fort.85 file
    if len(names)>1:
        time[nj-n1+1]=filetime[nj-n1+1]
    if time[nj-n1+1] not in time[0:nj-n1] and time[nj-n1+1]>0:
        seds.write(str(time[nj-n1+1])+'\n')
        #tbb.write(str(time[nj-n1+1])+'\n')
        
        eps85 , I85 =eps85[1::] , I85[1::]
        g89, N89 = g89[1::], N89[1::]  #electrons later sectionn
        
        """DATA REDUCTION"""    
        """Remove Extremely High Energy Gamma Rays/ Right Cut-off of Frequencies"""
        zz=np.where(eps85>nu_cutoff,eps85,0)
        zzz=zz.nonzero()
        for i in zzz:
            I85[i]=-100 #will be removed by THRESHOLD
        zz=np.where(g89>loggmax+3,g89,0)
        zzz=zz.nonzero()
        for i in zzz:
            N89[i]=-100  #will be removed by THRESHOLD
        
        """Sorting"""
        nuel=0
        for el in sorted(zip(eps85,I85)):
                eps85[nuel]=el[0]
                I85[nuel]=el[1]
                nuel=nuel+1
        
        """Cut off noise """
        bm=max(I85)-THRESHOLD
        eps85=eps85[I85>bm]
        I85=I85[I85>bm]
        
        if len(g89)==0:
            g89=g89n
            N89=g89n**-2
        
        Nm=max(N89)-THRESHOLD/2
        g89=g89[N89>Nm]
        N89=N89[N89>Nm]
        if len(eps85)>1 and len(I85)>1:
            ##Variable parameter
            if name55!='fakeTC.txt' or len(names)==1:
                fktime[nj-n1+1]=filetime[nj-1]
                fkvar[nj-n1+1]=filevar[nj-1]
            delta = (0.5*(fkvar[nj-n1+1]+fkvar[nj-n1])/np.mean(fkvar))**(POWER/deltapower)*delta0
            tcr=tcross*dilate/delta
            tobs[nj-n1+1]=round(time[nj-n1+1]*tcr,3)
            
            
            """Add extrapolated values at high energies """#for plotting reasons// avoid truncation
            eps85=np.append(eps85,(2*eps85[-1]-eps85[-2]))
            I85=np.append(I85,max(I85)-THRESHOLD)
            
            """Observed SED transformation"""
            x85=eps85+np.log10(delta*(me*c**2/h)/(1+zsh)) #observing (SED) units
            y85=I85+np.log10(delta**4/D**2*R*me*c**3/sigmaT/3)  #observing (SED) units
            ymodel=inpl.interp1d(10**x85,10**y85) #model as a function for the above arrays
            funcy85=inpl.interp1d(x85,y85,kind='slinear')
            if SED=='con':
                xa=xsed+np.log10(delta*(me*c**2/h)/(1+zsh)) #observing (SED) units
                xb=xa[xa>min(x85)+0.01]
                xc=xb[xb<max(x85)-0.01]
                for i,j in zip(xc,funcy85(xc)):
                    seds.write(str(i)+'\t'+str(j)+'\n')
            else:
                for i,j in zip(x85,y85):
                    seds.write(str(i)+'\t'+str(j)+'\n')
            
            """BLACK BODY ADDITION"""    
            if SSConly=='no':
                x81tbb=xbb[xbb<max(x85)][xbb[xbb<max(x85)]>min(x85)]
                y81tbb=np.log10(10**ybb[xbb<max(x85)][xbb[xbb<max(x85)]>min(x85)]\
                                        +10**funcy85(x81tbb))
                funcy81tbb=inpl.interp1d(x81tbb,y81tbb)
                if secondBB=='on':
                        funcybb2=inpl.interp1d(xbb2,ybb2)
                        y81tbb=np.log10(10**funcybb2(x81tbb)+10**y81tbb)
                #for i,j in zip(x81tbb,y81tbb):
                       # tbb.write(str(i)+'\t'+str(j)+'\n')
            else:
                    funcy81tbb=inpl.interp1d(x85,y85)
            
            
            ##Intergrate for SMARTS and FERMI bands
            lnm=np.linspace(range_smarts[0],range_smarts[1],Nint)
            lnx=np.linspace(range_xsoft[0],range_xsoft[1],Nint)
            lnxx=np.linspace(range_xhard[0],range_xhard[1],Nint)
            lng=np.linspace(range_fermi[0],min(range_fermi[1],10**max(x85)),Nint)
            if SMARTSmulti!='on':
                if min(lnm)>10**min(x85):
                    smarts[nj-n1+1]=intg.simps(ymodel(lnm)/lnm,lnm)
                else:
                    smarts[nj-n1+1]=smarts[nj-n1+1]
            else:
                for ssm in range(ism+1):
                    if globals()['sm'+str(ssm)]>10**min(x85):
                        globals()['smarts'+str(ssm)][nj-n1+1]=ymodel(globals()['sm'+str(ssm)])*np.log10(sm2/sm1)
                    else:
                        globals()['smarts'+str(ssm)][nj-n1+1]=globals()['smarts'+str(ssm)][nj-n1]
                        print('Error with soft O/IR at:\t'+str(nj))
            if min(lnx)>10**min(x85):
                xsoft[nj-n1+1]=intg.simps(ymodel(lnx)/lnx,lnx)
            else:
                xsoft[nj-n1+1]=xsoft[nj-n1]
                print('Error with soft X-rays at:\t'+str(nj))
            if min(lnxx)>10**min(x85):
                xhard[nj-n1+1]=intg.simps(ymodel(lnxx)/lnxx,lnxx)
            else:
                xhard[nj-n1+1]=xhard[nj-n1]
                print('Error with hard X-rays at step:\t'+str(nj))
            fermi[nj-n1+1]=intg.simps(ymodel(lng)/lng,lng)
            
            """PLOTTING"""
            if TIMELAPSE=='on':
                    plt.figure(1)
                    plt.xlabel(r'$'+xlab+'_{obs}\;\;'+' ['+xaxis+']$',fontsize=15)
                    plt.ylabel(r'$'+yaxis+'_{obs}\;\;['+yunits[yaxis]+']$',fontsize=15)
                    plt.plot(x85,y85,'k-',linewidth=1.5,label='t= {:.2f}'.format(round(tobs[nj-n1+1],2))+' days')
                    if SSConly!='on':
                        if namevar!='lext':
                            plt.plot(xbb,ybb,'b--',linewidth=0.5,label='BLR')
                        else:
                            plt.plot(xbb,ybb+np.log10(fkvar[nj-n1+1]/lext),'b--',linewidth=0.5,label='BLR')                   
                        if secondBB=='on':
                            plt.plot(xbb2,ybb2,'b-.',linewidth=0.5,label='Disk')
                    plt.legend(fontsize='large',framealpha=0.5)
                    plt.title(obj+'  SED',fontsize=18)
                    lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),\
                              round(min(vFv_pts)-HEIGHT),int(max(vFv_pts)+HEIGHT)]
                    plt.axis(lims)
                    if Multicolor=='on':
                        for i in range(len(dtpts)):
                            plt.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                                     globals()['errorv'+str(i)],globals()['errorvFv'+str(i)],\
                                     elinewidth=1.5,capsize=2.5,capthick=0.5,\
                                     fmt=form[i],color=colors[i],ms=3.5)
                    else:
                        plt.errorbar(v_pts,vFv_pts,yerr=errorvFv,elinewidth=1.5,\
                            capsize=2.5,capthick=0.5,fmt='r.',ecolor='red',ms=3.5)
                   # plt.legend(legend_pts)
                   #TRACEBACK CURVES
                    xpst0=x85
                    ypst0=y85                    
                    for pst in range(pastlength-1,int(pastlength/2),-1):
                          if len(globals()['xpst'+str(pst)])!=0:
                           plt.figure(1)
                           plt.plot(globals()['xpst'+str(pst)],globals()['ypst'+str(pst)],'k-',linewidth=lnwdth2)
    
                    for pst in range(int(pastlength/2),0,-1):
                          if len(globals()['xpst'+str(pst)])!=0:
                           plt.figure(1)
                           plt.plot(globals()['xpst'+str(pst)],globals()['ypst'+str(pst)],'k-',linewidth=lnwdth1)
                    xtemp=[xpst0]
                    ytemp=[ypst0]
                    for pst in range(1,pastlength+1):
                        xtemp.append(globals()['xpst'+str(pst-1)])
                        ytemp.append(globals()['ypst'+str(pst-1)])
                    for pst in range(pastlength):
                        globals()['xpst'+str(pst)]=xtemp[pst]
                        globals()['ypst'+str(pst)]=ytemp[pst]
                    #saveframe
                    if saveframes=='on':
                       plt.figure(1).savefig(routeimage+'vFv_'+str(nj)+'.jpg')
                    plt.pause(pausetime)
                    plt.clf()
            elif SED=='on':
                if realhist=='on':
                    ax10.plot(x85,y85,'-',color=onecolor,alpha=0.5,linewidth=0.6-0.1*(nmax//1000))
                else:
                    plt.figure(1)
                    if onecolor!='no':
                        plt.plot(x85,y85,'-',color=onecolor,alpha=0.5,linewidth=0.6-0.1*(nmax//1000))
                    else:
                            plt.plot(x85,y85,'-',linewidth=0.5)
    
            #MULTI WINDOW DIAGRAM      
#            if TIMELAPSE=='on':
#              figall=plt.figure(2)
#              figall.suptitle(obj+'  t='+str(round(tobs[nj-n1+1],2))+' days',fontsize=16)
#              if monochromVIR=='on':
#                
#                nu85=np.log10(delta*10**eps85*(me*c**2/h)/(1+zsh))
#                Fv85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3)-nu85
#                
#                #plt.figure(6)
#                plt.subplot(2,2,2)
#                plt.xlabel(r'$\nu_{obs}\;[Hz]$')
#                plt.ylabel(r'$ F\nu_{obs}\; [erg/cm^2/Hz/s]$')
#                plt.plot(nu85,Fv85,'y-')
#                plt.title('VIR')
#                if Multicolor=='on':
#                  for i in range(len(dtpts)):
#                    xp=globals()['v_pts'+str(i)]
#                    yp=globals()['vFv_pts'+str(i)]-globals()['v_pts'+str(i)]
#                    #erxp=globals()['errorv'+str(i)]
#                    eryp=globals()['errorvFv'+str(i)]
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
#                else:
#                    xp=v_pts
#                    yp=vFv_pts-v_pts
#                    #eryp=errorv
#                    eryp=errorvFv
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt='r.',ecolor='red',ms=3.5)
#                limlow=12.0
#                limhigh=16.5
#                lims=[limlow,limhigh,\
#                min(yp[xp<limhigh][xp[xp<limhigh]>limlow])-HEIGHT,\
#                max(yp[xp<limhigh][xp[xp<limhigh]>limlow])+HEIGHT+0.5]
#                plt.axis(lims)
#    #            if TIMELAPSE=='on':
#    #                legs=[]
#    #                plt.pause(0.5)
#    #                plt.clf()
#                
#              if Xrays=='on':
#                eps5=np.log10(delta*10**eps85/(1+zsh))
#                X85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3)-eps5
#                #plt.figure(8)
#                plt.subplot(2,2,3)
#                plt.xlabel(r'$E_{obs}\;[mec^2]$')
#                plt.ylabel(r'$F_X$    $counts\; [1/cm^2/s]$')
#                plt.plot(eps5,X85,'b-')
#                plt.title('X - rays')
#                if Multicolor=='on':
#                  for i in range(len(dtpts)):
#                    xp=globals()['v_pts'+str(i)]+np.log10(h/me/c**2)
#                    yp=globals()['vFv_pts'+str(i)]+np.log10(1/h)-globals()['v_pts'+str(i)]
#                    #erxp=globals()['errorv'+str(i)]
#                    eryp=globals()['errorvFv'+str(i)]
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
#                else:
#                    xp=v_pts+np.log10(h/me/c**2)
#                    yp=vFv_pts+np.log10(1/h/v_pts)
#                    #eryp=errorv
#                    eryp=errorvFv
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt='r.',ecolor='red',ms=3.5)
#                xp=v_pts+np.log10(h/me/c**2)
#                yp=vFv_pts-xp-np.log10(me*c**2)
#                limlow=-5.0
#                limhigh=1.0
#                lims=[limlow,limhigh,\
#                min(yp[xp<limhigh][xp[xp<limhigh]>limlow])-HEIGHT,\
#                max(yp[xp<limhigh][xp[xp<limhigh]>limlow])+HEIGHT+0.5]
#                plt.axis(lims)
#                
#              if radio=='on':        
#                nu5=np.log10(delta*10**eps85*(me*c**2/h)/(1+zsh))
#                R85=np.log10(10**I85*delta**4/D**2*R*me*c**3/sigmaT/3/Jy)-nu5
#                
#                #plt.figure(10)
#                plt.subplot(2,2,1)
#                plt.xlabel(r'$\nu_{obs}\;\;[Hz]$')
#                plt.ylabel(r'$F\nu\;\; [Jy]$')
#                plt.plot(nu5,R85,'r-')
#                plt.title('Radio')
#                if Multicolor=='on':
#                  for i in range(len(dtpts)):
#                    xp=globals()['v_pts'+str(i)]
#                    yp=globals()['vFv_pts'+str(i)]+np.log10(1/Jy)-globals()['v_pts'+str(i)]
#                    #erxp=globals()['errorv'+str(i)]
#                    eryp=globals()['errorvFv'+str(i)]
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
#                else:
#                    xp=v_pts
#                    yp=vFv_pts+np.log10(1/Jy)-v_pts
#                    #eryp=errorv
#                    eryp=errorvFv
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt='r.',ecolor='red',ms=3.5)   
#                limlow=8.0
#                limhigh=12.5
#                lims=[limlow,limhigh,-5,1]
#                plt.axis(lims)            
#            
#                
#              if gammarays=='on':
#                eps5=np.log10(delta*10**eps85/(1+zsh))
#                gev5=eps5+np.log10(1/eV/10**9*me*c**2)
#                G85=np.log10(10**I85*delta**4/D**2*R*c/sigmaT/3*3600)-eps5
#                
#                #plt.figure(12)
#                plt.subplot(2,2,4)
#                plt.xlabel(r'$E_{obs}\;[GeV]$')
#                plt.ylabel(r'$F_\gamma$   $events\,/h$')
#                plt.plot(gev5,G85,'k-')
#                plt.title('$\gamma$ - rays')
#                if Multicolor=='on':
#                  for i in range(len(dtpts)):
#                    xp=globals()['v_pts'+str(i)]+np.log10(h/eV/10**9)
#                    yp=globals()['vFv_pts'+str(i)]+np.log10(3600/h)-globals()['v_pts'+str(i)]
#                    #erxp=globals()['errorv'+str(i)]
#                    eryp=globals()['errorvFv'+str(i)]
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt=form[i],color=colors[i],ms=3.5)
#                else:
#                    xp=v_pts+np.log10(h/eV/10**9)
#                    yp=vFv_pts+np.log10(3600/h)-v_pts
#                    #erxp=errorv
#                    eryp=errorvFv
#                    plt.errorbar(xp,yp,eryp,elinewidth=1.5,capsize=2.5,\
#                                 capthick=0.5,fmt='r.',ecolor='red',ms=3.5)
#                xp=v_pts+np.log10(h/eV/10**9)
#                yp=vFv_pts+np.log10(3600/h)-v_pts
#                limlow=-3.0
#                limhigh=4.0
#                lims=[limlow,limhigh,\
#                min(yp[xp<limhigh][xp[xp<limhigh]>limlow])-HEIGHT,\
#                max(yp[xp<limhigh][xp[xp<limhigh]>limlow])+HEIGHT+3.5]
#                plt.axis(lims)
#    
#    
#            if TIMELAPSE=='on':
#                plt.pause(pausetime)
#                figall.subplots_adjust(hspace=0.6, wspace=0.5)
#                if saveframes=='on':
#                    figall.savefig(routeimage+'bands_'+str(nj)+'.jpg')
#                plt.clf()
#                
            if sampleframes=='on':
                nnw=nj+0 #offset to catch specific point
                if nnw%int(dtsmpl/tcr)==0:
            #    if nnw in [238,470,661,856,1119,1320]:
                    plt.figure(1)
                    xn85=np.linspace(min(x85),max(x85),300)
                    plt.plot(xn85,funcy85(xn85),'-',linewidth=1.3)
                    sflgd.append(str(round(nnw*tcr+initt,1))+' days')
                    print('Frame index to check in lightcurves/ statitstics:\t '+str(nj))
                    if SSConly=='no':
                        frms.write(str(time[nj-n1+1])+'\n')
                        for i,j in zip(x81tbb,y81tbb):
                            frms.write(str(i)+'\t'+str(j)+'\n')
    
                
            """ ELECTRONS"""     
    #        plt.figure(89)
    #        plt.plot(g89,N89+(p-2)*g89)
    #        plt.axis([0,loggmax+1.0,-12.0,0.0])
    #        plt.xlabel(r' $log\,\gamma$')
    #        plt.ylabel(r'$\gamma^{p-2}N(\gamma)$')
    #        print('\nlog(gmax)=\t'+str(max(g89)))
           
            """Cooling Breaking Calc"""
            if namevar in ['B','lext']:
                dif1=N89[1::]-N89[0:-1]
                dif1mod=inpl.interp1d(g89[1::],dif1)
                g89n=np.linspace(max(loggmin,g89[1]),g89[-1],100)
                dif1=dif1mod(g89n)
            ###Using 2nd derivative for finding breaking point (alter loop below though dif1->dif2)
               # dif2=np.array(dif1[0:-1]-dif1[1::])
               # dif2mod=inpl.interp1d(g89[1:-1],dif2)
               # dif2=dif2mod(g89n)
               # dec89=0.2 #declination of dN89/dg (change of slope) for which to assume the start of the breaking
               # c89=0.81 #offset to higher gamma for actual gbr value (value found from calibration)
                #c89 is found by the difference of the theoretical and analytical value:
                #np.mean(gbr[gbr>loggmin+1][gbr[gbr>loggmin+1]<5.4])-np.mean(gbrnum[gbr>loggmin+1][gbr[gbr>loggmin+1]<5.4])
                for idif in dif1:
                    if idif<(3-2*p)*0.1:
                        gbrnum[nj-n1+1]=g89n[list(dif1).index(idif)]#+c89
                        if gbrnum[nj-n1+1]>g89[-1]:
                            gbrnum[nj-n1+1]=g89[-1]
                        break
                    if fkvar[nj-n1+1]>globals()[namevar]:
                        gbrnum[nj-n1+1]=g89n[0]
                    else:
                        gbrnum[nj-n1+1]=g89n[-1]
            
    #        if TIMELAPSE=='on':
    #             plt.figure(89)
    #             plt.plot(g89,N89+(p-2)*g89)
    #             if namevar in ['B','lext']:
    #                 plt.plot(gbrnum[nj-n1+1]*np.ones(50),np.linspace(-10.0,0.0,50),'k--')
    #             mingobs=loggmin-1.0
    #             if mingobs<0.: mingobs=0.
    #             plt.axis([mingobs,loggmax+0.5,-10.,0.0])
    #             plt.xlabel(r' $log\,\gamma$')
    #             plt.ylabel(r'$\gamma^{p-2}N(\gamma)$')
    #             plt.legend([r' $t= '+str(round(tobs[nj-n1+1],1))+'\; days$',r'$\gamma^{num}_{br}$'])
    #             plt.pause(pausetime)
    #             if saveframes=='on':
    #                plt.figure(89).savefig(routeimage+'electrons_'+str(nj)+'.jpg')
    #             plt.clf()
    else:
        time[nj-n1+1]=0.0
    
#tbb.close()
seds.close()

#remove zeros
if SMARTSmulti=='on':
    for ssm in range(ism+1):
                globals()['smarts'+str(ssm)]=globals()['smarts'+str(ssm)][tobs>0]
else:
    smarts=smarts[tobs>0]
fermi=fermi[tobs>0]
xsoft=xsoft[tobs>0]
xhard=xhard[tobs>0]
if name55!='fakeTC.txt':
    fktime=fktime[tobs>0]+initt/tcr
    fkvar=fkvar[tobs>0]
if namevar=='B':
    lb=fkvar**2/8/np.pi/me/c**2*sigmaT*R
    gbr=np.log10(3/4/(lb+lext))
    gbrnum=gbrnum[tobs>0]
elif namevar=='lext':
    lb=B**2/8/np.pi/me/c**2*sigmaT*R
    gbr=np.log10(3/4/(fkvar+lb))
    gbrnum=gbrnum[tobs>0]

tobs=tobs[tobs>0]+initt
Time=max(time)+initt/tcr


"""X-rays selection"""
xband=xsoft
range_xband=range_xsoft
if ixband>0:
    xband=xhard
    range_xband=range_xhard

"""Smart filter selection"""
#SMARTS filter selection
if SMARTSmulti=='on':
    #smband=1 #from array of created smarts channels based on code bins
    smfq=str(round(np.log10(globals()['sm'+str(smband)]),1))
    smarts=globals()['smarts'+str(smband)]
else:
    smfq=''


"""Add Noise and BB radiation to SMARTS band"""
if addBB=='on' and SSConly=='no':
    noise=np.ones(len(smarts))
    if addnoise=='on':
            noise=np.random.poisson((100)**2,len(smarts))/(100)**2
    lnsm=np.linspace(range_smarts[0],range_smarts[1],Nint)
    smarts=smarts+noise*(intg.simps(bbmodel(lnsm)/lnsm,lnsm)+intg.simps(bb2model(lnsm)/lnsm,lnsm))
    lnx=np.linspace(range_xband[0],range_xband[1],Nint)
    xband=xband+intg.simps(bbmodel(lnx)/lnx,lnx)+intg.simps(bb2model(lnx)/lnx,lnx)
    if namevar=='lext':
            smarts=smarts+noise*(intg.simps(bbmodel(lnsm)/lnsm,lnsm)*fkvar/np.mean(fkvar)+intg.simps(bb2model(lnsm)/lnsm,lnsm))
            xband=xband+intg.simps(bbmodel(lnx)/lnx,lnx)*fkvar/np.mean(fkvar)+intg.simps(bb2model(lnx)/lnx,lnx)    
ci=[]
filterA, filterB= [], []
for t in range(len(tobs)):
    Aa,Bb=[],[]
    for ssm in range(ism+1):
        Aa.append(globals()['sm'+str(ssm)]) #nuFnu
        Bb.append(globals()['smarts'+str(ssm)][t]/np.log10(sm2/sm1)) #nuFnu
    smmodel=inpl.interp1d(Aa,Bb)
    af=smmodel(SMdic[cibnda])*FWHMdic[cibnda]
    bf=smmodel(SMdic[cibndb])*FWHMdic[cibndb]
    filterA.append(-2.5*np.log10(af/(SMdic[cibnda]*ZPdic[cibnda]*Jy*FWHMdic[cibnda])))
    filterB.append(-2.5*np.log10(bf/(SMdic[cibnda]*ZPdic[cibndb]*Jy*FWHMdic[cibnda])))
    ci.append(filterA[-1]-filterB[-1])
filterA=np.array(filterA).astype(np.float)
filterB=np.array(filterB).astype(np.float)
ci=np.array(ci).astype(np.float)
#label color
cilbl='('+cibnda+' - '+cibndb+') '
      

"""NEW BINNING for output LCs /  Default -1 for code's output binning"""
#NEW LCs save as %nane%2
#defaults if tbin = -1
arn=['xsoft','xhard','xband','fermi','smarts','fkvar','fktime']
if CI=='on':
    arn.append('ci')
    arn.append('filterA')
    arn.append('filterB')
if namevar in ['B','lext']:
    arn.append('gbr')
    arn.append('gbrnum')
arn.append('tobs')
for ar in arn:
    globals()[ar+'B']=globals()[ar]
#change if different binning selected
if tbin!=-1:
    dtold=tobs[2]-tobs[1]
    #slices = np.linspace(0,max(tobs), int((max(tobs)-1)/tbin)+1, True).astype(np.int)
    binsB=int(max(tobs)/tbin)+1
    pcntg=[] #for generating randomness in SMARTS obs while binning
    for ar in arn:
        globals()[ar+'B']=st.binned_statistic(tobs,globals()[ar],statistic='mean',bins=binsB)[0]
        #sample randomly from the observations of one day for SMARTS band to create smartsB , ciB
        if ar in ['smarts','ci']:
            stmin=st.binned_statistic(tobs,globals()[ar],statistic='min',bins=binsB)[0]
            stmax=st.binned_statistic(tobs,globals()[ar],statistic='max',bins=binsB)[0]
            if len(pcntg)==0:
                pcntg=np.random.random(binsB)
            globals()[ar+'B']=stmin+(stmax-stmin)*pcntg 
        #periodic condition for averaging, else we end up with huge last value and neglegible first value of timeseries
        #firlas=np.mean(globals()[ar])
        #globals()[ar][0]=firlas
        #globals()[ar][-1]=firlas
    taubin=tbin
    segment_size=int(tbin/dtold*segment_size)
    tbin=round(tobsB[1]-tobsB[0],2)


if SAVE=='on':
    """"Write file""" #does not consider new binning
    np.savetxt('FERMI.txt', np.c_[tobsB,fermiB])
    if sampleframes=='on':
        frms.write('STEADY\n')
        for i,j in zip(xsteady,ysteady):
            frms.write(str(i)+'\t'+str(j)+'\n')
        frms.close()
    np.savetxt('SMARTS.txt',np.c_[tobsB,globals()['smartsB']])
    np.savetxt('X_SOFT.txt', np.c_[tobsB,xsoftB])
    np.savetxt('X_HARD.txt', np.c_[tobsB,xhardB])
    if CI=='on':
        np.savetxt('CI.txt', np.c_[tobsB,filterAB,filterBB,ciB])
    np.savetxt('fakeTCB.txt', np.c_[fktimeB,fkvarB])  

    print('TIME CURVES CREATED AND SAVED')
#%% continue with plotting
if TIMELAPSE!='on' and (SED!='no' or sampleframes=='on'):
        """Contour Plot"""
        plt.figure(1)
        if SED=='con':
        # 2-D Histogram --> to be commented out
            seds=open('SEDs.txt','r')
            ss, xseds, yseds=[],[],[]
            for j in seds:
                ss.append(j)
            for j in ss:
                if len(j.split())>1:
                      xseds.append(j.split()[0])
                      yseds.append(j.split()[1])
            seds.close()
            xseds=np.array(xseds).astype(np.float)
            yseds=np.array(yseds).astype(np.float)
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
        if Multicolor=='on':                     
            for i in range(len(dtpts)):
                ax10.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                          yerr=globals()['errorvFv'+str(i)],\
                          elinewidth=1.5,capsize=2.5,capthick=0.5,ecolor=colors[i],\
                          fmt=form[i],color=colors[i],ms=5.0,fillstyle=fillstyle[i],\
                          label=legend_pts[i])
        else:
            ax10.errorbar(v_pts,vFv_pts,yerr=errorvFv,elinewidth=1.5,capsize=2.5,\
                         fmt='r.',ecolor='red',ms=3.5,label='Obs')
        if SSConly=='no':
            ax10.plot(xbb,ybb,'r--',linewidth=0.9,label='BLR')
        if secondBB=='on':
            ax10.plot(xbb2,ybb2,'r-.',linewidth=0.9,label='Disk')
        lg4=ax10.legend( loc='upper left',fontsize='large',framealpha=0.5)
        if sampleframes=='on'  and nmax/dtsmpl<10:
                ax10.add_artist(plt.legend(sflgd,fontsize='large',framealpha=0.5))
        if lgd2!=-1:
                ax10.annotate(lgd2, xy=(0.9,0.9), xycoords='axes fraction', size=14,weight='bold')
        ax10.add_artist(lg4)
        ax10.xaxis.set_minor_locator(AutoMinorLocator())
        ax10.yaxis.set_minor_locator(plt.LinearLocator(4*(lims[3]-lims[2])+1))
        if namevar=='lext':
            xTbb=np.log10(kb*T/me/c**2*np.array([0.3,9.0]))+np.log10(delta0*(me*c**2/h)/(1+zsh))
            ax10.fill_between(np.linspace(xTbb[0], xTbb[1]),-20,20, facecolor='grey', alpha=0.2, hatch='/',edgecolor='grey')
            
        if 'on' in SED:
                #plot scale errors
            wbd1, wbd2, wbdx, wbdf  =  FWHMdic[cibndb], FWHMdic[cibnda],\
                np.log10(range_xhard[1]/range_xhard[0])/2,\
                np.log10(range_fermi[1]/range_fermi[0])/2
            ax10.errorbar(np.log10(SMdic[cibndb]),lims[3]-.9,xerr=wbd1,\
                          color='red',ms=0.0,capsize=3.5)
            ax10.annotate(cibndb,(np.log10(SMdic[cibndb])-0.25,lims[3]-0.8),fontsize=13)
            ax10.errorbar(np.log10(SMdic[cibnda]),lims[3]-0.9,xerr=wbd2,\
                          color='c',ms=0.0,capsize=3.5)
            ax10.annotate(cibnda,(np.log10(SMdic[cibnda])-0.25,lims[3]-0.8),fontsize=13)
            ax10.errorbar(np.mean(np.log10(range_xsoft)),lims[3]-0.9,\
                          xerr=wbdx,color='grey',ms=0.0,capsize=3.5)
            ax10.annotate('2-10 keV',(np.log10(range_xsoft[0]),lims[3]-0.8),fontsize=13)
            ax10.errorbar(np.mean(np.log10(range_fermi)),lims[3]-0.9,\
                          xerr=wbdf,color='blue',ms=0.0,capsize=3.5)
            ax10.annotate('Fermi-LAT',(np.log10(range_fermi[0]),lims[3]-0.8),fontsize=13)


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
        
        if realhist=='on':
            #HISTOGRAM
            ax1[1].hist(np.log10(fermi),bins=50,orientation='horizontal',density=True,color=onecolor,alpha=0.5,label='sim')
            ax1[1].hist(obsf,bins=50,orientation='horizontal',density=True,color='tab:orange',alpha=0.5,label='obs')
            #ax1[1].set_title('Simulation/Observations',size=16)
            ax1[1].set_xlabel(r'Fermi-LAT PDF')
            ax1[1].legend(fontsize='large',framealpha=0.5)
            ax1[1].axis([0,max(max(np.histogram(np.log10(fermi),density=True)[0]),\
               max(np.histogram(obsf,density=True)[0]))*1.1,lims[2],lims[3]])

if SAVE=='on':
    titlenums='_'+str(n1)+'-'+str(nmax)
    if tbin!=-1:
        titlenums=titlenums+'_tbin'+str(tbin)
    if (nmax-n1)<limSEDs and SED!='no':
                    plt.figure(1).savefig(routeimage+SED.strip('on')+'vFv_'+namevar+titlenums+'.eps',bbox_inches='tight')
                    plt.figure(1).savefig(routeimage+SED.strip('on')+'vFv_'+namevar+'.png',bbox_inches='tight')
    elif (nmax-n1)>limSEDs and SED!='no':
        print('ERROR: Too many SEDs(over '+str(limSEDs)+') to introduce into one diagram')
        if realhist!='on':
                    plt.figure(1).savefig(routeimage+SED.strip('on')+'vFv_'+namevar+'.png',bbox_inches='tight')
        else:
                    plt.figure(1).savefig(routeimage+SED.strip('on')+'vFv_'+namevar+'.png',bbox_inches='tight')
plt.show()
#%% SEGMENTATION
for tt in ['','B']:
    #set Time Curves to segment initial and 1-day binned
    tx=globals()['tobs'+tt]
    sm, xs, xh , xb, fer ,fkvr, fkt = globals()['smarts'+tt], globals()['xsoft'+tt],\
                            globals()['xhard'+tt], globals()['xband'+tt],globals()['fermi'+tt], globals()['fkvar'+tt] , globals()['fktime'+tt]
    if CI=='on':
        cii, fiA, fiB = globals()['ci'+tt], globals()['filterA'+tt], globals()['filterB'+tt],  
    if namevar in ['B','lext']:
        ggb, ggbn = globals()['gbr'+tt],  globals()['gbrnum'+tt]
    seg, topr= 0 , 0
    lenseg=int(np.ceil(len(tx)/SEGnum))
    for to in range(1,len(tx)+1):
        if to%lenseg==0 or to==len(tx):
            seg+=1
            #initially 1 tcross resolved TCs
            globals()['TCs'+tt+str(seg)]=[tx[topr:to],sm[topr:to],xs[topr:to],xh[topr:to],\
                   xb[topr:to],fer[topr:to],fkvr[topr:to], fkt[topr:to]]
            if CI=='on':
              globals()['TCs'+tt+str(seg)].append(cii[topr:to])
              globals()['TCs'+tt+str(seg)].append(fiA[topr:to])
              globals()['TCs'+tt+str(seg)].append(fiB[topr:to])
    
            if namevar in ['B','lext']:
              globals()['TCs'+tt+str(seg)].append(ggb[topr:to])
              globals()['TCs'+tt+str(seg)].append(ggbn[topr:to])
            topr=to

for seg in range(1,SEGnum+1):
    tobs, smarts ,xsoft ,xhard ,xband ,fermi ,fkvar ,fktime = globals()['TCs'+str(seg)][0:8]
    tobsB,smartsB,xsoftB,xhardB,xbandB,fermiB,fkvarB , fktimeB = globals()['TCsB'+str(seg)][0:8]
    if CI=='on':
        ci ,filterA ,filterB  =globals()['TCs'+str(seg)][8:11]
        ciB,filterAB,filterBB =globals()['TCsB'+str(seg)][8:11]
    if namevar in ['B','lext']:
        gbr , gbrnum = globals()['TCs'+str(seg)][11::]
        gbrB, gbrnumB= globals()['TCsB'+str(seg)][11::]
    if SEGnum>1:
          pr=imtyp
          imtyp='_seg'+str(seg)+pr
          imtyp=pr
          print('SEGMENT\t\t'+str(seg))
#%%"""PLOTTING TIMECURVES"""
    nmv=namevar
    fktm, fkvr , nunit = fktimeB, fkvarB , units[nmv]
    ftmin, fvarmin =(fktm[list(fkvr).index(min(fkvr))]+addtcross/np.log(2))*tcr,\
                  np.log10(min(fkvr)*np.mean(fermi*2)/np.mean(fkvr))
    ftmax, fvarmax =(fktm[list(fkvr).index(max(fkvr))]+addtcross/np.log(2))*tcr,\
                  np.log10(max(fkvr)*np.mean(fermi*2)/np.mean(fkvr))
                  
    if TIMECURVES=='on':  
        #TIMECURVES
        
        fig20,ax20=plt.subplots(num=20,figsize=(15,5))
        ax20.plot(tobsB,np.log10(smartsB),'r-')
        ax20.plot(tobsB,np.log10(fermiB),'c-',tobsB,np.log10(xbandB),'y-')
        ax20.plot((fktm+addtcross/np.log(2))*tcr,\
                      np.log10(fkvr*(np.mean(fermiB*2)/np.mean(fkvr))))
        ax20.plot(ftmax,fvarmax,'r*',ftmin,fvarmin,'k*',ms=7.5)
        plt.title(obj,fontsize=16)
        plt.xlabel(r'$t_{obs}$[ $days$ ]',fontsize=14)
        plt.ylabel(r'$F_{obs}$[ $erg\,/cm^2\,/sec$ ]',fontsize=14)
        tclgd=['O/IR '+smfq,'gamma-rays','X-rays',r'$\log  '+nmv+'$ offset ',\
               r'max ${'+nmv+'} = '+str(round(max(fkvr),int(-np.log10(max(fkvr)))+1))+nunit+' $',r'min $ '+nmv+' = '+\
               str(round(min(fkvr),int(-np.log10(min(fkvr)))+1))+nunit+' $' ]
        plt.legend(tclgd,fontsize=legsize,markerscale=2,framealpha=0.5,loc='lower left')
        limya20=np.floor(min(min(np.log10(fermiB)),min(np.log10(smarts)),\
                             min(np.log10(xband)),np.log10(max(fkvr)*\
                                 np.mean(fermiB*2)/np.mean(fkvr))))
        limyb20=np.ceil(max(max(np.log10(fermiB)),max(np.log10(smarts)),\
                            max(np.log10(xband)),np.log10(max(fkvr)*\
                                np.mean(fermiB*2)/np.mean(fkvr))))
        ax20.axis([initt,int(Time*tcr)+1,limya20,limyb20])
        #ax20.xaxis.set_minor_locator(plt.LinearLocator(5*int(Time+1)))
        ax20.xaxis.set_minor_locator(AutoMinorLocator())
        ax20.yaxis.set_minor_locator(AutoMinorLocator())
        
        #ZOOM IN 
        #Zoom in Time-LightCurve 50days
        if timespan>2*(max(tobs)-min(tobs)):
            timespan=max(tobs)-min(tobs)
            rndtime=int(n1*tcr*dilate/delta0)
        else:
            rndtime=int(np.random.random()*max(tobs-2*timespan))\
                     +min(tobs)+timespan+initt
        irnd=list(tobs).index(min(tobs[tobs>rndtime]))
        trnd=list(tobs).index(max(tobs[tobs<rndtime+50]))
    #    irndB=list(tobsB).index(min(tobsB[tobsB>rndtime]))
    #    trndB=list(tobsB).index(max(tobsB[tobsB<rndtime+50]))
        
        
        if max(tobs)>150:
            # location for the zoomed portion  in fig 20
            ax24= plt.axes([0.65,0.4, .3, .5]) # units unity for axis scales
            ax24.tick_params(axis='y',labelright='on',\
                      labelleft='off',right='on',direction='in',labelsize=10)
            ax24.tick_params(axis='x',labelsize=9,pad=-15)
            plt.title(str(timespan)+' days')
            # insert the zoomed figure
            plt.setp(ax24)    
        else:
            fig24,ax24=plt.subplots(num=24,figsize=(10,5))
            plt.title(obj,fontsize=16)
            plt.xlabel(r'$t_{ obs} $[ $days$ ]',fontsize=14)
            plt.ylabel(r'$F_{ obs}$ [ $erg\,/cm^2\,/sec$ ]',fontsize=14)
            tclgd2=['O/IR '+smfq,'Gamma-rays','X-rays',r'$\log '+nmv+' $+ offset ',\
            r'max $'+nmv+' = '+str(round(max(fkvr),2))+nunit+' $',r'min $ '+nmv+' = '+\
            str(round(min(fkvr),int(-np.log10(min(fkvr))+1)))+nunit+' $' ]
            if CI=='on':
                ax24.plot(tobs,ci-np.mean(ci)+np.mean(np.log10(smarts)),'g') #oldci
                tclgd2.append(cilbl+' (offset -<CI>+'+str(round(np.mean(np.log10(smarts))))+')')
            plt.legend(tclgd2,fontsize=legsize,markerscale=2,framealpha=0.5)
        ax24.plot(tobsB,np.log10(smartsB),'r-')
        ax24.plot(tobsB,np.log10(fermiB),'c-',tobsB,np.log10(xbandB),'y-')
        ax24.plot((fktm+addtcross/np.log(2))*tcr,\
                          np.log10(fkvr*(np.mean(fermi*2)/np.mean(fkvr))))
        limya24=round(min(min(np.log10(fermi[irnd:trnd])),min(np.log10(smarts[irnd:trnd])),\
                          min(np.log10(xband[irnd:trnd])),np.log10(min(fkvar[irnd:trnd])*\
                              np.mean(fermi*2)/np.mean(fkvar))),1)-0.05
        limyb24=round(max(max(np.log10(fermi[irnd:trnd])),max(np.log10(smarts[irnd:trnd])),\
                          max(np.log10(xband[irnd:trnd])),np.log10(max(fkvar[irnd:trnd])*\
                              np.mean(fermi*2)/np.mean(fkvar))),1)+0.05
        ax24.axis([rndtime,rndtime+timespan,limya24,limyb24])
        #ax24.xaxis.set_minor_locator(AutoMinorLocator())
        #ax24.yaxis.set_minor_locator(AutoMinorLocator())
        
        """COLOR""" #CI + SMARTS +var  Zoomed in
        if CI=='on':
            fig19,ax19=plt.subplots(num=19,figsize=(15,5))
            ax19.plot(tobs,np.log10(smarts),color=colors[0],alpha=0.8)
            ax19.plot(tobs,ci-np.mean(ci)+np.mean(np.log10(smarts)),color='g',lw=2.0) #oldci
            #ax19.plot((fktm+addtcross/np.log(2))*tcr,\
            #              np.log10(fkvar*(np.mean(smarts*2)/np.mean(fkvar))),alpha=0.5)
            plt.xlabel(r'$t_{ obs}$ [ $days$ ]',fontsize=14)
            plt.ylabel(r'$F_{ obs}$ [ $erg\,/cm^2\,/sec$ ]',fontsize=14)
            #plt.title(obj+'    Color '+cilbl,size=16)
            tclgd=['O/IR '+smfq,cilbl+' (offset'+str(round(np.mean(np.log10(smarts))))+')']  #,'log '+namevar+' (offset'+str(round(np.mean(np.log10(smarts))))+')']
            plt.legend(tclgd,fontsize=legsize,markerscale=2,framealpha=0.5)
            limya19=min(ci-np.mean(ci))+np.mean(np.log10(smarts))-0.1 #oldci
            limyb19=max(ci-np.mean(ci))+np.mean(np.log10(smarts))+0.1 #oldci
            ax19.axis([rndtime,rndtime+timespan,limya19-0.1,limyb19+0.1])
            ax19.xaxis.set_minor_locator(AutoMinorLocator())
            ax19.yaxis.set_minor_locator(AutoMinorLocator())
        
        """Gamma Break""" #numerical , theoretical (varying process) + color
        if namevar in ['B','lext']:
            fig90,ax90=plt.subplots(num=90,figsize=(15,5))
            ax90.plot(tobs,gbrnum,label=r'numerical $\gamma_{br}$')
            ax90.plot(tobs,gbr,label=r'analytical $\gamma_{br}$')
            ax90.plot(tobs,np.ones(len(tobs))*loggmin,'k--',linewidth=0.5,\
                      label=r'$\gamma_{min}='+str(round(loggmin,2))+'$')
            ax90.plot(tobs,np.ones(len(tobs))*loggmax,'k-.',linewidth=0.5,\
                      label=r'$\gamma_{max}='+str(round(loggmax,2))+'$')
            ax90.plot(tobs,ci+np.mean(gbrnum-ci),label=cilbl)
            ax90.axis([0,max(tobs),loggmin-0.2,loggmax+0.2])
            plt.xlabel(r'$t_{ obs}$ [ $days$ ]',fontsize=14)
            plt.ylabel(r'$log\,\gamma_{br}\; / C.I. $',fontsize=14)
            plt.title(r'Cooling Break $\gamma_{break}$    '+obj,fontsize=16)
            plt.legend(tclgd,fontsize=legsize,markerscale=2,framealpha=0.5)
    
    
        if SAVE=='on':                
            plt.figure(20).savefig(routeimage+'LCs_'+namevar\
                      +titlenums+imtyp,bbox_inches='tight')
            #plt.figure(24).savefig(routeimage+'LCs_'+namevar+'_'+str(round(timespan))+'d_'+str(tbin)+imtyp,bbox_inches='tight')
            if CI=='on':
                plt.figure(19).savefig(routeimage+'CI_TC_'\
                          +namevar+titlenums+imtyp,bbox_inches='tight')
            if namevar in ['B','lext']:
                plt.figure(90).savefig(routeimage+'gbr_TC_'\
                          +namevar+titlenums+imtyp,bbox_inches='tight')    
    #%% HISTOGRAMS / PDFs
    if TIMECURVES=='on':
        plt.figure(6)
        plt.hist(np.log10(fermiB),bins=50,density=True,color='b',alpha=0.5,label='gamma-rays')
        plt.hist(np.log10(xsoftB),bins=50,density=True,color='tab:orange',alpha=0.5,label='2-10 keV')
        plt.hist(np.log10(xhardB),bins=50,density=True,color='grey',alpha=0.5,label='10-80 keV')
        plt.title('Probability Density Functions',size=16)
        plt.ylabel(r'PDFs',size=14)
        plt.xlabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]',size=14)
        plt.legend()
        if SAVE=='on':
            plt.figure(6).savefig('XraysPDFs'+imtyp,bbox_inches='tight')
    #%% CORRELATIONS
    #SPOTS='no'
    #CONTOUR='on'
    if CORR=='on':
        #FIT POWER LAW
        def line(x,a):
            return a*x
        incl=optmz.curve_fit(line,np.log10(smartsB/np.mean(smartsB)),np.log10(fermiB/np.mean(fermiB)))[0][0]
        fit21='no' # 'on' /'no'
        
        fig21,ax21=plt.subplots(num=21)
        if SPOTS=='on':
        ##PLOT MARKERS
            ax21.plot(np.log10(smartsB/np.mean(smartsB)),np.log10(fermiB/np.mean(fermiB)),\
                      '.',color='lightcoral',ms=mscorr,zorder=0)
            if CONTOUR=='on':
                X,Y,Z=surface_density(np.log10(smartsB/np.mean(smartsB)),np.log10(fermiB/np.mean(fermiB)),[])
                ax21.contour(X,Y,Z,colors=[colcont]*len(lvls),linestyles='dashed',\
                             linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
            spots=1
        else:
        ##PLOT LINES
            ax21.plot(np.log10(smartsB/np.mean(smartsB)),np.log10(fermiB/np.mean(fermiB)),'r-',linewidth=0.6)
            spots=0
        aa,bb,cc,dd=[ax21.get_xlim()[0],ax21.get_xlim()[1],ax21.get_ylim()[0],ax21.get_ylim()[1]]
        ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
        axcorr= [ ae, be ,ce ,de]
        axlim1, axlim2 =min(ae,ce) , max(be,de)
        ax21.plot(np.linspace(axlim1,axlim2),np.linspace(axlim1,axlim2),'k-')
        ax21.plot(np.linspace(axlim1,axlim2),2*np.linspace(axlim1,axlim2),'k-.')
        if fit21=='on':
            ax21.plot(np.log10(smartsB/np.mean(smartsB)),incl*np.log10(smartsB/np.mean(smartsB)),'m-')
        plt.xlabel(r'$log\,F_{'+smfq+'}/<F_{opt}>$',fontsize=corrfont)
        plt.ylabel(r'$log\,F_{\gamma}/<F_\gamma>$',fontsize=corrfont)
        plt.title(obj+r'   O/IR${'+smfq+'}$ vs $\gamma$-rays',fontsize=corrfont+2)
        ax21.legend([r'$\tau_{lag}=0\;$ days',r'$F_\gamma \propto F_{'+smfq+'}$', r'$F_\gamma \propto (F_{'+smfq+'})^2$',\
                     r'$F_{opt} \propto F_\gamma^{'+str(round(1/incl,2))+'}$'],fontsize=legsize,markerscale=2)
        ax21.axis(axcorr)
        ax21.xaxis.set_minor_locator(AutoMinorLocator())
        ax21.yaxis.set_minor_locator(AutoMinorLocator())
        ax21.tick_params(axis='y',labelsize=corrfont-3)
        ax21.tick_params(axis='x',labelsize=corrfont-3)
        plt.tight_layout()
    #    #plot scale errors
    #    xerr, yerr =  erro, errg
    #    plt.errorbar(axcorr[0]+0.3*(axcorr[1]-axcorr[0]),axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #    plt.annotate('obs error',(axcorr[0]+0.1*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr))
    
        #
        #
        #
        #CORRELATION SMARTS - Xrays
        fig22,ax22=plt.subplots(num=22)
        ##PLOT MARKERS
        if SPOTS=='on':
            ax22.plot(np.log10(smartsB/np.mean(smartsB)),np.log10(xbandB/np.mean(xbandB)),\
                      '.',color='lightsteelblue',ms=mscorr,zorder=0)
            if CONTOUR=='on':
                X,Y,Z=surface_density(np.log10(smartsB/np.mean(smartsB)),np.log10(xbandB/np.mean(xbandB)),[])
                ax22.contour(X,Y,Z,colors=[colcont]*len(lvls),linestyles='dashed',\
                             linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
        ##PLOT LINES
        else:
            ax22.plot(np.log10(smartsB/np.mean(smartsB)),np.log10(xbandB/np.mean(xbandB)),\
                      marker='-',color='lightsteelblue',linewidth=0.6)
        aa,bb,cc,dd=[ax22.get_xlim()[0],ax22.get_xlim()[1],ax22.get_ylim()[0],ax22.get_ylim()[1]]
        ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
        axcorr= [ ae, be ,ce ,de]
        axlim1, axlim2 =min(ae,ce) , max(be,de)
        ax22.plot(np.linspace(axlim1,axlim2),np.linspace(axlim1,axlim2),'k-')
        ax22.plot(np.linspace(axlim1,axlim2),2*np.linspace(axlim1,axlim2),'k-.')
        plt.xlabel(r'$log\,F_{'+smfq+'}/<F_{opt}>$',fontsize=corrfont)
        plt.ylabel(r'$log\,F_{X}/<F_X>$',fontsize=corrfont)
        plt.title(obj+r'    O/IR ${'+smfq+'}$ vs X-rays',fontsize=corrfont+2)
        ax22.legend([r'$\tau_{lag}=0\;$ days',r'$ F_X \propto F_{'+smfq+'}$',\
                     r'$ F_X \propto (F_{'+smfq+'})^2$'],fontsize=legsize,markerscale=2)
        ax22.axis(axcorr)
        ax22.xaxis.set_minor_locator(AutoMinorLocator())
        ax22.yaxis.set_minor_locator(AutoMinorLocator())
        ax22.tick_params(axis='y',labelsize=corrfont-3)
        ax22.tick_params(axis='x',labelsize=corrfont-3)
        plt.tight_layout()
    #    #plot scale errors
    #    xerr, yerr =  erro, errx
    #    plt.errorbar(axcorr[0]+0.3*(axcorr[1]-axcorr[0]),axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #    plt.annotate('obs error',(axcorr[0]+0.1*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr))
    
        #
        #
        #
        #
        #CORRELATION X-rays Fermi
        fig23,ax23=plt.subplots(num=23)
        ##PLOT MARKERS
        if SPOTS=='on':
            ax23.plot(np.log10(xbandB/np.mean(xbandB)),np.log10(fermiB/np.mean(fermiB)),'.',color='cornflowerblue',ms=mscorr,zorder=0)
            if CONTOUR=='on':
                X,Y,Z=surface_density(np.log10(xbandB/np.mean(xbandB)),np.log10(fermiB/np.mean(fermiB)),[])
                ax23.contour(X,Y,Z,colors=[colcont]*len(lvls),linestyles='dashed',linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
        else:
        ##PLOT LINES
            ax23.plot(np.log10(xbandB/np.mean(xbandB)),np.log10(fermiB/np.mean(fermiB)),'b-',linewidth=0.6)
        aa,bb,cc,dd=[ax23.get_xlim()[0],ax23.get_xlim()[1],ax23.get_ylim()[0],ax23.get_ylim()[1]]
        ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
        axcorr= [ ae, be ,ce ,de]
        axlim1, axlim2 =min(ae,ce) , max(be,de)
        ax23.plot(np.linspace(axlim1,axlim2),np.linspace(axlim1,axlim2),'k-')
        ax23.plot(np.linspace(axlim1,axlim2),2*np.linspace(axlim1,axlim2),'k-.')
        plt.xlabel(r'$log\,F_{X}/<F_{X}>$',fontsize=corrfont)
        plt.ylabel(r'$log\,F_{\gamma}/<F_\gamma>$',fontsize=corrfont)
        plt.title(obj+r'    X-rays vs $\gamma$-rays',fontsize=corrfont+2)
        ax23.legend([r'$\tau_{lag}=0\;$ days',r'$F_{\gamma} \propto F_X$',\
                     r'$F_{\gamma} \propto (F_X)^2$'],fontsize=legsize,markerscale=2)
        ax23.axis(axcorr)
        ax23.xaxis.set_minor_locator(AutoMinorLocator())
        ax23.yaxis.set_minor_locator(AutoMinorLocator())
        ax23.tick_params(axis='y',labelsize=corrfont-3)
        ax23.tick_params(axis='x',labelsize=corrfont-3)
        plt.tight_layout()
    #    #plot scale errors
    #    xerr, yerr =  errx, errg
    #    plt.errorbar(axcorr[0]+0.3*(axcorr[1]-axcorr[0]),axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #    plt.annotate('obs error',(axcorr[0]+0.1*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr))
       
    
        if CI=='on':
            #CORRELATION smarts - color intensity
            #smvci=str(round(np.log10(sm2),2)) #same as below (frequency of the bin used )
            smvci=cibndb
            smci=10**filterBB #select x axis from SMARTS bands and make it color
            fig17,ax17=plt.subplots(num=17)
            if SPOTS=='on':
            ##PLOT MARKERS
                #ax17.plot(np.log10(smci)-np.mean(np.log10(smci)),ci-np.mean(ci),'c.',ms=mscorr,zorder=0) #oldci
                ax17.plot(np.log10(smci),ciB,'c.',ms=mscorr,zorder=0)
    
                if CONTOUR=='on':
                    #X,Y,Z=surface_density(np.log10(smci)-np.mean(np.log10(smci)),ci-np.mean(ci),[]) #oldci
                    X,Y,Z=surface_density(np.log10(smci),ciB,lims=[]) #oldci
                    ax17.contour(X,Y,Z,colors=[colcont]*len(lvls),linestyles='dashed',\
                                 linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
            else:
            ##PLOT LINES
                ax17.plot(np.log10(smci)-np.mean(np.log10(smci)),ciB,'c-',linewidth=0.6)
            #plt.xlabel(r'$'+smvci+'\,-<'+smvci+'>$',size=corrfont) #oldci
            #plt.ylabel(r'$'+cilbl+'\, -<'+cilbl+'>$',labelpad=-0.5,size=corrfont) #oldci
            plt.xlabel(r'$'+smvci+'$',size=corrfont) #oldci
            plt.ylabel(r'$'+cilbl+'$',labelpad=-0.5,size=corrfont) #oldci
            plt.title(obj+r'    $'+smvci+'$ vs '+cilbl,size=corrfont+2)
            ax17.legend([r'$\tau_{lag}=0\;$ days'],fontsize=legsize,markerscale=2)
            aa,bb,cc,dd=[ax17.get_xlim()[0],ax17.get_xlim()[1],ax17.get_ylim()[0],ax17.get_ylim()[1]]
            ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
            axcorr= [ ae, be ,ce ,de]
            ax17.axis(axcorr)
            ax17.invert_xaxis()
            ax17.xaxis.set_minor_locator(AutoMinorLocator())
            ax17.yaxis.set_minor_locator(AutoMinorLocator())
            ax17.tick_params(axis='y',labelsize=corrfont-3)
            ax17.tick_params(axis='x',labelsize=corrfont-3)
            #plot scale errors
    #        xerr, yerr =  erro*2.5, 1.4*erro*2.5
    #        plt.errorbar(axcorr[0]+0.2*(axcorr[1]-axcorr[0])+xerr,axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #        plt.annotate('obs error',(axcorr[0]+0.03*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr),fontsize=corrfont-3)
    #        plt.tight_layout()
            
            
            
        #CORRELATION GAMMA RAYS  COLOR INTENSITY
            fig18,ax18=plt.subplots(num=18)
            if SPOTS=='on':
                ##PLOT MARKERS
                #ax18.plot(np.log10(fermi/np.mean(fermi)),ci-np.mean(ci),'.',color='limegreen',ms=mscorr,zorder=0) #oldci
                ax18.plot(np.log10(fermiB/np.mean(fermiB)),ciB,'.',color='limegreen',ms=mscorr,zorder=0) #oldci
                if CONTOUR=='on':
                    #X,Y,Z=surface_density(np.log10(fermi/np.mean(fermi)),ci-np.mean(ci),[]) #oldci
                    X,Y,Z=surface_density(np.log10(fermiB/np.mean(fermiB)),ciB,[]) #oldci
                    ax18.contour(X,Y,Z,colors=[colcont]*len(lvls),linestyles='dashed',\
                                 linewidths=2.0,levels=list(np.array(lvls)*np.amax(Z)))
            else:
                ##PLOT LINES
                #ax18.plot(np.log10(fermi/np.mean(fermi)),ci-np.mean(ci),'g-',linewidth=0.6) #oldci
                ax18.plot(np.log10(fermiB/np.mean(fermiB)),ciB,'g-',linewidth=0.6) 
            plt.xlabel(r'$log\,F_{\gamma}/<F_{\gamma}>$',size=corrfont)
            #plt.ylabel(r'$'+cilbl+'\, -<'+cilbl+'>$',labelpad=-0.5,size=corrfont) #oldci
            plt.ylabel(r'$'+cilbl+'$',labelpad=-0.5,size=corrfont)
            plt.title(obj+r'     $\gamma$-rays vs '+cilbl,size=corrfont+2)
            ax18.legend([r'$\tau_{lag}=0\;$ days'],fontsize=legsize,markerscale=2)
            aa,bb,cc,dd=[ax18.get_xlim()[0],ax18.get_xlim()[1],ax18.get_ylim()[0],ax18.get_ylim()[1]]
            ae,be,ce,de = (1-np.sign(aa)*0.05)*aa, (1+np.sign(bb)*0.05)*bb, (1-np.sign(cc)*0.05)*cc, (1+np.sign(dd)*0.05)*dd
            axcorr= [ ae, be ,ce ,de]
            plt.axis(axcorr)
            #ax18.axis([8.5,14.5,0.0,2.5])
            ax18.xaxis.set_minor_locator(AutoMinorLocator())
            ax18.yaxis.set_minor_locator(AutoMinorLocator())
            ax18.tick_params(axis='y',labelsize=corrfont-3)
            ax18.tick_params(axis='x',labelsize=corrfont-3)
            #plot scale errors
    #        xerr, yerr =  errg, 1.4*erro*2.5
    #        plt.errorbar(axcorr[0]+0.2*(axcorr[1]-axcorr[0])+xerr,axcorr[2]+0.05*(axcorr[3]-axcorr[2])+yerr,yerr=yerr,xerr=xerr,color='k',ms=0.0,capsize=1)
    #        plt.annotate('obs error',(axcorr[0]+0.03*(axcorr[1]-axcorr[0]),axcorr[2]+0.04*(axcorr[3]-axcorr[2])+yerr),fontsize=corrfont-3)
    #        plt.tight_layout()
    
    
    ##CORRELATION CI - gamma break
    #        if namevar in ['B','lext']: 
    #            fig92,ax92=plt.subplots(num=92)
    #            axlim1=min(gbrnum)-0.1
    #            axlim2=max(gbrnum)+0.1
    #            axcorr=[axlim1,axlim2,min(ci-np.mean(ci))*0.9,max(ci-np.mean(ci))*1.1]
    #            ##PLOT MARKERS
    #            if SPOTS=='on':
    #                ax92.plot(gbrnum,ci -np.mean(ci),'b.',ms=mscorr)
    #            else:
    #                ##PLOT LINES
    #                ax92.plot(gbrnum,ci-np.mean(ci),'b-',linewidth=0.6)
    #            ax92.plot(np.linspace(axlim1,axlim2),np.linspace(axlim1,axlim2),'k-')
    #            plt.ylabel(r'$'+cilbl+'\, -<'+cilbl+'>$',fontsize=corrfont)
    #            plt.xlabel(r'$log\,\gamma_{br,cool}$',fontsize=corrfont)
    #            plt.title(obj+'   '+cilbl+r' vs $\gamma_{br,cool}$',fontsize=corrfont+2)
    #            ax92.legend([r'$CI \propto \gamma_{br,cool}$','correlation'],fontsize=legsize,markerscale=2)
    #            ax92.axis(axcorr)
    #            ax92.xaxis.set_minor_locator(AutoMinorLocator())
    #            ax92.yaxis.set_minor_locator(AutoMinorLocator())
    #            ax92.tick_params(axis='y',labelsize=corrfont-3)
    #            ax92.tick_params(axis='x',labelsize=corrfont-3)
    #            plt.tight_layout()
    #    
        
        if SAVE=='on':
            plt.figure(21).savefig(routeimage+'corr_'+namevar+'s-g'\
                      +spots*'_SP'+titlenums+imtyp,bbox_inches='tight')
            plt.figure(22).savefig(routeimage+'corr_'+namevar+'s-x'\
                      +spots*'_SP'+titlenums+imtyp,bbox_inches='tight')
            plt.figure(23).savefig(routeimage+'corr_'+namevar+'x-g'\
                      +spots*'_SP'+titlenums+imtyp,bbox_inches='tight')
            if CI=='on':
                plt.figure(18).savefig(routeimage+'corr_'+namevar+'g-ci'\
                          +spots*'_SP'+titlenums+imtyp,bbox_inches='tight')
                plt.figure(17).savefig(routeimage+'corr_'+namevar+'s'\
                          +smvci+'-ci'+spots*'_SP'+titlenums+imtyp,bbox_inches='tight')
    #        if namevar in ['B','lext']:
    #            plt.figure(92).savefig(routeimage+'corr_'+namevar+'s'+smvci+'-gbr'+spots*'_SP'+titlenums+imtyp)
    
    #%%    
    """TIMELAGS"""
    #"""Find timelag from peaks or minima// QUICK METHOD"""
    #if TIMELAGS=='on':
    #    
    #    smod=inpl.interp1d(tobs,smarts,kind='quadratic')
    #    fmod=inpl.interp1d(tobs,fermi,kind='quadratic')
    #    xmod=inpl.interp1d(tobs,xband,kind='quadratic')
    #    tobsm=np.linspace(min(tobs),max(tobs),5*len(tobs))
    #    
    #    if len(tobs)>num_peaks*lendel*2:
    #        sm=smod(tobsm)
    #        fm=fmod(tobsm)
    #        xm=fmod(tobsm)
    #        inds=[]
    #        indf=[]
    #        indx=[]
    #        for i in range(num_peaks):
    #            if peaktype=='maxima':
    #                x=list(sm).index(max(sm))
    #                y=list(fm).index(max(fm))
    #                z=list(xm).index(max(xm))
    #            elif peaktype=='minima':
    #                x=list(sm).index(min(sm))
    #                y=list(fm).index(min(fm))
    #                z=list(xm).index(min(xm))
    #            else: break
    #            inds.append(x)
    #            indf.append(y)
    #            indx.append(z)
    #            sm=np.concatenate([sm[lendel:max(x-lendel,lendel)],sm[min(x+lendel,len(sm)-1-lendel):len(sm)-1-lendel]])
    #            fm=np.concatenate([fm[lendel:max(y-lendel,lendel)],fm[min(y+lendel,len(sm)-1-lendel):len(sm)-1-lendel]])
    #            xm=np.concatenate([xm[lendel:max(z-lendel,lendel)],xm[min(z+lendel,len(sm)-1-lendel):len(sm)-1-lendel]])
    #        
    #        timelags=[[],[],[]] # select below the delayed band in loop
    #        for i in range(len(inds)) :
    #            if abs(inds[i]-indf[i])<50 and abs(indf[i]-indx[i])<50 and abs(inds[i]-indx[i])<50 and peaktype in ['minima','maxima']:
    #                timelags[0].append(round(tobsm[indf[i]]-tobsm[inds[i]],3))
    #                timelags[1].append(round(tobsm[indx[i]]-tobsm[inds[i]],3))
    #                timelags[2].append(round(tobsm[indf[i]]-tobsm[indx[i]],3))
    #                
    #        print('\n Calculated timelags from 25 flares [delayed - firstly observed]:\n',timelags[0],'\n',timelags[1],'\n',timelags[2],\
    #              '\n Gamma Rays after SMARTS (Av.+/- std):\t '+str(round(np.mean((timelags[0])),3))+'+/-'+str(round(np.std(timelags[0]),3))+' days',\
    #              '\n Xrays after SMARTS (Av. +/- std):\t '+str(round(np.mean((timelags[1])),3))+'+/-'+str(round(np.std(timelags[1]),3))+' days',\
    #              '\n Xrays after Gammarays (Av. +/- std):\t '+str(round(np.mean((timelags[2])),3))+'+/-'+str(round(np.std(timelags[2]),3))+' days')
        
    """ DCF """
    if DCFcalc=='on':
            #CALCULATING
            dtau=tobs[3]-tobs[2]
            if taubin<dtau:
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
                tlags.append(np.array(tlag).astype(np.float))
                DCFs.append(np.array(DCF).astype(np.float))
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
            pk=DCFs[0][pospeak] #peak of SMARTS fermi DCF
            if abs(tpk)<2.0:
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
                    ax34.axis([tpk-2.5,tpk+2.5,pk*0.85,min(1.0,pk*1.05)])
                else:
                    ax34.axis([tpk-2.5,tpk+2.5,max(pk*1.05,-1.0),pk*0.85])
                from mpl_toolkits.axes_grid1.inset_locator import mark_inset
                mark_inset(ax25, ax34, loc1=1, loc2=2, fc="none", ec="0.5")
            plt.show()
    
            if SAVE=='on':
                fig25.savefig(routeimage+'tlags_'+namevar+titlenums+imtyp,bbox_inches='tight')
                np.savetxt('DCFs.txt', np.c_[np.concatenate(tlags, axis=0 ),np.concatenate(DCFs, axis=0 )])

                
                
    """ SF """
    if SFcalc=='on':
            logtau=np.logspace(np.log10(tbin/2.),np.log10(tobsB[-1]/2.))#SF
            fsf=logtau[2]/logtau[1]
            sflags, SFs= [], []
            fig27=plt.figure(27)
            lnbd=len(bands)
            errbands=[erro,errx,errg]
            minSF=1.0
            maxSF=1.0
            for bnd in range(lnbd):
                a=-2.5*np.log10(globals()[bands[bnd-1]+'B']) #mags
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
                sflags.append(np.array(sflag).astype(np.float))
                SFs.append(np.array(SF).astype(np.float))
                plt.loglog(sflags[bnd],SFs[bnd],'-',label=bands[bnd])
                minSF=min(min(SF),minSF)
                maxSF=max(max(SF),maxSF)
            plt.legend()
            plt.xlabel(r'$\tau_{ obs}\;\;[days]$',size=14)
            plt.ylabel(r'$SF$($\tau_{ obs}$)',size=14)
            plt.title(obj+'     Structure Function',size=16)
            ll=len(tlags[bnd])
            plt.axis([min(logtau),max(logtau),0.95*minSF,1.05*maxSF])
            plt.show()
            
            if SAVE=='on':
                fig27.savefig(routeimage+'SF_'+namevar+titlenums+imtyp,bbox_inches='tight')
    #%%
    """Stingray/ Timelags per frequency range"""
    if STINGRAY=='on':
        long_times=tobsB
        nn=0
        figstg,axstg=plt.subplots(3,1,num=26,sharex='col')
        axstg[0].set_title(obj+'    Frequency-dependent lags',size=16)
        for bnd in range(len(bands)):
            long_signal_1=globals()[bands[bnd-1]+'B']  #+B--> if new binning 
            long_signal_2=globals()[bands[bnd]+'B']
            
    
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
                
        if SAVE=='on':
            figstg.savefig(routeimage+'tlags_freq_'+namevar+titlenums+imtyp,bbox_inches='tight')    
    #%%
    """TIMELAGS with Extrapolation (sort of) using neighbour points of DCF"""
    if 'DCFs.txt' in nmlist:
        alldcf=open(route+'DCFs.txt','r')
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
            

        
    if DCFcalc=='on':
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
        globals()['tlag'+str(seg)] = tlag
        globals()['DCFs'+str(seg)] = DCFs
            
    """Fractional Variability"""
    #turn to CVs
    erro, errx, errg = 0., 0. , 0. 
    FVin=np.sqrt(np.var(fkvar)/np.mean(fkvar)**2)
    FVo=np.sqrt(np.var(smartsB)/np.mean(smartsB)**2-(10**erro-1)**2)
    FVx=np.sqrt(np.var(xbandB)/np.mean(xbandB)**2-(10**errx-1)**2)
    FVg=np.sqrt(np.var(fermiB)/np.mean(fermiB)**2-(10**errg-1)**2)
    Nn=len(fermiB)
    dFVo=np.sqrt(FVo**2+(((2/Nn)**0.5*(10**erro-1))+(((10**erro-1)/Nn)**0.5*2*FVo)**2)**0.5)-FVo
    dFVx=np.sqrt(FVx**2+(((2/Nn)**0.5*(10**errx-1))+(((10**errx-1)/Nn)**0.5*2*FVx)**2)**0.5)-FVx
    dFVg=np.sqrt(FVg**2+(((2/Nn)**0.5*(10**errg-1))+(((10**errg-1)/Nn)**0.5*2*FVg)**2)**0.5)-FVg
    
    print( 'Fractional Variabilities (daily binned LCs) \n SMARTS band:\t'+str(FVo)+'\t+/-'+str(dFVo)\
          +'\n X-rays band:\t'+str(FVx)+'\t+/-'+str(dFVx)\
          +'\n FERMI band:\t'+str(FVg)+'\t+/-'+str(dFVg))
    
    #SAVE Segment statistics and properties
    globals()['FVin'+str(seg)] =FVin
    globals()['FVo'+str(seg)], globals()['dFVo'+str(seg)]=FVo, dFVo
    globals()['FVx'+str(seg)], globals()['dFVx'+str(seg)]=FVx, dFVx
    globals()['FVg'+str(seg)], globals()['dFVg'+str(seg)]=FVg, dFVg
    if SEGnum>1:
        TIMECURVES='no'
        plt.pause(1)
        plt.close('all')
#%%COMPARE SEGMENTS
if SEGnum>1:
    fig71, ax71 = plt.subplots(num=71)
    for seg in range(1,SEGnum):
        FVin=globals()['FVin'+str(seg)]
        ax71.plot(FVin, globals()['FVo'+str(seg)],'r.')
        ax71.plot(FVin, globals()['FVx'+str(seg)],'gx')
        ax71.plot(FVin, globals()['FVg'+str(seg)],'b+')    
    ax71.set_title('Segment Size: '+str(lenseg)+' observational days',size=corrfont+2)
    ax71.set_xlabel(r'$CV_{in}$',size=corrfont)
    ax71.set_ylabel(r'$CV_{sim}$',size=corrfont)
    ax71.legend(['O/IR',r'$X-$rays', r'$\gamma-$rays'])
    fig71.savefig('CVg_CVin_'+obj+imtyp,bbox_inches='tight')
    
    
    #BOOTSTRAP DCFs
    if DCFcalc=='on':
        sDCF , xDCF, fDCF = np.zeros(len(taurange)) ,np.zeros(len(taurange)) ,np.zeros(len(taurange))
        for seg in range(1,SEGnum):
            for tl in range(0,len(taurange)):
                ttl=taurange[tl]
                ind=list(globals()['tlag'+str(seg)]).index(ttl)
                sDCF[tl]+=globals()['DCFs'+str(seg)][0][ind]
                xDCF[tl]+=globals()['DCFs'+str(seg)][1][ind]
                fDCF[tl]+=globals()['DCFs'+str(seg)][2][ind]
        
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
        fig72.savefig('DCFsboot_'+str(SEGnum)+imtyp,bbox_inches='tight')
    
#%%
INDEXES='no'
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
    
    gsm=0.5*np.log10(sm1/2.8/10**6/delta/B)
    gx=0.5*np.log10(range_xsoft[0]/2.8/10**6/delta/B*(1+zsh))
    #
    tcsm=tcsyn/10**gsm
    tcx=tcsyn/10**gx
#%
#%
#%
#%
#"""
if INDEXES=='on':
##actual unrenormalized histogram
    fakevar=inpl.interp1d(fktime*tcr,fkvar)
    try:
        vec=fakevar(tobs)
    except ValueError:
        vec=fkvar
    Nbin=50
    vect=[np.log10(vec),np.log10(smarts),np.log10(xband),np.log10(fermi)]
    vecnames=[namevar,'smarts','xband','fermi']
    
    Nb2=25
    bvec,hvec=np.histogram(vect[0],bins=Nb2)
    bsm,hsm=np.histogram(vect[1],bins=Nb2)
    bx,hx=np.histogram(vect[2],bins=Nb2)
    bf,hf=np.histogram(vect[3],bins=Nb2)
        
    
    import scipy.signal as sg
    pvec=hvec[list(bvec).index(max(bvec))+1]
    psm=hsm[list(bsm).index(max(bsm))+1]
    px=hx[list(bx).index(max(bx))+1]
    pf=hf[list(bf).index(max(bf))+1]
    
    wvec=sg.peak_widths(bvec,[list(bvec).index(max(bvec))])[0][0]*(hvec[1]-hvec[0])
    wsm=sg.peak_widths(bsm,[list(bsm).index(max(bsm))])[0][0]*(hsm[1]-hsm[0])
    wx=sg.peak_widths(bx,[list(bx).index(max(bx))])[0][0]*(hx[1]-hx[0])
    wf=sg.peak_widths(bf,[list(bf).index(max(bf))])[0][0]*(hf[1]-hf[0])
    
    #indexes of variability F ~ (B)**s_i
    ssm=round(np.log10(10**wsm)/np.log10(10**wvec),2)
    sx=round(np.log10(10**wx)/np.log10(10**wvec),2)
    sf=round(np.log10(10**wf)/np.log10(10**wvec),2)
    print('Indexes of variability:\n SMARTS:',ssm,'\t X:',sx,' \t FERMI:',sf)
    
    plt.figure(56)
    plt.plot((hvec[1::]-pvec),bvec,label=vecnames[0])
    plt.plot((hsm[1::]-psm),bsm,label=vecnames[1])
    plt.plot((hx[1::]-px),bx,label=vecnames[2])
    plt.plot((hf[1::]-pf),bf,label=vecnames[3])
    plt.xlabel('mag',fontsize=14)
    plt.ylabel('N',fontsize=14)
    plt.title(obj+'   PDFs',size=16)
    plt.legend()
    
#    #Histogram
#    plt.figure(57)
#    plt.hist(vect[0]-pB,alpha=0.8,label=vecnames[0])
#    plt.hist(vect[1]-psm,alpha=0.65,label=vecnames[1])
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
    nos1=10**np.mean(np.log10((smarts)))
    nog1=10**np.mean(np.log10((fermi)))
    nox1=10**np.mean(np.log10((xband)))
    noB=10**np.mean(np.log10(B))
    B=np.array(B).astype(np.float)
    plt.plot(np.log10(B[B>indxx*np.mean(B)]),np.log10(smarts[B>indxx*np.mean(B)]/nos1)/np.log10(B[B>indxx*np.mean(B)]/noB),'r.',label='smarts')
    plt.plot(np.log10(B[B>indxx*np.mean(B)]),np.log10(fermi[B>indxx*np.mean(B)]/nog1)/np.log10(B[B>indxx*np.mean(B)]/noB),'b.',label='fermi')
    plt.axis([min(vect[0]),max(vect[0]),-0.5,2.0])
    plt.legend()
    
    if namevar=='B':
        vi=0.3*(range_smarts[0]+range_smarts[1])
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
        plt.ylabel(r'$s_{smarts}=\frac{dlog\,F_{smarts}}{dlog\,B}$',fontsize=14)
        plt.title(r''+obj+'\t ($p='+str(p)+'$)',fontsize=16)
        
        
        plt.plot(np.linspace(min(np.log10(B))-0.5,max(np.log10(B))+0.5),np.ones(50)*ssm,'k-',label='numerical')
        plt.axis([min(np.log10(B)-0.1),max(np.log10(B)+0.1),min(pvar)-0.1,max(pvar)+0.1])
        fig101.legend()
        if SAVE=='on':
            fig101.savefig('ssm_'+obj+imtyp,bbox_inches='tight')
