#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Obs4fit
Created on Tue Sep 17 15:15:49 2019
@author: Markos Polkas
"""

import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import interpolate as inpl
import scipy.stats as stats
import os


"""SETTINGS"""
objects=['3C273','3C279','PKS2155304']
print(objects)
oo= 0 #int(input('Choose index .. from list above ( 0, 1 .. )')) #from above list 0 or 1 or 2 for modelling
obj=objects[oo]

save_code='on' #save code.inp file with parameters
RUN='on' #Run a steady state with the suggested parameters  #!! DEPENDENCE MK95 code executable!!
routerun='./'
Interactive='no' #Makes plot to choose the input parameters of peaks and slopes from diagram
Bibext='on' #for entering the values of bibliography for the external photons component (single BB)
Suggest_Gamma='no' 
Manual_params='no' #Overwrites suggested values with typed-input

"""INPUT"""
#Manual values FLAG -123 if you want the suggested value to be used
pman= 3.31
Bman= 0.04
Rman= 3.89e16
logleman= -4.01
loggminman= 3.81
loggmaxman= 5.7
Tman= 1e-06
loglextman= -12 #np.log10(lextb)-2*np.log10(Gamma)+0.5
deltaman= 30.0
Gammaman= 30.0
htaman= -123
factman=-123
# LISTED PARAMETERS EACH COLUMN IS A VALUE ACQUIRED FOR A DIFFERENCE OBJECT
#Bibliographical values
redsh=[0.158,0.536,0.0]
Dists=[749,3115,543]
tvar=[1*86400,2*86400,200] #sec
Lext=[0.1*2.*10**46,0.01*2.*10**45,'None']
Text=[1.0e3,9e2,'None']
Rext=[4.4e17,2.5e17,'None']
#DIAGRAM PARAMETERS ------> program FINDS code parameters for the fit
#Here you provide the "observables" from the SED plot
aa=[    0.45,    0.45,     0.05] #Uncooled  Slope #0.16 PKS
vbr=[   14.5,     14.6,    15.0] #log Syn peak frequency
vss=[14.5   , 14.6  ,    18.5] #log Syn maximum emitted energy
vcc=[   23.0,   24.0,    29.] #log IC max
vbb=[   12,     12,      12] #log breaking from escape
FF2=[   -9.7,  -9.8,    -9.5] #log peak value of vFv diagram
syncpeaked=[ 'on',  'no',   'on'] #'ys' for higher synchrotron peak/ 'no' otherwise
hh=[    3,     3.7,     0.25] #lext/lb compare the two humps (IC/synchrotron)
    
ff=[    10.,    100.0,     1.0] #Estimate of le/lb if no breaking point given
dd=[    10.,      21.0,     30.] #delta
gmin=[  0.0,  0.0,	 3.61] #log gamma_min
G2d=[   1.5,    0.5,   1.0 ] # Gamma to delta ratio

delta_man=dd[oo]
SSConly='no' #on tries SSC ONLY while no uses an appropriate external photon field
if oo==2:
    SSConly='on'
    Bibext='no' #use bibliographical values for external field
""" End of INPUT"""

"""Plotting Options"""
xoptions=['mec2','Hz','GeV']
yoptions=['le','vLv','vFv','radio','Fv','Xrays','gamma'] #radio -> Jy / Xrays -> counts /sec / Gamma -> counts/hour
yunits={'le':' ','vLv':'erg/s','vFv':'erg/s/cm^2','radio':'Jy',\
       'Fv':'erg/s/cm^2/Hz','Xrays':'counts/s/cm^2','gamma':'counts/hour/cm^2'}

xaxis=xoptions[1]# CAUTION! indexing of the list begins from zero
yaxis=yoptions[2]# CAUTION! indexing of the list begins from zero
Chi2test='on'
save='on'
Multicolor='on' #Plot data points with different color for every interval of 2 years

#delta_man=36 #-1 for None
#Gamma_man=delta_man #-1 for None
delta=delta_man

"""Second BB"""
secondBB='no' #add an additional BB on spectrum (usually DISK) but that is not involved in external photon fields of code
Gamma2=1
factbb=500. #(L_disk-L_BLR)/L_BLR
Tbb2=20000
Rbb2=10**15.62
#
#Bibliograpgical Data for objects
Rbb=[4.47e17,2e17,1.26e16] #Radii of the BB radiation zone
#
#
#Optional Plotting Defaults
nu_cutoff=13 #energy in mec^2 to cutoff plotting 
HEIGHT=1 #height of y-axis (threshold) in logscale measured below the minimum value of flux
FACTOR=1#Shift to the y-axis values (so adding up a constant to the log values)/Default 1 
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""


#"""Constants"""
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
Lsun=3.9*10**33
Ledd=3.2*10**4*Lsun*10**8 #M/10^8 Msun

route= './'+obj+'/'
filename=route+obj+".ascii"
points=ascii.read(filename)
v_pts=np.array([points["col1"]])[0]
vFv_pts=np.array([points["col3"]])[0]
errorv=abs(points["col4"])
errorvFv=abs(points["col4"])

if Interactive=='on':
    plt.figure(6)
    plt.plot(v_pts,vFv_pts,'r.')
    plt.errorbar(v_pts,vFv_pts,errorvFv,xerr=None,fmt='r.',ecolor='black')
    plt.pause(0.05)
    aa[oo]=float(input('\n slope (without cooling) =\t'))
    vss[oo]=float(input('\n log v_synchr,max (Hz) =\t'))
    vcc[oo]=float(input('\n log v_compt,max (Hz)=\t'))
    vbb[oo]=float(input('\n log v_nreak_esc (Hz) (if NON right 10) =\t'))
    FF2[oo]=float(input('\n nu*Flux peak=\t'))
    hh[oo]=10**float(input('\n log hta = log(lext/lb)=\t')) #lext/lb
    ff[oo]=float(input('\n fact (le/lb)=\t'))
    gmin[oo]=float(input('\n =log gmin\t'))
    dd[oo]=float(input('\n delta =\t'))
    SSC=str(input('\ SSC only, with NO external photons (on\no) =\t'))#on tries SSC ONLY while no uses an appropriate external photon field 
    plt.close()
    
a=aa[oo]
zsh=redsh[oo]
v0=10**vss[oo]*(1+zsh)
vs=10**(vbr[oo])*(1+zsh)
vc=10**(vcc[oo])*(1+zsh)
vb=10**(vbb[oo])*(1+zsh)
D=Dists[oo]*10**6*pc
FD2=10**(FF2[oo])*D**2     #carefull for BB Interferencesyncpeaked
if syncpeaked[oo]=='no':
    FD2=FD2/hh[oo]

hta=hh[oo]  #carefull for BB Interference
fact=ff[oo]
delta=dd[oo]
Gamma=G2d[oo]*delta
loggmin=gmin[oo]
Rvar=c*tvar[oo]*delta
if SSConly=='no':
    uextb=Lext[oo]/4/np.pi/c/Rext[oo]**2*Gamma**2
    Textb=Text[oo]*Gamma
    opacity=uextb/aBB/Textb**4
    Bibext='on'

"""Suggested Values for Simulation"""
p=3-2*a #index
KN=0 #Klein -Nishina regime flag

#for simplifying calculations
normgbr=3/4/(1/8/np.pi/me/c**2*sigmaT)/(1+hta)
normR=(1.2*10**18*(FD2/6/10**44)**0.5*delta**-2*(normgbr/10**loggmin)**(p/2-1))**(2/p)

if SSConly=='on':
    vs=vs/delta
    loggmax=0.5*np.log10(0.75*vc/v0)
    B=v0/(3*10**6*(10**loggmax)**2*delta)
    R=1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*B**-1*delta**-2
    T=0.000001
    lext=10**-15
    if a<0.5:
            R=(1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*(3e6*normgbr**2/vs)**(-1/3)*delta**-2)**(3)
            B=(3e6*normgbr**2/vs/(R**2))**(1/3)
            loggmax=0.25*np.log10(vc/3/10**6/B/delta)
            gbr=normgbr*B**-2*R**-1
            fact=((gbr**(2-p)-10**(loggmax*(2-p)))+(p-2)/(3-p)/gbr*(gbr**(3-p)-10**(loggmin*(3-p))))/\
                  (10**(loggmin*(2-p))-10**(loggmax*(2-p)))
    if loggmax>4.53-0.333*np.log10(B)>1:
        KN=1
        loggmax=np.log10(0.5*np.sqrt(3)*h*vc/delta/me/c**2*(1+np.sqrt(1+(me*c**2/h*delta)**2/vc/v0)))
        B=vs/(3*10**6*(10**loggmax)**2*delta)
        R=1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*B**-1*delta**-2
        if a<0.5:
            R=(1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*(3e6*normgbr**2/vs)**(-1/3)*delta**-2)**(3)
            B=(3e6*normgbr**2/vs/(R**2))**(1/3)
            loggmax=0.25*np.log10(vc/3/10**6/B/delta)
            gbr=normgbr*B**-2*R**-1
            fact=((gbr**(2-p)-10**(loggmax*(2-p)))+(p-2)/(3-p)/gbr*(gbr**(3-p)-10**(loggmin*(3-p))))/\
                  (10**(loggmin*(2-p))-10**(loggmax*(2-p)))
else:
    T=(kb/(10**6*h*(6*np.pi*aBB)**0.5)) *hta**0.5*v0/vc
    B=(8*np.pi*opacity*aBB/3)**0.5 *hta**-0.5*T**2
    v0=2.7*kb*T*delta/h
    loggmax=np.log10(0.5*np.sqrt(3)*h*vc/delta/me/c**2*(1+np.sqrt(1+(me*c**2/h*delta)**2/vc/v0)))
    R=1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*B**-1*delta**-2
    lext=sigmaT*R*(opacity*aBB*T**4)/me/c**2
    if loggmax>5.4+np.log10(20000/T):
        KN=1
        loggmax=np.log10(0.5*np.sqrt(3)*h*vc/delta/me/c**2*(1+np.sqrt(1+(me*c**2/h*delta)**2/vc/v0)))
        B=vs/(3*10**6*(10**loggmax)**2*delta)
        T=np.sqrt(B/(8*np.pi*opacity*aBB/3)**0.5 *hta**0.5)
        R=1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*B**-1*delta**-2
        lext=sigmaT*R*(opacity*aBB*T**4)/me/c**2

    if Bibext=='on':
        T=Textb
        v0=2.7*kb*T*delta/h
        loggmax=np.log10(0.5*np.sqrt(3)*h*vc/delta/me/c**2*(1+np.sqrt(1+(me*c**2/h*delta)**2/vc/v0)))
        B=vs/delta/3/10**6/10**(2*loggmax)
        R=1.2*10**18*(FD2/6/10**44)**0.5**hta**-0.5*B**-1*delta**-2
        lext=uextb/me/c**2*sigmaT*R
        if loggmax>5.4+np.log10(20000/T):
            KN=1
            loggmax=np.log10(0.5*np.sqrt(3)*h*vc/delta/me/c**2*(1+np.sqrt(1+(me*c**2/h*delta)**2/vc/v0)))
            B=vs/(3*10**6*(10**loggmax)**2*delta)
            R=1.2*10**18*(FD2/6/10**44)**0.5*hta**-0.5*B**-1*delta**-2
            lext=sigmaT*R*(opacity*aBB*T**4)/me/c**2  

le=3*sigmaT*FD2/delta**4/R/me/c**3
lb=B**2/8/np.pi/me/c**2*sigmaT*R
if a<0.5:
    le=(hta*lb)/fact #fact matches lgamma^sync/le

logle=np.log10(le)
loglb=np.log10(lb)
loglext=np.log10(lext)
#Escape
tesc=19924618*delta**0.5*B**-1.5*(R)**-1*vb**-0.5
#Coolings
gbreakcoolsyn=3/4/lb
tcoolsyn1=3/4/lb/10**loggmin
vbreakcsyn=3*10**6*B*gbreakcoolsyn**2*delta/(1+zsh)
gbreakcoolc=3/4/lext
tcoolext1=3/4/lext/10**loggmin
vbreakcc=2/h*kb*T*gbreakcoolc**2*delta/(1+zsh)
#percentage of coolage
percsyn=round( (1-(np.log10(gbreakcoolsyn)-loggmin)/(loggmax-loggmin))*100)
percic=round( (1-(np.log10(gbreakcoolc)-loggmin)/(loggmax-loggmin))*100)

gbreak=0.5*np.sqrt(gbreakcoolsyn*gbreakcoolc)
tbreak=1/(1/tcoolext1+1/tcoolsyn1)

if SSConly=='on':
    gbreak=gbreakcoolsyn
    tbreak=tcoolsyn1
    mingbreak=gbreak

#Corrections
if np.log10(vc/vs)<2*loggmax and SSConly=='no':
    loggmin=loggmax-0.5*np.log10(vc/vs)
mingbreak=min(gbreakcoolc,gbreakcoolsyn)
#if p>2:
#    logle=logle+(p-2)*(np.log10(mingbreak)-loggmin)     
#else:
#    logle=logle-(p-2)*(np.log10(mingbreak)-loggmin)
 
#le=10**logle

vssa=(((p+1)/32/np.pi**2)**2*(3*fortio/2/np.pi/me**3/c)*delta**(-1)*vs**(p-3)*\
(4*np.pi*FD2)**2*B/R**4)**(1/(p+4))/(1+zsh)

if Manual_params!='on':
    print(' delta=',delta,'\n Gamma=',Gamma,'\n hta (lext/lb)=',hta,'\n fact (le/lb)= ',fact,'\n\n index (p)=',p,'\n log gmin= ',loggmin,\
      '\n log gmin=',loggmin,'\n log gmax=',loggmax,'\t Klein Nishina Regime:',KN*'YES','\n B=',B,'\n R=',R,'   (logR= ',round(np.log10(R),1),')','   (logRvar= ',round(np.log10(Rvar),1),')'\
      '\n log le=',round(logle,2),'\n log  lb=',round(np.log10(lb),2),'\n log lext=',round(loglext,2),'\n T=',T,'\n tesc=',tesc,\
      '\n besc=',1/tesc,'\n t_cool,syn(1)=',tcoolsyn1,'\t',int(1-tcoolsyn1//1)*'COOLED',\
      '\n t_cool,ext(1)=',tcoolext1,'\t',int(1-tcoolext1//1)*'COOLED','\n log gamma_break,cool,syn=',\
      np.log10(gbreakcoolsyn),'\t',percsyn,'%','\n log v_break,cool,syn=',np.log10(vbreakcsyn),'\n log gamma_break,cool,ext=',\
      np.log10(gbreakcoolc),'\t',percic,'%','\n log v_break,cool,ext=',np.log10(vbreakcc),'\n log v_ssa=',np.log10(vssa),'\n le/lb/fact= ',le/lb/fact)

""" MANUAL CHANGE / !!!OVERWRITES THE DISPLAYED VALUES !!!"""
#Manual_params='on'
if Manual_params=='on':
    params4man=['p','B','R','logle','loggmin','loggmax','T','loglext','delta','Gamma','hta','fact']
    for par in params4man:
        if globals()[par+'man']!=-123:
            globals()[par]=globals()[par+'man']
    
    le=10**logle
    lb=B**2/8/np.pi/me/c**2*sigmaT*R
    loglb=np.log10(lb)
    if Bibext=='on':
        lext=uextb/me/c**2*sigmaT*R
    #Escape
    tesc=19924618*delta**0.5*B**-1.5*(R)**-1*vb**-0.5
    #Coolings
    gbreakcoolsyn=3/4/lb
    tcoolsyn1=3/4/lb/10**loggmin
    vbreakcsyn=3*10**6*B*gbreakcoolsyn**2*delta/(1+zsh)
    gbreakcoolc=3/4/lext
    tcoolext1=3/4/lext/10**loggmin
    vbreakcc=2/h*kb*T*gbreakcoolc**2*delta/(1+zsh)
    #percentage of coolage
    percsyn=round( (1-(np.log10(gbreakcoolsyn)-loggmin)/(loggmax-loggmin))*100)
    percic=round( (1-(np.log10(gbreakcoolc)-loggmin)/(loggmax-loggmin))*100)

    print('\n MANUAL OUTPUT\n delta=',delta,'\n Gamma=',Gamma,'\n hta (lext/lb)=',hta,'\n fact (le/lb)= ',fact,'\n\n index (p)=',p,\
      '\n log gmin=',loggmin,'\n log gmax=',loggmax,'\t Klein Nishina Regime:',KN*'YES','\n B=',B,'\n R=',R,'   (logR= ',round(np.log10(R),1),')', '   (logRvar= ',round(np.log10(Rvar),1),')'\
      '\n log le=',round(logle,2),'\n log  lb=',round(np.log10(lb),2),'\n log lext=',round(loglext,2),'\n T=',T,'\ntesc=',tesc,\
      '\n besc=',1/tesc,'\n t_cool,syn(1)=',tcoolsyn1,'\t',int(1-tcoolsyn1//1)*'COOLED',\
      '\n t_cool,ext(1)=',tcoolext1,'\t',int(1-tcoolext1//1)*'COOLED','\n log gamma_break,cool,syn=',\
      np.log10(gbreakcoolsyn),'\t',percsyn,'%','\n log v_break,cool,syn=',np.log10(vbreakcsyn),'\n log gamma_break,cool,ext=',\
      np.log10(gbreakcoolc),'\t',percic,'%','\n log v_break,cool,ext=',np.log10(vbreakcc),'\n log v_ssa=',np.log10(vssa),'\n le/lb/fact= ',le/lb/fact)   

 
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""


"""FIND SUGGESTED GAMMA FOR BB DEBOOSTING"""
if Suggest_Gamma=='on':
    theta=np.logspace(-2-round(np.log10(delta)),1-round(np.log10(delta)),500)
    breakgam=0
    while  breakgam==0 and len(theta)>1:
        gam=(1-np.sin(theta)*np.sqrt(1+delta**2*np.cos(theta)**2))/delta/np.sin(theta)**2
        minGam=min(gam[gam>0])
        beta=np.sqrt(1-1/minGam**2)
        thetamin=np.arccos((1-1/(delta*minGam))/beta)
        if str(thetamin)=='nan' or minGam<0:
            TH=list(theta)
            TH.pop(list(gam).index(minGam))
            theta=np.array(TH)
            GMM=list(gam)
            GMM.pop(list(gam).index(minGam))
            gam=np.array(GMM)
        else:breakgam=1
        #minGam=1
    print("\n\tMinimum value of Bulk Lorentz Factor, given random ",\
          "values of angle thetad\n theta=",thetamin,"\n Gamma=",minGam,\
          "\n Degradation of BB (logs)(delta^4+Gamma^2)): ",round(4*np.log10(delta),2),
          " + ",round(2*np.log10(minGam),2)," = ",round(np.log10(delta**4*minGam**2),2),"\t Must Be ~le-lext=",round(logle-loglext,2),\
          "\n You can print dictionary 'thetaGamma',for a different set of 'theta-Gamma' parameters\n\n")
    thetaGamma=dict(zip(theta,gam))
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""


"""SAVING OF PARAMETERS TO FILES"""
if save_code=='on':
    i=0
    LL=[]
    cn4=open('code.inp','r')
    for line in cn4:
        if i==2:
            LL.append(str(round(R,int(-np.log10(R)+2)))+'\t'+str(round(B,int(-np.log10(B)+2)))+'\t1.\t1.\n')
        elif i==4:
            if tesc<1:
                tesc=1.0
            LL.append('1\t'+str(round(loggmin,3))+'\t'+str(round(loggmax,3))+'\t'+str(round(p+0.01,3))+'\t'\
                  +str(round(logle,3))+'\t'+str(round(1/tesc,3))+' 0. 0\n')
        elif i==5:
            LL.append('1\t'+str(round(T,int(-np.log10(T)+2)))+'\t'+str(round(10**loglext,int(-loglext+2)))+'\n')
            lines=cn4.readlines()
            LL=LL+lines
            break
        else:
            LL.append(line)
        i=i+1
        
    cn4.close()
    newcn4=open('code.inp','w')
    newcn4.writelines(["%s" % item  for item in LL])
    newcn4.close()

if RUN=='on':
    os.system(routerun+'code_clean')
#%%
"""
%
%
%
%
%
%
%
"""

obj=objects[oo]
data_name=obj+".ascii"
dtpts=[]
if Multicolor=='on':
   dtpts.append(obj+'_rest.ascii')
   dtpts.append(obj+'_08-10.ascii')
   dtpts.append(obj+'_10-12.ascii')
   dtpts.append(obj+'_12-14.ascii')
   dtpts.append(obj+'_14-16.ascii')
   #dtpts.append(obj+'_16-18.ascii')
legend_pts=['<2008','2008-2010','2010-2012','2012-2014','2014-2018']
"""Observational Data"""
points=ascii.read(route+data_name)
Cunitx=0 #4 Hz 
Cunity=0
v_pts=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
vFv_pts=np.array([points["col3"]])[0]-np.ones(1)*Cunity
errorv=abs(points["col4"])
errorvFv=abs(points["col4"])
if Multicolor=='on':
    for i in range(len(dtpts)):
        points=ascii.read(route+dtpts[i])
        globals()['v_pts'+str(i)]=(np.array([points["col1"]])[0]-np.ones(1)*Cunitx)
        globals()['vFv_pts'+str(i)]=(([points["col3"]])[0]-np.ones(1)*Cunity)
        globals()['errorv'+str(i)]=(abs(points["col4"])[0])
        globals()['errorvFv'+str(i)]=(abs(points["col4"])[0])
        if 'rest' in dtpts[i]: 
            globals()['v_pts'+str(i)] = np.log10(globals()['v_pts'+str(i)])
            globals()['vFv_pts'+str(i)] = np.log10(globals()['vFv_pts'+str(i)])
            globals()['errorvFv'+str(i)] = 0*np.log10(globals()['errorvFv'+str(i)])


func_obs=inpl.interp1d(v_pts,vFv_pts)
#IF MULTIPLE TIMESTEPS DENOTE LENGTH
Names=['steady']

for n in Names:
    """Reading Files"""
    n5=n+'.85'
    file85=open(n5,'r')
    f85=file85.readlines()
    file85.close()
    eps85=[]
    I85=[]

    LEN=int(len(f85)/3-1)
    for y in f85:
        x=y.split()
        z=x[0].strip()
        w=x[1].replace('\n','').strip()
        if w=='' or w=='NAN':
            z='-100'
            w='-100'
        eps85.append(z)
        I85.append(w)
    eps85=np.array(eps85[-LEN::]).astype(np.float)
    I85=np.array(I85[-LEN::]).astype(np.float)
    


    """REMOVE UNWANTED POINTS"""    
    """Remove Extremely High Energy Gamma Rays/ Right Cut-off of Frequencies""" 
    zz=np.where(eps85>nu_cutoff,eps85,0)
    zzz=zz.nonzero()
    for i in zzz:
        eps85[i]=-100
        I85[i]=-100
    
    """Remove Extremely High Outlier Points"""
    THRESHOLD=10 #orders of magnitude of allowed difference of consecutive points 
    eps85=eps85[I85>max(I85)-THRESHOLD]
    I85=I85[I85>max(I85)-THRESHOLD]
    I85=I85+np.log10(FACTOR)
        
    
    """BLACK BODY ADDITION"""
    eps85bb=eps85
    Ibb=eps85
    if SSConly=='no':
        if len(eps85bb)<200:
            eps85bb=np.linspace(min(eps85),max(eps85),200)
            Rext=Rbb[oo]
        def BB(x,T,lext,Rbb):
            Ivbb=3*x+np.log10((2*me**3*c**4/h**2)/(np.exp(10**x*me*c**2/kb/T)-1)+10**-100)
            opacity=lext*me*c**2/sigmaT/R/aBB/T**4 #a factor occuring by the fact that opacity must be computed in the jets frame of reference
            #make it dimensionelss / scaling
            Ibb=x+np.log10(me*c**2/h)+Ivbb+np.log10(opacity*4*np.pi*Rbb**2)+np.log10(sigmaT/(4*np.pi*R*me*c**3))
            return Ibb
        Ibb=BB(eps85bb,T/Gamma,lext/Gamma**2,Rext)+np.log10(Gamma) #due to \nu boosting that wasn't used while computing nu *Fnu
        if secondBB=='on':
            lext2=1/Gamma**2*factbb*(Rext/Rbb2)**2*lext
            Ibb2=BB(eps85bb,Tbb2,lext2,Rbb2)+np.log10(Gamma)
            xbb2=eps85bb+np.log10((me*c**2/h))
            ybb2=np.log10(10**Ibb2/D**2*R*me*c**3/sigmaT)
    
    """TRANSFORMATION to the F.O.F of the Observer"""
    if xaxis=='mec2':
        x85=eps85+np.log10(delta/(1+zsh))
        xbb=eps85bb
    if xaxis=='Hz':
        x85=eps85+np.log10(delta*(me*c**2/h)/(1+zsh))
        xbb=eps85bb+np.log10((me*c**2/h))
    if xaxis=='GeV':
        x85=np.log10(delta*10**eps85/eV/10**9*me*c**2/(1+zsh))
        xbb=np.log10(10**eps85bb/eV/10**9*me*c**2)
        #
        #   
        #
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
    if SSConly=='no':
        x85nn=x85[y85>(max(y85)-12)] #without noise and unwanted points for interpolation to work
        y85nn=y85[y85>(max(y85)-12)]
        funcy85=inpl.interp1d(x85nn,y85nn)
        #x81tbb=np.concatenate([x85,xbb]) #overpose NOT addition
        #y81tbb=np.concatenate([y85,ybb]) #overpose NOT addition
        x81tbb=xbb[xbb<max(x85nn)][xbb[xbb<max(x85nn)]>min(x85nn)]
        y81tbb=np.log10(10**ybb[xbb<max(x85nn)][xbb[xbb<max(x85nn)]>min(x85nn)]+10**funcy85(x81tbb))
        funcy81tbb=inpl.interp1d(x81tbb,y81tbb)
        if secondBB=='on':
                funcybb2=inpl.interp1d(xbb2,ybb2)
                y81tbb=np.log10(10**funcybb2(x81tbb)+10**y81tbb)
        tbb=open('fort.81tbb','w')
        for i,j in zip(x81tbb,y81tbb):
                tbb.write(str(i)+'\t'+str(j)+'\n')
        tbb.close()
    else:
            x81tbb=x85
            y81tbb=y85
            funcy81tbb=inpl.interp1d(x81tbb,y81tbb)
    """PLOTTING"""
    xlab='v'
    legl=[r'SED, $\delta=$'+str(delta)]
    if SSConly=='no':
        legl.append('Black Body')
        if secondBB=='on':
            legl.append('Disk')
    if xaxis=='mec2':
        xlab='E'
    if Multicolor=='on':
        colors =['grey','tomato','gold','limegreen','purple']
        form=['.','o','h','s','+']
           
    if SSConly=='no':
            fig2, ax2 = plt.subplots(num=4)
            plt.xlabel(r'$'+xlab+'\;\;'+' ['+xaxis+']$',fontsize=15)
            plt.ylabel(r'$'+yaxis+'\;\;['+yunits[yaxis]+']$',fontsize=15) 
            if Multicolor=='on':
                for i in range(len(dtpts)):
                    alpha=1.0
                    if colors[i]=='grey': alpha =0.3
                    ax2.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                                 globals()['errorv'+str(i)],globals()['errorvFv'+str(i)],marker=form[i], color = colors[i],ms=3.5, linewidth=0.0,alpha=alpha)
            else:
                    ax2.errorbar(v_pts,vFv_pts,errorv,errorvFv,fmt='r.',ecolor='red',ms=3.5)
            leg1=ax2.legend(legend_pts)
            ax2.legend(legend_pts)
            ax2.plot(x85,y85,'b-')
            plt.title(obj+'  SED 85',fontsize=18)
            lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),round(min(vFv_pts)-HEIGHT),int(max(max(vFv_pts),max(y85))+0.1)]
            ax2.axis(lims)
            leg2=ax2.legend([legl[0]])
            ax2.add_artist(leg1)
            ax2.xaxis.set_minor_locator(plt.LinearLocator(lims[1]-lims[0]+1))
            ax2.yaxis.set_minor_locator(plt.LinearLocator(10*(lims[3]-lims[2])+1))
        
    #FIGURE FORT85+CORRECTED BLACKBODY (a.k.a. CORRECTED 81)
    fig3, ax3 = plt.subplots(num=5)
    plt.xlabel(r'$'+xlab+'\;\;'+' ['+xaxis+']$',fontsize=15)
    plt.ylabel(r'$'+yaxis+'\;\;['+yunits[yaxis]+']$',fontsize=15) 
    if Multicolor=='on':
        for i in range(len(dtpts)):
            alpha=1.0
            if colors[i]=='grey': alpha =0.3
            ax3.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                         globals()['errorv'+str(i)],globals()['errorvFv'+str(i)],marker=form[i], color = colors[i],ms=3.5, linewidth=0.0,alpha=alpha)
    else:
            ax3.errorbar(v_pts,vFv_pts,errorv,errorvFv,fmt='r.',ecolor='red',ms=3.5)
    leg1=ax3.legend(legend_pts,loc='upper right',bbox_to_anchor=(0.55, 0.05, 0.35, 0.35))
    ax3.plot(x81tbb,y81tbb,'b-')
    if SSConly=='no':
        ax3.plot(xbb,ybb,'b--',linewidth=0.5)
    if secondBB=='on':
        ax3.plot(xbb2,ybb2,'b-.',linewidth=0.5)
    
    plt.title(obj,fontsize=18)
    lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),round(min(vFv_pts)-HEIGHT),int((max(max(vFv_pts),max(y81tbb)))+0.1)]
    ax3.axis(lims)
    leg2=ax3.legend(legl)
    ax3.add_artist(leg1)        
    ax3.xaxis.set_minor_locator(plt.LinearLocator(lims[1]-lims[0]+1))
    ax3.yaxis.set_minor_locator(plt.LinearLocator(10*(lims[3]-lims[2])+1))
    plt.show()
       

    if save=='on':
        if SSConly!='on':
            plt.figure(4).savefig(yaxis+'_'+n5+'.ps')
            plt.figure(4).savefig(yaxis+'_'+n5+'.png')
        plt.figure(5).savefig(yaxis+'_fort.81tbb.ps')
        plt.figure(5).savefig(yaxis+'_fort.81tbb.png')
        import datetime
        dth=str(datetime.datetime.now()).split(':')[0:2]
        savename=obj+'_'+dth[0].split('-')[2].replace(' ','-')+'-'+dth[1]
        copyfile=open('run_copy.sh','w')
        copyfile.write('#'+savename+'\n')
        copyfile.write('cp vFv_fort.81tbb.ps ./Steady_States/'+savename+'_81tbb.ps \n')
        copyfile.write('cp code.inp ./Steady_States/'+'code_'+savename+'.inp \n')
        copyfile.close()
        print('Graphs Saved in Folder')

"""Chi^2 Test"""
if Chi2test=='on':
        v_pts.sort()
        x81tbb.sort()
        nucomp=v_pts[v_pts>max(v_pts[0],x81tbb[0])]
        nucomp=nucomp[nucomp<min(v_pts[-1],x81tbb[-1])]
        num_params=7
        DF= len(nucomp)-num_params #Degrees of freedom
        if SSConly=='on': #decrease for absence of T, lext
            DF=DF+2
        if float(lines[4].split()[1])!=0.01: #if gmin different than unity
            DF=DF-1
        crit=stats.chi2.ppf(q=0.16, df=DF) #minimum acceptable value
        #Calc chi2 and p-value
        expected=10**func_obs(nucomp)
        observed=10**funcy81tbb(nucomp)
        chi2s=(((observed/observed.sum()-expected/expected.sum())**2)/(expected/expected.sum())).sum()
#        rchi2_single=round(chi2s/crit,int(2-np.log10(chi2s/crit)))
#        #pvalue_single[l]= 1 - stats.chi2.cdf(x=chi2,df=DF)
#        print('One-zone, leptonic model fit chi^2 divided by the criterion value of 1-sigma confidence: :\n')
#        print('\t x^2/χ^2(5sigma)=\t{} '.format(rchi2_single))
        print('\n x^2=\t {} (value for real flux, not logscale)'.format(chi2s))
    
