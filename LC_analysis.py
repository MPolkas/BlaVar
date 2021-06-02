#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 18:52:18 2020

@author: Markos Polkas

"""
from matplotlib import use
use('Qt5Agg')

from DELCgen import *
import scipy.stats as st
import csv
import numpy as np
from scipy import interpolate as inpl
import matplotlib.pyplot as plt 
#from astropy.io import fits
from astropy.table import Table
import os
import datetime

obj=os.getcwd().split('/')[-1]
objects=['3C273','3C279','PKS2155304']
if obj not in objects:
    print(objects)
    oo=int(input("\n Give Number of object [0 , 1 ..] from following list:\t"))
    obj=objects[oo]
oo=objects.index(obj)
route = os.getcwd()+'/'
SAVE='on' #save plots and time curve info

""" READ INPUT FILE """
cn4=open(route+'fkTC.inp','r')
lines=cn4.readlines()
cn4.close()
#Settings for usage of created fake lightcurve
namevar=str(lines[9].split(' ')[0]) #variability on these parameter
lendays=int(str(lines[9].split(' ')[1])) #days of observations to recreate
NEW_delc=str(lines[9].split(' ')[2]) #create new fake timecurve with 'on'
POWER=float(lines[9].replace('\n','').split(' ')[3]) #sqrt of final fake lightcurve

#multi option (not for delta), not for newly created LC
MULTI='no'
remt=lendays
if 'multi' in NEW_delc:
    remt=float(NEW_delc.split('multi')[1])
    MULTI='on'
    NEW_delc='no'


"""EMMANOUELOPOULOS FAKE TIME CURVE CREATION SETTINGS"""
CSVfile='no'
FITSfile='on'
plotLC='no'
#PSD option
oneindex_model='no'
no_background = 'on'

tbin=1
initfile ='lc_daily_'+obj+'.txt'  #test file: initfile = "pks1510-089_gamma_lc.txt"

#redshifts and bibliography
Dist=[749,3115,543,1937][oo] #Mpc Luminosity Distances
zsh=[0.158,0.536,0.0,0.36][oo] #redshifts
#units
eV=1.602192e-12 #erg
MeV=10**6*eV #erg
range_fermi=[0.1*10**9*(2.4184*10**14),300*10**9*(2.4184*10**14)] #eV/h = 2.41e14
units=1. #set to one if default units/ needed for correct normalization  check IF line 

init=-1#index of initial point from fake LC delc1.time[ ] to segment //  set -1 for segment around maximum
widthgamma=1.5 #for 'gmax' namevar / forced width of variation in orders of magnitude

#Data treatment for NO MASK METHOD
DATAcut='no' #cuts point below the threshold flux
DATAfill='no' #interpolates between cut points
DATAreplace='no' #replaces with threshold value , fluxes below threshold
thresflx=1e-13*units # ph/cm^2/s Check for every specific object

#Interpolation
imethod='quadratic'
nint=4 #number of points between t_cross intervals 
##


#VARIABLES
tend=float(lines[0].split()[-1].replace('\n',''))
nsteps=float(lines[0].split()[2])
R=float(lines[2].split()[0])
B=float(lines[2].split()[1])
loggmin=loggmax=float(lines[4].split()[1])
loggmax=float(lines[4].split()[2])
p=float(lines[4].split()[3])
logle=float(lines[4].split()[4])
besc=float(lines[4].split()[5])
tesc=1/besc
delta=float(lines[7].split()[0])
Gamma=float(lines[8].split()[0])
T=float(lines[5].split()[1])
lext=float(lines[5].split()[2])

#
le=10**logle
gmin=10**loggmin
gmax=10**loggmax

#calculations
norm=globals()[namevar]
lenfk=lendays+1 #for fine array creation
dilate=1-((((1+zsh)**2-1)/((1+zsh)**2+1))**2)**2 #dilation factor, https://www.scirp.org/journal/paperinformation.aspx?paperid=65598
tcross=R/3/10**10/3600./24. #t_cross in days
tcr=tcross*dilate/delta

eflux, eflux_err, MJD, flux, flux_err = [[],]*5


if CSVfile=='on':
    csvdata = np.genfromtxt(route+initfile, delimiter=',',dtype=float)
    for row in csvdata:
        MJD.append(row[0])
        eflux.append(row[1])
        eflux_err.append(row[2])
    #use real dates
    np.savetxt(route+'myfile.txt', np.c_[MJD,eflux,eflux_err])
    t_dummy = [0] * len(t1)
    t_dummy = [i for i in xrange(len(t1))]
    #save LC with some dummy delta step of 1 day 
    np.savetxt(route+'myfile2.txt', np.c_[t_dummy,eflux,eflux_err])
    datfile='myfile.txt'
    
elif FITSfile=='on':
#       ##OLD METHOD
#        datfile=fits.open(route+'lc_daily_'+obj+'.fits')
#        tbdata=datfile[1].data
#        for tb in tbdata:
#            MJD.append(tb[0])
#            flux.append(tb[7]) #photon flux
#            flux_err.append(tb[6])
#            eflux.append(tb[11])  #energy flux
#            eflux_err.append(tb[12])
#        MJD=np.array(MJD)
#        eflux=np.array(eflux)*MeV
#        flux=np.array(flux)
#        norma=np.mean(eflux*MeV)
#        eflux_err=np.array(eflux_err)*MeV
#        flux_err=np.array(flux_err)
#        np.savetxt(route+'myfile.txt', np.c_[MJD,eflux,eflux_err])
#        datfile='myfile.txt'
        t = Table.read(route+'lc_daily_'+obj+'.fits')
        # mask to get fluxes
        m = (t['ts'] > 4)  & (t['eflux']>t['eflux_err'])
        #m = (t['ts']>-100) #keep all
        
        # central bin times for times with detection
        t_cen = np.array(0.5 * (t['tmin'] + t['tmax'])[m])  # in MJD
        
        # flux and error for times with detection
        eflux = np.array(t['eflux'][m]*eV*10**6)   # in photons / cm^2 / s
        eflux_err = np.array(t['eflux_err'][m]*eV*10**6)   # in photons / cm^2 / s
        #save edited file as .txt
        np.savetxt(route+'myfile.txt', np.c_[t_cen, eflux, eflux_err])
        datfile='myfile.txt'
        
        #upper limits 
        t_cen_ul = np.array(0.5 * (t['tmin'] + t['tmax'])[~m])  # in MJD
        eflux_ul = np.array(t['eflux'][~m]*eV*10**6)   # in photons / cm^2 / s
        #eflux_err_ul = np.array(t['eflux_err'][~m]*eV*10**6)   # in photons / cm^2 / s
        UL = zip(t_cen_ul, eflux_ul)
else:
     datfile=initfile
     #change all units to erg/cm^2/s , if erg/cm^2/h use 1/3600
#%%
datalc = Load_Lightcurve(route+datfile,tbin)
norma = np.mean(datalc.flux)*units
#normalize to the mean flux for the Fit_PDF to work
datalc.errors=datalc.errors/np.mean(datalc.flux)
datalc.flux = datalc.flux/np.mean(datalc.flux)
thresflx= thresflx/norma


if DATAreplace=='on':
    for dp in range(0,len(datalc.time)):
        if datalc.flux[dp]<thresflx:
            datalc.flux[dp]=thresflx
elif DATAcut=='on':
    cp=len(datalc.time[datalc.flux<thresflx])
    print str(cp)+' points are below the threshold, they are deleted as fake\n'
    datalc.time=np.linspace(datalc.time[0],datalc.time[0]+tbin*len(datalc.time[datalc.flux>thresflx]), len(datalc.time[datalc.flux>thresflx]))
    datalc.errors=datalc.errors[datalc.flux>thresflx]
    datalc.flux = datalc.flux[datalc.flux>thresflx]
elif DATAfill=='on':
        fluxmodel=inpl.interp1d(datalc.time[datalc.flux>thresflx],datalc.flux[datalc.flux>thresflx])
        errormodel=inpl.interp1d(datalc.time[datalc.flux>thresflx],datalc.errors[datalc.flux>thresflx])
        tnmin=min(datalc.time[datalc.flux>thresflx])
        tnmax=max(datalc.time[datalc.flux>thresflx])
        tntime=datalc.time[datalc.time>tnmin]
        datalc.time=tntime[tntime<tnmax]
        datalc.flux = fluxmodel(datalc.time)
        datalc.errors= errormodel(datalc.time)
        

np.savetxt(route+'myfile_new.txt', np.c_[datalc.time, datalc.flux,datalc.errors])
newfile='myfile_new.txt'
datalc = Load_Lightcurve(route+newfile,tbin)
#%%   
if NEW_delc=='on': 
        # plot the light curve
        if plotLC=='on':
            datalc.Plot_Lightcurve()
            
        # Bending power law params
        datalc.psdFit=None
        no_background =  'on'
        A,v_bend,a_low,a_high,c = 0.93, 4e-2, 0.21, 0.23, 0.00
        #Simple power law params/  psd(v)=cc0+cc*v^mm
        oneindex_model = 'no'
        mm,cc,cc0=-1, 1, 0.1
        
        if datalc.psdFit==None:
            if oneindex_model=='on':
                def func_powerlaw(v, m, c, c0):
                    return c0 + v**m * c
                datalc.Fit_PSD(initial_params= [mm,cc,cc0],model=func_powerlaw)
                print('Used single power law with background')
            if no_background =='on':
                def func_doublepowerlaw(v, A, v_bend, a_low, a_high):
                        numer = v**-a_low
                        denom = 1 + (v/v_bend)**(a_high-a_low)
                        return A * (numer/denom)
                datalc.Fit_PSD(initial_params= [A,v_bend,a_low,a_high], model=func_doublepowerlaw)
                print('Used double power law with no background')
            else:
                datalc.Fit_PSD(initial_params=[A,v_bend,a_low,a_high,c])
                print('Used default BendingPL expression')
#%%                
        # create mixture distribution to fit to PDF
        plt.close()
        datalc.pdfFit=None
        
        
        #new parameters of mixture model sub-models
        # Probability density function params
        
                
        #model1
        model1 =  st.lognorm
        skew = 0.0
        par1 , par2 = 0.51 , 0.18
        frozen1 = [[2],[skew]] #skew is position 2 for lognorm
        
        #model2
        model2=  st.lognorm
        cen, skew = 0.33, 0.0
        par3,par4 = 4.05, 0.68
        #array position of variables    array value of variables frozen
        frozen2 = [ [1,2], [cen,skew]]
        
        weight=  0.95 
        
        if model2!=None:
            mix_model = Mixture_Dist([model1,model2],[3,3], [frozen1, frozen2])
            inp_pars = [par1,par2,par3,par4,weight]
        else:
            mix_model  = Mixture_Dist([model1],[3], [frozen1])
            inp_pars  = [par1, par2]

        # fit the PDF of the light curve using the default mixture model
        if datalc.pdfFit==None:
            datalc.Fit_PDF(initial_params=inp_pars, model=mix_model)
        # simulate artificial light curve with Emmanoulopoulos method, using the PSD
        # and PDF of the data light curve, with default parameters (bending power law)
        # for PSD and mixture distribution of gamma and lognormal distribution for PDF)
        if datalc.pdfFit and datalc.psdFit:
            # 1st simulated LC
            delc1= Simulate_DE_Lightcurve(datalc.psdModel, datalc.psdFit.x, datalc.pdfModel, datalc.pdfFit.x, tbin=1, LClength=lendays)
            #delc1.time = delc1.time[delc1.flux<1000] 
            #delc1.flux = delc1.flux[delc1.flux<1000]
            delc1.errors = 0.0*delc1.flux   
            delc1.pdfFit = datalc.pdfFit
            delc1.psdFit = datalc.psdFit
            delc1.psdModel = datalc.psdModel
            delc1.pdfModel = datalc.pdfModel
            Comparison_Plots([datalc,delc1],names=["OBSERVATIONS","SYNTHETIC"],bins=50, norma= (norma/1e-9), add_logLC=True, add_logPDF=False) #Comparison_Plots([datalc,delc1,delc2],names=["Data LC","fake 1","fake 2"],bins=25)
            plt.show()
            plt.pause(1)
 #%%
            SATISFIED = 4 #float( input('Satisfied (Yes: >0 / No <=0):\t'))
            if SATISFIED>0:
                # Plot Compare Diagram and save it
                dttm = datetime.datetime.now()
                datenum=str(dttm.date()).split('-',1)[1]+'_'+str(dttm.time()).split('.')[0]
                datenum = datenum.replace(':','_')
                afile=open(route+'info_DELCfit_'+datenum+'.txt','w')
                afile.write('##PDF fit##')
                afile.write('\n'+str(type(model1)).split('.')[-1])
                afile.write('\n'+str(type(model2)).split('.')[-1])
                afile.write('\nFrozen 1:\t'+str(frozen1)) 
                afile.write('\nFrozen 2:\t'+str(frozen2)) 
                afile.write('\nFitted_parameters:\t '+str(datalc.pdfFit.x))
                afile.write('\n\n ##PSD fit##')
                if oneindex_model=='on':
                    afile.write('\n One Index Power Law +Bacground used')
                elif no_background=='on':
                    afile.write('\n Brokern Power Law used')
                else:
                    afile.write('\n Default Bending PL used')
                afile.write('\nFitted_parameters:\t '+str(datalc.psdFit.x))
                afile.close()
                plt.savefig(route+'DELC_'+obj+'_'+datenum+'.eps',bbox_inches='tight')
                plt.savefig(route+'DELC_'+obj+'_'+datenum+'.png',bbox_inches='tight')
                #save fake timecurves in output file
                np.savetxt(route+'fakeLC'+'_'+datenum, np.c_[delc1.time,delc1.flux,delc1.errors])
            else:
                while SATISFIED<0:
                    # 2nd simulated LC
                    delc2 = Simulate_DE_Lightcurve(datalc.psdModel, datalc.psdFit.x, datalc.pdfModel, datalc.pdfFit.x, tbin=1, LClength=lendays)
                    delc2.errors = 0*delc2.flux
                    Comparison_Plots([datalc,delc1,delc2],names=["OBS","PREV FAKE","NEW FAKE"],bins=50,norma=1.0)
                    plt.show()
                    SATISFIED =float( input('Satisfied (Yes: >0 / No <=0):\t'))
                    delc1=delc2
            flx = delc1.flux
            ft = delc1.time*delta/dilate
    ##OLD METHOD: cut a part of the new LC
    #        if init==-1:
    #            init=list(delc1.flux).index(max(delc1.flux))-int(lendays/tbin/2)
    #            if init<0:
    #                init=0
    #            if init>len(delc1.time)*tbin-lendays:
    #                init=len(delc1.time)*tbin-lendays-1
    #        flx=delc1.flux[init:init+lenfk]
    #        ft=delc1.time[init:init+lenfk]*delta/dilate

if NEW_delc!='on':
     delc1 = Load_Lightcurve(route+'fakeLC',tbin)
     flx = delc1.flux
     ft = delc1.time*delta/dilate     
#     if lendays/tbin < len(delc1.flux):
#         if init==-1:
#            init=list(delc1.flux).index(max(delc1.flux))-int(lendays/tbin/2)
#            if init<0:
#                    init=0
#            if init>len(delc1.time)*tbin-lendays:
#                    init=len(delc1.time)*tbin-lendays-1
#     ft=delc1.time[init:init+lenfk]*delta/dilate
#     flx=delc1.flux[init:init+lenfk]
#%%
""" Transform to Timecurve the cut lightcurve [ft, flx]"""
flx[0]=1.0 #first point from steady state value
fpoints=norm*flx
#use the sqruare root of variations (for SSC/ le parameter)
if POWER!='no':
    fpoints=norm*(flx)**POWER #/np.mean((flx)**POWER)
minnon0=min(fpoints)

if namevar=='gmax':
    s1=np.log10(max(flx)/min(flx))/widthgamma
    fpoints=norm*flx**(1/s1)
    
#Interpolation
funcfake=inpl.interp1d(ft,np.log(fpoints),kind=imethod)
if tcr<tbin:
    N=(int(tbin/tcr)+1)*nint
else:
    N=nint+1
if namevar=='delta':
    N=1
tbin_new=tcross/tcr*tbin/N
print 'New bin interval\n tbin_new=\t {} ={} t_cross'.format(tbin_new,round(tbin_new/tcross,5))
lenfake=N*(lenfk-1)

ftime=np.linspace(ft[0],ft[lenfk-2],lenfake)
fflux=np.exp(funcfake(ftime))

TIME=float(int(lenfake/tcross*tbin_new))
STEPS=int(lenfake/tcross*tbin_new)
#set first date t=0
t0=ftime[0]
ftime=ftime-np.ones(len(ftime))*t0

"""Create fakeTC to conjoin with crashed run"""
conshort = 1 #num of tcrosses of steady state before entering new time-varying period
logf=open(route+'multi.log','a')
if MULTI=='on' and (lendays/tcr-remt>10.0):
    tprev=TIME-(remt+conshort)
    Tcr0=ftime[ftime>max(ftime)-(remt+conshort)*tcross][0]     #set first date t=0
    fflux=fflux[ftime>max(ftime)-(remt+conshort)*tcross]
    ftime=ftime[ftime>max(ftime)-(remt+conshort)*tcross]-Tcr0
    #steady state for sometime
    fflux[ftime<conshort*tcross]=fflux[ftime<conshort*tcross][-1]
    TIME=float(remt)
    STEPS=int(remt)
    logf.write('\n'+str(round(tprev,3)))
    SAVE='no'
else:
    Tcr0=0.000
    logf.write('0.000')
logf.close()

#correct for initial steady state
ftime[0]=10**-15 #non-zero first time
#correct for steady state initial conditions
fflux[1]=fflux[0]
#correct for negative values
for i in range(0,len(fflux)-1):
    if fflux[i]<minnon0:
        fflux[i]=minnon0
        print('Corrected some values')
#save timcurve
np.savetxt(route+'fakeTC'+str(lendays)+'_'+namevar+'_'+obj+'.txt', np.c_[ftime/tcross,fflux])

cn4=open(route+'code.inp','r')
LL=[]
i=0
for line in cn4:
                if i==0:
                    LL.append('0  5  '+str(STEPS)+'  '+str(TIME)+'\n')
                elif i==2:
                    LL.append(str(round(R,int(-np.log10(R)+2)))+'   '+str(round(B,int(-np.log10(B)+2)))+' 1. 1.\n')
                elif i==4:
                        LL.append('1   '+str(round(loggmin,3))+'   '+str(round(loggmax,3))+'   '+str(round(p+0.01,3))+'   '\
                          +str(round(logle,3))+'\t'+str(round(besc,3))+'\t 0.\t0' +'\n')
                elif i==5:
                    LL.append('1   '+str(round(T,int(-np.log10(T)+2)))+'   '+str(round(lext,int(-np.log10(lext)+2)))+'\n')
                else:
                    LL.append(line)
                i=i+1                
cn4.close()
newcn4=open('code.inp','w')
newcn4.writelines(["%s" % item  for item in LL])
newcn4.close()


## Save info of the time-curve
if SAVE=='on':
    import datetime
    dttm = datetime.datetime.now()
    datenum=str(dttm.date()).split('-',1)[1]+str(dttm.time()).split('.')[0].replace(':','_')
    if NEW_delc=='on':
        np.savetxt(route+'fakeTC'+str(lendays)+namevar+'_'+datenum+'.txt', np.c_[ftime/tcross,fflux])
    savename='Fake Timecurve of'+namevar+' for object'+obj
    cfile=open(route+'info_fakecurve_'+datenum+'.txt','w')
    cfile.write('#'+savename+'\n')
    cfile.write('Steady state value of '+namevar+':\t'+str(norm)+'\t log('+namevar+')='+str(np.log10(norm))+' \n')
    cfile.write('Initial Time in rest frame:\t'+str(t0)+' MJD [Position in 10yr array+'+str(init)+'/ '+str(len(delc1.time))+']')
    cfile.write('Time interval (rest frame):\t'+str(lenfake/N)+' days \n')
    cfile.write('delta :\t'+str(delta)+'\n')
    cfile.write('Contaction of Cosmic Time for Small Time Interval:\t'+str(dilate))
    cfile.write('Time Interval jet frame:\t'+str(lenfake*tbin_new)+' days \n')
    cfile.write('Radius R\':\t'+str(R)+' cm \n')
    cfile.write('Initial t_bin on jet frame:\t'+str(tbin*delta/dilate)+'\n')
    cfile.write('Selected t_bin_new < t_cross:\t'+str(tbin_new)+'days \t ='+str(tbin_new/tcross)+' t_cross \n')
    cfile.write('New number of points:\t'+str(lenfake)+'\n')
    cfile.write('Method of interpolation:\t'+imethod+'\n')
    cfile.write('Power Index of observed Variability:\t'+str(POWER)+'\n')
    cfile.write('Maximum value / Mean value:\t'+str(round(max(fflux)/np.mean(fflux),3)))
    cfile.write('\nInitial Normalization (erg/cm^2/s):\t'+str(norma))
    if NEW_delc=='on':
        cfile.write('\n DATAcut:\t'+DATAcut)
        cfile.write('\n DATAreplace:\t'+DATAreplace)
        cfile.write('\n\n FIT PARAMETERS\n')
        cfile.write('PSD one index model:\t'+oneindex_model+'\n')
        cfile.write('PSD:\t'+str(datalc.psdFit['x'])+'\n')
        model2name=None
        if model2!=None: model2name= model2.name
        cfile.write('PDF models:\t'+str(model1.name)+'\t'+str(model2name)+'\n')
        cfile.write('PDF:\t'+str(datalc.pdfFit['x'])+'\n')
    else:
        cfile.write('\n\n Previous fake LC used, check previous info files')
    cfile.close()
    print('Info and Graph Saved into /Lightcurve_Simulation/obj Folder')
    


"""PLOTS"""
plt.figure(4,[15,5])
tused = ftime-(ftime[0]-Tcr0)*np.ones(len(ftime))
plt.plot(tused,fflux,'b-',markersize=1.5)
plt.plot(ft-ft[0]*np.ones(len(ft)),fpoints,'r.')
plt.plot(ft-ft[0]*np.ones(len(ft)),norm*np.ones(len(fpoints)),'k--',markersize=0.3)
plt.plot(ft-ft[0]*np.ones(len(ft)),np.mean(fflux)*np.ones(len(fpoints)),'g--',markersize=0.3)
plt.legend(['interpolation','initial points','normalization of fake LC: '+namevar+'_0','<'+namevar+'> ='+str(round(np.mean(fflux),3))+'= '+str(round(np.mean(fflux)/norm,2))+namevar+'_0'])
plt.xlabel(r'$t\;(days)$ jet frame')
plt.ylabel(r'log $'+namevar+'$')
plt.yscale('log')
plt.xlim(min(tused),max(tused))
plt.title(obj)
if SAVE=='on':
    plt.figure(4).savefig(route+'fake_lightcurve_'+namevar+'.png')


plt.figure(22,[15,5])
plt.plot(ftime[0:-1]/tcross,np.log10(fflux[0:-1]/fflux[1::])/(ftime[1]-ftime[0])*tcross,'r-',markersize=1.5)
plt.plot(ftime[0:-1]/tcross,1.*np.ones(len(ftime)-1),'k--')
plt.plot(ftime[0:-1]/tcross,-1.*np.ones(len(ftime)-1),'k--')
plt.xlabel(r'$t\;(t_{cross})$ jet frame')
plt.ylabel(r'$d(log\,'+namevar+')/dt$')
plt.title(obj)
if SAVE=='on':
    plt.figure(22).savefig(route+'Steepness.png')
    
#%% EXTRA PLOTS
""" Compare Real and Fake Timecurves"""
fig5=plt.figure(5)
plt.hist(np.log10(fflux/np.mean(fflux)),bins=int(len(fflux)/100),alpha=0.5)
plt.hist(np.log10(datalc.flux/np.mean(datalc.flux)),bins=int(len(datalc.flux)/100),alpha=0.5)
plt.hist(np.log10(delc1.flux/np.mean(delc1.flux)),bins=int(len(delc1.flux)/100),alpha=0.5)
plt.legend(['Fake TC + Interpolation','Real','Fake'])
plt.ylabel(r'$N$',size=14)
plt.xlabel(r'$log\;(\,Flux)$',size=14)
#plt.xlim(1.1*min(np.log10(datalc.flux/np.mean(datalc.flux))),0.9*max(max(np.log10(datalc.flux/np.mean(datalc.flux/np.mean(datalc.flux)))),\
#                  np.log10(max(fflux)/np.mean(fflux))))
if SAVE=='on':
    fig5.savefig(route+'PDF_comp.png')
#plt.pause(2)
#plt.close()
if lendays/tcr-remt<10.00:
    plt.pause(5)
    plt.close('all')

fig6=plt.figure(6)
plt.hist(np.log10(fflux/np.mean(fflux)),bins=50,density=True,alpha=0.5)
plt.hist(np.log10(eflux/np.mean(eflux)),bins=50,density=True,alpha=0.5)
plt.hist(np.log10(datalc.flux/np.mean(datalc.flux)),bins=50,density=True,alpha=0.5)
plt.legend(['Fake+Interpolation','Unedited Real','Edited Real'])
plt.ylabel(r'$PDF$',size=14)
plt.xlabel(r'$log\;(\,Flux)$',size=14)
if SAVE=='on':
    fig6.savefig(route+'PDF_density.png')
    
#%%
#%%
    
"""PLOTS"""
i0  = list(flx).index(max(flx))+250
ft0  = ft[i0] #days jet frame
fig , ax  = plt.subplots(3,1, figsize=  [8,5], sharex= True)
lend=150


ftime=np.linspace(ft[0],ft[lenfk-2],lenfake)
fflux=np.exp(funcfake(ftime))

tused = ft0-lend + np.linspace(0,lend*2,lend*2)
ax[0].set_title('Synthetic Time Curves',size=16)
ax[0].plot(ft-ft0-lend, norma/1e-9*flx, 'bD',markersize=4.0)
ax[0].set_ylabel(r'$F_{\gamma,-9} \;[{\rm erg\, s^{-1} cm^{-2}}]$',fontsize=12) 
ax[0].set_yscale('log')

aa = ax[1]
aa.plot(ftime-ft0-lend,fflux,'-',color='grey',linewidth=0.5)
aa.plot(ft-ft0-lend,fpoints,'mo',markersize=4.0)
#aa.plot(ftime-ft0-50,fflux,'ms',markersize=2.5)
plt.legend(['Interpolation','Generated points'])
aa.set_ylabel(r'$ l_e(t)  [G] $',size=14)
aa.set_yscale('log')

ax24= plt.axes([0.6,0.35, .35, .35]) # units unity for axis scales
ax24.tick_params(axis='y',labelright='on',labelleft='off',right='on',direction='in',labelsize=10)
ax24.tick_params(axis='x',labelsize=9,pad=-15)
#ax24.set_title('50 days')
ax24.set_xlim(0,30)
ax24.plot(ftime-ft0-lend,fflux,'.',color='grey',markersize=1.0)
ax24.plot(ft-ft0-lend,fpoints,'mo',markersize=4.0)
ax24.set_yscale('log')
# insert the zoomed figure
plt.setp(ax24)    
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(aa, ax24, loc1=1, loc2=3, fc="none", ec="0.5")


POWER_B=-1.0
fpointsB = B*(fpoints/np.mean(fpoints))**(POWER_B/POWER) #e.g. if fflux is le or lext
ffluxB= B*(fflux/np.mean(fflux))**(POWER_B/POWER) #e.g. if fflux is le or lext

aa = ax[2]
aa.plot(ftime-ft0-lend,ffluxB,'-',color='grey',linewidth=0.5)
aa.plot(ft-ft0-lend,fpointsB,'co',markersize=4.0)
#aa.plot(ftime-ft0-50,ffluxB,'cs',markersize=2.5)

aa.legend(['Interpolation','Generated points'])
aa.set_ylabel(r'$ B(t)  [G] $',size=14)
aa.set_yscale('log')
ax[-1].set_xlabel(r'$t_{jet}$ [days]',size=14)
ax[-1].tick_params(axis='x', labelsize=12)
ax[-1].tick_params(axis='y', labelsize=12)
ax[0].tick_params(axis='y', labelsize=12)
ax[1].tick_params(axis='y', labelsize=12)
aa.set_xlim(0,2*lend)
if SAVE=='on':
    plt.savefig(route+'fake_curves_'+namevar+'.png')
    plt.savefig(route+'fake_curves_'+namevar+'.eps')