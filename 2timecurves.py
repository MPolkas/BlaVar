
# -*- coding: utf-8 -*-
"""
Obs4fit
Created on Tue Sep 17 15:15:49 2019
@author: Markos Polkas
"""

from matplotlib import use
use('Qt5Agg')

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import  AutoMinorLocator
from astropy.io import fits
try:
    from stingray import Lightcurve, AveragedCrossspectrum #timelags
    STINGRAY='on'
except:
    STINGRAY='no'


obj=        '3C273'
nmlc1 = 'fakeLC' #name of LC 1
nmlc2=  'fakeLC' # name of  LC 2
filts = ['fermi'] # Please, go into script to define in which columns the flux and errors are 

route , routeimage =      './' , './' 
namevar=     'fermi'
imtyp=       '.png'
#OPTIONS / functions
save='on'
segment_size=  160

#DCF main settings
lentau=50 # days|| integer: -before/ + after zero timelag to calculate
taubin= 1.0 #days
verbose = True

diccol={'fermi':'b','B':'g','V':'y','I':'orange','B-I':'m'} #color of errobars
names=[nmlc1,nmlc2]
for ifilter,filt in enumerate(filts):
    ii=1
    for nm in names:
        ttemp, rtemp, ertemp,rtemp2,ertemp2,col = [], [], [],[],[],[]
        if '.csv' in nm:
            with open(nm) as csvfile:
             # dtl= list(csv.reader(csvfile))        
                 globals()['data'+str(ii)]= np.genfromtxt(csvfile, delimiter=',',dtype=float)
            for row in globals()['data'+str(ii)]:
                if str(row[3])!='nan':
                    ttemp.append(int(row[3]-2400000.5)) #MJD
                    rtemp.append(row[5])
                    ertemp.append(row[6])
                    rtemp2.append(row[7])
                    ertemp2.append(row[8])
            globals()['t'+str(ii)]=np.array(ttemp)
            globals()['f'+str(ii)]=0.5*(np.array(rtemp)+np.array(rtemp2))
            globals()['err'+str(ii)]=np.sqrt(np.array(ertemp)**2+np.array(ertemp2)**2)
        elif '.fits' in nm:
            datfile=fits.open(nm)
            globals()['data'+str(ii)]=datfile[1].data
            for tb in  globals()['data'+str(ii)]:
                ttemp.append(tb[0])
                rtemp.append(tb[1])
                ertemp.append(tb[2])
            globals()['t'+str(ii)]=np.array(ttemp)
            globals()['f'+str(ii)]=np.log10(rtemp)
            globals()['err'+str(ii)]=np.array(ertemp)/np.array(rtemp)/np.log(10)    
        elif '.dat' in nm or '.txt' in nm:
            globals()['data'+str(ii)] = [i.strip().split() for i in open(nm).readlines()]
            for tb in  globals()['data'+str(ii)][1::]:
                ttemp.append(float(tb[0]))
                # if clauses because we use the same file for both LC
                rtemp.append(float(tb[1]))
                ertemp.append(float(tb[2]))
            globals()['t'+str(ii)]=np.array(ttemp)
            globals()['f'+str(ii)]=np.array(rtemp)
            globals()['err'+str(ii)]=np.array(ertemp)
        else:
            try:
                globals()['data'+str(ii)] = [i.strip().split() for i in open(nm).readlines()]
                for tb in  globals()['data'+str(ii)][1::]:
                    ttemp.append(float(tb[0]))
                    # if clauses because we use the same file for both LC
                    rtemp.append(float(tb[1]))
                    ertemp.append(float(tb[2]))
                globals()['t'+str(ii)]=np.array(ttemp)
                globals()['f'+str(ii)]=np.array(rtemp)
                globals()['err'+str(ii)]=np.array(ertemp)
            except:
              print( 'Please, make *object*.*filter*raw.csv/ .fits /.csv/.dat input files')
        ii=ii+1
    #%%
    """CORRELATION DIAGRAM"""
    #Find time window to compare (times are rounded to 1.0 MJD)
    tt,lc1,lc2,erlc1,erlc2=[],[],[],[],[]
    for i in t1:
        for j in t2:
            if abs(i-j)<1 and (err1[list(t1).index(i)]<10 and err2[list(t2).index(j)]<10):
                tt.append(i)
                lc1.append(f1[list(t1).index(i)])
                lc2.append(f2[list(t2).index(j)])
                erlc1.append(err1[list(t1).index(i)])
                erlc2.append(err2[list(t2).index(j)])
    lc1 , lc2 , err1 ,err2 =np.array(lc1) , np.array(lc2) ,np.array(erlc1), np.array(erlc2)
    #average flux for better data manipulation
    nlc1=lc1+2.5*np.log10(np.mean(10**(-0.4*lc1)))  
    nlc2=lc2+2.5*np.log10(np.mean(10**(-0.4*lc1)))
    #CORRELATION DIAGRAM
    plt.figure(1+ifilter)
    plt.errorbar(nlc1,nlc2,yerr=erlc2,xerr=erlc1,fmt='r.',ecolor=diccol[filt],capsize=2.4)
    limx,limy=max(abs(min(nlc1)),max(nlc1)),max(abs(min(nlc2)),max(nlc2))
    plt.plot(np.linspace(-max(limx,limy),max(limx,limy)),np.linspace(-max(limx,limy),max(limx,limy)),'k-')
    #plt.axis([-limx*0.95,limx*1.05,-limy*0.95,limy*1.05])
    #plt.axis([-2.5,2.5,-2.5,2.5]) #BOX
    plt.xlabel(nmlc1,fontsize=14)
    plt.ylabel(nmlc2,fontsize=14)
    plt.title(obj+'   '+nmlc1+' vs '+nmlc2,fontsize=16)
    if save=='on': plt.figure(1).savefig('Corrs_'+filt+imtyp)
    #%% STINGRAY
    """Stingray/ Timelags per frequency range"""
    if STINGRAY=='on':
        figstg,ax=plt.subplots(1,1,num=10+ifilter,sharex='col')
        
        long_times=tt
        long_signal_1=np.array(lc1) #mag to flux 
        long_signal_2=np.array(lc2) #mag to flux
    
        long_lc1 = Lightcurve(long_times, long_signal_1)
        long_lc2 = Lightcurve(long_times, long_signal_2)
        
        avg_cs = AveragedCrossspectrum(long_lc1, long_lc2,segment_size)
        freq_lags, freq_lags_err = avg_cs.time_lag()
        #PLOT
        ax.hlines(0, avg_cs.freq[0], avg_cs.freq[-1], color='black', linestyle='dashed', lw=2)
        ax.errorbar(avg_cs.freq, freq_lags, yerr=freq_lags_err,fmt=".", lw=1, color='r',label=filt+' - B')
        ax.legend()
     #   axax=ax.axis()
     #   ax.axis([axax[0],avg_cs.freq[-1]/2,axax[2],axax[3]]) 
        ax.set_ylabel("Time lag (days)")
        plt.title('Observations')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.tick_params(which='major', width=1.5, length=7)
        ax.tick_params(which='minor', width=1.5, length=4)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)
            ax.set_xlabel("Frequency (days^-1)")
                
        if save=='on':  figstg.savefig(routeimage+'tlags_freq'+'_'+filt+imtyp)
#%%  DCFS
    #taubin=0.15
    tobs=tt
    dtau=tobs[3]-tobs[2]
    taurange=np.linspace(-lentau,lentau,2*int(lentau/taubin)+1)
    tlags , DCFs =[] , []
    a, b = nlc1 , nlc2 #mag to flux scale 
    erra , errb =  np.log(10)*a*err1, np.log(10)*b*err2 # 0.0*np.ones(len(a)), 0.0*np.ones(len(a)) # can be set to zero if outliers not removed 
    amean, bmean =np.mean(a) , np.mean(b) #averages
    astd , bstd = np.std(a) , np.std(b) #standard deviation
    tlag,  DCF, step =[] , [] , 0 #initialize variables
    fig25=plt.figure(20+ifilter)
    for tau in taurange:
        udcf, M = 0 , 0
        for i in range(len(a)-1):
         aa=a[i]
         era = erra[i]
         for j in range(len(b)-1):
            bb=b[j]
            erb = errb[j]
            dt=tobs[i]-tobs[j]
            if dt<tau+taubin/2 and dt>tau-taubin/2:
                udcf=udcf+(aa-amean)*(bb-bmean)/np.sqrt(astd**2-era**2)/np.sqrt(bstd**2-erb**2) #/astd/bstd
                
                M=M+1
        tlag.append(tau)
        if M!=0:
            DCF.append(udcf/M)
        else:
            DCF.append(0.)
        if verbose: print('Step  '+str(step+1)+' / '+str(len(taurange))+'     bands pair compared:  '+nmlc1+'-'+nmlc2)    
        step=step+1 #count steps
    tlags.append(np.array(tlag).astype(np.float))
    DCFs.append(np.array(DCF).astype(np.float))
    np.savetxt('DCFs.txt', np.c_[np.concatenate(tlags, axis=0 ),np.concatenate(DCFs, axis=0 )])
    plt.plot(tlags[0],DCFs[0],'r-',label='fermi vs  fermi')
    plt.legend()
    plt.xlabel(r'$\tau\;\;days$')
    plt.ylabel(r'$DCF$($\tau$)')
    plt.title('DCF    '+str(lentau))
    ll=len(tlags)
    plt.plot(np.linspace(-lentau,lentau,50),0*np.linspace(-1,1,50),'k--',linewidth=0.5)
    plt.plot(0*np.linspace(-1,1,50),np.linspace(-1.,1,50),'k--',linewidth=0.5)
    plt.axis([-lentau,lentau,-1,1])
    plt.show()
    
    if save=='on': fig25.savefig(routeimage+'tlags_'+str(taubin)+'_'+obj+imtyp)
