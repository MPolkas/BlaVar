#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 13:54:06 2019

Converter: Code Output to Observable Diagram for "name".81 and "name".85 files stimutaniously
@author: markos
"""
import numpy as np
from scipy import interpolate as inpl
from matplotlib import pyplot as plt
from astropy.io import ascii
import scipy.stats as stats


"""MODES: Plotting Options"""
xoptions=['mec2','Hz','GeV']
yoptions=['le','vLv','vFv','radio','Fv','Xrays','gamma'] #radio -> Jy / Xrays -> counts /sec / Gamma -> counts/hour
yunits={'le':' ','vLv':'erg/s','vFv':'erg/s/cm^2','radio':'Jy',\
       'Fv':'erg/s/cm^2/Hz','Xrays':'counts/s/cm^2','gamma':'counts/hour/cm^2'}
""" OPTIONS """
xaxis=xoptions[1]# CAUTION! indexing of the list begins from zero
yaxis=yoptions[2]# CAUTION! indexing of the list begins from zero
Chi2test='on'
save='on'
Multicolor='on' #Plot data points with different color for every interval of 2 years
manualBoost='on'
delta_man=40.
Gamma_man=36. #delta_man
SSConly='on'
secondBB='no' #add an additional BB on spectrum (usually DISK) but that is not involved in external photon fields of code
Gamma2=1
factbb=10. #(L_disk-L_BLR)/L_BLR
Tbb2=11000
Rbb2=10**15.89
#
#
Names=['fort'] # Name of code files .81 and .85
objects=['3C273','3C279','PKS2155304']
oo=2 #from above list 0 or 1 or 2
#Bibliograpgical Data for objects
Dists=[749,3115,543] #Mpc
Rbb=[4.47e17,5e16,1.26e16] #Radii of the BB radiation zone
#
#
#Optional Plotting Defaults
nu_cutoff=13 #energy in mec^2 to cutoff plotting 
HEIGHT=1 #height of y-axis (threshold) in logscale measured below the minimum value of flux
FACTOR=1#Shift to the y-axis values (so adding up a constant to the log values)/Default 1 
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
#
#
#
#
#
#
obj=objects[oo]
route='./'+obj+'/'
data_name=obj+".ascii"
dtpts=[]
if Multicolor=='on':
   dtpts.append(obj+'_08-10.ascii')
   dtpts.append(obj+'_10-12.ascii')
   dtpts.append(obj+'_12-14.ascii')
   dtpts.append(obj+'_14-16.ascii')
   #dtpts.append(obj+'_16-18.ascii')
#
#
#
#

"""Constants"""
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


"""INPUT Observable Parameters"""
cn4=open('./code.inp','r')
lines=cn4.readlines()
cn4.close()

tend=float(lines[0].split()[-1].replace('\n',''))
nsteps=float(lines[0].split()[-2])
R=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
B=float(lines[2].split()[1]) #dimension of source in cm
p=float(lines[4].split()[3])
loggmin=float(lines[4].split()[1])
loggmax=float(lines[4].split()[2])
T=float(lines[5].split()[1])
lext=float(lines[5].split()[2])
D=Dists[oo]*10**6*pc #"Distance of Source in Mpc"
tcross=R/c/3600/24 #days


if manualBoost=='on':
    delta=delta_man
    Gamma=Gamma_man

"""
%
%
%
%
%
%
%
"""


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

func_obs=inpl.interp1d(v_pts,vFv_pts)

for n in Names:
    """Reading Files"""
    n1=n+'.81'
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
        x85=eps85+np.log10(delta)
        xbb=eps85bb
    if xaxis=='Hz':
        x85=eps85+np.log10(delta*(me*c**2/h))
        xbb=eps85bb+np.log10((me*c**2/h))
    if xaxis=='GeV':
        x85=np.log10(delta*10**eps85/eV/10**9*me*c**2)
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
    legend_pts=['2008-2010','2010-2012','2012-2014','2014-2018']
    legl=[r'SED, $\delta=$'+str(delta)]
    if SSConly=='no':
        legl.append('Black Body')
        if secondBB=='on':
            legl.append('Disk')
    if xaxis=='mec2':
        xlab='E'
    if Multicolor=='on':
        colors =['r','k','b','g']
        form=['o','.','+','s']
           
    if SSConly=='no':
            fig2, ax2 = plt.subplots(num=4)
            plt.xlabel(r'$'+xlab+'\;\;'+' ['+xaxis+']$',fontsize=15)
            plt.ylabel(r'$'+yaxis+'\;\;['+yunits[yaxis]+']$',fontsize=15) 
            if Multicolor=='on':
                for i in range(len(dtpts)):
                    ax2.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                                 globals()['errorv'+str(i)],globals()['errorvFv'+str(i)],fmt=form[i]+colors[i],ms=3.5)
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
            ax3.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
                         globals()['errorv'+str(i)],globals()['errorvFv'+str(i)],fmt=form[i]+colors[i],ms=3.5)
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
        plt.figure(5).savefig(yaxis+'_fort.81tbb.ps')
        import datetime
        dth=str(datetime.datetime.now()).split(':')[0:2]
        savename=obj+'_'+dth[0].split('-')[2].replace(' ','-')+'-'+dth[1]
        copyfile=open('run_copy.sh','w')
        copyfile.write('#'+savename+'\n')
        #copyfile.write('cp vFv_fort.81.ps /home/markos/Desktop/Blavar/Results/Steady_States/'+savename+'_81.ps \n')
        #copyfile.write('cp vFv_fort.85.ps /home/markos/Desktop/Blavar/Results/Steady_States/'+savename+'_85.ps \n')
        copyfile.write('cp vFv_fort.81tbb.ps /home/markos/Desktop/BlaVar/Results/Steady_States/'+savename+'_81tbb.ps \n')
        copyfile.write('cp code.inp /home/markos/Desktop/BlaVar/Results/Steady_States/'+'code_'+savename+'.inp \n')
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
        



