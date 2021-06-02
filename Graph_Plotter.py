# -*- coding: utf-8 -*-
"""
Thesis Plotter
This is a temporary script file.
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
Ledd=3.2*10**4*Lsun*10**8 #M/10^8 Ms    un
#save plots at
route = os.getcwd()
routeimage= route +'/'

#Observ Soldi (s)/ Chetalier(c)
#calculated via averaging of points in the mid-(log)range frequency
##xs=np.array([ 9.4169,  9.6762,  9.8828,  9.952 , 10.156 , 10.348 , 10.536 ,\
#       11.148 , 10.943 , 11.335 , 11.439 , 11.523 , 11.783 , 11.881 ,\
#      13.153 , 13.438 , 13.799 ,
xs=np.array([ 13.929 , 14.105 , 14.243 , 14.365 ,\
       14.699 , 14.786 , 14.89  , 15.029 , 15.028 , 15.097 , 15.149 ,\
       15.236 , 15.322 , 16.35  , 16.668 , 17.067 , 17.362 , 17.656 ,\
       18.076 , 18.355 , 18.668 , 19.067 , 19.029 , 19.345 , 19.692 ,\
       19.759 , 20.022 ])#  21.83  , 22.367 , 22.828 , 23.346 ])

#ys=np.array([0.04016 , 0.093522, 0.15521 ,  0.16689 , 0.29024 , 0.2836  ,\
 #      0.37195 , 0.34535 , 0.63533 , 0.44871 , 0.46205 , 0.53039 ,\
#       0.58209 , 0.73377 , 0.35557 , 0.23227 , 0.34065 ,
ys=np.array([0.15733 ,\
       0.11735 , 0.15736 , 0.16238 , 0.10241 , 0.11576 , 0.1291  ,\
       0.14078 , 0.15245 , 0.17079 , 0.17913 , 0.19247 , 0.20915 ,\
       0.2426  , 0.15763 , 0.20101 , 0.22604 , 0.27774 , 0.24946 ,\
       0.26449 , 0.27785 , 0.3329  , 0.40789 , 0.34626 , 0.4013  ,\
       0.46964 , 0.43134 ]) # 0.57654 , 0.6466  , 0.48165 , 0.19671 ])

rm=17
eys=np.ones(len(xs))*0.05
eys[10-rm] ,eys[11-rm] ,eys[12-rm], eys[13-rm]  ,eys [14-3], eys[16-rm] , eys[29-rm] ,eys[30-rm]  =0.1 , 0.075, 0.1 , 0.15 , 0.1 , 0.1 , 0.1, 0.075
eys[39-rm] ,eys[41-rm], eys[42-rm] , eys[43-rm] = 0.075 , 0.1 , 0.075, 0.1
#eys[44] ,eys[45] , eys[46] ,eys[47] =  0.2 , 0.1 , 0.2 , 0.3


xc=[14.837,14.74, 14.67, 14.383, np.log10(1e3*eV/h),np.log10( 3e3*eV/h), np.log10(5e3*eV/h),\
    np.log10(0.3e9*eV/h),0.5*(np.log10(300e9*eV/h)+np.log10(0.1e9*eV/h)),\
    np.log10(2e12*eV/h)]     
dxc=[0.08,0.08,0.08,0.08, np.log10(2/0.3), np.log10(2)/2, np.log10(7/4)/2,  np.log10(10)/2, np.log10(300)/2,np.log10(50)/2]
#yc=np.array([0.39,0.8,0.42])
#eyc=np.array([0.02,0.05,0.025])
yc = [ 0.378, 0.371 , 0.369 , 0.383, 0.591, 0.716 ,  0.796, 0.36 , 0.43   , 0.657 ]
eyc=  [ 0.004 , 0.004 , 0.003 , 0.005, 0.004 , 0.003 , 0.004 , 0.075 , 0.04, 0.008]


#average errorbars from LCs
erro=[0.00307/2.5, 0.00307/2.5, 0.005036082474226804/2.5] #average magnitude error from SMARTS LCs
errx=[0.003695, 0.003695 , 0.003149]
errg=[0.1456, 0.1456 , 0.26115] # average relevant error for error < norms cases

#Fractional Variabilities and errors from Sims
range_smarts0=np.log10([10**13.8,10**14.2])
range_smarts=np.log10([10**14.2,10**14.6])
range_smarts1=np.log10([10**14.6,10**15.0])
range_fermi=np.log10([0.1*10**9*(2.4184*10**14),300*10**9*(2.4184*10**14)]) #eV/h = 2.41e14
range_xsoft=np.log10([2e3*(2.4184*10**14),10e3*(2.4184*10**14)]) #eV/h = 2.41e14
range_xhard=np.log10([1e4*(2.4184*10**14),8e4*(2.4184*10**14)]) #eV/h = 2.41e14
range_EUV=np.log10([10*(2.4184*10**14) ,124*(2.4184*10**14)]) #eV/h = 2.41e14
range_fermiB=np.log10([0.1*10**9*(2.4184*10**14),1*10**9*(2.4184*10**14)]) #eV/h = 2.41e14
range_fermiC=np.log10([100*10**9*(2.4184*10**14),4*10**12*(2.4184*10**14)])
energies=np.array([(range_smarts0[0]+range_smarts0[1])/2,
         (range_smarts[1]+range_smarts[0])/2,
         (range_smarts1[1]+range_smarts1[0])/2,
          (range_EUV[1]+range_EUV[0])/2,
          (range_xsoft[1]+range_xsoft[0])/2,
          (range_xhard[1]+range_xhard[0])/2,
          (range_fermiB[1]+range_fermiB[0])/2,
          (range_fermi[1]+range_fermi[0])/2,
          (range_fermiC[1]+range_fermiC[0])/2])
denergies=np.array([(range_smarts0[1]-range_smarts0[0])/2,
          (range_smarts[1]-range_smarts[0])/2,
          (range_smarts1[1]-range_smarts1[0])/2,
          (range_EUV[1]-range_EUV[0])/2,
          (range_xsoft[1]-range_xsoft[0])/2,
          (range_xhard[1]-range_xhard[0])/2,
          (range_fermiB[1]-range_fermiB[0])/2,
          (range_fermi[1]-range_fermi[0])/2,
          (range_fermiC[1]-range_fermiC[0])/2])
##3C273
#sim1=   [ 0.53  , 0.35 , 0.28  ] #B , B=10G  365 days
#er1 =   [  0.06 , 0.04 ,0.17   ] 
#sim2=   [ 0.18  , 0.383 ,1.008   ] # B, B=26 G 365 days
#er2 =   [  0.1 ,0.03  ,0.05   ]
#sim3=   [  0.42 , 1.57 , 0.84  ]  #le B=26 G 365 days
#er3 =   [  0.05 , 0.001 ,0.06   ]
#sim4=   [  0.2 ,0.65  , 0.62  ] #lext ,B=26G 365 days B=26G 
#er4 =   [  0.2 ,0.02,   0.07]
#sim5=  [  0.33 , 0.78 ,   1.10] #delta B=26G ind=1/5  delta=(fermi)^ind
#er5= [  0.06 , 0.01 ,   0.04]  
#
#sim6=   [  0.45 , 3.088 , 1.97  ] #le,B=26G for 10yrs days of obs add as le2
#er6 =   [  0.019 , 0.002 ,0.014   ]
#sim7=   [  0.225 , 2.57 , 1.37  ] #le,B=10G for 10yrs days of obs add as le2
#er7 =   [  0.05 , 0.003 ,0.019   ]
#sim8=   [  0.04 , 0.51, 0.856  ] #lext ,B=26G for 10yrs days of obs 
#er8 =   [  0.10 , 0.01 ,0.03  ]
#sim9 = [ 0.055 ,  0.470,  0.461]  #lext ,B=10G for 10yrs days of obs 
#er9 =   [ 0.11 ,  0.015,  0.053]
#sim10=   [  0.12 , 0.453 , 1.29 ] #B  ,B=26G for 10yrs days of obs
#er10 =   [  0.07 , 0.0155 ,0.02 ]
#sim11=   [  0.157 , 0.30 , 0.66 ] #B  ,B=10G for 10yrs days of obs
#er11 =   [  0.07 , 0.02 ,0.03 ]
#sim12=  [  0.21 , 0.88 ,   1.26] #delta B=26G ind=1/5  delta=(fermi)^ind 
#er12= [  0.05 , 0.01 ,   0.02] 
#sim13=   [  0.382 , 0.906 , 1.286 ] #delta ,B=10G for 10yrs days of obs add 
#er13 =   [  0.03 , 0.008 ,0.02 ]
#
#
#
##sim5=   [  0.2 , 1.63 ,0.16   ] #delta 0.25
##er5 =   [   0.11, 0.011 ,0.21   ] #delta 0.25
##sim7= [  0.25, 0.61 ,   0.83]  #delta ind=1/6
##er7= [  0.07 , 0.01 ,   0.05] #delta ind=1/6
#
#
##sim=   [    0.12   ,  0.376   ,   1.26    ] #B^-1 ~ same as B
##er =   [   0.13   ,  0.03     ,   0.04    ] #B^-1


##PKS2155304
#simB1=  [  0.33 ,0.31  ,0.57   ] #le 3yrs
#erB1 =  [  0.05 ,0.03  ,0.05   ]
#simB2=  [  0.68 ,0.72  ,0.41   ]# B 1yr
#erB2 =  [  0.03 , 0.016 , 0.1  ]
#simB3=  [  0.66 , 0.95 ,  0.78 ] #delta 0.25 1000days
#erB3 =  [  0.03 , 0.01 , 0.05  ]
##simB4=  [  0.07 , 0.37 , 0.07  ] #gmax 0.25
##erB4=   [  0.3 , 0.05 ,  0.3 ]
#simB4=  [  0.43 , 0.413 , 0.834  ] #le 3590 days
#erB4=   [  0.026 , 0.015 ,  0.03 ]
#simB5=  [  1.11 , 1.16 , 0.54  ] #B 3590 days
#erB5=   [  0.01 , 0.006 ,  0.04]
#simB6=  [  1.027 , 1.777 , 1.344  ] #delta 3590 days
#erB6=   [  0.013, 0.004 ,  0.017]
#simB3= [0.48716667770597566, 0.7142049902429926, 0.5803481540605783] #delta 0.2
#erB3= [0.04620883804872089, 0.01664457330043878, 0.06849240731082773] #delta 0.2
#simB3= [0.41591764431139866, 0.6025758427463513, 0.4893311701207448] #delta 1/6
#erB3=[0.05304755144749307, 0.019441600468463327, 0.07904891531683922] #delta 1/6



#3C273
sim1=   [-1.0, 0.53 ,-1.0,-1.0, 0.35 ,-1.0,-1.0, 0.28  ,-1.0] #B , B=10G  365 days
er1 =   [ -1.0, 0.06,-1.0,-1.0, 0.04 ,-1.0,-1.0,0.17 ,-1.0   ] 
sim2=   [ -1.0,0.18 ,-1.0,-1.0, 0.383 ,-1.0,-1.0,1.008 ,-1.0   ] # B, B=26 G 365 days
er2 =   [-1.0,  0.1 ,-1.0,-1.0, 0.03  , -1.0,-1.0,0.05  ,-1.0  ]
sim3=   [-1.0,  0.42 ,-1.0,-1.0, 1.57 ,-1.0,-1.0, 0.84  ,-1.0 ]  #le B=26 G 365 days
er3 =   [-1.0,  0.05 , -1.0,-1.0, 0.001 ,-1.0,-1.0,0.06 ,-1.0  ]
sim4=   [-1.0,  0.2 , -1.0,-1.0, 0.65  , -1.0,-1.0,0.62 ,-1.0 ] #lext ,B=26G 365 days B=26G 
er4 =   [-1.0,  0.2 , -1.0,-1.0, 0.02,   -1.0,-1.0,0.07 ,-1.0]
sim5=  [ -1.0, 0.33 , -1.0,-1.0, 0.78 , -1.0,-1.0,  1.10 ,-1.0] #delta B=26G ind=1/5  delta=(fermi)^ind
er5= [ -1.0 , 0.06 , -1.0,-1.0, 0.01 ,  -1.0,-1.0, 0.04 ,-1.0]  

#sim6=   [  0.439, 0.4469 , 0.448 , 2.11, 3.088 , 2.87, 2.02,  1.932, -1.0] #le,B=26G for 10yrs days of obs add as le2
#er6 =   [  0.009, 0.019 ,0.009,  0.002 ,0.003, 0.003 , 0.027, 0.027, -1.0]
#sim7=   [ 0.219, 0.228 , 0.234, 2.31, 2.57, 2.297 , 1.39, 1.31 ,-1.0] #le,B=10G for 10yrs days of obs add as le2
#er7 =   [  0.017, 0.017 , 0.016, 0.003, 0.003, 0.003 ,0.036 ,0.038, -1.0 ]
#sim8=   [  0.0528 , 0.053, 0.057, 0.13 , -1.0,  0.51, 0.683, 0.693, -1.0 ] #lext ,B=26G for 10yrs days of obs 
#er8 =   [  0.0518, 0.051 , 0.05, 0.04, -1.0 , 0.01, 0.07, 0.067, -1.0 ]
#sim9 = [ 0.06, 0.0622 , 0.065, 0.3195,  -1.0, 0.50, 0.24, 0.2589, -1.0]  #lext ,B=10G for 10yrs days of obs 
#er9 =   [ 0.048, 0.11 , 0.046, 0.018, -1.0, 0.014, 0.15,  0.053, -1.0]
#sim10=  [ 0.132,  0.125 ,0.12, 0.53, 0.453, 0.56, 1.14, 1.17,-1.0] #Β, B=26G ind=1/5  delta=(fermi)^ind 
#er10 =   [  0.03, 0.03, 0.03, 0.015, 0.015, 0.012 , 0.043 ,0.04, -1.0]
#sim11=   [  0.17, 0.159 ,0.149, 0.819, 0.30 , 0.346, 0.517, 0.54,-1.0 ] #B  ,B=10G for 10yrs days of obs
#er11 =   [  0.02, 0.02 ,0.025, 0.007, 0.007 , 0.007, 0.09, 0.08 ,-1.0]
#sim12=  [ 0.40,  0.383, 0.359, 1.285, 0.453, 1.03, 1.16, 1.22, 0.95] #delta B=26G ind=1/5  delta=(fermi)^ind  10yrs
#er12= [ 0.015, 0.015 , 0.011, 0.004, 0.007 ,0.007, 0.04, 0.04, 0.05] 
#sim13=   [ 0.226, 0.216 ,0.202, 1.236,  0.889 , 0.974, 1.137 , 1.286, 0.966 ] #delta ,B=10G  ind=1/5 10yrs days 
#er13 =   [ 0.017, 0.018 , 0.02, 0.005, 0.008 , 0.007, 0.04, 0.04 ,0.05 ]

sim6=   [  0.439, 0.4469 , 0.448 , -1.0, 3.088 , 2.87, -1.0,  1.964, -1.0] #le,B=26G for 10yrs days of obs add as le2
er6 =   [  0.009, 0.019 ,0.009,  0.002 ,0.003, 0.003 , 0.027, 0.027, -1.0]
sim7=   [ 0.219, 0.228 , 0.234,-1.0, 2.57, 2.297 , -1.0, 1.31 ,-1.0] #le,B=10G for 10yrs days of obs add as le2
er7 =   [  0.017, 0.017 , 0.016, 0.003, 0.003, 0.003 ,0.036 ,0.038, -1.0 ]
sim8=  [ 0.132,  0.125 ,0.115, -1.0, 0.453, 0.56, -1.0, 1.23,-1.0] #Β, B=26G ind=1/5  delta=(fermi)^ind 
er8 =   [  0.03, 0.03, 0.03, 0.015, 0.015, 0.012 , 0.043 ,0.04, -1.0]
sim9=   [  0.17, 0.159 ,0.149, -1.0, 0.30 , 0.346, -1.0, 0.54,-1.0 ] #B  ,B=10G for 10yrs days of obs
er9 =   [  0.02, 0.02 ,0.025, 0.007, 0.007 , 0.007, 0.09, 0.08 ,-1.0]
sim10=  [ 0.40,  0.383, 0.359, -1.0, 0.453, 1.03, -1.0, 1.265, 0.95] #delta B=26G ind=1/5  delta=(fermi)^ind  10yrs
er10= [ 0.015, 0.015 , 0.011, 0.004, 0.007 ,0.007, 0.04, 0.04, 0.05] 
sim11=   [ 0.226, 0.216 ,0.202, -1.0,  0.889 , 0.974, -1.0 , 1.286, 0.966 ] #delta ,B=10G  ind=1/5 10yrs days 
er11 =   [ 0.017, 0.018 , 0.02, 0.005, 0.008 , 0.007, 0.04, 0.04 ,0.05 ]
sim12=   [  0.0528 , 0.053, 0.057, -1.0 , -1.0,  0.51, -1.0, 0.801, -1.0 ] #lext ,B=26G for 10yrs days of obs 
er12 =   [  0.0518, 0.051 , 0.05, 0.04, -1.0 , 0.01, 0.07, 0.067, -1.0 ]
sim13 = [ 0.06, 0.0622 , 0.065, -1.0,  -1.0, 0.50, -1.0, 0.2589, -1.0]  #lext ,B=10G for 10yrs days of obs 
er13 =   [ 0.048, 0.11 , 0.046, 0.018, -1.0, 0.014, 0.15,  0.053, -1.0]


#TS >4 & eflux > its error
sim14=   [  0.4056, 0.4125 , 0.4171 , -1.0, 3.18, 2.94, -1.0,  1.904, -1.0] #le,B=26G for 10yrs days of obs add as le2
er14 =   [  0.009, 0.019 ,0.009,  0.002 ,0.003, 0.003 , 0.027, 0.027, -1.0]
sim15=   [ 0.2001, 0.2017 , 0.2125,-1.0, 2.61, 2.304 , -1.0, 1.2922 ,-1.0] #le,B=10G for 10yrs days of obs add as le2
er15 =   [  0.017, 0.017 , 0.016, 0.003, 0.003, 0.003 ,0.036 ,0.038, -1.0 ]

sim16=  [ 0.1094,  0.1031 ,0.0986, -1.0, 0.3901, 0.4891, 1.2146, 1.23,-1.0] #Β, B=26G ind=1/5  delta=(fermi)^ind 
er16 =   [  0.03, 0.03, 0.03, 0.015, 0.015, 0.012 , 0.043 ,0.04, -1.0]
sim17=   [  0.1353, 0.1248, 0.1169,  -1.0, 0.1969 , 0.2241, -1.0, 0.5758,-1.0 ] #B  ,B=10G for 10yrs days of obs
er17 =   [  0.02, 0.02 ,0.025, 0.007, 0.007 , 0.007, 0.09, 0.08 ,-1.0]

sim18=  [ 0.3873,  0.3687, 0.3463, -1.0, 0.8481, 0.9807, -1.0, 1.2729, 0.9837] #delta B=26G ind=1/5  delta=(fermi)^ind  10yrs
er18=   [ 0.015, 0.015 , 0.011, 0.004, 0.007 ,0.007, 0.04, 0.04, 0.05] 
sim19=   [ 0.217, 0.2082 ,0.1944, -1.0,  0.8304 , 0.9223, -1.0 , 1.252, 0.9994] #delta ,B=10G  ind=1/5 10yrs days 
er19 =   [ 0.017, 0.018 , 0.02, 0.005, 0.008 , 0.007, 0.04, 0.04 ,0.05 ]

sim20=   [  0.051 , 0.0563, 0.0587, -1.0 , -1.0,  0.63, -1.0, 0.7255, -1.0 ] #lext ,B=26G for 10yrs days of obs 
er20 =   [  0.0518, 0.051 , 0.05, 0.04, -1.0 , 0.01, 0.07, 0.067, -1.0 ]
sim21 = [ 0.046, 0.0622 , 0.0606, -1.0,  -1.0, 0.6324, -1.0, 0.3938, -1.0]  #lext ,B=10G for 10yrs days of obs 
er21 =   [ 0.048, 0.11 , 0.046, 0.018, -1.0, 0.014, 0.15,  0.053, -1.0]

#Raw Unedited CVs
# le 3c corr 3600 b26 
# SMARTS band:	0.3945566633851785	+/-0.0
# X-rays band:	3.187645848460143	+/-0.0
# FERMI band:	1.905061676721904	+/-0.0
# SMARTS 0:	0.40569283738745227	+/-0.0
# X-rays hard:	2.9425132283727913	+/-0.0
# FERMI 0.1-10:	1.9051417285234924	+/-0.0
# SMARTS 2:	0.4171905789391356	+/-0.0
# EUV:	2.0797428506020923	+/-0.0

 
#le 3c corr 3600 b10
# SMARTS band:	0.201763717203663	+/-0.0
# X-rays band:	2.6119847992048375	+/-0.0
# FERMI band:	1.292180606210388	+/-0.0
#  SMARTS 0:	0.20016498648056322	+/-0.0
# X-rays hard:	2.3043151535690547	+/-0.0
# FERMI 0.1-10:	1.2921462335856728	+/-0.0
# SMARTS 2:	0.21255807598756876	+/-0.0
# EUV:	2.322046013177166	+/-0.0



#B 3c corr 3600 b26
#SMARTS band:	0.1031211906051966	+/-0.0
#X-rays band:	0.3901312583953679	+/-0.0
#FERMI band:	1.2146182011735462	+/-0.0
#  SMARTS 0:	0.10942995314770988	+/-0.0
#  X-rays hard:	0.4891214312563056	+/-0.0
#  FERMI 0.1-10:  1.2144786113949277 +/-0.0
# SMARTS 2:	0.09859703918241222	+/-0.0
# EUV:	0.4054392628166188	+/-0.0

#B 3c corr 3600 b10
# SMARTS band:	0.12486420600081542	+/-0.0
# X-rays band:	0.19693721912626264	+/-0.0
# FERMI band:	0.5758658378293686	+/-0.0
# SMARTS 0:	0.13530734912439105	+/-0.0
# X-rays hard:	0.22409380086572767	+/-0.0
# FERMI 0.1-10:	0.5759134468761115	+/-0.0
# SMARTS 2:	0.11691025470181898	+/-0.0
# EUV:	0.5971235473783574	+/-0.0
 
##delta 3c corr 3600 b26 
# SMARTS band:	0.3687812964157029	+/-0.0
# X-rays band:	0.8481272745414699	+/-0.0
# FERMI band:	1.2729048390547655	+/-0.0
# SMARTS 0:	0.38729827021975355	+/-0.0
# X-rays hard:	0.9807327508879665	+/-0.0
# FERMI 0.1-4TeV:	0.9837099626174286	+/-0.0
#  SMARTS 2:	0.34633086214322273	+/-0.0
# EUV:	1.2675462572664873	+/-0.0

# delta 3c corr 3600 b10
# SMARTS band:	0.20824202294817099	+/-0.0
# X-rays band:	0.8304030607620317	+/-0.0
# FERMI band:	1.252144410398455	+/-0.0
# SMARTS 0:	0.2174639968863494	+/-0.0
# X-rays hard:	0.922338662150839	+/-0.0
# FERMI 0.1-4 TeV:	0.9994288909905029	+/-0.0
#  SMARTS 2:	0.1944318148652205	+/-0.0
# EUV	1.2185393369047848	+/-0.0


##lext 3c corr 3600 b26
#SMARTS band:	0.056310817025487056	+/-0.0
#X-rays band:	0.5413582903363195	+/-0.0
#FERMI band:	0.7255514195581553	+/-0.0
# SMARTS 0:	0.05103467904612822	+/-0.0
# X-rays hard:	0.6324107456340355	+/-0.0
# FERMI 0.1-10:	0.7278487071899953	+/-0.0
#  SMARTS 2:	0.058681274590510236	+/-0.0
# EUV:	0.09799913736828832	+/-0.0


##lext 3c corr 3600  b10
# SMARTS band:	0.06216366928681941	+/-0.0
# X-rays band:	0.453095524112025	+/-0.0
# FERMI band:	0.3938326424432775	+/-0.0
#  SMARTS 0:	0.046240343570929175	+/-0.0
# X-rays hard:	0.6324107456340355	+/-0.0
# FERMI 0.1-10:	0.7278487071899953	+/-0.0
#  SMARTS 2:	0.060557391967956675	+/-0.0
# EUV:	0.26007303148138505	+/-0.0


#sim5=   [  0.2 , 1.63 ,0.16   ] #delta 0.25
#er5 =   [   0.11, 0.011 ,0.21   ] #delta 0.25
#sim7= [  0.25, 0.61 ,   0.83]  #delta ind=1/6
#er7= [  0.07 , 0.01 ,   0.05] #delta ind=1/6


#sim=   [    0.12   ,  0.376   ,   1.26    ] #B^-1 ~ same as B
#er =   [   0.13   ,  0.03     ,   0.04    ] #B^-1
#PKS2155304
simB1=  [  -1.0, 0.33 ,-1.0, -1.0 ,0.31  ,-1.0, -1.0, 0.57, -1.0   ] #le 3yrs
erB1 =  [ -1.0,  0.05 ,-1.0, -1.0 ,0.03  ,-1.0, -1.0, 0.05 , -1.0   ]
simB2=  [ -1.0,  0.68 ,-1.0, -1.0 ,0.72  ,-1.0, -1.0, 0.41  , -1.0 ]# B 1yr
erB2 =  [ -1.0,  0.03,-1.0, -1.0  , 0.016 , -1.0, -1.0, 0.1 , -1.0 ]
simB3=  [ -1.0,  0.66 ,-1.0, -1.0 , 0.95 ,  -1.0, -1.0, 0.78 , -1.0] #delta 0.25 1000days
erB3 =  [ -1.0, 0.03 ,-1.0, -1.0 , 0.01 ,-1.0, -1.0,  0.05 , -1.0 ]
#simB4=  [  0.07 , 0.37 , 0.07  ] #gmax 0.25
#erB4=   [  0.3 , 0.05 ,  0.3 ]
simB4=  [  0.4349, 0.4183, 0.41326, -1.0, 0.40306, 0.40300 ,-1.0, 0.8520 , 0.8112] #le 3590 days
erB4=   [  0.012, 0.015 ,0.012, 0.011,  0.015 , 0.015, 0.25, 0.19, 0.37]
simB5=  [ 1.0255, 1.0455 ,1.15, -1.0, 1.18, 2.02 , -1.0, 0.806,-1.0] #B 3590 days
erB5=   [ 0.005, 0.005 , 0.004, 0.004, 0.006, 0.006 , 0.37,  0.37, -1.0]
simB6=  [ 0.94,  1.027 ,1.155, -1.0, 1.777 , 2.53 , -1.0, 1.326 ,1.673] #delta 3590 days
erB6=   [0.005,  0.005, 0.0046, 0.0037, 0.004 , 0.08, 0.06, 0.05,0.05]

simB7=  [  0.3699, 0.3658, 0.3616, -1.0, 0.3408, 0.3414 ,-1.0, 0.7773 , 0.7302] #le 3600 days
erB7=   [  0.012, 0.015 ,0.012, 0.011,  0.015 , 0.015, 0.25, 0.19, 0.37]
simB8=  [ 0.6257, 0.5955 ,0.5979,  -1.0, 0.62, 1.088 , -1.0, 0.2733, 0.326] #B 3600 days
erB8=   [ 0.005, 0.005 , 0.004, 0.004, 0.006, 0.006 , 0.37,  0.2733, 0.326]
simB9=  [ 0.934,  1.03648 , 1.1986, -1.0, 2.0947, 3.4033 , -1.0, 1.454 ,1.9696] #delta 3600 days
erB9=   [0.005,  0.005, 0.0046, 0.0037, 0.004 , 0.08, 0.06, 0.05,0.05]

#le pks corr 3600
# SMARTS band:	0.3658132048204554	+/-0.0
# X-rays band:	0.3408202706909622	+/-0.0
# FERMI band:	0.7772682555541374	+/-0.0
# SMARTS 0:  	0.36994426118720125	+/-0.0
# X-rays hard:	0.34143231096184795	+/-0.0
# FERMI 0.1-10:	0.7969515683903577	+/-0.0
# EUV:	0.34763677426299106	+/-0.0


#B pks corr 3600
# SMARTS band:	0.5955692138912473	+/-0.0
# X-rays band:	0.6288258462355305	+/-0.0
# FERMI band:	0.2733323970823095	+/-0.0
# SMARTS 0:	0.6257119777476021	+/-0.0
# X-rays hard:	1.0884514926655893	+/-0.0
# FERMI 0.1-10:	0.3026267051678513	+/-0.0
# SMARTS 2:	0.5979655470116699	+/-0.0
# EUV:	0.6340649633099391	+/-0.0
# TeV: 0.3206


#delta pks corr 3600
# SMARTS band:	1.0364801328855786	+/-0.0
# X-rays band:	2.094745449661599	+/-0.0
# FERMI band:	1.4540947060152187	+/-0.0
# SMARTS 0:	0.9340878342419926	+/-0.0
# X-rays hard:	3.4033460811962284	+/-0.0
# FERMI 0.1-4 TeV:	1.9696806525966741	+/-0.0
# SMARTS 2:	1.1986423766518366	+/-0.0
# EUV:	1.508964165993394	+/-0.0


""" SETTINGS"""

legsize='x-large' #medium, large x-large etc
fntsz=17 #fontsize
coldic={'l_{ext}':'tab:orange','B':'c','l_e':'m','g_{max}':'r','delta':'g'}
clprev=''
marks=['o','s','^','>','<'] #for same kind of simulations
lnstyles=[(0, (3, 1, 1, 1)),'dashed', 'dashdot', 'dotted','solid',(0, (3, 2, 1, 1)), (0, (3, 2, 2, 1))]
markersize=10

listsimA=[r'$B$, $B_0=10$G 365d',\
          r'$B$, $B_0=26$G 280d',\
          r'l_e, $B_0=26$G 1000d',\
          r'$l_{ext}$ ($B_0=26G$)',\
          r'$\delta$ ($B_0=26$G)',\
          r'$l_e$ ($B_0=26$G)',\
          r'$l_e$ ($B_0=10$G)',\
          r'$B$ ($B_0=26$G)',\
          r'$B$ ($B_0=10$G)',\
          r'$\delta$ ($B_0=26$G)',\
          r'$\delta$ ($B_0=10$G)',\
          r'$l_{ext}$ ($B_0=26$G)',\
          r'$l_{ext}$ ($B_0=10$G)',\
          r'$l_e$ ($B_0=26$G)',\
          r'$l_e$ ($B_0=10$G)',\
          r'$B$ ($B_0=26$G)',\
          r'$B$ ($B_0=10$G)',\
          r'$\delta$ ($B_0=26$G)',\
          r'$\delta$ ($B_0=10$G)',\
          r'$l_{ext}$ ($B_0=26$G)',\
          r'$l_{ext}$ ($B_0=10$G)']

listsimB=[r'$l_e$ 1,000d',\
          r'$B$ 365d',\
          r'$\delta$ 1,000d',\
          r'$l_e$',\
          r'$B$',\
          r'$\delta$',\
          r'$l_e$',\
          r'$B$',\
          r'$\delta$']


axxy=[13.5,27.5,0,3.20] #watch fpr maximum FV
numA=21 #len(listsimA)-1 #for all sims of 3C
numB=9 #len(listsimB)-1 #for all sims of PKS
#Arange=np.array([6,8,10,12])#range(6,numA+1)
#Brange=range(4,numB+1)
Arange=np.array([14,16,18,20]) #B=10G
Brange=np.array([7,8,9])

fig=plt.subplots(1,figsize=[10,11],sharey='all',num=1)
ax1=plt.subplot(2,1,2)
ax2=plt.subplot(2,1,1)

#3C273
ax1.errorbar(xs,ys,yerr=eys,fmt='D',color='grey',fillstyle='full', ms=7.0,elinewidth=1.5,capsize=4.5,\
            label='obs FV')
#CV
ax1.errorbar(np.log10(np.sqrt(30)*1e9*eV/h), 1.1625,\
             yerr=0.00,xerr=(range_fermi[1]-range_fermi[0])/2,marker='s',color='black',ms=markersize,\
             label='Fermi-LAT LC')
#PKS2155-304
ax2.errorbar(xc+np.log10(h/eV),yc,yerr=eyc,xerr=dxc,fmt='D',color='grey',fillstyle='full',ms=7.0,elinewidth=1.5,capsize=4.5,\
             label='obs FV')
ax2.errorbar(np.log10(np.sqrt(30)*1e9),1.55,yerr=0.00,xerr=(range_fermi[1]-range_fermi[0])/2,marker='s',color='black',ms=markersize,\
             label='Fermi-LAT LC')
#FV
#ax1.errorbar(np.log10(np.sqrt(30)*1e9*eV/h),1.27,yerr=0.0055,xerr=(range_fermi[1]-range_fermi[0])/2,marker='s',color='black',ms=6.5,\
#             label='Fermi-LAT LC')
#ax2.errorbar(np.log10(np.sqrt(30)*1e9),0.3459,yerr=0.0055,xerr=(range_fermi[1]-range_fermi[0])/2,marker='s',color='black',ms=6.5,\
 #            label='Fermi-LAT LC')
 #CV


#3C 273 CVs SIMS
marks.reverse()
mk, clprev, lns = len(marks)-1, [], len(lnstyles)-1
for a in Arange:
    lbl=listsimA[a-1]
    for cl in list(coldic):
        if cl in lbl.split()[0]:
            col=coldic[cl]
    if col==clprev:
        mk=mk+1
    fillstyle='full'
    if '10$G' in lbl:
        fillstyle='none'
    mk=mk-1
    marker=marks[mk]
    clprev=col
    lnst=lnstyles[lns]
    lns=lns-1
    #VAR
    y=np.array(globals()['sim'+str(a)])
    ya=np.sqrt(y**2 + (10**np.array([erro[0], erro[0], erro[0], errx[0], errx[0], errx[0], errg[0], errg[0], errg[0]])-1)**2)
    #era=np.array(globals()['er'+str(a)])
    era=np.zeros(len(y[y>0])) #era[ya>0]
    x1=energies[y>0]
    dx1=denergies[y>0]
    ya=ya[y>0]
    if 'l_{ext}' not in lbl:
        ax1.errorbar(x1,ya,yerr=0,xerr=dx1,marker=marker,\
                 color=col,markersize=markersize,fillstyle=fillstyle,linestyle=lnst,label=lbl)
    else:
        ax1.errorbar(x1,ya,yerr=0,xerr=dx1,marker=marker,\
                 color=col,markersize=markersize,fillstyle=fillstyle,linestyle=lnst,label=lbl)
#ax1.set_title('3C 273',fontsize=fntsz+4)
axxy[3]=axxy[3]-0.01
ax1.axis(axxy)
axxy[3]=axxy[3]+0.01
#leg=ax1.legend(fontsize=legsize,bbox_to_anchor=(-0.2, 1.0),title='3C 273')
leg=ax1.legend(fontsize=legsize,framealpha=0.5,loc='upper left',title='3C 273')

leg.get_title().set_fontsize(fntsz+2)
ax1.set_xlabel(r'$\log({\nu_{obs}})$ [Hz]',fontsize=fntsz+2)
#ax1.set_ylabel(r'FV',fontsize=fntsz+2)
#ax2.set_ylabel(r'FV',fontsize=fntsz+2)
ax1.set_ylabel(r'$CV$',fontsize=fntsz+2)
ax2.set_ylabel(r'$CV$',fontsize=fntsz+2)
ax1.tick_params(axis='x',labelsize=fntsz-2)
ax1.tick_params(axis='y',labelsize=fntsz-2)
ax1.tick_params(which='major', width=1.2, length=7)
ax1.tick_params(which='minor', width=1.1, length=4)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.1)
        
#PKS2155-304 CVs SIMS
mk, clprev, lns = len(marks)-1, [], len(lnstyles)-1
for b in Brange:
    lbl=listsimB[b-1]
    for cl in list(coldic):
        if cl in lbl.split()[0]:
            col=coldic[cl]
    if col==clprev:
        mk=mk+1
    mk=mk-1
    marker=marks[mk]
    lns=lns-1
    lnst=lnstyles[lns]
#FV
#    yb=np.array(globals()['simB'+str(b)])
#    erb=np.array(globals()['erB'+str(b)])
#    erb=erb[yb>0]
#    x1=energies[yb>0]
#    dx1=denergies[yb>0]
#    yb=yb[yb>0]
#VAR
    y=np.array(globals()['simB'+str(b)])
    yb=np.sqrt(y**2 + (10**np.array([erro[2], erro[2], erro[2], errx[2], errx[2], errx[2], errg[2], errg[2], errg[2]])-1)**2)
    erb=np.zeros(len(y[y>0])) #era[ya>0]
    x1=energies[y>0]
    dx1=denergies[y>0]
    yb=yb[y>0]
    
    ax2.errorbar(x1+np.log10(h/eV),yb,yerr=0,xerr=dx1,marker=marker,color=col,\
                 markersize=markersize,linestyle=lnst,label=lbl)
ax2.set_title(r'$\epsilon_{obs}$  [eV]',fontsize=fntsz,pad=28.0)
#eV axis
axxy[0]=axxy[0]+np.log10(h/eV)
axxy[1]=axxy[1]+np.log10(h/eV)
ax2.axis(axxy)
#leg=ax2.legend(fontsize=legsize,bbox_to_anchor=(-0.2,1.0),title='PKS 2155-304')
leg=ax2.legend(fontsize=legsize,framealpha=0.5,loc='upper center',title='PKS 2155-304')
leg.get_title().set_fontsize(str(fntsz+2))
ax2.tick_params(axis='x',labelbottom='off', labeltop='on', bottom='off',top='on',labelsize=fntsz-2)
ax2.tick_params(axis='y', labelsize=fntsz-2)
ax2.tick_params(which='major',bottom='on',top='on',width=1.2, length=7)
ax2.tick_params(which='minor',bottom='on',top='on', width=1.1, length=4)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(1.1)
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0)
#plt.figure(1).savefig('FVs.eps',bbox_inches='tight')
#plt.figure(1).savefig('FVs.png',bbox_inches='tight')
plt.figure(1).savefig('CVs.eps',bbox_inches='tight')
plt.figure(1).savefig('CVs.png',bbox_inches='tight')
#%%
""" Fermi PDFs"""
BINS=50
#arrays exist labeled fermiAB, obsfAB etc 
dimvec, dimhor= 2, 4
fig,ax=plt.subplots(dimvec,dimhor,figsize=[15,5],sharex=True,num=99)
index=['AA','AB','BA','BB','AD','BD','BE'] # plotted PDF diagrammes
locc=['PKS2155304/le_^0.5_3600d','PKS2155304/B_^1.0_3600d','3C273/le_^1.0_3600d_B=26G','3C273/B_^-1.0_3600d_B=26G',\
      'PKS2155304/delta_^0.25_3600d','3C273/delta_^0.2_3600d_B=26G','3C273/lext_^1.0_3600d_B=26G']
index=index+['BC','BF','BG','BH'] #additional Fermi LC for overplot
locc=locc+['3C273/B_^-1.0_3600d_B=10G','3C273/lext_^1.0_3600d_B=10G',\
       '3C273/le_^1.0_3600d_B=10G','3C273/delta_^0.2_3600d_B=10G']
lbs={'A':r'$l_e(t)$', 'B':r'$B(t)$', 'E':r'$l_{ext}(t)$','D':r'$\delta(t_{obs})$'}

#OBSRVATIONS
obs=[]
ff=open(route+'/Lightcurve_Simulation/PKS2155304/myfile_new.txt')
for f in ff:
    obs.append(f.split()[1])
globals()['obsA']=np.array(obs).astype(np.float)
obsA = obsA*3.3e-10 #correct for using wrong distribution for PKS2155-304
fermi=obsA
#srt=np.array(sorted(np.log10(fermi)))
#mn=list(srt).index(srt[srt>np.mean(srt)][0])
#n1=int(np.floor(len(fermi)*0.34))
#low, high = srt[max(mn-n1,0)] , srt[mn+n1]
cnt=0
for i in sorted(np.log10(fermi)):
    cnt=cnt+1
    if cnt>(1-0.16)*len(fermi):
        high=sorted(np.log10(fermi))[cnt]
        break
cnt=0
for i in sorted(np.log10(fermi)):
    cnt=cnt+1
    if cnt>0.16*len(fermi):
        low=sorted(np.log10(fermi))[cnt]
        break
#print(r'PKS obs '+lbs[inn[1]]+' mean= '+str(round(np.mean(np.log10(fermi)),2))+'\t median='+str(round(np.median(np.log10(fermi)),2))+\
#       ' \t range '+str(round(high-low,2)))

obs=[]
ff=open(route+'/Lightcurve_Simulation/3C273/myfile_new.txt')
for f in ff:
    obs.append(f.split()[1])
globals()['obsB']=np.array(obs).astype(np.float)
obsB = obsB*1.214e-10 #correct for using wrong distribution for PKS2155-304
fermi=obsB
#srt=np.array(sorted(np.log10(fermi)))
#mn=list(srt).index(srt[srt>np.mean(srt)][0])
#n1=int(np.floor(len(fermi)*0.34))
#low, high = srt[max(mn-n1,0)] , srt[mn+n1]
cnt=0
for i in sorted(np.log10(fermi)):
    cnt=cnt+1
    if cnt>(1-0.16)*len(fermi):
        high=sorted(np.log10(fermi))[cnt]
        break
cnt=0
for i in sorted(np.log10(fermi)):
    cnt=cnt+1
    if cnt>0.16*len(fermi):
        low=sorted(np.log10(fermi))[cnt]
        break
#print(r'3C obs '+lbs[inn[1]]+' mean='+str(round(np.mean(np.log10(fermi)),2))+'\t median='+str(round(np.median(np.log10(fermi)),2))+\
#       ' \t range \n'+str(round(high-low,2)))

for inn in index:
    fermi=[]
    lcc=locc[list(index).index(inn)]
    ff=open(route+'/Var/Results/'+lcc+'/FERMI.txt')
    for f in ff:
        fermi.append(f.split()[1])
    globals()['fermi'+inn]=np.array(fermi).astype(np.float)
    ff.close()
    
for inn in index[0:int(dimvec*dimhor-1)]:
    ab=[0,0]    
    for i in range(len(inn)):
        if inn[i]=='B':
            ab[i]+=1
        if inn[i]=='D':
            ab[i]+=2 
        if inn[i]=='E':
            ab[i]+=3
    aa=ax[ab[0],ab[1]]
    fermi=globals()['fermi'+inn]
    obsf=np.log10(globals()['obs'+inn[0]])
    color='m'
    MAX=max(max(np.histogram(np.log10(fermi),bins=BINS,density=True)[0]),\
                   max(np.histogram(obsf,bins=BINS,density=True)[0]))*1.1
    if inn[1]=='B':
        color='c'
    if inn[1]=='D':
        color='g'
    if inn[1]=='E':
        color='y'
    if inn[0]=='A':
        aa.set_title(r'$l_e-$varying',size=fntsz)
        if inn[1]!='A':
           aa.set_title(r'$B-$varying',size=fntsz)
        else:
           aa.set_ylabel('PKS 2155-304',size=fntsz)

    if inn[1]=='D':
        if inn[0]=='A':
           aa.set_title(r'$\delta-$varying',size=fntsz)

    if inn[0]=='B':
        aa.set_xlabel(r'$F_{obs}\; [erg/cm^2/s]$',size=fntsz)
        if inn[1]=='A':
            aa.set_ylabel(r'3C 273',size=fntsz)
        if inn[1]=='E':
            aa.set_title(r'$l_{ext}-$varying',size=fntsz)
    label='sim'
    
    if inn=='BB':
        fermiB=fermiBC
        MAX=1.7
    if inn=='BE':
        fermiB=fermiBF
        MAX=1.7
    if inn=='BA':
        fermiB=fermiBG
        label=r'$B_0=26~$G'
    if inn=='BD':
        fermiB=fermiBH
        label=r'$B_0=26~$G'
    if inn[0]=='B':
        aa.hist(np.log10(fermiB),bins=BINS,density=True,color='grey',alpha=0.5,label=r'$B_0=10~$G')
        label=r'$B_0=26~$G'
        cnt=0
        for i in sorted(np.log10(fermiB)):
         cnt=cnt+1
         if cnt>(1-0.16)*len(fermiB):
            high=sorted(np.log10(fermiB))[cnt]
            break
        cnt=0
        for i in sorted(np.log10(fermiB)):
         cnt=cnt+1
         if cnt>0.16*len(fermiB):
            low=sorted(np.log10(fermiB))[cnt]
            break
        objs={'A':'PKS','B':'3C'}
        print(objs[inn[0]]+r' B=26G '+lbs[inn[1]]+' mean='+str(round(np.mean(np.log10(fermiB)),2))+'\t median='+str(round(np.median(np.log10(fermiB)),2))+\
           ' \t range '+str(round(high-low,2)))
    aa.hist(np.log10(fermi),bins=BINS,density=True,color=color,alpha=0.5,label=label)
    aa.hist(obsf,bins=BINS,density=True,color='tab:orange',alpha=0.5,label='obs')
    aa.legend(handlelength=0.5,framealpha=0.5)
    aa.axis([-13.,-8.0, 0.,MAX])
#    srt=np.array(sorted(np.log10(fermi)))
#    mn=list(srt).index(srt[srt>np.mean(srt)][0])
#    n1=int(np.floor(len(fermi)*0.34))
#    low, high = srt[mn-n1] , srt[mn+n1]
    cnt=0
    for i in sorted(np.log10(fermi)):
     cnt=cnt+1
     if cnt>(1-0.16)*len(fermi):
        high=sorted(np.log10(fermi))[cnt]
        break
    cnt=0
    for i in sorted(np.log10(fermi)):
     cnt=cnt+1
     if cnt>0.16*len(fermi):
        low=sorted(np.log10(fermi))[cnt]
        break
    objs={'A':'PKS','B':'3C'}
    print(objs[inn[0]]+r'sim '+lbs[inn[1]]+' mean='+str(round(np.mean(np.log10(fermi)),2))+'\t median='+str(round(np.median(np.log10(fermi)),2))+\
           ' \t range '+str(round(high-low,2)))


ax[0,3].remove()
fig.suptitle('All PDFs', fontsize=16,y=1.00)
fig.savefig('AllPDFs.png',bbox_inches='tight')
fig.savefig('AllPDFs.eps',bbox_inches='tight')
#%%
""" Multiple DCFs"""

#arrays exist labeled fermiAB, obsfAB etc 
dimvec, dimhor= 2, 4
fig,ax=plt.subplots(dimvec,dimhor,figsize=[15,5],sharex=True,num=98)
bands=[r'$\gamma-$rays',r'$X-$rays','$J-$band'] #more to less energetic
lnbd=len(bands)
interest=['AA','BE']

index=['AA','AB','BA','BB','AD','BD','BE'] # plotted PDF diagrammes
locc=['PKS2155304/le_^0.5_3600d','PKS2155304/B_^1.0_3600d','3C273/le_^1.0_3600d_B=10G','3C273/B_^-1.0_3600d_B=10G',\
      'PKS2155304/delta_^0.25_3600d','3C273/delta_^0.2_3600d_B=10G','3C273/lext_^1.0_3600d_B=10G']


for bnd in range(lnbd):
    if bnd==0:
        band1=bands[-1]
        band2=bands[0]
        bands[-1]=bands[0] #change sequence of comparison for visual reasons  less -> more energetic
        bands[0]=band1
        
for inn in index:
    tlags, DCFs = [], []
    lcc=locc[list(index).index(inn)]
    ff=open(route+'/Var/Results/'+lcc+'/DCFs.txt')
    for f in ff:
        tlags.append(f.split()[0])
        DCFs.append(f.split()[1])
    globals()['tlags'+inn]=np.array(tlags).astype(np.float)
    globals()['DCFs'+inn]=np.array(DCFs).astype(np.float)
    ff.close()
    
for inn in index:
    ab=[0,0]
    tlags=[]
    DCFs=[]
    for i in range(len(inn)):
        if inn[i]=='B':
            ab[i]+=1
        if inn[i]=='D':
            ab[i]+=2 
        if inn[i]=='E':
            ab[i]+=3
    aa=ax[ab[0],ab[1]]
   
    if inn[0]=='A':
        aa.set_title(r'$l_e-$varying',size=fntsz)
        if inn[1]!='A':
           aa.set_title(r'$B-$varying',size=fntsz)
        else:
           aa.set_ylabel('PKS 2155-304',size=fntsz)


    if inn[1]=='D':
        if inn[0]=='A':
           aa.set_title(r'$\delta-$varying',size=fntsz)    
    if inn[0]=='B':
        aa.set_xlabel(r'$\tau_{obs}\; $ [days]',size=fntsz)
        if inn[1]=='A':
            aa.set_ylabel(r'3C 273',size=fntsz)
        if inn[1]=='E':
            aa.set_title(r'$l_{ext}-$varying',size=fntsz)    
    
    if inn[1]!='A':
        aa.tick_params(axis='y',labelleft='off',left='off',direction='in',labelsize=10,pad=3)

    #PLOTTING
    yy=globals()['DCFs'+inn]
    xx=globals()['tlags'+inn]
    ll=int(len(xx)/3)
    lentau=int((ll-1)/2)
    DCFs=[yy[0:ll],yy[ll:2*ll],yy[2*ll::]]
    tlags=[xx[0:ll],xx[ll:2*ll],xx[2*ll::]]
    
    for bnd in range(lnbd):
            aa.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+' - '+bands[bnd-1])
    aa.plot(np.linspace(-lentau,lentau,50),0*np.linspace(-1,1,50),'k--',linewidth=0.5)
    aa.plot(0*np.linspace(-1,1,50),np.linspace(-1.,1,50),'k--',linewidth=0.5)
    #Plot ACF
    tacf, ACF = [] , []
    if inn[0]=='B':
        acf=open(route+'/Lightcurve_Simulation/3C273/ACFs_3C273.txt')
    else:
        acf=open(route+'/Lightcurve_Simulation/PKS2155304/ACFs_PKS2155304.txt')
    for f in acf:
            tacf.append(f.split()[0])
            ACF.append(f.split()[1])
    acf.close()
    tacf = np.array(tacf).astype(np.float)
    ACF = np.array(ACF).astype(np.float)
    aa.plot(tacf,ACF,'ks',label='$\gamma-$rays - $\gamma-$rays',markersize=1.5)    
    aa.axis([-lentau,lentau,-1,1.05])
        #calculate timelags for gamma-rays - O/IR
#    if inn in interest:
    bnd=0  #gamma  - O/IR
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
    print('Most strongly detected timelag (peak) at:\t'+str(tlag_extra*24)+' hours')
#    aa.legend([r'$\tau^{o\gamma}_{obs}='+str(round(tlag_extra*24,1))+'$ h'],fontsize='large',loc='lower left') 
    
    #pospeak=list(abs(DCFs[0])).index(max(abs(DCFs[0])))
    #tpk=tlags[0][pospeak]
    #pk=DCFs[0][pospeak] #peak of SMARTS fermi DCF
    #if abs(tpk)<2.0:
    #        # location for the zoomed portion  in fig 20
    #        ax34= plt.axes([0.65,0.16, 0.25, 0.25])# units unity for axis scales
    #        ax34.set_xlabel(r'$\tau_{ obs}$ [days]',size=10)
    #        ax34.tick_params(axis='y',labelright='on',labelleft='off',right='on',direction='in',labelsize=10,pad=3)
    #        ax34.tick_params(axis='x',labelbottom='on',bottom='on',direction='in',labelsize=10,pad=-12)
    #        # insert the zoomed figure
    #        plt.setp(ax34)
    #        for bnd in range(lnbd):
    #            plt.plot(tlags[bnd],DCFs[bnd],'-',label=bands[bnd]+'-'+bands[bnd-1])
    #        if pk>0:
    #            ax34.axis([tpk-2.5,tpk+2.5,pk*0.85,min(1.0,pk*1.05)])
    #        else:
    #            ax34.axis([tpk-2.5,tpk+2.5,max(pk*1.05,-1.0),pk*0.85])
    #        from mpl_toolkits.axes_grid1.inset_locator import mark_inset
    #        mark_inset(aa, ax34, loc1=1, loc2=2, fc="none", ec="0.5")
ax[1,3].add_artist(ax[1,3].legend(fontsize=legsize,bbox_to_anchor=(1.0,2.0)))
#plt.legend([r'$\tau_{obs}^{o\gamma}='+str(round(tlag_extra*24,1))+'$ h'],fontsize='large',loc='lower left') 
ax[0,3].remove()
plt.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12,hspace=0.07,wspace=0.02)
plt.show()
fig.suptitle('All DCFs', fontsize=16,y=1.05)
fig.savefig('AllDCFs.png',bbox_inches='tight')
fig.savefig('AllDCFs.eps',bbox_inches='tight')
#%%
""" Parameter PDFs"""
fntsz=14
#arrays exist labeled fermiAB, obsfAB etc 

dimvec, dimhor= 2,4
fig,ax=plt.subplots(dimvec,dimhor,figsize=[15,6],num=97)
index=['AA','AB','BA','BB','AD','BD', 'BE'] # plotted PDF diagrammes
locc=['PKS2155304/le_^0.5_3600d','PKS2155304/B_^1.0_3600d','3C273/le_^1.0_3600d_B=26G','3C273/B_^-1.0_3600d_B=26G',\
      'PKS2155304/delta_^0.25_3600d','3C273/delta_^0.2_3600d_B=26G','3C273/lext_^1.0_3600d_B=26G']
index=index+['BG','BC','BH','BF'] #additional Fermi LC for overplot
locc=locc+['3C273/le_^1.0_3600d_B=10G','3C273/B_^-1.0_3600d_B=10G',\
       '3C273/delta_^0.2_3600d_B=10G','3C273/lext_^1.0_3600d_B=10G']
lbs={'A':'l_e', 'B':'B', 'E':'l_{ext}','D':r'\theta_{obs}',
     'G':'l_e', 'C':'B', 'F':'l_{ext}','H':r'\theta_{obs}'}

for inn in index:
    filevar,filetime, globals()['var'+inn]=[] , []  , []
    if inn[1]!='D' and inn[1]!='H':
        lcc=locc[list(index).index(inn)]
        root=route+'/Var/Results/'+lcc+'/'
        nmlist=os.listdir(root)
        names=[]
        cnt=0
        for nm in nmlist:
            if 'fort_' in nm and '.85' in nm and 'steady' not in nm:
                cnt=cnt+1 #.85, .89 files
        for ia in range(cnt):
            names.append('fort_'+str(ia))
        if len(names)==0:
            names.append('fort')
        
        for name in names:
            #read fort.55
            name55=name+'.55' #name55='fakeTC.txt' #read analytically all the input fake variables
            varfile=open(root+name55,'r')
            for y in varfile:
                x=y.split()
                z=x[0]
                w=x[1].replace('\n','')
                filetime.append(z)
                filevar.append(w)
            varfile.close()
        filetime=np.array(filetime).astype(np.float) 
        addt=0
        #    for it in range(len(filetime)):
        #        if filetime[it]<1e-5 and it>1:
        #            addt=filetime[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
        #        # the time of steady state (conjoint part of two TC)
        #        filetime[it]=filetime[it]+addt
        filevar=np.array(filevar).astype(np.float)
        
        for tt in range(0,len(filetime)-1):
            tm=filetime[tt]
            if tm!=filetime[tt-1] and tm>0.0001:
                globals()['var'+inn].append(filevar[tt])        
        globals()['var'+inn]=np.array(globals()['var'+inn]).astype(np.float) 
        

    else:
        var=[]
        ff=open(route+'/Var/Results/'+locc[index.index(inn)]+'/fakeTC.txt')
        for f in ff:
            globals()['var'+inn].append(f.split()[1])
        globals()['var'+inn]=np.array(globals()['var'+inn]).astype(np.float)
        
        
        
for inn in index[0:int(dimvec*dimhor-1)]:
    ab=[0,0]
    for i in range(len(inn)):
        if inn[i]=='B':
            ab[i]+=1
        if inn[i]=='D':
            ab[i]+=2
        if inn[i]=='E':
            ab[i]+=3
    aa=ax[ab[0],ab[1]]
    var=globals()['var'+inn]        
    
    if inn[0]=='A':
        ledd=1.2*Ledd
        R=4.9e16
        Gamma=30.0
        ldisk=1.0
        B0=0.08
        le0= 10**-4.931
        delta0= 30
    else:
        ledd=8.8*Ledd
        R=6.64e15
        Gamma=35.4
        Rext=6.3e17
        ldisk=2e46
        le0=  10**-1.95
        B0=26
        if '10' in locc[list(index).index(inn)]:
            B0 = 10
        delta0 = 10
        lext0 = 0.033
        
        
    color='m'
    tr=4*np.pi*me*c**3*R/3/sigmaT*2/3*Gamma**2/ledd
    lum=tr*var
    title=r'$L_{e,jet} /L_{edd}$'
    var_st  =  le0
    if inn[1]=='B':
        var=10**var
        tr=c*R**2*1/4*Gamma**2/ledd
        lum=tr*var**2
        color='c'
        title=r'$L_{B,jet} /L_{edd}$'
        var_st = B0
    if inn[1]=='E':
        tr=1/(sigmaT*Gamma**2*R/4/np.pi/me/c**3/Rext**2)/ledd
        lum=tr*var
        color='y'
        title=r'$L_{ext}/L_{edd}$'
        var_st  = lext0
    if inn[1]=='D':
        color='g'
        label=r'$\theta$'
        beta=np.sqrt(1-1/Gamma**2)
        gd=1/beta*(1-1/Gamma/var)
        theta=np.zeros(len(var))
        for i in range(0,len(theta)-1):
            if gd[i]>1.0:
                theta[i]=np.arccos(1/beta*(1-2/Gamma**2))
            else:
                theta[i]=np.arccos(gd[i])/np.pi*180
        delta=var
        var=theta[theta>0]   
        gd0=1/beta*(1-1/Gamma/delta0)
        var_st=np.arccos(gd0)/np.pi*180
        
    if inn[0]=='A':
        if inn[1]=='A':
           aa.set_ylabel('PKS 2155-304',size=fntsz+2)

    #if inn[1]=='D':
    #    if inn[0]=='A':
    #       aa.set_title(r'$\delta-$varying',size=fntsz)    
    if inn[1]=='A':
            aa.set_xlabel(r'$\log{l_e}$',size=fntsz)
    if inn[1]=='B':
            aa.set_xlabel(r'$\log{B}$ [G]',size=fntsz)
    if inn[1]=='D':
            aa.set_xlabel(r'$\theta_{obs}$ [deg]',size=fntsz)
    if inn[1]=='E':
           aa.set_xlabel(r'$\log{l_{ext}}$',size=fntsz)
    if inn=='BA':
            aa.set_ylabel(r'3C 273',size=fntsz+2)
            
    if inn[1]!='D':
        #remove spike
        spk=np.histogram(np.log10(var),bins=50)[1][1]
        print('\n Spike bin neglected:\t',spk)
        var=var[np.log10(var)>spk]
        
        #if inn=='BB':
           #varB=10**varBC[varBC>np.histogram(varBC,bins=50)[1][1]-0.01]
           #lumBC=np.log10((varB)**2*c*R**2*1/4*Gamma**2/ledd)
           #aa.hist(np.log10(varB),bins=50,density=True,color='limegreen',alpha=0.5,\
           #        label=r'$\log<'+lbs[inn[1]]+'>='+str(round(np.log10(np.mean(varB)),1))+\
           #        '$ \n $1\sigma$ range $'+str(round(2*np.std(np.log10(varB)),1))+'$')
        alpha=0.8
        if inn[1]=='A':
            alpha=0.5
        aa.hist(np.log10(var),bins=50,density=True,color=color,alpha=alpha,\
                #label=r'$\log<'+lbs[inn[1]]+'>='+str(round(np.log10(np.mean(var)),1))+\
                label=r'$<\log{'+lbs[inn[1]]+'}>='+str(round(np.mean(np.log10(var)),1))+\
           '$ \n $1\sigma$ range $'+str(round(2*np.std(np.log10(var)),1))+' $')
        aa.vlines(np.log10(var_st),0.0, max(np.histogram(np.log10(var),bins=50,density=True)[0]),linestyles='dashed', color ='grey', alpha=0.5,label='steady')
        if inn[1]=='BB':
            aa.hist(np.log10(varBC),bins=50,density=True,color=color,alpha=alpha,\
                    #label=r'$\log<'+lbs[inn[1]]+'>='+str(round(np.log10(np.mean(var)),1))+\
                    label=r'$<\log{'+lbs[inn[1]]+'}>='+str(round(np.mean(np.log10(varBC)),1))+\
               '$ \n $1\sigma$ range $'+str(round(2*np.std(np.log10(varBC)),1))+' $')
        ax34= aa.twiny()
        ax34.set_xlabel(title,size=fntsz,labelpad=10)
        if inn[1]!='B':
            ax34.axis([ aa.axis()[0]+np.log10(tr), aa.axis()[1]+np.log10(tr), aa.axis()[2], aa.axis()[3] ])
        else:
            ax34.axis( [2*aa.axis()[0]+np.log10(tr), 2*aa.axis()[1]+np.log10(tr),aa.axis()[2],aa.axis()[3]])
            ax34.xaxis.set_minor_locator(AutoMinorLocator()) 
    else:
        aa.hist(var,bins=50,density=True,color=color,alpha=alpha,\
                label=r'$<'+lbs[inn[1]]+'>='+str(round(np.mean(var),1))+\
           '$ \n $1\sigma$ range $'+str(round(2*np.std(var),1))+' $')
        aa.vlines(var_st,0.0, max(np.histogram(var,bins=50,density=True)[0]), linestyles='dashed', color ='grey', alpha=0.5)


    ax2=list(ax34.axis())
    
#    if np.log10(ldisk/ledd)>ax2[0] and np.log10(ldisk/ledd)<ax2[1] and inn[0]=='B':
#        lnsp=np.linspace(-0.1,2.0)
#        ax34.plot(np.log10(ldisk/ledd)*np.ones(50),lnsp,'--',color='grey',lw=2.0,label=r'$L_{disk}$')
    #if inn not in ['AA','AB','BE']:
    aa.legend(framealpha=0.25,fontsize='large',loc='upper right',handlelength=0.0)
    #else:
    #    aa.legend(framealpha=0.25,fontsize='large',loc='upper right',handlelength=0.0)#markerfirst=False)
    aa.tick_params(axis='x', labelsize=fntsz-2)
    ax34.tick_params(axis='x', labelsize=fntsz-2)
    aa.tick_params(axis='y',labelsize=fntsz-2,pad=-0.5)
    aa.xaxis.set_minor_locator(AutoMinorLocator())


ax[0,3].remove()
plt.subplots_adjust(left=0.1,right=0.95,top=0.85,bottom=0.12,hspace=0.85,wspace=0.2)
#ax[2,0].remove()
fig.savefig('ParamPDFs.png',bbox_inches='tight')
fig.savefig('ParamPDFs.eps',bbox_inches='tight')

#%%
"""  CIs """
fntsz=12
#arrays exist labeled fermiAB, obsfAB etc 
dimvec, dimhor= 2,3
fig17,ax17=plt.subplots(dimvec,dimhor,figsize=[15,5],num=96)
index=['AA','AB','AD','BA','BB','BE'] # plotted PDF diagrammes
locc=['PKS2155304/le_^0.5_3600d','PKS2155304/B_^1.0_3600d','PKS2155304/delta_^0.25_3600d',\
      '3C273/le_^1.0_3600d_B=26G','3C273/B_^-1.0_3600d_B=26G','3C273/lext_^1.0_3600d_B=26G']
#index.append('BC') #additional Fermi LC for overplot
#locc.append('3C273/B_^-1_3620d_B=10G')
lbs={'A':r'$l_e(t)$', 'B':r'$B(t)$', 'E':r'$l_{ext}(t)$','D':r'$\delta(t_{obs})$'}
CONTOUR = 'on'

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
    return X,Y,Z


for inn in index:
    lcc=locc[list(index).index(inn)]
    root=route+'/Var/Results/'+lcc+'/'
    filex,filey,globals()['filterB_'+inn], globals()['ciB_'+inn] = [] , []  , [] , []
    varfile=open(root+'CI.txt','r')
    for y in varfile:
        x=y.split()
        t=x[0]
        fa=x[1]
        fb=x[2]
        ci=x[3].replace('\n','')
        filex.append(fb)
        filey.append(ci)
    varfile.close()
    globals()['filterB_'+inn]=np.array(filex).astype(np.float) 
    globals()['ciB_'+inn]=np.array(filey).astype(np.float) 
    filterBB=globals()['filterB_'+inn]
    ciB=globals()['ciB_'+inn]     
    smvci, cilbl='J' , 'B-J'
    
    #CORRECTIONS TO SCALES
    if inn[1]=='E':
        filterBB=filterBB[ciB>0.5]
        ciB=ciB[ciB>0.5]
    if inn[1]=='B':
        ciB=ciB[filterBB<20]
        filterBB=filterBB[filterBB<20]
    
        
    smci=10**filterBB #select x axis from SMARTS bands and make it color
    ab=[0,0]
    for i in range(len(inn)):
        if inn[i]=='B':
            ab[i]+=1
        if inn[i]=='E' or inn[i]=='D':
            ab[i]+=2
    aa=ax17[ab[0],ab[1]]
    
    
    color='m'
    if inn[1]=='B':
        color='c'
    if inn[1]=='E':
        color='y'
    if inn[1]=='D':
        color='g'
    #TITLE
    if inn[0]=='A' and inn[1]=='A':
        aa.set_title(r'PKS 2155-304',size=fntsz+2)
    if inn[0]=='B' and inn[1]=='A':
        aa.set_title(r'3C 273',size=fntsz+2)
        
    #LABELS
    if inn[1]=='A':   
        aa.set_ylabel(r'$'+cilbl+'$',labelpad=-0.5,size=fntsz+1) #oldci
    
    aa.set_xlabel(r'$'+smvci+'$',size=fntsz+1) #oldci
    aa.plot(np.log10(smci),ciB,'.',color=color,ms=1.5,zorder=0)
    if CONTOUR=='on':
            lims=[]
            if inn[1]=='E':
                lims=[10,15,min(ciB)-0.5,0.9]
            if inn[1]=='D':
                colcont='orange'
            else:
                colcont='black'
            X,Y,Z=surface_density(np.log10(smci),ciB,lims=lims) #oldci
            aa.contour(X,Y,Z,colors=colcont,linestyles=['-.','--'],linewidths=1.5,levels=[(1- 0.956)*np.max(Z),(1-0.682)*np.max(Z)])
    #aa.plot(np.log10(smci)-np.mean(np.log10(smci)),ciB,'c-',linewidth=0.6) #LINES
    aa.legend([lbs[inn[0]]],handlelength=0.5,framealpha=0.5)
    a1,a2,b1,b2=[aa.get_xlim()[0],aa.get_xlim()[1],aa.get_ylim()[0],aa.get_ylim()[1]]
    ae,be,ce,de = (1-np.sign(a1)*0.05)*a1, (1+np.sign(a2)*0.05)*a2, (1-np.sign(b1)*0.05)*b1, (1+np.sign(b2)*0.05)*b2
    axcorr= [ ae, be ,ce ,de]
    aa.axis(axcorr)
    aa.invert_xaxis()
    aa.xaxis.set_minor_locator(AutoMinorLocator())
    aa.yaxis.set_minor_locator(AutoMinorLocator())
    aa.tick_params(axis='y',labelsize=fntsz-3)
    aa.tick_params(axis='x',labelsize=fntsz-3)

#ax[2,0].remove()
plt.subplots_adjust(left=0.3,right=0.95,top=0.95,bottom=0.12,hspace=0.3,wspace=0.35)
fig17.savefig('allCIs.png',bbox_inches='tight')
fig17.savefig('allCIs.eps',bbox_inches='tight')
#%%


#%%
#""" Luminosities """
#fntsz=15
##arrays exist labeled fermiAB, obsfAB etc 
#
#dimvec, dimhor= 3,2
#fig,ax=plt.subplots(dimvec,dimhor,figsize=[7.5,6],num=97)
#index=['AA','AB','BA','BB','AC','BC'] # plotted PDF diagrammes
#locc=['PKS2155304/le_^0.5_3600d','PKS2155304/B_^1.0_3600d','3C273/le_^1.0_3600d_B=10G','3C273/B_^-1.0_3600d_B=10G',\
#      '3C273/lext_^1.0_3600d_B=10G']
##index=index+['BC','BF','BG','BH'] #additional Fermi LC for overplot
##locc=locc+['3C273/B_^-1_3620d_B=10G','3C273/lext_^1_3620d_B=10G',\
##       '3C273/le_^1_3620d_B=10G','3C273/delta_^0.2_3620d_B=10G']
#lbs={'A':'l_e', 'B':'B', 'E':'l_{ext}',
#     'G':'l_e', 'C':'B', 'F':'l_{ext}'}    
#
#for inn in index:
#    filevar,filetime, time, var = [], [] , [],[]
#    lcc=locc[list(index).index(inn)]
#    root=route+'/Var/Results/'+lcc+'/'
#    nmlist=os.listdir(root)
#    names=[]
#    cnt=0
#    for nm in nmlist:
#        if 'fort_' in nm and '.85' in nm and 'steady' not in nm:
#            cnt=cnt+1 #.85, .89 files
#    for ia in range(cnt):
#        names.append('fort_'+str(ia))
#    if len(names)==0:
#        names.append('fort')
#    
#    for name in names:
#        #read fort.55
#        name55=name+'.55' #name55='fakeTC.txt' #read analytically all the input fake variables
#        varfile=open(root+name55,'r')
#        for y in varfile:
#            x=y.split()
#            z=x[0]
#            w=x[1].replace('\n','')
#            filetime.append(z)
#            filevar.append(w)
#        varfile.close()
#    filetime=np.array(filetime).astype(np.float) 
#    addt=0
#    for it in range(len(filetime)):
#            if filetime[it]<1e-5 and it>1:
#                addt=filetime[it-1]-1#-1 comes from 'conshort' parameter in my_analysis_v2.py  and matches to 
#            # the time of steady state (conjoint part of two TC)
#            filetime[it]=filetime[it]+addt
#    filevar=np.array(filevar).astype(np.float)
#    
#    for tt in range(0,len(filetime)-1):
#        tm=filetime[tt]
#        if tm!=filetime[tt-1] and tm>0.0001:
#            time.append(filetime[tt])        
#            var.append(filevar[tt])     
#    time = np.array(time).astype(np.float) 
#    var = np.array(var).astype(np.float) 
#    lgamma,le,Le,Lgamma,tobs0,tobs = [], [],[],[],[], []
#    lcc=locc[list(index).index(inn)]
#    ff=open(route+'/Var/Results/'+lcc+'/LUMS.txt')
#    for f in ff:
#        tobs0.append(f.split()[0])
#        lgamma.append(f.split()[1])
#        le.append(f.split()[2])
#    tobs.append(float(tobs0[0]))
#    Le.append(float(le[0]))
#    Lgamma.append(float(lgamma[0]))
#    tobs0=np.array(tobs0).astype(np.float)
#    
#    for it in range(0,len(tobs0)):
#        if tobs0[it]-tobs0[it-1]>0.001:
#            tobs.append(tobs0[it])
#            Lgamma.append(lgamma[it])
#            Le.append(le[it])
#    tobs = np.array(tobs).astype(np.float)
#    globals()['tobs'+inn], globals()['Lgamma'+inn] ,  globals()['Le'+inn] = tobs ,np.array(Lgamma).astype(np.float) ,np.array(Le).astype(np.float)
#    #interpolate for input values since some may miss for the resulted arrays (the last ones) or t=0.0
#    modvar = inpl.interp1d(time,var) 
#    tcr = 0.63058
#    if inn[0]=='B': tcr=0.25635
#    newtime = tobs/tcr+0.005
#    globals()['time'+inn] = newtime
#    globals()['var'+inn] = modvar(newtime)
#    ff.close()
#        
#for inn in index:
#    ab=[0,0]
#    for i in range(len(inn)):
#        if inn[i]=='B':
#            ab[i]+=1
#        if inn[1]=='C':
#            ab[1]=2
#    aa=ax[ab[1],ab[0]]
#    
#    color='m'
#    
#    if inn[0]=='A':
#        le0=10**-4.93
#        B0=0.08
#        R=4.9e16
#        Gamma=30.0
#        delta=30.0
#        ldisk=1.0
#        gmin, gmax , p = 10**3.61, 10**5.55, 2.71
#        redsh,  Dist=  0.116 , 543*pc*10**6
#    else:
#        le0=10**-1.95
#        B0=26.0
#        R=6.64e15
#        Gamma=35.4
#        delta=10.0
#        Rext=6.3e17
#        ldisk=2e46
#        gmin,gmax,p =  10**0.01, 10**3.475, 2.15
#        redsh, Dist= 0.158 , 749*pc*10**6
#    tcross=R/c
#
#    var=globals()['var'+inn]      
#    Lgamma = globals()['Lgamma'+inn]/delta**4*2/3*Gamma**2
#    Le = globals()['Le'+inn]*2/3*Gamma**2
#    
#    Le0=4*np.pi*me*c**3*R/sigmaT*le0*2/3*Gamma**2
#    Lb0=3/2*B0**2*R**2*c*1/2*Gamma**2
#    fact = (-p+1)/(-p+2)*(gmax**(-p+2)-gmin**(-p+2))/ (gmax**(-p+1)-gmin**(-p+1))
#    Lp0=3/4*(mp/me)/fact*Le0
#    Lj0=Le+Lb0+Lgamma+Lp0
#    if inn[1]=='A':       
#        norm = 4*np.pi*me*c**3*R/3/sigmaT*2/3*Gamma**2
#        lelum = norm*var
#        Lj = Lj0-Le0+lelum
#        aa.plot(np.log10(Lj),np.log10(lelum/Lj),'.',color=color,markersize=4)
#        title=r'$L_{\rm e,jet}/\,L_{\rm j}$'
#    if inn[1]=='B':
#        color='c'
#        var=10**var
#        norm =c*R**2*1/4*Gamma**2
#        Lb=norm*var**2
#        Lj = Lj0-Lb0+Lb
#        aa.plot(np.log10(Lj), np.log10(Lb/Lj),'.',color=color,markersize=4)
#        title=r'$L_{\rm B,jet}/\,L_{\rm j}$'
#    if inn[1]=='C':
#        norm = 4*np.pi*me*c**3*R/3/sigmaT*2/3*Gamma**2
#        lelum = norm*var
#        Lj = Lj0-Le0+lelum
#        aa.plot(np.log10(Lj),np.log10(Lgamma/Le),'.',color=color,markersize=4)
#        title=r'$L_{\rm \gamma,jet}/\,L_{\rm e,jet}$'
#        aa.set_xlabel(r'$\log{L_{\rm j}}$',size=fntsz)
#        
#    if inn[0]=='A':
#        aa.set_ylabel(title,size=fntsz)
#        
#    if inn=='AA':
#           aa.set_title('PKS 2155-304',size=fntsz+2)
#    if inn=='BA':
#           aa.set_title('3C 273',size=fntsz+2)
#
#    #remove spike
#    spk=np.histogram(np.log10(var),bins=50)[1][1]
#    var=var[np.log10(var)>spk]
#    
#    aa.tick_params(axis='x', labelsize=fntsz-2)
#    aa.tick_params(axis='y',labelsize=fntsz-2,pad=-0.5)
#    aa.xaxis.set_minor_locator(AutoMinorLocator())
#    print(inn)
#
#plt.subplots_adjust(left=0.1,right=0.95,top=0.85,bottom=0.12,hspace=0.50,wspace=0.2)
#fig.savefig('All_lums.png',bbox_inches='tight')
#fig.savefig('All_lums.eps',bbox_inches='tight')
#%%
#""" Steady state plotter"""
#
##Choose le simulation folders for getting the steady state
#routess= [ '/home/markos/Desktop/BlaVar/Var/Results/PKS2155304'+'/B_365_new/',\
#          '/home/markos/Desktop/BlaVar/Var/Results/3C273'+'/le_1000_B=10G/']
#for rt in routess :
#    route=rt
#    name='fort'
#    """Selection of Object and Features"""
#    objects=['3C273','3C279','PKS2155304']
#    namevar='le'
#    obj=rt.split('/')[-3]
#    oo=objects.index(obj)
#    
#    
#    units={'le':'','B':'G','lext':'','gmax':'','theta':'^o','Gamma':''}
#    nunits=units[namevar]
#    
#    #Bibliograpgical Data for objects
#    Dists=[749,3115,543] #Mpc Luminosity Distances
#    redsh=[0.158,0.536,0.0] #redshifts
#    Rbb=[6.3e17,5e17,10e12] #Radii of the BB radiation zone
#    SSC=['no','no','on']
#    BBextra=['on','on','no']
#    #Second BB into the model
#    secondBB=BBextra[oo] #add an additional BB on spectrum (usually DISK) but that is not involved in external photon fields of code
#    Gamma2=[1,1,1][oo] #Deboosting/Boosting for observation
#    factbb=[10,500,'-'][oo] #(L_disk-L_BLR)/L_BLR
#    Tbb2=[11000,20000,'-'][oo] #T_disk
#    Rbb2=[10**15.89,10**15.62,'-'][oo] #disk radii
#    
#    
#    """Data Reduction"""
#    nu_cutoff=8 #energy in mec^2 to cutoff plotting 
#    HEIGHT=1 #height of y-axis (threshold) in logscale measured below the minimum value of flux
#    THRESHOLD=20 #orders of magnituded for data reduction from maximum value
#    
#    limSEDs=500 #select nmax-n1 less than this limit to plot multiple SEDs
#    Multicolor='on' #Plot data points with different color for every interval of 2 years
#    diconec={'le':'m','B':'c','lext':'y','gmax':'m'}
#    onecolor=diconec[namevar] #'no' for mixed color and 1.0 transparancy
#    SED='on' #'no', 'on'(overplot), 'con' (contour)  vFv Diagram [CAUTION!]only up to *limSEDs* timepoints overlayed
#   
#    #x,y labels/units
#    xlab='v_{obs}'
#    xoptions=['mec2','Hz','GeV']
#    yoptions=['le','vLv','vFv','radio','Fv','Xrays','gamma'] #radio -> Jy / Xrays -> counts /sec / Gamma -> counts/hour
#    yunits={'le':' ','vLv':'erg/s','vFv':'erg/s/cm^2','radio':'Jy',\
#       'Fv':'erg/s/cm^2/Hz','Xrays':'counts/s/cm^2','gamma':'counts/hour/cm^2'}
#    #Select here from above
#    xaxis=xoptions[1]# CAUTION! indexing of the list begins from zero
#    yaxis=yoptions[2]# CAUTION! indexing of the list begins from zero
#    plt.close('all')
#    
#    #Multicolor Points 
#    routef='/home/markos/Desktop/BlaVar/Var/'
#    data_name=routef+obj+'/'+obj+".ascii"
#    dtpts=[]
#    if Multicolor=='on':
#       dtpts.append(routef+obj+'/'+obj+'_08-10.ascii')
#       dtpts.append(routef+obj+'/'+obj+'_10-12.ascii')
#       dtpts.append(routef+obj+'/'+obj+'_12-14.ascii')
#       dtpts.append(routef+obj+'/'+obj+'_14-16.ascii')
#       legend_pts=['2008-2010','2010-2012','2012-2014','2014-2018']
#       colors =['r','c','tab:orange','limegreen']
#       form=['o','.','+','s']
#       #dtpts.append(obj+'_16-18.ascii')
#    
#    
#    """%%%%%%%%%%%%%%%%%%    END OF SETTINGS   %%%%%%%%%%%%%%"""
#    #
#    #
#    #
#    

#    
#    
#    """INPUT Observable Parameters"""
#    cn4=open(routef+obj+'/code_new4.inp','r')
#    lines=cn4.readlines()
#    cn4.close()
#    
#    tend=float(lines[0].split()[-1].replace('\n',''))
#    nsteps=float(lines[0].split()[-2])
#    R=float(lines[2].split()[0].replace('d','e')) #dimension of source in cm
#    B=float(lines[2].split()[1]) #dimension of source in cm
#    p=float(lines[4].split()[3])
#    loggmin=float(lines[4].split()[1])
#    loggmax=float(lines[4].split()[2])
#    logle=float(lines[4].split()[4])
#    delta=float(lines[7].split()[0]) #set zero if used the above for computation
#    Gamma=float(lines[8].split()[0])
#    T=float(lines[5].split()[1])
#    lext=float(lines[5].split()[2])
#    D=Dists[oo]*10**6*pc #"Distance of Source in Mpc"
#    zsh=redsh[oo] #redshift
#    dilate=1-((((1+zsh)**2-1)/((1+zsh)**2+1))**2)**2 #dilation factor, https://www.scirp.org/journal/paperinformation.aspx?paperid=65598
#    SSConly=SSC[oo] #if SSC or EC modelling is used
#    tcross=R/c/3600/24 #days jet frame
#    tcr=tcross*dilate/delta
#    
#    if SSConly=='no':
#        def BB(x,T,lext,Rbb):
#            Ivbb=3*x+np.log10((2*me**3*c**4/h**2)/(np.exp(10**x*me*c**2/kb/T)-1)+10**-100)
#            opacity=(lext*me*c**2/sigmaT/R)/aBB/T**4
#            #make it dimensionelss / scaling
#            Ibb=x+np.log10(me*c**2/h)+Ivbb+np.log10(opacity*4*np.pi*Rbb**2)+np.log10(sigmaT/(4*np.pi*R*me*c**3))
#            return Ibb
#    """
#    %
#    %
#    %
#    """
#    
#    """Observational Data"""
#    points=ascii.read(data_name)
#    Cunitx=0 #4 Hz 
#    Cunity=0 #for erg/cm^2/s
#    v_pts=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
#    vFv_pts=np.array([points["col3"]])[0]-np.ones(1)*Cunity
#    #errorv=abs(points["col4"])
#    errorvFv=np.array(abs(points["col4"]))
#    #average measurement errors
#    erro=np.mean(errorvFv[v_pts<np.log10(10**15.2)][v_pts[v_pts<np.log10(10**15.2)]>np.log10(10**14.0)])
#    errx=np.mean(errorvFv[v_pts<np.log10(range_xhard[1])][v_pts[v_pts<np.log10(range_xhard[1])]>np.log10(range_xsoft[0])])
#    errg=np.mean(errorvFv[v_pts<np.log10(range_fermi[1])][v_pts[v_pts<np.log10(range_fermi[1])]>np.log10(range_fermi[0])])
#    if Multicolor=='on':
#      for i in range(len(dtpts)):
#        points=ascii.read(dtpts[i])
#        globals()['v_pts'+str(i)]=np.array([points["col1"]])[0]-np.ones(1)*Cunitx
#        globals()['vFv_pts'+str(i)]=np.array([points["col3"]])[0]-np.ones(1)*Cunity
#        globals()['errorv'+str(i)]=0.009
#        globals()['errorvFv'+str(i)]=np.array([points["col4"]])[0]
#    
#    func_obs=inpl.interp1d(v_pts,vFv_pts)
#    """
#    %
#    %
#    %
#    %
#    """
#    """READING FILES"""
#     #read fort.81 or .85
#    fphotons=open(route+name+'.85','r')
#    f85=fphotons.readlines()
#    fphotons.close()
#    epsfile=[]
#    Ifile=[]
#    
#    LEN=0
#    for ln in f85:
#       if float(ln.split("  ")[1])==0. and LEN!=0:
#           break
#       LEN+=1
#    
#    stcnt=0
#    LENst=0
#    for ln in f85[2*LEN+1::]:
#       if stcnt>1:
#           LENst+=1
#       if float(ln.split("  ")[1])==5.:
#           stcnt+=1
#       
#    n0=int(len(f85)/LEN)
#    for y in f85:
#        x=y.split("  ")
#        z=x[1].strip()
#        w=x[2].replace('\n','').strip()
#        if w=='':
#            w='-100'
#        epsfile.append(z)
#        Ifile.append(w)
#    epsfile=np.array(epsfile).astype(np.float)
#    Ifile=np.array(Ifile).astype(np.float)
#    
#    timebin=tend/nsteps
#    """ STEADY STATE"""
#    #STEADY STATE
#    epsst=epsfile[-LENst::]
#    Ist=Ifile[-LENst::]
#    epsst=epsst[Ist>max(Ist)-THRESHOLD]
#    Ist=Ist[Ist>max(Ist)-THRESHOLD]
#    xsteady=np.log10(delta*10**epsst*(me*c**2/h)/(1+zsh))  #observing (SED) units
#    ysteady=np.log10(10**Ist*delta**4/D**2*R*me*c**3/sigmaT/3)  #observing (SED) units
#    
#
#    lims=[round(min(v_pts)-0.5),round(max(v_pts)+0.5),\
#                          round(min(vFv_pts)-HEIGHT),int(max(vFv_pts)+HEIGHT+0.5)]
#    
#    """BLACK BODIES (BB)"""
#    #Correct for 3c273 (lessen the BB from bibliographical values)
#    if oo==0:
#        Gamma=2.36*Gamma
#    
#    eps85bb=epsst
#    Ibb=epsst
#    
#    if SSConly=='no':
#     if len(eps85bb)<200:
#        eps85bb=np.linspace(min(epsst),max(epsst),200)
#        Rext=Rbb[oo]
#        Ibb=BB(eps85bb,T/Gamma,lext/Gamma**2,Rext)+np.log10(Gamma)
#        xbb=eps85bb+np.log10((me*c**2/h))
#        ybb=np.log10(10**Ibb/D**2*R*me*c**3/sigmaT)  #observing (SED) units
#        bbmodel=inpl.interp1d(10**xbb,10**ybb)  #model BB
#        if secondBB=='on':
#                Ibb2=BB(eps85bb,Tbb2*Gamma2,(Gamma2)**2/Gamma**2*factbb*(Rext/Rbb2)**2*lext,Rbb2)+np.log10(Gamma)
#                xbb2=eps85bb+np.log10((me*c**2/h)/Gamma2)
#                ybb2=np.log10(10**Ibb2*Gamma2**-2/D**2*R*me*c**3/sigmaT)
#                bb2model=inpl.interp1d(10**xbb2,10**ybb2)
#    
#    fig1,ax1=plt.subplots(num=11)
#    ax10=ax1
#    plt.axis(lims)
#    plt.xlabel(r'$log\,('+xlab+')\;\;'+' ['+xaxis+']$',fontsize=15,labelpad=-1)
#    plt.ylabel(r'$log\,('+yaxis+')\;\;['+yunits[yaxis]+']$',fontsize=15)
#    plt.title(obj+'  SED',fontsize=18)
#    
#    plt.figure(11)
#    ax10.plot(xsteady,ysteady,'k-',linewidth=2.0,alpha=1.0,label='Steady')
#    if Multicolor=='on':                     
#        for i in range(len(dtpts)):
#             ax10.errorbar(globals()['v_pts'+str(i)],globals()['vFv_pts'+str(i)],\
#                         yerr=globals()['errorvFv'+str(i)],\
#                         elinewidth=1.5,capsize=2.5,capthick=0.5,ecolor='k',\
#                         fmt=form[i],color=colors[i],ms=2.5,\
#                         label=legend_pts[i])
#    else:
#        ax10.errorbar(v_pts,vFv_pts,yerr=errorvFv,elinewidth=1.5,capsize=2.5,\
#                 fmt='r.',ecolor='red',ms=3.5,label='Obs')
#    if SSConly=='no':
#        ax10.plot(xbb,ybb,'r--',linewidth=0.9,label='BLR')
#    if secondBB=='on':
#        ax10.plot(xbb2,ybb2,'r-.',linewidth=0.9,label='Disk')
#    lg4=ax10.legend( loc='upper left')
#    ax10.add_artist(lg4)
#    ax10.xaxis.set_minor_locator(AutoMinorLocator())
#    ax10.yaxis.set_minor_locator(plt.LinearLocator(4*(lims[3]-lims[2])+1))
#    ax10.fill_between(np.linspace( np.log10(10**range_smarts[0]), np.log10(10**range_smarts[1])),-20,20,
#               facecolor='red', alpha=0.2)
#    ax10.fill_between(np.linspace(np.log10(10**range_xsoft[0]),np.log10(10**range_xsoft[1])),-20,20,
#               facecolor='grey', alpha=0.2)
#    ax10.fill_between(np.linspace(np.log10(10**range_fermi[0]),np.log10(10**range_fermi[1])),-20,20,
#               facecolor='blue', alpha=0.2)
#    ax10.tick_params(axis='x', labelsize=12)
#    ax10.tick_params(axis='y', labelsize=12)
#    ax10.tick_params(which='major', width=1.2, length=7)
#    ax10.tick_params(which='minor', width=1.1, length=4)
#    for axis in ['top', 'bottom', 'left', 'right']:
#        ax10.spines[axis].set_linewidth(1.1)
#    plt.figure(11).savefig(routeimage+'SS_'+obj+'.eps',bbox_inches='tight')
#    plt.figure(11).savefig(routeimage+'SS_'+obj+'.png',bbox_inches='tight')
#    
#    
#    
#    
