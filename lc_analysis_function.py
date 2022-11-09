from matplotlib import use
use('Qt5Agg')

from DELCgen import *
import scipy.stats as st
import csv
import numpy as np
from scipy import interpolate as inpl
import matplotlib.pyplot as plt 
from astropy.io import fits
import os


def write_input_file(name_inp,route='./',START=-300, RESdex=-300,STEPout=-300, TIME=-300, 
                         R=-300, B=-300,
                         loggmin= -300, loggmax=-300, p=-300, logle=-300, tesc=-300, 
                         T=-300, loglext=-300, Gamma=-300, delta=-300,
                         LOG_NOISE_E=-300, LOG_NOISE_G = -300, TOL=-300
                         ):
        LL=[]
        if route[-1]!='/':
            route+='/'
        cn5=open(route+name_inp+'.inp','r')
        i=0
        for line in cn5:
            if i==0:
                START_prev, RESdex_prev, STEPout_prev, TIME_prev=  line.split()[0],  line.split()[1],  line.split()[2],  line.split()[3]
                if START==-300: START= int(START_prev)
                if RESdex==-300: RESdex=int(RESdex_prev)
                if STEPout==-300: STEPout=int(STEPout_prev)
                if TIME==-300: TIME=float(TIME_prev)
                LL.append('0\t'+str(RESdex)+'\t'+str(STEPout)+'\t'+str(TIME)+'\n')
            elif i==1:
                    #LL.append('7.   -10.   -34.  -53.  1.e-4\n')
                    LOG_NOISE_E_prev, LOG_NOISE_G_prev, TOL_prev = line.split()[2], line.split()[3], line.split()[4]
                    if LOG_NOISE_E==-300: LOG_NOISE_E==LOG_NOISE_E_prev
                    if LOG_NOISE_G==-300: LOG_NOISE_G=LOG_NOISE_G_prev
                    if TOL==-300: TOL=TOL_prev
                    LL.append('7.   -10.   {}  {}  {}\n'.format(LOG_NOISE_E,LOG_NOISE_G,TOL))
            elif i==2:
                R_prev = line.split()[0]
                B_prev = line.split()[1]
                if R==-300: R=float(R_prev)
                if B==-300: B=float(B_prev)
                LL.append(str(round(R,int(-np.log10(R)+2)))+'\t'+str(round(B,int(-np.log10(B)+2)))+'\t1.\t1.\n')
            elif i==4:
                loggmin_prev = line.split()[1]
                loggmax_prev = line.split()[2]
                p_prev = line.split()[3]
                logle_prev = line.split()[4]
                besc_prev = line.split()[5]
                if tesc<1.0: tesc=1.0
                if loggmin==-300: loggmin=float(loggmin_prev)
                if loggmax==-300: loggmax=float(loggmax_prev)
                if p==-300: p=float(p_prev)
                if logle==-300: logle=float(logle_prev)
                if tesc==-300: tesc=1/float(besc_prev)
                LL.append('1\t'+str(round(loggmin,3))+'\t'+str(round(loggmax,3))+'\t'+str(round(p+0.01,3))+'\t'\
                      +str(round(logle,3))+'\t'+str(round(1/tesc,3))+' 0. 0\n')
            elif i==5:
                T_prev = line.split()[1]
                lext_prev = line.split()[2]
                if T==-300: T=float(T_prev)
                if loglext==-300: loglext=np.log10(float(lext_prev))
                LL.append('1\t'+str(round(T,int(-np.log10(T)+2)))+'\t'+str(round(10**loglext,int(-loglext+2)))+'\n')
            elif i==7 and Gamma!=-300 and 'new4' in name_inp:
                LL.append(str(float(Gamma))+'\n')
            elif i==8 and 'new4' in name_inp:
                delta_prev = line.split()[0]
                if delta==-300: delta = float(delta_prev)
                LL.append(str(int(delta))+'\n')
                lines=cn5.readlines()
                LL=LL+lines
                break
            elif i==9  and 'new4' not in name_inp:
                LL.append('1  0  1  1  1  1  1 \n')
                lines=cn5.readlines()
                LL=LL+lines
                break
            else:
                LL.append(line)
            i=i+1
        cn5.close()
        newcn5=open(route+name_inp+'.inp','w')
        newcn5.writelines(["%s" % item  for item in LL])
        newcn5.close()


def lc_analysis(route_main='./', lendays=None, NEW_delc=None,POWER=None,namevar=None, nameday=None,
                 SAVE_LOG=True , PLOT= False, obj='blazar', ALLOW_MULTI=False,
                 DATA_PROCESS=   False , SEGMENTS=False,
                 datfile_follow_up = 'current_residuals.txt', DOUBLE_VAR= False #allow all the following transformations if activated + SEDextract
                 ):
    #route = os.getcwd()+'/'
    #route_main = route.split('Lightcurve_Simulation')[0]
    if route_main[-1]!='/':
        route_main = route_main+'/'
    route_steady = route_main+'Steady_States/'
    route_data  = route_main +'data/' 
    route_lc = route_main+'/Lightcurve_Simulation/'
    route_timelapse = route_main+'/data/seds_all_flares/'
    
    #for mrk421 x-ray 56393 and MAGIC >56393 
    
    #units, use for different types of input light curves
    eV=     1.602192e-12 #erg
    MeV=        10**6*eV #erg
    units=            1. #set to one if default units/ needed for correct normalization  check IF line 
    
    objects=['3C273','3C279','PKS2155304','PKS1510-089', 'Mrk421']
    if obj not in objects:
        print(objects)
        oo= int(input("\n Give Number of object [0 , 1 ..] from following list:\t"))
        obj=objects[oo]
    oo=objects.index(obj)
    
    POWER_daysplit = {'56393half':0.5,   '56394':1.5,  '56395':0.5,   '56396':0.5,  '56397':0.5,
                      '56398':1.0,   '56399':1.5,   '56400':1.5, '56393':0.5, '56392':0.5}
    
    initfile_daysplit = {'56393half': 'NUSTAR_LC_7-30keV_3cols.txt',   
                         '56394':'MAGIC_LC_3cols.txt',  '56395':'MAGIC_LC_3cols.txt',
                         '56396':'MAGIC_LC_3cols.txt',  '56397':'MAGIC_LC_3cols.txt',
                            '56398':'MAGIC_LC_3cols.txt',   '56399':'MAGIC_LC_3cols.txt',  
                            '56400':'MAGIC_LC_3cols.txt'}
    
    len_daysplit = {'56393half':121, '56392':60, '56393':60, '56394':96,  '56395':96,   '56396':100,  '56397':93,
                    '56398':91,   '56399':99,   '56400':102} #for full filling of gaps this many datapoints (automatically changes)
    
    #for t_prior = 0.0
    tfirst_daysplit = {'56394': 56393.882575, '56395': 56394.896345,'56396':56395.872605,'56397': 56396.898165,'56398':
           56397.930555, '56399':56398.928735,'56400': 56399.940195,'56401': 56401.03893,
           '56393half':56392.88998, '56393':56393.882575, '56392':56392.921865}
        
    tlast_daysplit = {'56394': 56394.896345,'56395':56395.872605,'56396': 56396.898165,'56397':
           56397.930555, '56398':56398.928735,'56399': 56399.940195,'56400': 56401.03893,
           '56393half':56394.896345, '56393':56394.896345, '56392':56393.15387}
    
    #Settings for usage of created fake lightcurve
    """ READ INPUT FILE """
    cn4=open(route_lc+'fkTC.inp','r')
    lines=cn4.readlines()
    cn4.close()
    namevar_txt=str(lines[9].split(' ')[0]) #variability on these parameter
    lendays_txt=int(str(lines[9].split(' ')[1])) #days of observatSions to recreate
    NEW_delc_txt=str(lines[9].split(' ')[2]) #create new fake timecurve with 'on'
    POWER_txt=float(lines[9].replace('\n','').split(' ')[3]) #sqrt of final fake lightcurve
    if SEGMENTS:
        nameday_txt= lines[9].split()[4] #sqrt of final fake lightcurve
        lendays_txt =len_daysplit[nameday]  #lendays and remt input are useless for Mrk421; all LC is simulated
        POWER_txt = POWER_daysplit[nameday]
    else:
        nameday_txt=None
    #end for mrk421
    
    if nameday==None:
        nameday = nameday_txt
        print('nameday not provided. Using fkTC.inp value')
    if namevar==None:
        namevar= namevar_txt
        print('namevar not provided. Using fkTC.inp value')
    if POWER==None:
        POWER= POWER_txt
        print('POWER not provided. Using fkTC.inp value')
    if NEW_delc==None:
        NEW_delc= False
        print('NEW_delc not provided. Set to False')
    if lendays==None:
        lendays= lendays_txt
        print('NEW_delc not provided. Set to False')
    #multi option (not for delta), not for newly created LC
    MULTI ='no'
    remt=lendays
    if ALLOW_MULTI:
        if 'multi' in NEW_delc:
            MULTI='on'
            remt=float(NEW_delc.split('multi')[1])
            if remt>lendays  and nameday!=None:remt=lendays
            NEW_delc=False
    else:
        NEW_delc=False
    
    ###HANDY (comment out for run_var.sh use)
    #namevar='B'
    #NEW_delc=False
    #POWER=0.5 #sqrt of final fake lightcurve
    #MULTI='no'
    #lendays=724 #maximum value days
    #nameday='56394'int(lines[9].replace('\n','').split(' ')[4]) #sqrt of final fake lightcurve
    #print('Manual input from my_analysis_vesrsion.py script')
    
    """EMMANOUELOPOULOS FAKE TIME CURVE CREATION SETTINGS"""
    #input
    initfile= 'NUSTAR_LC_7-30keV_3cols.txt'
    #initfile= 'MAGIC_LC_3cols.txt'
    #initfile = initfile_daysplit[nameday]
    tbin=           0.01  # 0.5 for fermi-LAT
    #Mrk421
    
    #how to read input, used one time
    SED_EXTRACT , logband = False , 18.38 #will create initfile so name it properly
    CSVfile=       False
    FITSfile=      False
    #other options
    plotLC=         True   #plot the output of Emanoulopoulos code
    oneindex_model= True   #use a simple model for PSD 
    GENERATE_FAKE=       False  #don't run Emanoulopoulos code overall use the imported data
    #initfile ='Fermi_DLgen'+'.txt'  #test file: initfile = "pks1510-089_gamma_lc.txt"
    
    #Data treatment
    FILL_WITH_STEADY_STATE , split_dt , time_prior = True , 0.5, 0.01 #will create initfile so name it properly
    DATAcut=        False #cuts point below the threshold flux
    DATAfill=       False #interpolates between cut points
    DATAreplace=    False #replaces with threshold value , fluxes below threshold
    thresflx=       1e-20*units # ph/cm^2/s Check for every specific object
    norma_stored = -9.6433  #stored value to be used if no processing is taking place
    
    #Interpolation
    imethod=    'slinear'
    nint=       4 #number of points between t_cross intervals 
    ##
    #redshifts and bibliography
    Dist=[749,3115,543,1937,142][oo] #Mpc Luminosity Distances
    zsh=[0.158,0.536,0.0,0.36, 0.031][oo] #redshifts
    
    #special correction
    lobjun=['PKS2155304']
    if obj in lobjun and FITSfile!='on':
            units=MeV*1680. #counts per day to energy flux with average energy  per count ~1.7 GeV
            
    init=-1#index of initial point from fake LC delc1.time[ ] to segment //  set -1 for segment around maximum
    widthgamma=1.5 #for 'gmax' namevar / forced width of variation in orders of magnitude
    
    
    #VARIABLES
    # tlast=float(lines[0].split()[-1].replace('\n',''))
    # nsteps=float(lines[0].split()[2])
    cn4=open(route_steady+nameday+'_new4.inp','r')
    lines=cn4.readlines()
    cn4.close()
    
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
    # linear quants
    le , gmin, gmax =   10**logle ,  10**loggmin ,   10**loggmax
    
    
    # if initfile=='NUSTAR_LC_7-30keV_3cols.txt' and namevar='B':
    #     POWER = 1/2 #2/(p+1) #for X-rays explicit relation
        
    if namevar=='gmax' or namevar=='p':
        POWER=None
        
    #calculations
    norm = 10**logle #default in case of not using norm
    if not DOUBLE_VAR:
        norm=locals()[namevar] #choose steady state parameter for normalization
    lenfk=lendays+1 #for fine array creation
    dilate=1-((((1+zsh)**2-1)/((1+zsh)**2+1))**2)**2 #dilation factor, https://www.scirp.org/journal/paperinformation.aspx?paperid=65598
    tcross=R/3/10**10/3600./24. #t_cross in days
    tcr=tcross*dilate/delta
    
    if CSVfile:
        with open(route_lc+initfile) as csvfile:
              data = list(csv.reader(csvfile))
    
        data2 = np.genfromtxt(route_lc+initfile, delimiter=',',dtype=float)
        t1, norms, errors = [], [], []
        
        for row in data2:
            t1.append(row[0])
            norms.append(row[1])
            errors.append(row[2])
        x = [0] * len(t1)
        x = [i for i in xrange(len(t1))]
        #use real dates
        np.savetxt(route_lc+'myfile.txt', np.c_[t1,norms,errors])
        #save LC with some dummy delta step of 1 day 
        np.savetxt(route_lc+'myfile2.txt', np.c_[x,norms,errors])
        datfile='myfile.txt'
    elif FITSfile:
            datfile=fits.open(route_lc+'lc_daily_'+obj+'.fits')
            tbdata=datfile[1].data
            norms,errors,MJD=[],[],[]
            norms2,errors2= [],[]
            for tb in tbdata:
                MJD.append(tb[0])
                norms2.append(tb[7]) #photon flux
                errors2.append(tb[6])
                norms.append(tb[11])  #energy flux
                errors.append(tb[12])
            MJD=np.array(MJD)
            norms=np.array(norms)*MeV
            norms2=np.array(norms2)
            norma=np.mean(norms*MeV)
            errors=np.array(errors)*MeV
            errors2=np.array(errors2)
            np.savetxt(route_lc+'myfile.txt', np.c_[MJD,norms,errors])
            datfile='myfile.txt'
    elif SED_EXTRACT:
        from os import listdir
        initfile =     'SEDextract.txt'  #test file: initfile = "pks1510-089_gamma_lc.txt"
        llf = [[float(i.split('_')[1]),float(i.split('_')[2].replace('.txt',''))] for i in listdir(route_timelapse) if 'sed_' in i]
        llf.sort()
        nlf  = np.array(llf)-llf[0][0]
        store=[]
        f = open(initfile,'w')
        f.write('#time #flux #error')
        for i in range(len(nlf)):
                namedata  = route_timelapse+'sed_'+str(llf[i][0])+'_'+str(llf[i][1])+'.txt'
                tdata = np.loadtxt(namedata)
                v_pts  = tdata[:,0]
                norms  = tdata[:,1] #vfv
                errors = tdata[:,2]
                zo = min(v_pts, key=lambda x:abs(x-10**logband))
                iz = list(v_pts).index(zo)
                store.append([(llf[i][0]+llf[i][1])/2, norms[iz] , errors[iz]])
                f.write('\n {} {} {}'.format((llf[i][0]+llf[i][1])/2, norms[iz], errors[iz]))
        f.close()
        store = np.array(store)
        np.savetxt(initfile,store)
        datfile=initfile
    else:
        datfile=initfile
    
    #load as an object
    datalc = Load_Lightcurve(route_lc+datfile,tbin)
    norma  = norma_stored
    #SPLIT in subflares and introduce initial time
    norma=np.mean(datalc.flux)*units
    #normalize to the mean flux for the Fit_PDF to work
    datalc.errors=datalc.errors/np.mean(datalc.flux)
    datalc.flux = datalc.flux/np.mean(datalc.flux)
    thresflx=thresflx/norma
    norm_init = 0.8 #units of norm_prev
    if namevar=='p':
        cnS = np.loadtxt(route_data+'powerlaw_fit_nustar.dat')
        cnS = cnS[cnS[:,1].argsort()]
        p_mjd    =( cnS[:,1]+cnS[:,2]) /2.0
        nustar_slopes = cnS[:,4]
        p_lc = (cnS[:,4]-2)*2+3
        p_lc_error = (cnS[:,5]*2 -cnS[:,3]*2)/2
        datalc.time = p_mjd 
        datalc.flux = p_lc
        datalc.errors  = p_lc_error
        norm_init = 1.0
    
    if namevar=='le' and DOUBLE_VAR:
        namevar_sec = 'p'
        datfile = datfile_follow_up
        nameday_p = str(nameday)
        if int(nameday) in [56393, 56394]: nameday_p = '56393half'
        fakeTC_first = np.loadtxt(route_lc+'fakeTC_p_Mrk421_'+nameday_p+'.txt')
        diffs = np.loadtxt(route_lc+datfile_follow_up)
        datalc.time = diffs[:,0]
        datalc.flux = (1+diffs[:,1])
        datalc.errors = (1+diffs[:,1])/10.0
    
    if DATA_PROCESS:
        if FILL_WITH_STEADY_STATE:
            diff = datalc.time[1::]-datalc.time[0:-1]
            times_start = datalc.time[1::][diff>split_dt] # 0.11 days found from histogram of differences last great subday difference is 0.104
            times_orig = datalc.time
            fluxes_orig = datalc.flux
            i_start =  [list(datalc.time).index(t) for t in times_start] +[len(times_orig)-1]
            renorm =    1.0 #units norm_prev, if 1.0 regulates what fraction of the first flux value to be the steady state 
            norm_prev = 1.0 #np.mean(datalc.flux[0:i_start[0]])
           
            err_orig = datalc.errors
            fullfill =  False
            if time_prior==None: 
                time_prior_0 = 5.0
            else:
                time_prior_0 = time_prior
            Nadded =int(np.floor(time_prior/tbin))
            er = np.array([err_orig[0]*renorm/norm_prev,]*Nadded + list(datalc.errors[0:i_start[0]]*renorm/norm_prev))
            splits =  [[np.array(list(np.linspace(times_orig[0]-time_prior+tbin, times_orig[0]-tbin,Nadded))+ list(datalc.time[0:i_start[0]])) , 
                       np.array([fluxes_orig[0]*renorm/norm_prev,]*Nadded + list(datalc.flux[0:i_start[0]]*renorm/norm_prev)), 
                       er]]
            totaladded=0
            for j,i in enumerate(i_start[0:-1]):
                #if j!=len(i_start)-1: 
                    #renorm = np.mean(fluxes_orig[i_start[j]:i_start[j+1]])
                if time_prior==None: 
                    time_prior = times_orig[i]-times_orig[i-1]
                    fullfill=True
                Nadded=int(np.floor(time_prior/tbin))
                tt, ff, eff = np.linspace(times_orig[i]-time_prior+tbin, times_orig[i]-tbin,Nadded)  ,  [fluxes_orig[i]*norm_init/norm_prev,]*Nadded , [err_orig[i]*norm_init/norm_prev,]*Nadded 
                ff = list(np.array(ff)[tt>times_orig[i-1]+tbin])
                eff= list(np.array(eff)[tt>times_orig[i-1]+tbin])
                if len(ff)!=len(tt):
                    print('subtracting:'+str( len(ff)- len(tt)))
                    Nadded+= len(ff)- len(tt) #subtract the point removed
                print( '\n Section '+str(j)+' t in ['+str(tt[0])+'-'+str(times_orig[i_start[j+1]])+']\n'+ 'Added '+str(Nadded))
                print('Total points in this time inteval: \t'+str(Nadded+i_start[j+1]-i))
                tt = list(tt[tt>times_orig[i-1]+tbin])
                datalc.time= np.array(list(datalc.time[0:i+totaladded]) +tt+ list(datalc.time[i+totaladded::]))
                datalc.flux= np.array(list(datalc.flux[0:i+totaladded]) +ff+ list(datalc.flux[i+totaladded::]*renorm/norm_prev))
                datalc.errors = np.array(list(datalc.errors[0:i+totaladded]) +eff+ list(datalc.errors[i+totaladded::]*renorm/norm_prev))
                totaladded+=Nadded
                tx =datalc.time[i+totaladded-Nadded:i_start[j+1]+totaladded]
                fx = datalc.flux[i+totaladded-Nadded:i_start[j+1]+totaladded]
                errr = datalc.errors[i+totaladded-Nadded:i_start[j+1]+totaladded]
                splits.append([tx, fx, errr])
                print( 'Renorm factor '+str(renorm/norm_prev))
                norm_prev = renorm
                if fullfill: time_prior = None
            print('Extrapolated for gaps with steady state')
            sp, t= splits[0], '56393half' #int(np.median(times_orig[prev_t_index:i_start[i-1]]))
            if 'MAGIC' in datfile: t = '56392'
            print( '\n Segment initial segment'+' t in ['+str(times_orig[0])+'-'+str(times_orig[i_start[1-1*('MAGIC' in datfile)]-1])+']')
            print('Name given:\t '+str(t)+'\n'+'Length segment:\t'+str(len(sp[0])))
            m_sp=np.array(list(sp[0][1::] == sp[0][0:-1])+[False]) #identicals
            np.savetxt(route_lc+'myfile_new_MJD_{}.txt'.format(t), np.c_[sp[0][~m_sp],sp[1][~m_sp],sp[2][~m_sp]])
            len_daysplit[str(t)] = len(sp[0][~m_sp])
            tfirst_daysplit[str(t)] = sp[0][~m_sp][0]
            tlast_daysplit[str(t)] = sp[0][~m_sp][-1]
            for j in range(1-1*('MAGIC' in datfile),len(i_start)-1):
                sp = splits[j]
                if len(sp[0])>1:
                    t = 56394+j -1*('MAGIC' in datfile)#int(times_orig[i_start[j]]) #int(np.median(times_orig[prev_t_index:i_start[i-1]]))
                    print('\n Segment '+str(j)+' t in ['+str(times_orig[i_start[j]])+'-'+str(times_orig[i_start[j+1]-1])+']')
                    print('Name given:\t '+str(t)+'\n Length segment:\t'+str(len(sp[0])))
                    m_sp=np.array(list(sp[0][1::] == sp[0][0:-1])+[False]) #identicals
                    np.savetxt(route_lc+'myfile_new_MJD_{}.txt'.format(t), np.c_[sp[0][~m_sp],sp[1][~m_sp],sp[2][~m_sp]])
                    len_daysplit[str(t)] = len(sp[0][~m_sp])
                    tfirst_daysplit[str(t)] = sp[0][~m_sp][0]
                    tlast_daysplit[str(t)] = sp[0][~m_sp][-1]
                else:
                    print('Oops, something empty was produced at segment '+str(j))
            np.savetxt(route_lc+'myfile_new.txt', np.c_[datalc.time,datalc.flux,datalc.errors])
            if PLOT:
                plt.figure(figsize=[16,9])
                for sp in splits:
                    x= sp[0]
                    y=sp[1]
                    if len(x)>0:
                        plt.plot(x-56390,np.log10(y),'.')
                plt.plot(times_orig-56390,np.log10(fluxes_orig*1.5),'kv',ms=2) 
        elif DATAcut:
            cp=len(datalc.time[datalc.flux<thresflx])
            print( str(cp)+' points are below the threshold, they are deleted as fake\n')
            datalc.time=np.linspace(datalc.time[0],datalc.time[0]+tbin*len(datalc.time[datalc.flux>thresflx]), len(datalc.time[datalc.flux>thresflx]))
            datalc.errors=datalc.errors[datalc.flux>thresflx]
            datalc.flux = datalc.flux[datalc.flux>thresflx]
        elif DATAfill:
            fluxmodel=inpl.interp1d(datalc.time[datalc.flux>thresflx],datalc.flux[datalc.flux>thresflx])
            errormodel=inpl.interp1d(datalc.time[datalc.flux>thresflx],datalc.errors[datalc.flux>thresflx])
            tnmin=min(datalc.time[datalc.flux>thresflx])
            tnmax=max(datalc.time[datalc.flux>thresflx])
            tntime=datalc.time[datalc.time>tnmin]
            datalc.time=tntime[tntime<tnmax]
            datalc.flux = fluxmodel(datalc.time)
            datalc.errors=errormodel(datalc.time)
        if DATAreplace:
            for dp in range(0,len(datalc.time)):
                if datalc.flux[dp]<thresflx:
                    datalc.flux[dp]=thresflx
    #save for retrieving them in var.py file
    np.save(route_lc+'tfirst.npy', tfirst_daysplit, allow_pickle=True)
    np.save(route_lc+'tlast.npy', tlast_daysplit, allow_pickle=True)
    #datalc = Load_Lightcurve(route_lc+'/myfile_new.txt',tbin)
    datalc_day = Load_Lightcurve(route_lc+'/myfile_new_MJD_'+nameday+'.txt',tbin)
    if not GENERATE_FAKE: 
        init = 0
        delc1=datalc_day
        if namevar!='p' and not DOUBLE_VAR:
            delc1.flux = delc1.flux/np.mean(delc1.flux)
        #%%   
    if NEW_delc: 
            # plot the light curve
            if plotLC:
                datalc.Plot_Lightcurve()
                
            # Bending power law params
            A,v_bend,a_low,a_high,c = 0.04, 2e-2, 2.0, 0.5, 0.0
            #Simple power law params/  psd(v)=cc0+cc*v^mm
            mm,cc,cc0=-1,1,0.1
            # Probability density function params
            kappa,theta,lnmu,lnsig,weight = 0.01,1.0, 0.01, 1.0,0.5
            
            # create mixture distribution to fit to PDF
            model1=st.lognorm
            model2=st.gamma
            #new parameters of mixture model sub-models
            #par1,par2,par3,par4,weight=kappa,theta,lnmu,lnsig,weight
            par1,par2,par3,par4,weight= 1.0,1.2, 3.1, 0.85, 0.5
            mix_model = Mixture_Dist([model1,model2],[3,3],[[[2],[0]],[[2],[0]]])
    
            if datalc.psdFit==None and GENERATE_FAKE:
                if oneindex_model:
                    def func_powerlaw(v, m, c, c0):
                        return c0 + v**m * c
                    datalc.Fit_PSD(initial_params= [mm,cc,cc0],model=func_powerlaw)
                else:
                    datalc.Fit_PSD(initial_params=[A,v_bend,a_low,a_high,c])
            
            # fit the PDF of the light curve using the default mixture model
            if datalc.pdfFit==None and GENERATE_FAKE:
                datalc.Fit_PDF(initial_params=[par1,par2,par3,par4,weight],model=mix_model)
            
            # simulate artificial light curve with Emmanoulopoulos method, using the PSD
            # and PDF of the data light curve, with default parameters (bending power law)
            # for PSD and mixture distribution of gamma and lognormal distribution for PDF)
            if datalc.pdfFit and datalc.psdFit and GENERATE_FAKE:
                # 1st simulated LC
                delc1= datalc.Simulate_DE_Lightcurve()
                # 2nd simulated LC
                delc2 = datalc.Simulate_DE_Lightcurve()
                
                #renormalize to units
                delc1.errors=0*datalc.errors
                # Plot Compare Diagram and save it
                Comparison_Plots([datalc,delc1],names=["OBS","FAKE"],bins=50,norma=1.0)#Comparison_Plots([datalc,delc1,delc2],names=["Data LC","fake 1","fake 2"],bins=25)
                plt.savefig(route_lc+'DELC.pdf',bbox_inches='tight')
                plt.savefig(route_lc+'DELC.png',bbox_inches='tight')
                #save fake timecurves in output file
                np.savetxt(route_lc+'fakeLC', np.c_[delc1.time,delc1.flux,delc1.errors])
            
            if init==-1:
                init=list(delc1.flux).index(max(delc1.flux))-int(lendays/tbin/2)
                if init<0:
                    init=0
                if init>len(delc1.time)*tbin-lendays:
                    init=len(delc1.time)*tbin-lendays-1
    elif GENERATE_FAKE:
         delc1 = Load_Lightcurve(route_lc+'fakeLC',tbin)
         if init==-1:
            init=list(delc1.flux).index(max(delc1.flux))-int(lendays/tbin/2)
            if init<0:
                    init=0
            if init>len(delc1.time)-lendays:
                    init=len(delc1.time)*tbin-lendays-1
    ft=delc1.time[init:init+lenfk]*delta/dilate
    flx=delc1.flux[init:init+lenfk]
    #%%
    """ Transform to Timecurve the cut lightcurve [ft, flx]"""
    if namevar=='p': 
        origlc = datalc
        norm= np.mean(datalc_day.flux)
    elif (namevar=='le'and namevar_sec=='p'): 
        origlc = datalc
        norm  = le
    else:
        origlc = Load_Lightcurve(route+datfile,tbin)
    msk_delc1 =( origlc.time>min(delc1.time) )&(origlc.time<max(delc1.time))
    origft=origlc.time[msk_delc1]*delta/dilate
    orignorm = np.mean(origlc.flux[msk_delc1])
    origflx = origlc.flux[msk_delc1] /orignorm
    origfpoints =  norm * origflx
    
    flx[0]=flx[1] #first point from steady state value
    fpoints=norm*flx
    lenfk = len(flx)#redefine if more than available is requested
    #use the sqruare root of variations (for SSC/ le parameter)
    if POWER!=None:
        fpoints=norm*(flx)**POWER #/np.mean((flx)**POWER)
        origfpoints = norm*(origflx)**POWER
    minnon0=min(fpoints)
    maxnon0=max(fpoints)
    if namevar=='gmax':
        s1=np.log10(max(flx)/min(flx))/widthgamma
        fpoints=norm*flx**(1/s1)
        
    if namevar=='p':
        fpoints = flx #np.log10(flx)*(-2.0) + norm
        minnon0= 2.0
        maxnon0= 5.0
    #Interpolation
    funcfake=inpl.interp1d(ft,fpoints,kind=imethod)
    tbin_avg = np.mean(delc1.time[1::]-delc1.time[0:-1])
    if tcr<tbin:
        N=(int(tbin_avg/tcr)+1)*nint
    else:
        N=nint+1
    if namevar=='delta':  N=1
    tbin_new=tcross/tcr*tbin_avg/N
    print( 'New bin interval\n tbin_new=\t {} ={} t_cross'.format(tbin_new,round(tbin_new/tcross,5)))
    lenfake=N*(lenfk-1)
    
    ftime=np.linspace(ft[0],ft[lenfk-1],lenfake)
    fflux=funcfake(ftime)
    
    TIME=float(int(lenfake/tcross*tbin_new))
    STEPS=int(lenfake/tcross*tbin_new)
    #set first date t=0
    t0=ftime[0]
    ftime=ftime-np.ones(len(ftime))*t0
    
    """Create fakeTC to conjoin with crashed run"""
    conshort = 1 #num of tcrosses of steady state before entering new time-varying period
    logf=open(route_lc+'multi.log','a')
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
        SAVE=False
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
        if fflux[i]>maxnon0:
            fflux[i] = maxnon0
    if namevar=='B': fflux[0] = B
    
    if namevar=='le' and namevar_sec=='p':
        fflux_first  =fakeTC_first[:,1]    
        if len(fflux) != len(fakeTC_first):
            if max(ftime)<max(fakeTC_first[:,0]*tcross):
                mod = inpl.interp1d(fakeTC_first[:,0], fakeTC_first[:,1])
                ftime[0]=1e-14
                fflux_first = mod(ftime/tcross)
            else:
                mod = inpl.interp1d(ftime/tcross, fflux)
                fakeTC_first[0,0]=1e-14
                fflux = mod(fakeTC_first[:,0])
                ftime = fakeTC_first[:,0]
        np.savetxt(route_lc+'fakeTC_'+namevar+'+'+namevar_sec+'_'+obj+'_'+nameday+'+'+namevar_sec+'.txt', np.c_[ftime/tcross,fflux_first,np.log10(fflux)])
        np.savetxt(route_main+'fakeTC_'+namevar+'+'+namevar_sec+'_'+obj+'_'+nameday+'.txt', np.c_[ftime/tcross,fflux_first,np.log10(fflux)])
        np.savetxt(route_main+'fakeTC.txt', np.c_[ftime/tcross,fflux_first,np.log10(fflux)])
        np.save(route_main+'tfirst.npy', tfirst_daysplit, allow_pickle=True)
        np.save(route_main+'tlast.npy', tlast_daysplit, allow_pickle=True)
    #save timcurve
    else: 
        np.savetxt(route_lc+'fakeTC_'+namevar+'_'+obj+'_'+nameday+'.txt', np.c_[ftime/tcross,fflux])
        np.savetxt(route_main+'fakeTC_'+namevar+'_'+obj+'_'+nameday+'.txt', np.c_[ftime/tcross,fflux])
        np.savetxt(route_main+'fakeTC.txt', np.c_[ftime/tcross,fflux])
        


    
    ##SAVE code.inp
    #if namevar=='B':
    #    B=min(fflux)
    
    # cn4=open(route+'code.inp','r')
    # LL=[]
    # i=0
    # for line in cn4:
    #     if i==0:
    #         LL.append('0  5  '+str(STEPS)+'  '+str(TIME)+'\n')
    #     elif i==1:
    #         #LL.append('7.   -10.   -34.  -53.  1.e-4\n')
    #         LL.append('7.   -10.   -37.  -53.  1.e-5\n')
    #     elif i==2:
    #         LL.append(str(round(R,int(-np.log10(R)+3)))+'   '+str(round(B,int(-np.log10(B)+3)))+' 1. 1.\n')
    #     elif i==4:
    #             LL.append('1   '+str(round(loggmin,3))+'   '+str(round(loggmax,3))+'   '+str(round(p+0.01,3))+'   '\
    #               +str(round(logle,3))+'\t'+str(round(besc,3))+'\t 0.\t0' +'\n')
    #     elif i==5:
    #         LL.append('1   '+str(round(T,int(-np.log10(T)+2)))+'   '+str(round(lext,int(-np.log10(lext)+2)))+'\n')
    #     elif i==9:
    #         LL.append('1  0  1  1  1  1  1 \n')
    #     elif i==10:
    #         LL.append('0  0  0  0  0  0  0 \n')
    #     else:
    #         LL.append(line)
    #     i=i+1                
    # cn4.close()
    # newcn4=open('code.inp','w')
    # newcn4.writelines(["%s" % item  for item in LL])
    # newcn4.close()
    
    """PLOTS"""
    if PLOT:
        plt.figure(4,[15,5])
        plt.plot(ftime-(ftime[0]-Tcr0)*np.ones(len(ftime)),fflux,'b.',markersize=1.5)
        plt.plot(origft-ft[0],origfpoints,'r*')
        plt.plot(ft-ft[0]*np.ones(len(ft)),norm*np.ones(len(fpoints)),'k--',markersize=0.3)
        plt.plot(ft-ft[0]*np.ones(len(ft)),np.mean(fflux)*np.ones(len(fpoints)),'g--',markersize=0.3)
        plt.legend(['interpolation','initial points','normalization of fake LC: '+namevar+'_0','<'+namevar+'> ='+str(round(np.mean(fflux),3))+'= '+str(round(np.mean(fflux)/norm,2))+namevar+'_0'])
        plt.xlabel(r'$t\;(days)$ jet frame')
        plt.ylabel(r'$'+namevar+'$')
        
        plt.figure(22,[15,5])
        plt.plot(ftime[0:-1]/tcross,np.log10(fflux[0:-1]/fflux[1::])/(ftime[1]-ftime[0])*tcross,'r-',markersize=1.5)
        plt.plot(ftime[0:-1]/tcross,1.*np.ones(len(ftime)-1),'k--')
        plt.plot(ftime[0:-1]/tcross,-1.*np.ones(len(ftime)-1),'k--')
        plt.xlabel(r'$t\;(t_{cross})$ jet frame')
        plt.ylabel(r'$d(log\,'+namevar+')/dt$')
        #plt.close('all')
    
    
    ## Save info of the time-curve
    if SAVE_LOG:
        if PLOT:
            plt.figure(4).savefig(route_lc+'fake_lightcurve_'+namevar+'.pdf')
            plt.figure(22).savefig(route_lc+'Steepness.png')
        import datetime
        dttm = datetime.datetime.now()
        datenum=str(dttm.date()).split('-',1)[1]+str(dttm.time()).split('.')[0]
        if NEW_delc:
            np.savetxt(route_lc+'fakeTC'+str(lendays)+namevar+'_'+datenum+'.txt', np.c_[ftime/tcross,fflux])
        savename='Fake Timecurve of'+namevar+' for object'+obj
        cfile=open('info_fakecurve_'+datenum+'.txt','w')
        cfile.write('#'+savename+'\n')
        cfile.write('Steady state value of '+namevar+':\t'+str(norm)+'\t log('+namevar+')='+str(np.log10(norm))+' \n')
        cfile.write('Initial Time in rest frame:\t'+str(t0)+' MJD [Position in 10yr array+'+str(init)+'/ '+str(len(delc1.time))+']')
        cfile.write('Time interval (rest frame):\t'+str(lenfake/N)+' days \n')
        cfile.write('delta :\t'+str(delta)+'\n')
        cfile.write('Contraction of Cosmic Time for Small Time Interval:\t'+str(dilate))
        cfile.write('Time Interval jet frame:\t'+str(lenfake*tbin_new)+' days \n')
        cfile.write('Radius R\':\t'+str(R)+' cm \n')
        cfile.write('Initial t_bin on jet frame:\t'+str(tbin*delta/dilate)+'\n')
        cfile.write('Selected t_bin_new < t_cross:\t'+str(tbin_new)+'days \t ='+str(tbin_new/tcross)+' t_cross \n')
        cfile.write('New number of points:\t'+str(lenfake)+'\n')
        cfile.write('Method of interpolation:\t'+imethod+'\n')
        cfile.write('Power Index of observed Variability:\t'+str(POWER)+'\n')
        cfile.write('Maximum value / Mean value:\t'+str(round(max(fflux)/np.mean(fflux),3)))
        cfile.write('\nInitial Normalization (erg/cm^2/s):\t'+str(norma))
        if NEW_delc:
            cfile.write('\n DATAcut:\t'+str(DATAcut))
            cfile.write('\n DATAreplace:\t'+str(DATAreplace))
            cfile.write('\n\n FIT PARAMETERS\n')
            cfile.write('PSD one index model:\t'+str(oneindex_model)+'\n')
            if datalc.psdFit:
                cfile.write('PSD:\t'+str(datalc.psdFit.x[0])+'\t'+str(datalc.psdFit.x[1])+'\t'+str(datalc.psdFit.x[2])+'\t'+\
                str(datalc.psdFit.x[3])+'\t'+str(datalc.psdFit.x[4])+'\n')
            if datalc.pdfFit:
                cfile.write('PDF models:\t'+str(model1.name)+'\t'+str(model2.name)+'\n')
                cfile.write('PDF:\t'+str(datalc.pdfFit.x[0])+'\t'+str(datalc.pdfFit.x[1])+'\t'+str(datalc.pdfFit.x[2])+'\t'+\
                str(datalc.pdfFit.x[3])+'\t'+str(datalc.pdfFit.x[4])+'\n')
        else:
            cfile.write('\n\n Previous fake LC used, check previous info files')
        cfile.close()
        print('Info and Graph Saved into /Lightcurve_Simulation/obj Folder')
    plt.close('all')
    
    # #%% EXTRA PLOTS
    # """ Compare Real and Fake Timecurves"""
    # fig5=plt.figure(5)
    # plt.hist(np.log10(fflux/np.mean(fflux)),bins=100,alpha=0.5)
    # plt.hist(np.log10(datalc.flux/np.mean(datalc.flux)),bins=100,alpha=0.5)
    # plt.hist(np.log10(delc1.flux/np.mean(delc1.flux)),bins=100,alpha=0.5)
    # plt.legend(['Fake TC+Interpolation','Real Whole','Fake Segment'])
    # plt.ylabel(r'$N$',size=14)
    # plt.xlabel(r'$log\;(\,Flux)$',size=14)
    # plt.axis([1.1*min(np.log10(datalc.flux/np.mean(datalc.flux))),0.9*max(max(np.log10(datalc.flux/np.mean(datalc.flux/np.mean(datalc.flux)))),\
    #                   np.log10(max(fflux)/np.mean(fflux))),0,160])
    # if SAVE_LOG:
    #     fig5.savefig(route_lc+'PDF_comp.png')
    # plt.pause(2)
    # plt.close('all')
    return TIME, STEPS, len(fflux)
