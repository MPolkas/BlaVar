!              The main driver for the solution of the kinetic 
!              equations describing particle injection and energy losses,
!              including photopair (Bethe-Heitler) and photomeson production.
!              The NAG library routine D02EJF is used (https://www.nag.com/numeric/fl/nagdoc_fl24/pdf/d02/d02ejf.pdf)

!              SUMMARY OF PHYSICAL PROCESSES
!              Electrons/positrons: 
!              external injection, physical escape, losses due synchrotron radiation, 
!              losses due to inverse Compton scattering, injection by photon-photon pair production,
!              injection from Bethe-Heitler pair production, injection from photomeson production.
!              Photons: 
!              production from synchrotron radiation, attenuation from synchrotron self-absorption, 
!              production from inverse Compton scattering, injection from photomeson production, 
!              attenuation from photon-photon pair production, physical escape.
!              Protons:
!              external injection, physical escape, losses due to Bethe-Heitler pair production,
!              losses due to photomeson production.
!              Neutrons:
!              physical escape, decay, injection from photomeson production, losses due to photomeson production
!              Neutrinos:
!              physical escape, injection from photomeson production, injection from neutron decay
!              Kaons, pions, muons: not treated with kinetic equations, but correction for syn cooling is included.
 

!               There are 13 FLAGS in code.inp controling these processes:
!               isyn=0/1   switches off/on electron  synchrotron cooling.
!               ipair=0/1  switches off/on proton-photon pair production.
!               icompt=0/1 switches off/on electron Compton cooling.
!               igg=0/1    switches off/on photon-photon pair production.
!               issa=0/1   switches off/on synchrotron self absorption.
!               iesc=0/1   switches off/on electron escape.
!               ipsc=0/1   switches off/on photon escape.
!               ipion=0/1  switches off/on photopion interactions.
!               ianni=0/1  switches off/on electron-positron annihilation.
!               imsyn=0/1  switches off/on muon synchrotron
!               ipsyn=0/1  switches off/on pion synchrotron
!               iksyn=0/1  switches off/on kaon synchrotron 
!               ineutron=0//1   switches off/on neutron-photon interactions.

!              USEFUL DIMENSIONLESS PHYSICAL QUANTITIES.  
!              TB: magnetic compactness defined as (B**2/8π)/σΤ*R, with R the size of the blob.
!              Q: the ratio of the magnetic field to the critical value Bcr=4.413*10**13 G

!              DIMENSIONS OF CODE QUANTITIES.
!              Time is measured in units of the photon crossing time R/c.
!              Particle number densities are in units of σT*R.
!              Photon energies are in units of me*c**2.


!              PARTICLE ARRAYS.
!              YP(NP): array of dimension NP containing the natural log of the dimensionless 
!                      proton number density in equally spaced intervals of deltap=d(logγ).    
!              GP(NP): array of dimension NP containing the log10 of proton Lorentz factor γ. 
!              YE(NE): array of dimension NE containing the natural log of the dimensionless 
!                      electron/positron number density in equally spaced intervals of deltap=d(logγ).   
!              GE(NE): array of dimension NE containing the log10 of electron/positron Lorentz factor γ. 
!              YG(NG): array of dimension NG containing the natural log of the dimensionless 
!                      photon number density in equally spaced intervals of deltax=d(logx)=2*deltap.  
!              XG(NG): array of dimension NG containing the log10 of the dimensionless photon energy x.
!                      The lower and upper limits of the array are: X(1)=XMIN=Q*GP(1)**2 and X(NG)=XMAX=GP(NP).
!                      XSCH=X(NP)=Q*GP(NP)**2 is the maximum energy of synchrotron photons, so there
!                      is no synchrotron contribution from XSCH to XMAX. The number of photon extra bins 
!                      NADD=NG-NP is given by NADD=LOG(XMAX/XSCH)/DELTAX. NADD must be integer, so 
!                      care must be taken when choosing GP(1), GP(NP), and NP.
!              YN(NP): array of dimension NP containing the natural log of the dimensionless 
!                      neutron number density in equally spaced intervals of deltap=d(logγ).    
!              YNT(NE):array of dimension NE containing the natural log of the dimensionless   
!                      neutrino (all flavor) number density in equally spaced intervals of deltap=d(logγ).
!              YM(NE): array of dimension NE containing the natural log of the dimensionless   
!                      muon number density in equally spaced intervals of deltap=d(logγ).
!              GM(NE): array of dimension NE containing the log10 of muon Lorentz factor γ.
!              YPI(NE):array of dimension NE containing the natural log of the dimensionless   
!                      pion number density in equally spaced intervals of deltap=d(logγ).
!              GPI(NE): array of dimension NE containing the log10 of pion Lorentz factor γ.
!              YK(NE): array of dimension NE containing the natural log of the dimensionless   
!                      kaon number density in equally spaced intervals of deltap=d(logγ).
!              GK(NE): array of dimension NE containing the log10 of kaon Lorentz factor γ.

      implicit real*8 (a-h,o-z)

      integer n,iw,ifail
      integer iprext,ielext,ielextbr, iprextbr
      integer iphotext, iphotext2
      integer isyn,ipair,icompt,igg,issa,iesc,ipsc,ipion 
      integer ianni,imsyn,ipsyn,iksyn,ineutron
      integer tarray
      real*8 npdec
      character*1 relabs
	
!     PARAMETERS USED IN d02ejf ROUTINE 
      PARAMETER(npmax=400,ntotal=800,nwork=440350,relabs='M')
	
!     PHYSICAL CONSTANTS
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
             
!     ARRAY DEFINITIONS
      DIMENSION tarray(3)
      DIMENSION y(ntotal),yp(npmax),ye(npmax),yg(npmax)
     $         ,gp(npmax),ge(npmax),x(npmax),work(nwork)
     $         ,ygamma(ntotal),ygbb(npmax)
     $	       ,yn(npmax), ynt(npmax) !N added neutron distributions SD-01/2010 !nt added neutrinos SD-01/2010
     $         ,ym(npmax), pms(npmax), gm(npmax) !m added muon distributions SD-05/2012
     $         ,ypi(npmax), ppis(npmax), gpi(npmax) !pk added pions and kaons SD-02/2014
     $         ,yk(npmax), pks(npmax), gk(npmax)    !pk
      DIMENSION fpbh(npmax,npmax,npmax)
      DIMENSION floss(40)
      DIMENSION Tpg0(3,101), Tpgm(300), Tpgpi(300), Tpgk(300)    !m !pk
      DIMENSION Rpg01(3,51,150)
      DIMENSION Rpg02(3,51,150)
      DIMENSION Rpg03(3,51,150)
      DIMENSION Rpg04(3,51,150)
      DIMENSION Rpg05(3,51,150)
      DIMENSION Rpg07(3,51,150)     !pk
      DIMENSION Rpg08(3,51,150)     !pk
      DIMENSION Rpg09(3,51,150)     !pk
      DIMENSION Rpg10(3,51,150)     !pk
      DIMENSION Rpg13(3,51,150)
      DIMENSION Rpg14(3,51,150)
      DIMENSION Rpg15(3,51,150)
      DIMENSION Rpg16(3,51,150)
      DIMENSION Rpg17(3,51,150)
      DIMENSION Rpg18(3,51,150)
      DIMENSION RatioM(300)                 !m
      DIMENSION RatioPi(300), RatioKph(300)   !pk
      DIMENSION RatioKe(300), RatioKpi(300), RatioKm(300), RatioKnt(300)   !pk
      DIMENSION Rpm02(130,150)
      DIMENSION Rpm15(130,150)
      DIMENSION Rpm18(130,150)
      DIMENSION Rpp04(130,150)
      DIMENSION Rpp17(130,150)
      DIMENSION Rpk01(130,150)
      DIMENSION Rpk03(130,150)
      DIMENSION Rpk05(130,150)
      DIMENSION Rpk07(130,150)
      DIMENSION Rpk08(130,150)
      DIMENSION Rpk16(130,150)
      DIMENSION Rpk18(130,150)
      DIMENSION fktime(100000),fakeB(100000) !FAKE TIMECURVE

!     DECLARATION OF SUBROUTINES      
      external out
      external deriv
      external d02ejf,d02ejy,d02ejw

      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/flags/ipair,isyn,icompt,igg,issa,iesc
     $  ,ipsc,ipion,ianni,ikn 
      common/ygm/ygamma
      common/nonlin/factor
      common/tfl/xnorm
      common/htm/h,tend,iout  
      common/fqq/nsteps
      common/phexter/gext,ygext,iphotext  
      common/ext/x1,x2,xbr,beta1,beta2,extph0,iphotext2 
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp 
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,iprext,ipexp  
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr 
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr 
      common/ebb/ygbb,nbbmx	
      common/tvv/iprcl 
      common/bet/fpbh
      common/prl/floss
      common/pgp/radius,xl10min,bfield
      common/pgr/Tpg0,Rpg01,Rpg02,Rpg03,Rpg04,Rpg05,Rpg13,
     $ Rpg14,Rpg15,Rpg16,Rpg17,Rpg18,Tpgm,Rpm02,Rpm15,Rpm18,pms,RatioM,
     $ Tpgpi,RatioPi,Tpgk,RatioKph,RatioKe,RatioKm,RatioKpi,RatioKnt, 
     $ Rpp04,Rpp17,Rpk01,Rpk03,Rpk05,Rpk07,Rpk08,Rpk16,Rpk18,
     $ Rpg07,Rpg08,Rpg09,Rpg10,ppis,pks
      common/pgm/gm,ym,gpi,ypi,gk,yk
      common/iflg/imsyn,ipsyn,iksyn,ineutron
      common/fake/tol,ifake
      common/tm/fktime,fakeB,dfk  !FAKE TIMECURVE
!      common/times/ttime,tprev 

      call cpu_time (t_start)
 
      call itime(tarray)
      tt1=tarray(1)*3600.+tarray(2)*60.+tarray(3)
 
      open (unit=13, file='code.inp', status='old')
      open (unit=98, file='code_r4.dat', status='unknown')      
      open (unit=99, file='dum_r4.dat', status='unknown') 

!     READING THE INPUT FILE
      read (13,*) ireadin,npdec,nsteps,tend
      read (13,*) gpexmx,ypmin,yemin,ygmin,tol	
      read (13,*) radius, bfield,belesc
      read (13,*) iprext,gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,ipexp
      read (13,*) ielext,geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ieexp
      read (13,*) iphotext,temperat,exlumth
      read (13,*) iphotext2,x1,x2,xbr,beta1,beta2,extph0
      read (13,*) ielextbr,geextbr,slelints1,slelints2,ae2	
      read (13,*) iprextbr,gpextbr,slprints1,slprints2,ap2
      read (13,*) isyn,ipair,icompt,igg,issa,iesc,ipsc
      read (13,*) ipion,ianni,imsyn,ipsyn,iksyn,ineutron,ikn

!      FAKE TIMECURVE
       ifake=57600
      open (unit=50, file='fakeTC.txt',
      1     status='unknown') 
      do i=1,ifake
              read (50,*) fktime(i), fakeB(i)
      enddo
      dfk=fktime(3)-fktime(2)
!     APPEND TO FILES INSTEAD OF RESETING / FAKE TIMECURVE      
      if (ireadin.gt.0) open(unit=81,file='fort.81',
     $ Access = 'append',Status='old')
      if (ireadin.gt.0) open(unit=85,file='fort.85',
     $ Access = 'append',Status='old')
      if (ireadin.gt.0) open(unit=89,file='fort.89',
     $ Access = 'append',Status='old')
      if (ireadin.gt.0) open(unit=55,file='fort.55',
     $ Access = 'append',Status='old')
      
      if (ielextbr.gt.0) ielext=0      	
      if (iprextbr.gt.0) iprext=0	
 
	geextbr=10.**geextbr
	gpextbr=10.**gpextbr		
	gpextmn=10.**gpextmn
	gpextmx=10.**gpextmx
	geextmn=10.**geextmn
	geextmx=10.**geextmx
	exlumel=10.**exlumel
	exlumpr=10.**exlumpr

!	bplus=belesc
!       dimensionless B field	
	qb=bfield/Bcr
!       magnetic energy density [erg/cm^3]	
	UB=bfield**2./8./pi
!       magnetic compactness	
	tb=sthom*UB*radius/elrms
!       co-moving isotropic-equivalent electron injection luminosity [erg/s]
        xLeinj=exlumel*4.*pi*radius*elrms*c/sthom
!       co-moving isotropic-equivalent proton injection luminosity [erg/s]
        xLpinj=exlumpr*4.*pi*radius*prms*c/sthom
          
        write(6,*) 'B [G], UB [erg/cm^3], R [cm]'
	write (6,1000) Bfield,UB,Radius
	write(6,*) 'Le,inj [erg/s], Lp,inj [erg/s]'
        write(6,1000) xLeinj, xLpinj 
        
!       For blackbody	
	xbb=temperat*boltz/elrms
	yglbb=exlumth
	bblum=4*pi*radius*elrms*c/sthom
	Ubb=bblum/4./pi/radius**2/c
	Ubbrl=astbo*temperat**4

!       Normalization factors         
        factor=tb*qb**(-5./3.)/1836.1
	xnorm=tb*qb**(-5./3.)
! 	write (6,1000) factor,xnorm

	

! bc-to be updated...
      if (ireadin.eq.0) then
      write(98,1002)ireadin,npdec,nsteps,tend
      write(98,1005)gpexmn,gpexmx,xexmn,tol 
      write(98,1000)bplus,tb,qb 
      write(98,1000)ypmin,yemin,ygmin,radius
      write(98,1003)isyn,ipair,icompt,igg,issa,iesc,ipsc,ipion,
     1 ianni,imsyn,ipsyn,iksyn,ineutron
!       write(98,1004)ineutesc,iphotext,xbb,yglbb
      end if
! ec


!     LIMITS FOR THE PARTICLE ENERGY ARRAYS -- NOT TO BE CHANGED
        gpexmn=1./(1.*npdec)
	iprcl=0

	gpmin=10.**gpexmn
	gpmax=10.**gpexmx	
	gemin=gpmin
	gemax=100.*gpmax
	gmmin=gpmin		!m
	gmmax=100.*gpmax	!m
	gpimin=gpmin		!pi
	gpimax=100.*gpmax	!pi
	gkmin=gpmin		!k
	gkmax=100.*gpmax	!k
	
!       total number of bins for the protons/electron arrays
	np=int(npdec*log10(1.01*gpmax))
	ne=int(npdec*log10(1.01*gemax))
	
!       logarithmic step for particle and photon arrays
        deltap=log(gpmax/gpmin)/real(np-1)
	deltax=2.*deltap
	
!       change in the binning of the photon array for high B.
        qbmin=1.E-4*qb
        if(bfield.GT.1.E2)then    
          abin = int(log10(bfield/1.E2))+1
          xmin=qbmin*gemin**2./(10**abin)
          else
	  xmin=qbmin*gemin**2.
        end if
	xlmin=log(xmin)
	xl10min=log10(xmin)
	
!       determine the dimension of the photon array NG
	do n=1,3*ne
	 x(n)=exp(xlmin+(n-1)*deltax)
	 if (x(n).gt.gemax) then
	 ng=n-1
	 goto 300
	 endif 
        enddo
300     continue
	xmax=x(ng)
        ntot=np+ne+ng+1+ne+np	!nt!N

      write (6,*) 'Np=',np
      write (6,*) 'Ne=',ne
      write (6,*) 'Ng=',ng
      write (6,*) 'Ntot=',ntot

      write(6,7000)gpmin,gpmax,gemax
      write(6,7001)xmin,xmax

7000  format(1x,'min/max limits for protons and electrons:',3(g12.5,1x))
7001  format(1x,'min/max limits for photons:',2(g12.5,1x))


!     INITIALIZATION OF THE PARTICLE ARRAYS
      yp(1)=ypmin*log(10.)
      y(1)=yp(1)
      ye(ne)=yemin*log(10.)
      y(np+ne)=ye(ne)
      yg(ng)=ygmin*log(10.)
      y(np+ne+ng)=yg(ng)
      yn(1)=yp(1)		!N
      y(np+ne+ng+1+ne+1)=yn(1)	!N
      ynt(ne)=yemin*log(10.)	!nt
      y(np+ne+ng+1+ne)=ynt(ne)	!nt
      gp(np)=gpmax
      ge(ne)=gemax
      gm(ne)=gmmax		!m
      ym(1)=ypmin*log(10.)      !m
      gpi(ne)=gpimax		!pi
      ypi(1)=ypmin*log(10.)     !pi
      gk(ne)=gkmax		!k
      yk(1)=ypmin*log(10.)      !k

!     background protons/neutrons
      do 1 n=1,np
         gp(n)=gpmin*(gpmax/gpmin)**(real(n-1)/real(np-1))
         yp(n)=(-2.)*log(gp(n)/gp(1))+ypmin*log(10.)-20.
         yn(n)=(-2.)* log(gp(n)/gp(1))+ypmin*log(10.)-40. !N
         y(n)=yp(n)
         y(np+ne+ng+1+ne+n)=yn(n) !N
         ygamma(n)=gp(n)
1     continue 

!     background electrons/muons/pions/kaons     
      do 2 n=ne-1,1,-1
         ge(n)=gemin*(gemax/gemin)**(real(n-1)/real(ne-1))
         gm(n)=gmmin*(gmmax/gmmin)**(real(n-1)/real(ne-1))	!m
         gpi(n)=gpimin*(gpimax/gpimin)**(real(n-1)/real(ne-1))	!pi
         gk(n)=gkmin*(gkmax/gkmin)**(real(n-1)/real(ne-1))	!k
         ye(n)=(ne-n)*2.5*deltap + ye(ne)
         y(n+np)=ye(n)
         ynt(n)=(ne-n)*slelints*deltap + ye(ne)	!nt (we use np for neutrons and ne for neutrinos)
         y(np+ne+ng+1+n)=ynt(n)                                !nt
         ym(n)=(ne-n)*slelints*deltap + ye(ne) -40.		!m
         ypi(n)=(ne-n)*slelints*deltap + ye(ne) -40.		!pi
         yk(n)=(ne-n)*slelints*deltap + ye(ne) -40.		!k
2     continue

         ym(ne)=ye(ne)-40.	        !m
         ypi(ne)=ye(ne)-40.	        !pi
         yk(ne)=ye(ne)-40.	        !k 
         
	 do n=1,ne
	 ygamma(np+n)=ge(n)
	 enddo 
      
!     background photons 
!     ygbb: the external photon field dimensionless number density
!     yglbb: external photon compactness for BB field (read from input file)
!     extph0: external photon compactness for BPL field (read from input file)
!     we set initially ygbb<<yg

!!!      external photons with blackbody (BB) energy spectrum 
      if(iphotext.eq.1)then 
      sum=0.
      do 3 n=1,ng
         yg(n)=-2.*log(x(n)/x(ng)) + yg(ng)
	 ygbb(n)=yg(n)-20.
            if (x(n).gt.30.*xbb) goto 301  
	bbnorm=45.*yglbb/pi**4./xbb**4./xnorm/3. !! bc - check the 3 - ec!!           
 	ygbb(n)=log(bbnorm*x(n)**2./(exp(x(n)/xbb)-1.))
	sum=sum+deltax*x(n)**2.*exp(ygbb(n))
301     continue  
	y(np+ne+n)=log(exp(yg(n))+exp(ygbb(n)))
3       continue
	endif
	
!!!      external photons with broken power law (BPL) spectrum 
        if(iphotext2.eq.1) then
        sum=0.
        do 4 n=1,ng
        yg(n)=-2.*log(x(n)/x(ng)) + yg(ng)
        ygbb(n)=yg(n)-20.
           if(x(n).ge.x1.and.x(n).le.xbr)then
           ygbb(n)=log(x(n)**(-beta1)/xnorm)
           elseif(x(n).gt.xbr.and.x(n).le.x2)then
           ygbb(n)=log(xbr**(beta2-beta1)*x(n)**(-beta2)/xnorm)
           end if
       sum=sum+deltax*x(n)**2.*exp(ygbb(n))
4      continue
!      find correct normalization for BPL 
       bbnorm=extph0/(xnorm*sum) 
       sum=0.
!      add to the photon bg the BPL field with correct normalization       
       do 5 n=1,ng
       yg(n)=-2.*log(x(n)/x(ng)) + yg(ng)
       ygbb(n)=yg(n)-20.
       if (x(n).gt.3.d0*x2) then
       nbbmx=n-1
       endif
           if(x(n).ge.x1.and.x(n).le.xbr)then
	   ygbb(n)=log(bbnorm*x(n)**(-beta1)/xnorm)
	   elseif(x(n).gt.xbr.and.x(n).le.x2)then
	   ygbb(n)=log(bbnorm*xbr**(beta2-beta1)*x(n)**(-beta2)/xnorm)
	   endif
       sum=sum+deltax*x(n)**2.*exp(ygbb(n))	   
       y(np+ne+n)=log(exp(yg(n))+ exp(ygbb(n)))
5      continue
       endif
	
	do  n=1,ng
        ygamma(np+ne+n)=x(n)

        enddo 
 
        
	do n=1,ng
	if (iphotext.eq.1)then
	   if (x(n).gt.30.*xbb) then
	   nbbmx=n-1
	   goto 302
	   endif
	elseif (iphotext2.eq.1)then
	   if (x(n).gt.3.d0*x2) then
	   nbbmx=n-1
	   goto 302
	   endif
	endif    
        enddo
302     continue  

! the last point is reserved for the cool electrons (γ<=1)
	ynecool=ye(1) 
	y(np+ne+ng+1)=ynecool	!N 
	y(np+ne+ng)=-200. 
	     
 
! 	t=0.
! 	iout=nsteps
	
	t=0.
	iout=nsteps

!      In case the program crashes:
!             Pick up where we left off...

         if(ireadin.eq.1)then
         read(99,*) iout,t
         do ijk=1,ntot
            read(99,*)y(ijk)
	    if (ijk.le.np) y(ijk)=y(ijk)-0.0*log(10.)
         enddo
         end if
         
	h=(tend-t)/real(iout)
	
!     INITIALIZATION OF ARRAYS FOR PHOTO-HADRONIC INTERACTIONS
! For photopair
        if (ipair.eq.1) then
         	call betheit(gp,ge,x,np,ne,ng,fpbh,floss)
        endif
        
! For photopion
	If (ipion.eq.1) then
          call PGPI_READER(Tpg0,Rpg01,Rpg02,Rpg03,Rpg04,Rpg05,Rpg13,
     $Rpg14,Rpg15,Rpg16,Rpg17,Rpg18,Rpm02,Rpm15,Rpm18,Tpgm,RatioM,
     $Tpgpi,RatioPi,Tpgk,RatioKph,RatioKe,RatioKm,RatioKpi,RatioKnt,
     $Rpp04,Rpp17,Rpk01,Rpk03,Rpk05,Rpk07,Rpk08,Rpk16,Rpk18,
     $Rpg07,Rpg08,Rpg09,Rpg10)

                  If (imsyn.eq.1.or.ipsyn.eq.1.or.iksyn.eq.1) then
	do i=1,ne
	pms(i) = (1./(6.*pi))*sthom*c*gm(i)**2*bfield**2.
     $ *Tpgm(i)*206.8**(-3.)*1.22E6*0.33 !m
 		if (pms(i).gt.gm(i)) then
		pms(i) = gm(i)
		end if
	  pms(i)=(1-pms(i)/gm(i))  
!pms is the ratio of energy NOT lost by the muon to synchrotron radiation during its lifetime.
!synchrotron losses take up (1-pms) of its energy.
!We divide pms by gm because Tpgm = t_decay*gm
C		write(*,*)i,pms(i)
                  If (ipsyn.eq.1.or.iksyn.eq.1) then
       ppis(i)= (1./(6.*pi))*sthom*c*gpi(i)**2*bfield**2.
     $ *Tpgpi(i)*273.9**(-3.)*1.22E6*0.33 !pk
 		if (ppis(i).gt.gpi(i)) then
		ppis(i) = gpi(i)
		end if
	  ppis(i)=(1-ppis(i)/gpi(i))
                  End If
                  
                  If (iksyn.eq.1) then
	pks(i) = (1./(6.*pi))*sthom*c*gk(i)**2*bfield**2.
     $ *Tpgk(i)*966.7**(-3.)*1.22E6*0.33 !pk
 		if (pks(i).gt.gk(i)) then
		pks(i) = gk(i)
		end if
	  pks(i)=(1-pks(i)/gk(i))
! Added the same for pions and kaons
                  End If
	end do					 
                  End If
  

                  If (ipsyn.eq.2.or.iksyn.eq.2) then
       !! Corrective extrapolation for lower energy tmuon, added SD-18/3/2014 !pk
	 do k=1,ne
	 tttt=0
	 i=ne+1-k
	 do j=1,ne
	 tttt=tttt+Rpp04(i,j-3)
	 end do
	If (tttt.EQ.0) then
          RatioPi(i)= 0.78778953
	do j=1,ne
	  Rpp04(i,j-3) = Rpp04(i+1,j-2)
	  Rpp17(i,j-3) = Rpp17(i+1,j-2)
	 end do
	End If
	end do
                  End If	

                  If (iksyn.eq.2) then
       !! Corrective extrapolation for lower energy tmuon from kaons, added SD-22/3/2014 !pk
	 do k=1,ne
	 tttt=0
	 i=ne+1-k
	 do j=1,ne
	 tttt=tttt+Rpk05(i,j-3)
	 end do
          RatioKe(i)= 0.015
          RatioKph(i)= 0.0145
          RatioKpi(i)= 0.17
          RatioKm(i)= 0.35
          RatioKnt(i)= 0.316
	do j=1,ne
	  Rpk01(i,j-3) = Rpk01(i+1,j-2)
	  Rpk03(i,j-3) = Rpk03(i+1,j-2)
	  Rpk05(i,j-3) = Rpk05(i+1,j-2)
	  Rpk07(i,j-3) = Rpk07(i+1,j-2)
	  Rpk08(i,j-3) = Rpk08(i+1,j-2)
	  Rpk16(i,j-3) = Rpk16(i+1,j-2)
	  Rpk18(i,j-3) = Rpk18(i+1,j-2)
	 end do
	end do
                  End If          

     	End If

!       CALL TO MAIN ROUTINE FOR SOLVING THE SYSTEM OF EQUATIONS	
	ifail=0
	iw=nwork	
      call  d02ejf(t,tend,ntot,y,deriv,d02ejy,tol,relabs,out,
     $             d02ejw,work,iw,ifail)
         
 
         if(tol.le.0)then
                      print*,' warning, no change in solution'
                      tol=abs(tol)
         end if
         
         if(ifail.ne.0) goto 303

         print*,' integration complete to t=',tend
         call cpu_time(t_stop)
         write(*,*) 'Elapsed CPU time [s] = ', t_stop-t_start
         call itime(tarray)
         tt2=tarray(1)*3600.+tarray(2)*60.+tarray(3)
         write(*,*) 'Elapsed wallclock time [s] = ', tt2-tt1
      stop

!                   diagnostics for failed D02... call
 
303   continue

      write(6,6200)tol,tend,t
      write(6,6201)ntot
      
      do n=1,ntot
         write(6,6202)yg(n),y(n)
      enddo

1000  format(1x,6(1pe12.4,1x))
1005  format(3x,4(1pd12.4),2x,i10)
1006  format(2x,i10,2x,1pd12.4)
1002  format(2x,i10,2x,d12.4,2x,i10,2x,d12.4)
1003  format(2x,13(i10))
1004  format(2x,2(i10),2x,2(1pd12.4))
1100  format(1x,g12.4,i5)
1200  format(1x,g12.4,i5,1x,1pe12.4,2x,e12.4)
6200  format(1x,'tol=',g12.4,' tend=',g12.4,' t=',g12.4)
6201  format(1x,'ntot=',i5)
6202  format(1x,e16.8,3x,e16.8)

!!!   END OF MAIN PROGRAM !!!      
      END
      
      
!!! SUBROUTINES !!! 
************************************************************************

      subroutine out(t,y)
!     OUTPUTS THE RESULTS      
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
     
      DIMENSION y(ntotal)
      DIMENSION yp(npmax),ye(npmax),yg(npmax),zd(npmax)
     $         ,gp(npmax),x(npmax),ygamma(ntotal),ge(npmax)
     $	       ,yn(npmax),ynt(npmax),ygbb(npmax)	!N!nt       
      common/param/bplus,tb,qb,gpmax,
     $  deltap,deltax,np,ne,ntot,ng
      common/flags/ipair,isyn,icompt,igg,issa,iesc,
     $ ipsc,ipion,ianni,ikn 
      common/ygm/ygamma
      common/htm/h,tend,iout  
      common/tfl/xnorm
      common/fqq/nsteps
      common/ebb/ygbb,nbbmx
      common/pyy/tauth 
      common/xuu/xen,gmax
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp
      common/pgp/radius,xl10min,bfield

      do n=1,np
	 gp(n)=ygamma(n) 
         yp(n)=y(n)
	 yn(n)=y(np+ne+ng+1+ne+n)	!N
      enddo

      do n=1,ne
	 ge(n)=ygamma(np+n) 
         ye(n)=y(np+n)
	 ynt(n)=y(np+ne+ng+1+n)	!nt
      enddo

      do n=1,ng
        x(n)=ygamma(np+ne+n) 
	if (x(n).lt.0.1) zd(n)=1.
	if (x(n).ge.0.1.and.x(n).le.1.) zd(n)=(1.-x(n))/.9
	if (x(n).gt.1.) zd(n)=0.
        yg(n)=y(np+ne+n)-log(1.+tauth*zd(n)/3.)
      enddo
	

!  the one before last bin is reserved for the cooled electrons
	ynecool=y(np+ne+ng+1)	!N
!  measure of the Thomson optical depth on cold electrons	
	tauth=exp(ynecool) 

         rewind 99
         write(99,*)iout+1,t
         do  ijk=1,ntot
            write(99,*)y(ijk)
         enddo

!  sump is a measure of the energy density in protons
!  sumnp is a measure of the density in protons
!  sume is a measure of the energy density in electrons
!  send is a measure of the density in electrons
!  sumx is a measure of the energy density in photons
!  sxnd is is a measure of the density in photons

	   sump=0.
	   sumnp=0.
	   sume=0.
	   send=0.
           sumx=0.
	   sxnd=0. 
	   
	do k=1,np
	   if (k.eq.1.or.k.eq.np) then
	      fact=.5
		else 
	      fact=1.
	   end if	  
	   sump=sump+fact*deltap*gp(k)**2.*exp(yp(k))
           sumnp=sumnp+fact*deltap*gp(k)*exp(yp(k))
        enddo

        do k=1,ne
            if (k.eq.1.or.k.eq.ne) then
	      fact=.5
		else 
	      fact=1.
	   end if 
	   sume=sume+fact*deltap*ge(k)**2.*exp(ye(k))
	   send=send+deltap*ge(k)*exp(ye(k))
        enddo 
        
	 do k=1,ng
	   if (k.eq.1.or.k.eq.ng) then
	      fact=.5
		else 
	      fact=1.
	   end if 
	   sumx=sumx+fact*deltax*x(k)**2.*exp(yg(k))
	   sxnd=sxnd+fact*deltax*x(k)*exp(yg(k))
        enddo 
           xen=sumx*xnorm

      write (6,1000) t,sump,xen,tauth
      write (90,1000) t,sump,xen,tauth 
      write (98,1000) t,sump,sume,xen,sumnp,send,sxnd
      write (72,1000) t,log10(xen)
      write (55,1000) t,log10(bfield)!FAKE TIMECURVE
      write (73,1000) t,log10(1836.1*sump) 
      
      write (81,1000) t,sump,xen,tauth	
      write (85,1000) t,sump,xen,tauth	
	do n=1,ng	
	if(exp(yg(n))-exp(ygbb(n)).gt.0.)then
	write (71,1000) log10(x(n)),log10(xnorm*(exp(yg(n))-exp(ygbb(n)))*x(n)**2),
     $	log10(xnorm*exp(ygbb(n))*x(n)**2)
	else
    	write (71,1000) log10(x(n)),log10(xnorm*(exp(yg(n)))*x(n)**2),
     $	log10(xnorm*exp(ygbb(n))*x(n)**2)
	endif
	write (81,1000) log10(x(n)),log10(xnorm*exp(yg(n))*x(n)**2)
	write (85,1000) log10(x(n)),log10(xnorm*(exp(yg(n))-exp(ygbb(n)))*x(n)**2)
        enddo

!      write (88,1000) t,sump,xen,tauth
!	do n=1,np
!        write (88,1000) log10(gp(n)), log10(gp(n)**2*exp(yp(n)))
!        enddo 
        
      write (89,1000) t,sump,xen,tauth
	do n=1,ne
        write (89,1000) log10(ge(n)), log10(ge(n)**2*exp(ye(n)))
        enddo 
        
!      write (82,1000) t,sump,xen,tauth		 
!	do n=1,ne
!	write (82,1000) log10(ge(n)),log10(ge(n)**2*exp(ynt(n))) !nt
!	enddo

!     write (87,1000) t,sump,xen,tauth
!	do n=1,np
!        write (87,1000) log10(gp(n)), log10(gp(n)**2*exp(yn(n))) !N
!        enddo

	do n=1,ntot
	if (n.lt.ntot) then
	  if (n.le.np+ne) then
	    tau=0.
	    else
	    tau=tauth*zd(n-np-ne)/3.
	  end if
	endif
        enddo

	t=tend-real(iout)*h
	iout=iout-1
	
1000  format(1x,6(1pe12.4,1x))
	return
	end
	
**************************************************************************

      subroutine deriv(time,y,yprime)
!! COMPUTES TIME DERIVATIVES OF ALL PHYSICAL PROCESSES
!! subroutine called by d02ejf

      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer j1,j2,k1
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,afine=7.3e-3,csrt=.2,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 

      DIMENSION y(ntotal),yprime(ntotal),ygamma(ntotal)
      DIMENSION yp(npmax),ye(npmax),yg(npmax),ge(npmax)
     $         ,gp(npmax),x(npmax),yn(npmax),ynt(npmax)
     $         ,ym(npmax),gm(npmax),ypi(npmax),gpi(npmax)
     $         ,yk(npmax),gk(npmax),pks(npmax),ygbb(npmax)
     $         ,yinit(npmax),tne(npmax),tnm(npmax),tna(npmax)
     $         ,gcsth(npmax),gcskn(npmax),gdens(npmax)
     $         ,cskn(npmax),syncel(npmax),syncph(npmax)
     $         ,csth(npmax),denkn(npmax),zd(npmax),relann(npmax)
     $         ,sst(npmax),extelin(npmax),extprinxp(npmax)
     $         ,sumsy(npmax),p(npmax,npmax)
     $         ,prloee(npmax),ypel(npmax),pgeera(npmax)
     $         ,pgray(npmax),fpgloss(40)
     $         ,protpg(npmax),photpg(npmax),elecpg(npmax)
     $         ,trinopg(npmax),tneutpg(npmax),tmuonpg(npmax) 
     $         ,tkaonpg(npmax),tkaonpgt(npmax),tmuonpgt(npmax)
     $         ,tpionpgt(npmax),trinopge(npmax)
     $         ,trinopgm(npmax),trinopgt(npmax)
      DIMENSION syncprph(npmax),sumpsy(npmax),syncpr(npmax)	 
     $   ,syncm(npmax),summsy(npmax)  
     $   ,syncmph(npmax),pms(npmax)             
     $   ,syncpi(npmax),sumpisy(npmax) 
     $   ,syncpiph(npmax),ppis(npmax)       
     $   ,synck(npmax),sumksy(npmax)   
     $   ,synckph(npmax)              
     
      DIMENSION  yggal(npmax),yggcral(npmax)
     $   ,anng(npmax),annel(npmax)
     $   ,ggabsr(npmax),gginje(npmax),pgeein(npmax)      
  

      common/param/bplus,tb,qb,gpmax,
     $ deltap,deltax,np,ne,ntot,ng
      common/flags/ipair,isyn,icompt,igg,issa,iesc,
     $ ipsc,ipion,ianni,ikn
      common/ygm/ygamma
      common/tfl/xnorm
      common/nonlin/factor      
      common/phexter/gext,ygext,iphotext 
      common/ext/x1,x2,xbr,beta1,beta2,extph0,iphotext2 
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp  
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ ap,iprext,ipexp  
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr  
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr       
      common/ebb/ygbb,nbbmx 
      common/pyy/tauth 
      common/xuu/xen,gmax
      common/tvv/iprcl  
      common/prl/fpgloss
      common/pgp/radius,xl10min,bfield
      
      
  !     FAKE TIMECURVE
       DIMENSION fktime(100000),fakeB(100000)
      common/fake/tol,ifake
      common/tm/fktime,fakeB,dfk
      if(time.gt.fktime(1))then
        do i=1,ifake
        if(time.ge.fktime(i).and.time.lt.fktime(i+1))then
        dfk=fktime(i+1)-fktime(i)
        Bint=(log10(fakeB(i+1))-log10(fakeB(i)))/dfk*
     &  (time-fktime(i))+log10(fakeB(i))
!         write(6,*) xlrint,i,fktime(i)
        endif 
        enddo 
        bfield=10.**Bint
!        write(55,1000)  time , bfield 
       endif

	bplus=belesc
!       dimensionless B field	
	qb=bfield/Bcr
!       magnetic energy density [erg/cm^3]	
	UB=bfield**2./8./pi
!       magnetic compactness	
	tb=sthom*UB*radius/elrms
!       Normalization factors         
        factor=tb*qb**(-5./3.)/1836.1
	xnorm=tb*qb**(-5./3.)
! 	write (6,1000) factor,xnorm       


	iprcl=iprcl+1
      
      do n=1,np
	 gp(n)=ygamma(n) 
         yp(n)=y(n)
	 yn(n)=y(np+ne+ng+1+ne+n)	!N
      enddo 

      do n=1,ne
	 ge(n)=ygamma(np+n) 
         ye(n)=y(np+n)
	 ynt(n)=y(np+ne+ng+1+n)		!nt
      enddo

	gemin=ge(1)
	gemax=ge(ne)

      do n=1,ng
	 x(n)=ygamma(np+ne+n) 
         yg(n)=y(np+ne+n)
	if (x(n).lt..1) zd(n)=1.
	if (x(n).ge..1.and.x(n).le.1.) zd(n)=(1.-x(n))/.9
 	if (x(n).gt.1.) zd(n)=0.
        enddo
	zd(ng)=0.
 
!  the one before last bin is reserved for the cooled electrons
	ynecool=y(np+ne+ng+1)	!N
!  measure of the Thomson optical depth on cold electrons	
	tauth=exp(ynecool) 
	if (iprcl.eq.1) then
	n=np+ne+ng+1+ne+np	!N
	yinit(n)=y(n)
	endif

!!! RATES OF PHYSICAL PROCESSES !!!

!! PROTON INJECTION 
        call prinj(gp,yp,np,deltap,extprinxp,sumprinj)
!! ELECTRON INJECTION
        call elinj(ge,ye,ne,deltap,extelin,sumelinj)        
!! SYNCHROTRON RADIATION
        if(isyn.eq.1) then
! i) electron loss   
        call synecool(ge,ye,syncel,sumsyncel)
! ii) photons from electron synchrotron
        call synephot(ge,ye,x,yg,syncph,sumsyncph)
! iii) proton losses
        call synpcool(gp,yp,syncpr,sumsyncpr)
! iv) photons from proton synchrotron        
        call synpphot(gp,yp,x,yg,syncprph,sumsyncprph)
        endif 
!! SYNCHROTRON SELF ABSORPTION
        if(issa.eq.1) then
! photon loss (attenuation)
        call ssaphot(ge,ye,x,yg,sst,sumphssa)
        endif 
!! INVERSE COMPTON SCATTERING (ICS)        
        if (icompt.eq.1) then
! electron loss
        call compe_sim(x,ge,ye,yg,ikn,csth,cskn,sumcsth,
     $ sumcskn,taukn) 
! photons from electron ICS using Bloumenthal-Gould emissivity
        call compph_bg(x,ge,yg,ye,gcsth,sumgcsth)
        endif
!! PHOTON-PHOTON PAIR PRODUCTION
        if (igg.eq.1) then 
! photon attenuation
        call ggphot(x,yg,ggabsr,sumggabs)       
! electron/positron injection 
        call ggelec(x,yg,ge,ye,gginje,sumggabs,sumgginj)  
        endif 
!! PROTON-PHOTON PAIR PRODUCTION
	if (ipair.eq.1) then
! proton loss 
        call bhloss(x,gp,yg,yp,prloee,sumprloee) 
! electron/positron injection
        call bhelec(x,gp,ge,yg,yp,ye,pgeera,sumelBH)
	endif 
!! PROTON-PHOTON PION INTERACTION
        if (ipion.eq.1) then
! proton/neutron loss and injection of secondaries
        call pgall(gp,yp,yn,ge,ye,ynt,x,yg
     $,protpg,elecpg,photpg
     $,tneutpg,trinopg,trinopge,trinopgm
     $,tmuonpg,syncmph,syncpiph, synckph)
        endif
        
!!! THE EQUATIONS !!!
!       derivatives computed at first bin

!       protons:  /FAKE TIMECURVE iesc=0
        yprime(1)=0
!       electrons:
        yprime(np+1)=extelin(1)-iesc*belesc+syncel(1)*isyn+
     $  (csth(1) - cskn(1))*icompt + igg*gginje(1) +
     $  ipair*pgeera(1) + ipion*elecpg(1) 
!       photons:
 	yprime(np+ne+1)=-ipsc*1./(1.+tauth*zd(1)/3.)+
     $	isyn*(syncph(1) +0*syncprph(1))+ issa*sst(1) +  !FAKE TIMECURVE
     $  icompt*gcsth(1) -igg*ggabsr(1)+ipion*photpg(1) +
     $  imsyn*syncmph(1) +ipsyn*syncpiph(1) +iksyn*synckph(1)       
!       cold electrons: 
	pycsth=4./3.*xnorm*gdens(1)*(gp(1)**2.-1.)*exp(ye(1))/tauth
	pysyn= 4./3.*tb*(gp(1)**2.-1.)*exp(ye(1))/tauth
        yprime(np+ne+ng+1)= pycsth+pysyn+taukn
     $    -3.*tauth/32.-belesc	     
!       neutrons:	!N 
        yprime(np+ne+ng+1+ne+1)=0.
!       neutrinos:	!nt 
        yprime(np+ne+ng+1+1)= -ipsc
CX       muons:		!m
CX     yprime(np+ne+ng+1+ne+np+1)=0.

!       derivatives computed at last bin

!       protons: 
	yprime(np)=0.
!       yprime(np)=extprinxp(np)-bpresc+syncpr(np)*isyn	
        if (y(np).lt.-150.) yprime(np)=0.
!       electrons:
	yprime(np+ne)=-.1
!       yprime(np+1)=extelin(ne)-belesc+syncel(ne)*isyn	
!        photons:
        yprime(np+ne+ng)=-ipsc*1./(1+tauth*zd(ng)/3.)+   
     $	isyn*(syncph(ng) +0*syncprph(ng)) +issa*sst(ng) +  !FAKE TIMECURVE
     $  icompt*gcsth(ng) -igg*ggabsr(ng) +ipion*photpg(ng) +
     $  imsyn*syncmph(ng) +ipsyn*syncpiph(ng) +iksyn*synckph(ng) 
!       neutrons:	!N
        yprime(np+ne+ng+1+ne+np)= .0
!       neutrinos:	!nt
        yprime(np+ne+ng+1+ne)= -ipsc
CX      muons:		!m
CX     yprime(np+ne+ng+1+ne+np+ne)=0.

***************************

!       derivatives computed at all other bins 

	sumprot=0.
	sumnepio=0.
	sumprpio=0.
	
!! PROTONS + NEUTRONS
	do n=2,np-1!    
! protons	/FAKE TIMECURVE iesc=0
	yprime(n)=0  
! neutrons     
 	yprime(np+ne+ng+1+ne+n) = ipion*tneutpg(n) -ipsc

	sumprpio = sumprpio+deltap*gp(n)**2.*exp(yp(n))*
     $  ratmpe*protpg(n)
        sumnepio = sumnepio+deltap*gp(n)**2.*exp(yn(n))*
     $  ratmpe*tneutpg(n)
        sumprot=sumprot+deltap*gp(n)**2.*exp(yp(n))*yprime(n)
        enddo 

        sumelec=0.
        sumelpio=0.
        sumntpio=0.
        summpio=0. 
     
!! electrons       
        do n=2, ne-1 
        yprime(np+n)=extelin(n)-iesc*belesc+
     $  isyn*syncel(n)+
     $  icompt*(csth(n) - cskn(n)) +
     $  igg*gginje(n) + ipair*pgeera(n) +
     $  ipion*elecpg(n)
!! neutrinos     
	yprime(np+ne+ng+1+n) = ipion*trinopg(n) -ipsc		!nt
        tna(n)=trinopg(n)*ipsc ! all flavor neutrinos
	tne(n)=trinopge(n)*ipsc ! electron neutrinos 
	tnm(n)=trinopgm(n)*ipsc ! muon neutrinos
!	write (86,1000) log10(ge(n)),log10(ge(n)**2.*exp(ynt(n))*tna(n)) !CHANGE_MARK 85->86
!	write (83,1000) log10(ge(n)),log10(ge(n)**2.*exp(ynt(n))*tne(n))      !FAKE TIMECURVE
!	write (84,1000) log10(ge(n)),log10(ge(n)**2.*exp(ynt(n))*tnm(n))  
CX	yprime(np+ne+ng+1+ne+np+n) = tmuonpg(n) -syncm(n)- belesc	!m
	if (y(np+n).lt.-200.) yprime(np+n)=0.
	if (y(np+ne+ng+1+n).lt.-200.) yprime(np+ne+ng+1+n)=0.	!nt

        sumelpio=sumelpio+deltap*ge(n)**2.*exp(ye(n))*elecpg(n)
	sumntpio=sumntpio+deltap*ge(n)**2.*exp(ynt(n))*trinopg(n)	        
 	summpio = summpio+deltap*gm(n)**2.*exp(ym(n))*
     $    206.8*tmuonpg(n)   
	sumelec=sumelec+deltap*ge(n)**2.*exp(ye(n))*yprime(np+n)
        enddo 
 
************************
!   PHOTONS
	
	sumphot=0.
	sumgrpio=0.
	sumsynckph=0.
	sumsyncpiph=0.
	sumsyncmph=0.
	
	do n=1,ng
   !FAKE TIMECURVE
	yprime(np+ne+n)=-ipsc*1./(1.+tauth*zd(n)/3.) +
     $	isyn*(syncph(n) + 0.*syncprph(n)) +issa*sst(n) +                                         
     $  icompt*gcsth(n) -igg*ggabsr(n) +ipion*photpg(n) +
     $  imsyn*syncmph(n) +ipsyn*syncpiph(n) +iksyn*synckph(n) 

     	sumsyncmph=sumsyncmph+xnorm*deltax*x(n)**2.*exp(yg(n))*
     $	syncmph(n)
        sumsyncpiph=sumsyncpiph+xnorm*deltax*x(n)**2.*exp(yg(n))*
     $  syncpiph(n)
	sumsynckph=sumsynckph+xnorm*deltax*x(n)**2.*exp(yg(n))*
     $  synckph(n)
	sumgrpio=sumgrpio+xnorm*deltax*x(n)**2.*exp(yg(n))*photpg(n)     
        sumphot=sumphot+deltax*x(n)**2.*xnorm*exp(yg(n))*yprime(np+ne+n)

! treatment of external photon field     
        if ((iphotext2.eq.1.and.x(n).le.x2.and.x(n).ge.x1).
     $  or.(n.le.nbbmx.and.iphotext.eq.1)) then
	if (ygbb(n).gt.yg(n)) then
           yprime(np+ne+n)=-1./(1.+tauth*zd(n)/3.) + exp(ygbb(n)-yg(n))
	else
           yprime(np+ne+n)=yprime(np+ne+n)+exp(ygbb(n)-yg(n))
        endif
        endif

        enddo
        

	yprime(np+ne+ng)=0.
     
	fend=99999.
	yprmax=abs(yprime(1))
	nmx=1

	do n=1,np+ne+ng+1+ne+np		!nt
	if (yprmax.lt.abs(yprime(n))) then
	yprmax=abs(yprime(n))
	nmx=n
	endif
        enddo

	sumprtot=0.
	sumeltot=0.
	sumphtot=0.
        sumptotl=0.


	do ik=1,500000
*%%
        if (ik*100.eq.iprcl) then
	sumptotl=sumelpio+sumgrpio+sumntpio+sumnepio
C +summpio
	sumpioloss=sumptotl/sumprpio 
C  Output on screen:
	write (6,1400) iprcl,time,nmx,y(nmx),yprime(nmx),tauth
	write (6,1000) bfield, log10(bfield),sst(1), tol 
!	write (6,1000) sumprot*ratmpe,sumelec,sumphot, sumnepio*ratmpe
	write (6,1000) sumsyncel,sumsyncph
!	write (6,1000) sumsyncpr*ratmpe, sumsyncprph
        write (6,1000) sumcsth,sumcskn,sumcsth+sumcskn, sumgcsth
     	write (6,1000) sumgginj,sumggabs
!	write (6,1000) sumprloee,sumelBH
!	write (6,1000) sumprpio,sumelpio,sumgrpio,sumntpio,sumnepio
!     $ ,summpio
! 	write (6,1000) sumsyncmph, sumsyncpiph
! 	write (6,1000) sumelann,sumanng
	write (6,*) 'loss %',sumelpio/sumptotl,sumgrpio/sumptotl,
     $sumntpio/sumptotl,sumnepio/sumptotl, sumptotl
            endif
           enddo 


!          
!        sumprtot= sumyprimalt+sumprloee+sumprpio+sumprprlo
!        sumeltot=sumsyncel+sumcsth+sumcskn+sumelpio+sumannel
!      1     +sumelupsc+sumelinpp+sumgginje+sumelbrms+sumpgeein+sumelss
!        sumphtot=sumsyncph+sumphesc+sumphssa+sumgcsth+sumgcskn
!      1	   +sumggabsr+sumgrpio+sumphinpp+sumphbrms
!      2     +sumanng

 
 
	
2001    format ('protlspi=',1pd12.4,2x,'neutlspi=',d12.4,2x,
     1    'empioinj=',d12.4)
2002    format ('protlspe=',1pd12.4,2x,'injpee=',d12.4
     1    ,2x,'inj/totls=',d12.4)

 
	
1400  format(1x,i10,1x,1pe12.4,i5,1x,1pe12.4,1x,e12.4,1x,e12.4,1x,e12.4
     1    ,1x,e12.4)
1300  format(1x,1pe12.4,i5,1x,i5,1x,i5,2x,e12.4,2x,e12.4)
 
!  
 
1000  format(1x,6(1pe12.4,1x))
1009    format(/)

	do 377 i=1,ntot
	   if( y(i).gt.1. .or. y(i).le.1.)then
              continue
           else
           print*,'overflow!!! i=',i,' y(i)=',y(i)
           stop
           end if
	   if( yprime(i).gt.1. .or. yprime(i).le.1.)then
              continue
           else
           print*,'overflow!!! i=',i,' yprime(i)=',yprime(i)
           stop
           end if
377      continue

	tprev=time

      return
      end
************************************************************************

!!    PROTON INJECTION   
      subroutine prinj(gp,yp,np,deltap,extprinxp,sumprinj)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer ipexp
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION gp(np),yp(np),extprinxp(np)
     
      common/prexter/gpextmn,gpextmx,slprints,exlumpr,bpresc,
     $ap,iprext,ipexp
      common/prexterbr/gpextbr,slprints1,slprints2,ap2,iprextbr 

! sumprinj: energy injected into protons   
      sumprinj=0.
      if (iprext.eq.1) then
	 slpr=slprints
         q2=(gpextmx**(2.-slpr)-gpextmn**(2.-slpr))/(2.-slpr)
         qextpr=3.*exlumpr/q2
         do n=1,np   
         if (gp(n).ge.gpextmn.and.gp(n).lt.gpextmx) then
! 	if (gp(n).ge.gpextmn) then 
        extprinxp(n)=qextpr/exp(yp(n))*gp(n)**(-slprints)
     $  *exp(-(gp(n)/gpextmx)**ap)**ipexp
         else
          extprinxp(n)=0.
         endif         
       sumprinj=sumprinj+deltap*gp(n)**2.*exp(yp(n))*extprinxp(n)
!        write(6,1000) gp(n), extprinxp(n), sumprinj, ap, ipexp
         enddo         
      endif
  
      if (iprextbr.eq.1) then
 	 slpr1=slprints1
	 slpr2=slprints2     
	 qp21=(gpextbr**(2.-slpr1)-gpextmn**(2.-slpr1))/(2.-slpr1)
         qp22=(gpextmx**(2.-slpr2)-gpextbr**(2.-slpr2))/(2.-slpr2)
	 fogp=gpextbr**(slpr1-slpr2)
	 qextpr=3.*exlumpr/(fogp*qp21+qp22)
	 do n=1,np
	 if (gp(n).ge.gpextmn.and.gp(n).le.gpextmx) then
	   if (gp(n).lt.gpextbr) then
	   extprinxp(n)=fogp*qextpr/exp(yp(n))*gp(n)**(-slpr1)
	   else
	 extprinxp(n)=qextpr/exp(yp(n))*gp(n)**(-slpr2)
     $   *exp(-(gp(n)/gpextmx)**ap2)**ipexp
	   end if  
         else
         extprinxp(n)=0.
         endif        
        sumprinj=sumprinj+deltap*gp(n)**2.*exp(yp(n))*extprinxp(n)
!         write(6,1000) gp(n), extprinxp(n), sumprinj
         enddo	
      endif 
      
      return 
1000  format(1x,6(1pe12.4,1x))      
      end
*************************************************************************

!!    ELECTRON INJECTION   
      subroutine elinj(ge,ye,ne,deltap,extelin,sumelinj)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      integer ieexp
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION ge(ne),ye(ne),extelin(ne)
      
      common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ae,ielext,ieexp  
      
      common/elexterbr/geextbr,slelints1,slelints2,ae2,ielextbr


! sumelinj: energy injected into electrons    
      sumelinj=0.

      if (ielext.eq.1) then
         slel=slelints
         q2=(geextmx**(2.-slel)-geextmn**(2.-slel))/(2.-slel)
         qextel=3.*exlumel/q2
         do n=1,ne  
         if(ge(n).ge.geextmn.and.ge(n).le.geextmx)then
! 	if (ge(n).ge.geextmn) then 
        extelin(n)=qextel/exp(ye(n))*ge(n)**(-slel)
     $  *exp(-(ge(n)/geextmx)**ae)**ieexp
         else
         extelin(n)=0.
         endif         
       sumelinj=sumelinj+deltap*ge(n)**2.*exp(ye(n))*extelin(n)
!        write(6,1000) ge(n), extelin(n)
         enddo         
      endif 
       
      if (ielextbr.eq.1) then
	 slel1=slelints1
	 slel2=slelints2    
         q21=(geextbr**(2.-slel1)-geextmn**(2.-slel1))/(2.-slel1)
         q22=(geextmx**(2.-slel2)-geextbr**(2.-slel2))/(2.-slel2)
         fog=geextbr**(slel1-slel2)
	 qextel=3.*exlumel/(fog*q21+q22)
	 do n=1,ne
	 if (ge(n).ge.geextmn.and.ge(n).le.geextmx) then
	   if (ge(n).lt.geextbr) then
	   extelin(n)=fog*qextel/exp(ye(n))*ge(n)**(-slel1)
	   else
	 extelin(n)=qextel/exp(ye(n))*ge(n)**(-slel2)
     $   *exp(-(ge(n)/geextmx)**ae2)**ieexp
	   end if  
         else
         extelin(n)=0.
         endif        
        sumelinj=sumelinj+deltap*ge(n)**2.*exp(ye(n))*extelin(n)
!        write(6,1000) ge(n), extelin(n), sumelinj
         enddo	
      endif 
      
      return 
1000  format(1x,6(1pe12.4,1x))   
      end       
    
*************************************************************************    

!!    ELECTRON SYNCHROTRON LOSSES
      subroutine synecool(ge,ye,syncel,sumsyncel)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION ge(ne),ye(ne),syncel(ne)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      
! sumsyncel: energy lost by electrons in synchrotron
	sumsyncel=0.

	do n=1,ne-1 
	syndf=(ge(n+1)**2.-1.)*exp(ye(n+1))-(ge(n)**2.-1.)*exp(ye(n))
	syncel(n)=4.*tb*syndf/deltap/ge(n)/exp(ye(n))/3.
        sumsyncel= sumsyncel+deltap*ge(n)**2.*syncel(n)*exp(ye(n))
        enddo
        return
      end         
*************************************************************************

!!    PROTON SYNCHROTRON LOSSES
      subroutine synpcool(gp,yp,syncpr,sumsyncpr)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
      DIMENSION gp(np),yp(np),syncpr(np)

      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng

!  sumsyncpr: energy lost by protons in synchrotron
	sumsyncpr=0.

	do n=1,np-1
	syndf=(gp(n+1)**2.-1.)*exp(yp(n+1))-(gp(n)**2.-1.)*exp(yp(n))
	syncpr(n)=4.*tb*syndf/deltap/gp(n)/exp(yp(n))/3./ratmpe**3.
        sumsyncpr=sumsyncpr+deltap*gp(n)**2.*syncpr(n)*exp(yp(n))
        enddo
        return
      end
************************************************************************* 

!!    PHOTON ATTENUATION DUE TO SYNCHROTRON SELF ABSORPTION
      subroutine ssaphot(ge,ye,x,yg,sst, sumphssa)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,afine=7.3e-3,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)
      DIMENSION ge(ne),ye(ne),x(ng),yg(ng),sst(ng)
      
       common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm
       common/elexter/geextmn,geextmx,slelints,exlumel,belesc,
     $ ae,ielext,ieexp  !FAKE TIMECURVE
!  uses delta-function approximation -- see equation (40) in Mastichiadis & Kirk (1995)
!  sumphssa: energy absorbed in synchrotron self absorption
        sumphssa=0.
!START FAKE TIMECURVE --------------------------------------------------
	xlminnew=log10(qb*(geextmn)**2.)
        xlmin=log10(qb*(ge(1))**2.)
        eldvthres=-200. !threshold for elderiv before set to zero
        sstthres=5.e11  !highest allowed calculated value for sst(n)
        ansst=0.0 !index of sst slop behind delta function of electrons 4/3-(ansst)=2

!       determine the dimension of the photon array NG
	do n=1,ng
	 xltmp=log10(x(1))+(n-1)*log10(exp(deltax))
	 if (xlminnew.lt.xltmp) then
	 nm=n-1
	 goto 322
	 endif 
        enddo
322     continue  

	do n=1,ng
	 xltmp=log10(x(1))+(n-1)*log10(exp(deltax))
	 if (xlmin.lt.xltmp) then
	 nmin=n-1
	 goto 367
	 endif 
        enddo
367     continue

       if (nm.lt.nmin) nm=nmin
      !SSA FOR MINIMUM GAMMA BIN
       elderiv=exp(ye(nm-nmin+2))/ge(nm-nmin+2)**(2.)/
     $ deltap-exp(ye(nm-nmin+1))/ge(nm-nmin+1)**(2.)/deltap
       if (elderiv.gt.0.) elderiv=-elderiv
       if (log10(-elderiv).lt.eldvthres) elderiv=0.
       if (1/elderiv.eq.0.) elderiv=0.
       sst(nm+1)=pi/6./afine/qb**(.5)/x(nm+1)**(.5)/ge(nm-nmin+1)*
     $ elderiv
       sstgmin=sst(nm+1)
       sumphssa=sumphssa+xnorm*deltax*x(nm+1)**2.*exp(yg(nm+1))*
     $ sst(nm+1)
	
       if (-sstgmin.gt.sstthres) then
        sstgmin=-sstthres
        sst(nm+1)=-sstthres
       endif

      !SSA FOR BINS BELOW MINIMUM GAMMA BIN LOOP 1
       do n=1,nm
       !Artificial power for photon bins that don't match electron's bins of power law
       sst(n)=sstgmin*(x(n)/x(nm+1))**(ansst) 
       sumphssa=sumphssa+xnorm*deltax*x(n)**2.*exp(yg(n))*sst(n)
       enddo

       !LOOP 2 SSA FOR DISTRIBUTION OF ELECTRONS
       do n=nm+2,ng-1  !FAKE TIMECURVE do n=1,ng-1
       if (ge(n-nmin).eq.0.) then
        sst(n)=0.
       else   
       elderiv=(exp(ye(n-nmin+1))/ge(n-nmin+1)**2.
     $ -exp(ye(n-nmin))/ge(n-nmin)**2.)/deltap
       if (elderiv.gt.0.) elderiv=-elderiv
       if (log10(-elderiv).lt.eldvthres) elderiv=0.
       if (1/elderiv.eq.0.) elderiv=0.
       sst(n)=pi/6./afine/qb**(.5)/x(n)**(.5)/ge(n-nmin)*elderiv
       endif
       if (x(n).gt.1.) sst(n)=0.
       if (-sst(n).gt.-sst(nm+1)) sst(n)=sst(nm+1)
       sumphssa=sumphssa+xnorm*deltax*x(n)**2.*exp(yg(n))*sst(n)
       enddo  

  !      do n=1,ng-1
  !      write(6,*) sst(n),ge(n-nmin),ye(n-nmin),x(n),yg(n)
   !     enddo
!END FAKE TIMECURVE ----------------------------------------------------
!ORIGINAL SCRIPT
!       do n=1,ng-1 
!       if (ge(n).eq.0.) then
!       sst(n)=0.
!       else         
!       elderiv=(exp(ye(n+1))/ge(n+1)**2.-exp(ye(n))/ge(n)**2.)/
!     $ deltap
!       if (elderiv.gt.0.) elderiv=-elderiv
!       sst(n)=pi/6./afine/qb**(.5)/x(n)**(.5)/ge(n)*elderiv
!       endif
!       if (x(n).gt.1.) sst(n)=0.
!       sumphssa=sumphssa+xnorm*deltax*x(n)**2.*exp(yg(n))*sst(n)
!       enddo    
      return
      end
*************************************************************************  

!!    ELECTRON INVERSE COMPTON LOSSES 
!     using 2 options for electron cooling
      subroutine compe_sim(x,ge,ye,yg,ikn,csth,cskn,sumcsth,
     $ sumcskn,taukn) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION ge(ne),ye(ne),x(ng),yg(ng), csth(ne), cskn(ne),
     $ gdens(ng),denkn(ng),p(npmax,npmax)
       
       common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 
       common/pyy/tauth 

       call densityold(x,ge,yg,ne,ng,deltax,gdens,denkn)
       
! sumcsth: electron energy losses due to ICS in Thomson regime
! sumcskn: electron energy losses due to ICS in KN regime
!  taukn: number of KN collisions which produce cool pairs

        sumcsth=0.
        sumcskn=0.
        taukn=0.
        if (ikn.eq.0) then
!  uses simplified scheme for KN losses -- see equation (45) in Mastichiadis & Kirk (1995)        
        do j=1, ne-1
	csdfth=(ge(j+1)**2.-1.)*exp(ye(j+1))-
     $	(ge(j)**2.-1.)*exp(ye(j))
	csth(j)=4.*tb*csdfth*qb**(-5./3.)*gdens(j)/
     $          deltap/ge(j)/exp(ye(j))/3. 
        cskn(j)=xnorm*denkn(j)/ge(j)
        sumcsth=sumcsth+deltap*ge(j)**2.*exp(ye(j))*csth(j)
	sumcskn=sumcskn-deltap*ge(j)**2.*exp(ye(j))*cskn(j)
	taukn=taukn+deltap*cskn(j)*ge(j)*exp(ye(j))/tauth	
!         taukn=taukn+deltap*xnorm*denkn(j)*exp(ye(j))/tauth	
        enddo 
        endif 
                
        if(ikn.eq.1)then
!  uses full expression from Bloumenthal & Gould (1970)   
        call bgel(x,ge,yg,ye,ne,ng,deltap,deltax,p) 
        
        do j=1, ne-1
       
	sume1=0.
	sume2=0.

	do m=1,j
	if (m.eq.1.or.m.eq.j) then
	fcr=.5
	else
	fcr=1.
	endif
	fp1=fcr*p(j,m)*ge(m)*deltap
	sume1=sume1+fp1 
        enddo 

	do m=j+1,ne-1
	if (m.eq.1.or.m.eq.j) then
	fcr=.5
	else
	fcr=1.
	endif
	fp2=fcr*exp(ye(m))*ge(m)*p(m,j)*deltap
        sume2=sume2+fp2
        enddo

	csdf=(sume1-sume2/exp(ye(j)))
        cskn(j)=csdf 
        csdfth=(ge(j+1)**2.-1.)*exp(ye(j+1))-
     $	(ge(j)**2.-1.)*exp(ye(j))
	csth(j)=4.*tb*csdfth*qb**(-5./3.)*gdens(j)/
     $          deltap/ge(j)/exp(ye(j))/3.  
        sumcsth=sumcsth+deltap*ge(j)**2.*exp(ye(j))*csth(j)
	sumcskn=sumcskn-deltap*ge(j)**2.*exp(ye(j))*cskn(j)
	taukn=taukn+deltap*cskn(j)*ge(j)*exp(ye(j))/tauth! to be checked
        enddo 
        endif 

        return
      end
************************************************************************* 
 
!!    PHOTON LOSSES DUE TO PHOTON-PHOTON PAIR PRODUCTION
      subroutine ggphot(x,yg,ggabsr,sumggabs) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION x(ng), yg(ng), ggabsr(ng), yggal(ng)

       common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 
       
        call ggabsal(deltax,x,yg,ng,yggal)
        
! sumggabs: energy lost due to photon-photon absortion
        sumggabs=0.
	do n=1,ng
        ggabsr(n)=tb*qb**(-5./3.)*yggal(n)
	sumggabs=sumggabs+deltax*x(n)**2.*ggabsr(n)*exp(yg(n))*xnorm
        enddo 
        return
      end
*************************************************************************    

!!    PROTON LOSSES DUE TO PROTON-PHOTON PAIR PRODUCTION
      subroutine bhloss(x,gp,yg,yp,prloee,sumprloee)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION x(ng), yg(ng), gp(np), yp(np),
     $ prloee(ne), ypel(np)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 

      call pgelos(x,gp,yg,np,ng,deltax,ypel)

! sumprloee: energy lost due to proton-photon pair production
	sumprloee=0.
	do n=2,np-1
	nplus1=n+1
	hcurr=exp(yp(n))*gp(n)*ypel(n)
	hplus1=exp(yp(nplus1))*gp(nplus1)*ypel(nplus1)
	dfpls=(hplus1-hcurr)/gp(n)/exp(yp(n))/(deltap)**2.	
! total proton losses 7.4.2002 (divide by deltap to normalise
! Ray's table)
	prloee(n)=dfpls
        sumprloee=sumprloee+deltap*exp(yp(n))*gp(n)**2.*dfpls*ratmpe
        enddo
        
        return
        end
************************************************************************* 

!!    PHOTON INJECTION FROM ELECTRON SYNCHROTRON
      subroutine synephot(ge,ye,x,yg,syncph,sumsyncph)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION ge(ne),ye(ne),x(ng),yg(ng),syncph(ng),sumsy(ng)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
	
! sumsyncph: energy injected in photons by electron synchrotron losses
	sumsyncph=0.

        do n=1,ng
	sumsy(n)=0.
! for each photon energy we form the ratio x/xcrit
            do m=1,ne
	xsycr=1.5*ge(m)**2.*qb
	ysycr=x(n)/xsycr
            if (ysycr.gt.30.) then
	fsyval=0.
	  else
	fsyval=fsynch(ysycr)/4./pi
            end if
	sumsy(n)=sumsy(n)+deltap*ge(m)*exp(ye(m))*
     $  (4.*sqrt(3.)*tb/qb/x(n))*fsyval/xnorm/exp(yg(n))
            enddo 
        enddo

	do n=1,ng
	syncph(n)=sumsy(n)
	sumsyncph = sumsyncph +
     $    deltax*xnorm*x(n)**2.*exp(yg(n))*syncph(n)	
        enddo
        return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************

!!    PHOTON INJECTION FROM PROTON SYNCHROTRON
      subroutine synpphot(gp,yp,x,yg,syncprph,sumsyncprph)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION gp(np),yp(np),x(ng),yg(ng),syncprph(ng),sumpsy(ng)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
 
!  sumsyncprph: energy injected in photons by proton synchrotron losses
	sumsyncprph=0.

	do n=1,ng
	sumpsy(n)=0.
! for each photon energy we form the ratio x/xcrit
            do m=1,np
	xsycr=1.5*gp(m)**2.*qb/ratmpe
	ysycr=x(n)/xsycr
	if (ysycr.gt.30.) then
	fsyval=0.
	  else
	fsyval=fsynch(ysycr)/4./pi
	end if
	sumpsy(n)=sumpsy(n)+deltap*gp(m)*exp(yp(m))*
     $      (4.*sqrt(3.)*tb/qb/x(n))*fsyval/xnorm/exp(yg(n))
	   enddo
        enddo 

	do n=1,ng
	syncprph(n)=0* sumpsy(n)/ratmpe !FAKE TIMECURVE
	sumsyncprph = sumsyncprph +
     $    deltax*xnorm*x(n)**2.*exp(yg(n))*syncprph(n)*0 !FAKE TIMECURVE
        enddo 
        
        return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************

!!    PHOTON INJECTION FROM ELECTRON INVERSE COMPTON SCATTERING
      subroutine compph_bg(x,ge,yg,ye,gcsth, sumgcsth)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION  ge(ne),ye(ne),x(ng),yg(ng),gcsth(ng)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
      
      call bg(x,ge,yg,ye,ne,ng,deltap,deltax,gcsth) 
! sumgcsth: energy injected in photons from electron inverse Compton scattering
       sumgcsth=0.
       do j=1,ng
       sumgcsth=sumgcsth+deltax*x(j)**2.*gcsth(j)*exp(yg(j))*xnorm
       enddo 

       return
1000  format(1x,6(1pe12.4,1x))          
        end 
*************************************************************************
 
!!    PAIR INJECTION DUE TO PHOTON-PHOTON PAIR PRODUCTION
      subroutine ggelec(x,yg,ge,ye,gginje,sumggabs,sumgginj) 
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 
       DIMENSION x(ng),yg(ng),ge(ne),ye(ne),gginje(ne),yggcral(ng)

       common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
       common/tfl/xnorm 

        call ggcreal(deltax,x,ge,yg,ne,ng,yggcral)
        
! sumgginj: energy injection into pairs due to photon-photon pair production
        sumgginj=0.
        tempgginj=0.
	do n=1,ne
        gginje(n)=4.*tb**2.*qb**(-10./3.)*yggcral(n)/exp(ye(n))
	tempgginj=tempgginj+deltap*ge(n)**2.*gginje(n)*exp(ye(n))
        enddo
                
! re-normalization based on photon energy losses 
        do n=1,ne
        gginje(n)=gginje(n)*sumggabs/tempgginj
        sumgginj=sumgginj+deltap*ge(n)**2.*gginje(n)*exp(ye(n))
        enddo 

        return        
       end 
*************************************************************************

!!    PAIR INJECTION DUE TO PROTON-PHOTON PAIR PRODUCTION
      subroutine bhelec(x,gp,ge,yg,yp,ye,pgeera,sumelBH)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)

      DIMENSION x(ng), yg(ng), ge(ne), ye(ne), gp(np), yp(np),
     $ pgeera(ne), pgray(ne)
      
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
      
      call pgeeray(x,gp,ge,yg,yp,np,ne,ng,deltax,pgray)
      
! sumelBH: energy injection into pairs due to proton-photon pair production
       sumelBH=0.
       do n=1,ne
       pgeera(n)=ratmpe*pgray(n)/exp(ye(n))/deltap
       sumelBH=sumelBH+deltap*ge(n)**2.*exp(ye(n))*pgeera(n)
       enddo 
       
      return
      end 
*************************************************************************    

!!    PROTON-PHOTON PION PRODUCTION
      subroutine pgall(gp,yp,yn,ge,ye,ynt,x,yg
     $,protpg,elecpg,photpg
     $,tneutpg,trinopg,trinopge,trinopgm
     $,tmuonpg,syncmph, syncpiph, synckph)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6,ratmpe=1836.1,
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13)
      
      DIMENSION gp(np), yp(np), ge(ne), ye(ne), x(ng), yg(ng)      
      DIMENSION yn(npmax), ynt(npmax)  
     $,ym(npmax), pms(npmax), gm(npmax)  
     $,ypi(npmax), ppis(npmax), gpi(npmax)  
     $,yk(npmax), pks(npmax), gk(npmax)    
     $,protpg(npmax), tneutpg(npmax), elecpg(npmax), photpg(npmax)
     $,tpionpg(npmax), tmuonpg(npmax), tkaonpg(npmax)
     $,trinopg(npmax), trinopge(npmax), trinopgm(npmax), trinopgt(npmax)
     $,syncmph(npmax), syncpiph(npmax), synckph(npmax)
     $,tne(npmax),tnm(npmax),tna(npmax) 
      DIMENSION Tpg0(3,101),rlum(18,150),Tpgm(300),Tpgpi(300),Tpgk(300)    !m !pk
      DIMENSION Rpg01(3,51,150)
      DIMENSION Rpg02(3,51,150)
      DIMENSION Rpg03(3,51,150)
      DIMENSION Rpg04(3,51,150)
      DIMENSION Rpg05(3,51,150)
      DIMENSION Rpg07(3,51,150)     !pk
      DIMENSION Rpg08(3,51,150)     !pk
      DIMENSION Rpg09(3,51,150)     !pk
      DIMENSION Rpg10(3,51,150)     !pk
      DIMENSION Rpg13(3,51,150)
      DIMENSION Rpg14(3,51,150)
      DIMENSION Rpg15(3,51,150)
      DIMENSION Rpg16(3,51,150)
      DIMENSION Rpg17(3,51,150)
      DIMENSION Rpg18(3,51,150)
      DIMENSION RatioM(300)
      DIMENSION RatioPi(300), RatioKph(300)   !pk
      DIMENSION RatioKe(300), RatioKpi(300), RatioKm(300), RatioKnt(300)   !pk
      DIMENSION Rpm02(130,150)
      DIMENSION Rpm15(130,150)
      DIMENSION Rpm18(130,150)
      DIMENSION Rpp04(130,150)
      DIMENSION Rpp17(130,150)
      DIMENSION Rpk01(130,150)
      DIMENSION Rpk03(130,150)
      DIMENSION Rpk05(130,150)
      DIMENSION Rpk07(130,150)
      DIMENSION Rpk08(130,150)
      DIMENSION Rpk16(130,150)
      DIMENSION Rpk18(130,150)
          
      common/pgr/Tpg0,Rpg01,Rpg02,Rpg03,Rpg04,Rpg05,Rpg13,
     $ Rpg14,Rpg15,Rpg16,Rpg17,Rpg18,Tpgm,Rpm02,Rpm15,Rpm18,pms,RatioM,
     $ Tpgpi,RatioPi,Tpgk,RatioKph,RatioKe,RatioKm,RatioKpi,RatioKnt, 
     $ Rpp04,Rpp17,Rpk01,Rpk03,Rpk05,Rpk07,Rpk08,Rpk16,Rpk18,
     $ Rpg07,Rpg08,Rpg09,Rpg10,ppis,pks
      common/pgm/gm,ym,gpi,ypi,gk,yk
      common/iflg/imsyn,ipsyn,iksyn,ineutron       
      common/param/bplus,tb,qb,gpmax
     $  ,deltap,deltax,np,ne,ntot,ng
      common/tfl/xnorm 
      common/nonlin/factor
      common/pgp/radius,xl10min,bfield  
      
        pcor=1 !pk pion binning correction factor
	pmcor=1 !m muon binning correction factor

        do i=1,ne
	syncmph(i)=0.	!m
	syncpiph(i)=0.	!pk
	synckph(i)=0.	!pk
	end do

	do i=1,np     ! calculation for all resulting protons
	protpg(i)=0
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
         ngtmp=int(k1+j2/2 -ij/2)
         
         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(i+30-(j2-ij)).LE.150.AND.
     $(i+30-(j2-ij)).GE.1) then

       protpg(i) = protpg(i)+ Rpg13(jjj,ngtmp,i+30-(j2-ij))*
     $ exp(yp(j)+yg(k)-yp(i))*factor*ratmpe*(gp(j)*x(k)/gp(i))*
     $(deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !protons 
	End If
           End do
         End do
        enddo

        do j=1,np        ! for all initial protons
         j1=j+30
         j2=j1-90
         do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
         ngtmp=int(k1+j2/2 -ij/2)
         
         If (Tpg0(jjj,ngtmp).GT.0.AND.ngtmp.LE.51.AND.ngtmp.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then
       protpg(j) = protpg(j)- exp(yp(j)+yg(k)-yp(j))*
     $ factor*ratmpe*(gp(j)*x(k)/gp(j))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)  !proton losses

! 	sumprpio = sumprpio+deltap*gp(j)**2.*exp(yp(j))*
!      $  ratmpe*protpg(j)
         End If
         End do
        enddo 

	do i=1,np     ! calculation for all resulting neutrons
	tneutpg(i)=0
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
!!!!!!!! produced neutrons
         ngtmp=int(k1+j2/2 -ij/2)
         
         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(i+30-(j2-ij)).LE.150.AND.
     $(i+30-(j2-ij)).GE.1) then
     
       tneutpg(i) = tneutpg(i)+ Rpg14(jjj,ngtmp,i+30-(j2-ij))*
     $ exp(yp(j)+yg(k)-yn(i))*factor*ratmpe*(gp(j)*x(k)/gp(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !neutrons
         End If
         End do
         End do
! 	sumnpio = sumnpio+deltap*gp(i)**2.*exp(yn(i))*
!      $    ratmpe*tneutpg(i)
        enddo
        
	do i=1,ne     
	elecpg(i)=0.
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
!!!!!!!! electrons
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=i-3-(j2-ij)

	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.
     $ngtmp.GE.1.AND.ngtmp.LE.51.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then

       elecpg(i) = elecpg(i)+ Rpg02(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg03(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !electrons

	 End If
          End do
         End do
! 	sumelpio = sumelpio+deltap*ge(i)**2.*exp(ye(i))*
!      $    elecpg(i)
        enddo 

                             IF (imsyn.EQ.1) THEN
	do i=1,ne     ! calculation for all resulting muons
	tmuonpg(i)=0.
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
!!!!!!!! produced muons							!m
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+21-(j2-ij))
         
         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tmuonpg(i) = tmuonpg(i)+ Rpg04(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gm(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pmcor
     $ + Rpg05(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gm(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pmcor     

	 End If
          End do
         End do
! 	summpio = summpio+deltap*gm(i)**2.*exp(ye(i))*
!      $    207.*tmuonpg(i)
        enddo 
                             END IF  
                             
                             IF (ipsyn.EQ.1) THEN
	do i=1,ne     ! calculation for all resulting pions
	tpionpg(i)=0.
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
!!!!!!!! produced pions							!pk
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+22-(j2-ij))

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tpionpg(i) = tpionpg(i)+ Rpg07(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gpi(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pcor
     $ + Rpg08(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gpi(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pcor
   
	 End If
          End do
         End do
! 	sumppio = sumppio+deltap*gpi(i)**2.*exp(ye(i))*
!      $    274.*tpionpg(i) 
        enddo 
                             END IF
                             
                             IF (iksyn.EQ.1) THEN
	do i=1,ne     ! calculation for all resulting pions
	tkaonpg(i)=0.
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! produced kaons							!pk
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+28-(j2-ij))

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tkaonpg(i) = tkaonpg(i)+ Rpg09(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gk(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg10(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gk(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)

	End If
          End do
         End do
! 	sumkpio = sumkpio+deltap*gk(i)**2*exp(ye(i))*
!      $    967.*tkaonpg(i)
        enddo         
                             END IF

	do i=1,ng    
	photpg(i)=0.
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! photons
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int((i+5.*xl10min-2.33)-(j2-ij)/2)

	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.ngtmp.LE.51.AND.
     $ngtmp.GE.1.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then
    
	If (i.LE.ng-1) then
       photpg(i) = photpg(i)+ 0.7*Rpg01(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $ (deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons
     $	                    + 0.3*Rpg01(jjj,ngtmp,ngtmp2+1)*
     $ exp(yp(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $ (deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons* interpolation
	Elseif (i.EQ.ng) then
       photpg(i) = photpg(i)+ 0.7*Rpg01(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $ (deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons	
	End If

	End If

          End do
         End do
! 	sumgpio = sumgpio+xnorm*deltax*x(i)**2.*exp(yg(i))*
!      $    photpg(i)
        enddo 

 	do i=1,ne     
	trinopg(i)=0
	trinopge(i)=0
	trinopgm(i)=0
	trinopgt(i)=0
         do j=1,np        ! for all possible initial protons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! neutrinos (ye(i) is used because there is no neutrino background)
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=i-3-(j2-ij)
         
	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.
     $ngtmp.GE.1.AND.ngtmp.LE.51.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then
     
       trinopg(i) = trinopg(i)+ Rpg15(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg16(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   
     $ + Rpg17(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   
     $ + Rpg18(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !neutrinos
     
        trinopge(i) = trinopge(i)+ Rpg15(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg16(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)               ! SD-06/2014

       trinopgm(i) = trinopgm(i)+ Rpg17(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg18(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)               ! SD-06/2014 
     
	End If

          End do
         End do         
! 	sumntpio = sumntpio+deltap*ge(i)**2.*exp(ynt(i))*
!      $    trinopg(i)
        enddo 
 
 
	sumsec=sumnpio+sumelpio+sumntpio+summpio+
     $  sumkpio+sumppio+sumgpio
	sumsecpi=summpio+
     $  sumkpio+sumppio


! neutron losses
                             IF (ineutron.EQ.1) THEN

	do i=1,np     ! calculation for all resulting neutrons
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! neutrons
         ngtmp=int(k1+j2/2 -ij/2)

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(i+30-(j2-ij)).LE.150.AND.
     $(i+30-(j2-ij)).GE.1) then

       tneutpg(i) = tneutpg(i)+ Rpg14(jjj,ngtmp,i+30-(j2-ij))*
     $ exp(yn(j)+yg(k)-yn(i))*factor*ratmpe*(gp(j)*x(k)/gp(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !new neutrons

	End If

         End do
         End do

        enddo 

        do j=1,np        ! for all initial neutrons
         j1=j+30
         j2=j1-90
         do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
         ngtmp=int(k1+j2/2 -ij/2)

         If (Tpg0(jjj,ngtmp).GT.0.AND.ngtmp.LE.51.AND.ngtmp.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then

       tneutpg(j) = tneutpg(j)- exp(yn(j)+yg(k)-yn(j))*
     $ factor*ratmpe*(gp(j)*x(k)/gp(j))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)  !neutron losses

         End If
         End do

        enddo


	do i=1,np     ! calculation for all resulting protons
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! produced protons
         ngtmp=int(k1+j2/2 -ij/2)

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(i+30-(j2-ij)).LE.150.AND.
     $(i+30-(j2-ij)).GE.1) then
     
       protpg(i) = protpg(i)+ Rpg13(jjj,ngtmp,i+30-(j2-ij))*
     $exp(yn(j)+yg(k)-yp(i))*factor*1836.1*(gp(j)*x(k)/gp(i))*
     $(deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !protons

	End If
         End do
         End do
        enddo

	do i=1,ne   
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! electrons
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=i-3-(j2-ij)

	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.
     $ngtmp.GE.1.AND.ngtmp.LE.51.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then

       elecpg(i) = elecpg(i)+ Rpg02(jjj,ngtmp,ngtmp2) *
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg03(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !electrons

	End If
          End do
         End do
        enddo 

                             IF (imsyn.EQ.1) THEN
	do i=1,ne     ! calculation for all resulting muons
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If
!!!!!!!! produced muons							!m
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+21-(j2-ij))
         
         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tmuonpg(i) = tmuonpg(i)+ Rpg04(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gm(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pmcor 
     $ + Rpg05(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gm(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pmcor     

	End If
         End do
         End do
        enddo 
                             END IF

                             IF (ipsyn.EQ.1) THEN
	do i=1,ne    
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! produced pions							!pk
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+22-(j2-ij))

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tpionpg(i) = tpionpg(i)+ Rpg07(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gpi(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pcor
     $ + Rpg08(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gpi(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)*pcor
        End If

         End do
         End do
      enddo
                             END IF
                             
                             IF (iksyn.EQ.1) THEN
	do i=1,ne     ! calculation for all resulting kaons
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! produced kaons							!pk
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=int(i+28-(j2-ij))

         If (ngtmp.GE.1.AND.ngtmp.LE.51.AND.
     $Tpg0(jjj,ngtmp).GT.0.AND.
     $j2.GE.-59.AND.j2.LE.60.AND.(ngtmp2).LE.150.AND.
     $(ngtmp2).GE.1) then
     
       tkaonpg(i) = tkaonpg(i)+ Rpg09(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gk(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg10(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ye(i))*factor*ratmpe*(gp(j)*x(k)/gk(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)                  !pk
	End If
         End do
         End do
        enddo 
                             END IF
                             
	do i=1,ng     
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! photons
         ngtmp=int(k1+j2/2 -ij/2)
        ngtmp2=int((i+5*xl10min-2.33)-(j2-ij)/2)

	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.ngtmp.LE.51.AND.
     $ngtmp.GE.1.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then
	If (i.LE.ng-1) then
	
       photpg(i) = photpg(i)+ 0.7*Rpg01(jjj,ngtmp,ngtmp2)*
     $exp(yn(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $(deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons
     $	+ 0.3*Rpg01(jjj,ngtmp,ngtmp2+1)*
     $exp(yn(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $(deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons* interpolation
	Elseif (i.EQ.ng) then
	
       photpg(i) = photpg(i)+ 0.7*Rpg01(jjj,ngtmp,ngtmp2)*
     $exp(yp(j)+yg(k)-yg(i))*factor*ratmpe*(gp(j)*x(k)/x(i))*
     $(deltax*deltap/deltax)/(Tpg0(jjj,ngtmp)*sthom*c*xnorm)   !photons	
	End If
	End If
          End do
         End do

        enddo


	do i=1,ne    
         do j=1,np        ! for all possible initial neutrons
         j1=j+30
         j2=j1-90
          do k = 1,ng     ! for all background photons
          k1=k+5*xl10min+48
               If (j2.LE.30.AND.k1.LE.25) then
               jjj=1
               ij=1
               ElseIf (j2.LE.30.AND.k1.GT.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.LE.25) then
               jjj=2
               ij=31
               ElseIf (j2.GT.30.AND.k1.GT.25) then
               jjj=3
               ij=60
               End If

!!!!!!!! neutrinos
         ngtmp=int(k1+j2/2 -ij/2)
         ngtmp2=i-3-(j2-ij)
	If (ngtmp2.LE.150.AND.Tpg0(jjj,ngtmp).GT.0.AND.
     $ngtmp.GE.1.AND.ngtmp.LE.51.AND.ngtmp2.GE.1.AND.
     $j2.GE.-59.AND.j2.LE.60) then
     
       trinopg(i) = trinopg(i)+ Rpg15(jjj,ngtmp,ngtmp2) *
     $ exp(yn(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg16(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   
     $ + Rpg17(jjj,ngtmp,ngtmp2)*
     $ exp(yn(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   
     $ + Rpg18(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)   !neutrinos
     
      trinopge(i) = trinopge(i)+ Rpg15(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg16(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)               ! SD-06/2014

      trinopgm(i) = trinopgm(i)+ Rpg17(jjj,ngtmp,ngtmp2) *
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)
     $ + Rpg18(jjj,ngtmp,ngtmp2)*
     $ exp(yp(j)+yg(k)-ynt(i))*factor*ratmpe*(gp(j)*x(k)/ge(i))*
     $ (deltax*deltap/deltap)/(Tpg0(jjj,ngtmp)*sthom*c)               ! SD-06/2014

	End If
          End do
         End do

        enddo 

                             END IF !(ineutron.EQ.1)


                             IF (iksyn.EQ.1) THEN

!!!!!!!! Kaons	
	do j=1,ne     

	Csyn=(1/(6.*pi))*sthom*c*966.7**(-3.)*1.22E6*bfield**2.
	gnew = gk(j)/(1.+gk(j)*Csyn*Tpgk(j))
	jk = nint(10*log10(gnew))

         do i=1,ng      
!!!!!!!! photons
         ik=int(jk+21)
	 ij = i+jk-j
         ie=int(i+5.*xl10min-2.33+15)

	If (ie.LE.150.AND.ie.GE.1) then

	      if (j-jk.EQ.1) then 
       photpg(ij) = photpg(ij)+ Rpk01(ik,ie)*tkaonpg(jk)
     $ *(x(ij)/ge(jk))**(-1)*1
     $ *(exp(yg(ij)-ye(jk)))**(-1)*RatioKph(ik)
     $ *deltap/deltax/xnorm
	      end if

       photpg(i) = photpg(i)+ Rpk01(ik,ie)*tkaonpg(j)
     $ *(x(i)/ge(j))**(-1.)
     $ *(exp(yg(i)-ye(j)))**(-1.)*RatioKph(ik)
     $ *pks(j)*deltap/deltax/xnorm

	End If
          End do

         do i=1,ne      
!!!!!!!! electrons
         ik=int(jk+3)
	 ij = i+jk-j
         ie=int(i-3)
	
	If (ie.LE.150.AND.ie.GE.1) then

	      if (j-jk.EQ.1) then
       elecpg(ij) = elecpg(ij)+ Rpk03(ik,ie)*tkaonpg(jk)
     $ *(ge(ij)/ge(jk))**(-1.)
     $ *(exp(ye(i)-ye(jk)))**(-1.)*RatioKe(ik)
	      end if 

       elecpg(i) = elecpg(i)+ Rpk03(ik,ie)*tkaonpg(j)
     $*(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKe(ik)
     $ *pks(j)

	End If
          End do

         do i=1,ne
!!!!!!!! muons
         ik=int(jk-20)
	 ij = i+jk-j
         ie=int(i-3)

	If (ie.LE.150.AND.ie.GE.1) then

        If (ik.GE.1) then

	      if (j-jk.EQ.1) then
       tmuonpg(ij) = tmuonpg(ij)+ Rpk05(ik,ie)*tkaonpg(jk)
     $ *(ge(ij)/ge(jk))**(-1.)
     $ *(exp(ye(ij)-ye(jk)))**(-1.)*RatioKm(ik)
	      end if 

       tmuonpg(i) = tmuonpg(i)+ Rpk05(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKm(ik)
     $ *pks(j)

	else

       tmuonpg(i) = tmuonpg(i)+ Rpk05(ik+20,ie+20)*tkaonpg(j)
     $*(ge(i)/ge(j))**(-1)
     $ *(exp(ye(i)-ye(j)))**(-1)*RatioKm(ik+20)
     $ *pks(j)

	End If
	End If
          End do


         do i=1,ne
!!!!!!!! pions
         ik=int(jk-20)
	 ij = i+jk-j
         ie=int(i-3)

	If (ie.LE.150.AND.ie.GE.1) then

        If (ik.GE.1) then

	      if (j-jk.EQ.1) then
       tpionpg(ij) = tpionpg(ij)+ Rpk07(ik,ie)*tkaonpg(jk)
     $ *(ge(ij)/ge(jk))**(-1.)
     $ *(exp(ye(ij)-ye(jk)))**(-1.)*RatioKpi(ik)
     $ + Rpk08(ik,ie)*tkaonpg(jk)
     $ *(ge(ij)/ge(jk))**(-1)
     $ *(exp(ye(ij)-ye(jk)))**(-1.)*RatioKpi(ik)
	      end if 

       tpionpg(i) = tpionpg(i)+ Rpk07(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKpi(ik)
     $ *pks(j)
     $ + Rpk08(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKpi(ik)
     $ *pks(j)

	else

       tpionpg(i) = tpionpg(i)+ Rpk07(ik+20,ie+20)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKpi(ik+20)
     $ *pks(j)
     $ + Rpk08(ik+20,ie+20)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioKpi(ik+20)
     $ *pks(j)

	End If
	End If
          End do

         do i=1,ne        ! for all resulting neutrinos
!!!!!!!! neutrinos
         ik=int(j+3)
	 ij = i+jk-j
         ie=int(i-3)


	If (ie.LE.150.AND.ie.GE.1) then

	      if (j-jk.EQ.1) then
       trinopg(ij) = trinopg(ij)+ Rpk16(ik,ie)*tkaonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *3.5
     $ + Rpk18(ik,ie)*tkaonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *3.5
     
        trinopge(ij) = trinopge(ij)+ Rpk16(ik,ie)*tkaonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *3.5
     
        trinopgm(ij) = trinopgm(ij)+ Rpk18(ik,ie)*tkaonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *3.5
	      end if        

       trinopg(i) = trinopg(i)+ Rpk16(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *pks(j)*3.5
     $ + Rpk18(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *pks(j)*3.5
     
        trinopge(i) = trinopge(i)+ Rpk16(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *pks(j)*3.5
     
         trinopgm(ij) = trinopgm(ij)+ Rpk18(ik,ie)*tkaonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*RatioKnt(ik)
     $ *pks(j)*3.5
     

	End If
           End do

        enddo ! end of kaon loop
 
                             END IF ! (iksyn.eq.1)
                             
                             IF (ipsyn.EQ.1.OR.iksyn.EQ.1) THEN

	do j=1,ne     

	Csyn=(1./(6.*pi))*sthom*c*273.9**(-3.)*1.22E6*bfield**2.
	gnew = gpi(j)/(1.+gpi(j)*Csyn*Tpgpi(j))
	jpi = nint(10*log10(gnew))

         do i=1,ne     
!!!!!!!! muons
         ipi=int(j-22)
	 ij = i+jpi-j
         ie=int(i-3)
 
	If (ie.LE.150.AND.ie.GE.1) then
        If (ipi.GE.1) then

	      if (j-jpi.EQ.1) then
       tmuonpg(ij) = tmuonpg(ij)+ Rpp04(ipi,ie)*tpionpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ye(ij)-ye(j)))**(-1.)*RatioPi(ipi)
              end if     

       tmuonpg(i) = tmuonpg(i)+ Rpp04(ipi,ie)*tpionpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioPi(ipi)
     $ *ppis(j)

        Else
       tmuonpg(i) = tmuonpg(i)+ Rpp04(ipi+22,ie+22)*tpionpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioPi(ipi+22)
     $ *ppis(j)	        
          
	End If
	End If
          End do

         do i=1,ne        ! for all resulting neutrinos
!!!!!!!! neutrinos
         ipi=int(j+1)
	 ij = i+jpi-j
         ie=int(i-3)

	If (ie.LE.150.AND.ie.GE.1) then

	      if (j-jpi.EQ.1) then
       trinopg(ij) = trinopg(ij)+ Rpp17(ipi,ie)*tpionpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*(1.-RatioPi(ipi))
     $ *3.  

              end if

       trinopg(i) = trinopg(i)+ Rpp17(ipi,ie)*tpionpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*(1.-RatioPi(ipi))
     $ *ppis(j)*4.  !neutrinos

	End If

           End do


        enddo 
                    END IF ! (ipsyn.EQ.1.OR.iksyn.EQ.1)

                        IF (imsyn.EQ.1.OR.ipsyn.EQ.1.OR.iksyn.EQ.1) THEN
!!!!!!! Muon decay							!m
 
	do j=1,ne     ! calculation for all possible initial muons

	Csyn=(1./(6.*pi))*sthom*c*206.8**(-3.)*1.22E6*bfield**2
	gnew = gm(j)/(1.+gm(j)*Csyn*Tpgm(j))
	jm = nint(10*log10(gnew))

         do i=1,ne        ! for all resulting electrons

!!!!!!!! electrons
         im=int(j+1)
	 ij = i+jm-j
         ie=int(i-3)
 
	If (ie.LE.150.AND.ie.GE.1) then

	      if (j-jm.EQ.1) then
       elecpg(ij) = elecpg(ij)+ Rpm02(im,ie)*tmuonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ye(ij)-ye(j)))**(-1.)*RatioM(im)
     $ *3. 
	      end if 

       elecpg(i) = elecpg(i)+ Rpm02(im,ie)*tmuonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ye(i)-ye(j)))**(-1.)*RatioM(im)
     $ *pms(j)*3.

	End If

          End do


         do i=1,ne        ! for all resulting neutrinos
!!!!!!!! neutrinos
         im=int(j+1)
	 ij = i+jm-j
         ie=int(i-3)

	If (im.LE.150.AND.im.GE.1) then

	      if (j-jm.EQ.1) then
       trinopg(ij) = trinopg(ij)+ Rpm15(im,ie)*tmuonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*(1.-RatioM(im))
     $  
     $ + Rpm18(im,ie)*tmuonpg(j)
     $ *(ge(ij)/ge(j))**(-1.)
     $ *(exp(ynt(ij)-ye(j)))**(-1.)*(1.-RatioM(im))

	      end if 

       trinopg(i) = trinopg(i)+ Rpm15(im,ie)*tmuonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*(1.-RatioM(im))
     $ *pms(j)*1.5 
     $ + Rpm18(im,ie)*tmuonpg(j)
     $ *(ge(i)/ge(j))**(-1.)
     $ *(exp(ynt(i)-ye(j)))**(-1.)*(1.-RatioM(im))
     $ *pms(j)*1.5

     
	End If
 	
           End do

 
        enddo 

                             END IF !(imsyn.EQ.1.OR.ipsyn.EQ.1.OR.iksyn.EQ.1)
 
                        IF (imsyn.EQ.1.OR.ipsyn.EQ.1.OR.iksyn.EQ.1) THEN
! What's left of muons for synchrotron				!m
	Do i=1,ne		! calculation for all possible initial muons

	j=int(1+(log(qb*ge(i)**2./206.8)-log(x(1)))/deltax)

		if (j.ge.1.AND.j.LE.ng) then

	syncmph(j) = syncmph(j)+ tmuonpg(i)*(1-pms(i+20))*
     $  (ge(i)/x(j))**2.*deltap/deltax*exp(ye(i)-yg(j))/xnorm*206.8

		end if
	End Do

                             END IF


                             IF (ipsyn.EQ.1.OR.iksyn.EQ.1) THEN
! What's left of pions for synchrotron				!pk
	Do i=1,ne		! calculation for all possible initial pions

	j=int(1+(log(qb*ge(i)**2./274.)-log(x(1)))/deltax)
		if (j.ge.1.AND.j.LE.ng) then

	syncpiph(j) = syncpiph(j)+ tpionpg(i)*(1-ppis(i+21))*
     $  (ge(i)/x(j))**2.*deltap/deltax*exp(ye(i)-yg(j))/xnorm*274.

		end if
	End Do

                             END IF
                             
                             IF (iksyn.EQ.1) THEN
! What's left of kaons for synchrotron				!pk

	Do i=1,ne		! calculation for all possible initial kaons

	j=int(1+(log(qb*ge(i)**2./966.7)-log(x(1)))/deltax)
		if (j.ge.1.AND.j.LE.ng) then

	synckph(j) = synckph(j)+ tkaonpg(i)*(1-pks(i+20))*
     $  (ge(i)/x(j))**2*deltap/deltax*exp(ye(i)-yg(j))/xnorm*966.7

		end if
	End Do

                             END IF 

        return
        end 
************************************************************************* 

      subroutine betheit(gp,ge,x,np,ne,ng,fbh,floss)
      implicit real*8(a-h,o-z)
      parameter(npmax=400,ntotal=800,nwork=440350)
      
!     PHYSICAL CONSTANTS        
      PARAMETER(sthom=6.65e-25,c=3.e10,elrms=1.6e-12*.511e6,
     $ prms=1.6e-12*938.257e6, 
     $ pi=3.14159265,boltz=1.38e-16,astbo=7.56e-15,Bcr=4.4e13) 

        dimension f(100,100)
        dimension irow(40),gp(npmax),time(40),floss(40),sray(40)
	dimension gerfmn(40),gerfmx(40),ge(npmax),x(npmax)
	dimension fbh(npmax,npmax,npmax)
	dimension nemin(40),nemax(40),npin(40)

	common/sgm/sray

        open(unit=19,file='ray.dat',status='unknown')

        do 1 j=1,40
        read (19,*) nevents,x0,dens,gp0,time(j),floss(j)
        read (19,*) nemin(j),nemax(j),npin(j)
********
        csnew=time(j)/(sthom*c)
        sray(j)=csnew 
********* 
        res=1.*(nemax(j)+1-nemin(j))/5.
        irow(j)=int(1.*(nemax(j)+1-nemin(j))/5.)
	gerfmn(j)=x0*10.**(.1*nemin(j))
        gerfmx(j)=x0*10.**(.1*nemax(j)) 
        
        do i=1,irow(j)
	k=5*(i-1)+1
        read (19,*) f(j,k),f(j,k+1),f(j,k+2),f(j,k+3),f(j,k+4)
        enddo
        
1      continue

	
	deltap=log(10.**(.1))

	do 2 j=1,40
	sume=0.
	  do k=1,nemax(j)-nemin(j)+1
	  l=k+nemin(j)-1
	  sume=sume+deltap*10.**(.1*l)*f(j,k)
	  enddo 
	ga=sume/10.**(.1*(npin(j)-1))/1836.1
        gb=floss(j)
	floss(j)=ga
2	continue	
	
	do 3 i=1,np
	if (gp(i).lt.1.e2) goto 3
	  do 21 j=1,ng
	  if (x(j).gt.1.e0) goto 3
	  yu=log10(x(j)*gp(i))
	  if (yu.lt.0.05.and.yu.gt.4.05) then
	  ipair=0
	  else
            do 22 k=1,40
	    yumn=0.1*k-.05
	    yumx=0.1*(k+1)-0.05
	    if (yumn.le.yu.and.yu.lt.yumx) then
	    yuint=yumn+.05
	    kint=int(10.02*yuint)
	    xint=10.**yuint/gp(i)
	    gemn=gerfmn(kint)/xint
	    gemx=gerfmx(kint)/xint
	    gelmn=log10(gemn)
	    gelmx=log10(gemx) 
	    nemn=int(10.01*gelmn)
	    nemx=int(10.01*gelmx)
                do 23 l=1,ne
	        if (l.lt.nemn.or.l.gt.nemx) then
	        fbh(i,l,j)=0.
	        else
	        fbh(i,l,j)=f(kint,l-nemn+1)
                endif
23	        continue
	   endif
22	   continue
	endif
21	continue
3	continue

	end
     
**************************************************************************      

      subroutine densityold(x,ge,yg,ne,ng,deltax,gdens,denkn)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),gdens(ng),denkn(ng)

      do 1 n=1,ne,2
      xtrns=3./4./ge(n)
        
        do m=1,ng
	  if (x(m).gt.xtrns) then
	  mtr=m-1
	  goto 3
	  endif
        enddo

	
3	sumth=.0
	if (mtr.le.1) goto 8

	do m=1,mtr
	gx=x(m)*ge(n)
! 	fkn=1.-gx
	fkn=1.	
	  if (m.eq.1.or.m.eq.mtr) then
	     fac=.5
 	  else
	     fac=1.
	  endif
	sumth=sumth+fac*fkn*deltax*x(m)**2.*exp(yg(m))
        enddo
   
8	sumkn=0.
	do m=mtr+1,ng
	  if (m.eq.mtr+1.or.m.eq.ng) then
	     fac=.5
 	  else
	     fac=1.
	  endif
	ylog=log(4.*x(m)*ge(n)-11./6.)
	sumkn=sumkn+fac*deltax*exp(yg(m))*ylog
        enddo 

        gdens(n)=sumth
        if (gdens(n).eq.0.) gdens(n)=1.e-40
	denkn(n)=sumkn
!   	write (6,1000) n,ge(n),gdens(n)
1       continue

	do n=2,ne-1,2

	if (gdens(n-1).gt.0.) then 
	gdlm=log(gdens(n-1))
	gdlp=log(gdens(n+1))
	gdl=gdlm+.5*(gdlp-gdlm)
        gdens(n)=exp(gdl)
	endif

	dknm=log(denkn(n-1))
	dknp=log(denkn(n+1))
	dkn=dknm+.5*(dknp-dknm)
	denkn(n)=exp(dkn)
        enddo

1000    format (2x,i10,2x,4(1pd12.4,2x))
	end
*************************************************************************

       subroutine bg(x,ge,yg,ye,ne,ng,deltap,deltax,gcsth)
       implicit real*8(a-h,o-z)
       dimension x(ng),ge(ne),yg(ng),ye(ne),gcsth(ng)

       do l=1,ng
        xsc=x(l)

	sumg=0.
	do 15 m=1,ne
	if (ge(m).lt.xsc) then
	goto 15
	endif

	sumx=0.
	do 20 n=1,ng
	xsoft=x(n)	

	if (xsoft.ge.xsc) then
	goto 20
	endif

	if (xsoft.gt.ge(m)) then
	goto 20
	endif

	e0mn=xsc/4./(ge(m)-xsc)/ge(m)
	if (xsoft.le.e0mn) goto 20
	
	Esc=xsc/ge(m)
	Gebg=4.*ge(m)*xsoft

        qs = Esc/(Gebg*(1.-Esc))
	
        sa=2.*qs*log(qs)+(1.+2.*qs)*(1.-qs)
        sb=(gebg*qs)**2.*(1.-qs)/(2.*(1+Gebg*qs))
	Fg=sa+sb
	
	if (Fg.lt.0.) then
	goto 20
	endif

	fux=deltax*exp(yg(n))*Fg
	sumx=sumx+fux
	
	if (sumx.gt.0.) then
	if (fux/sumx.lt.1.e-4) goto 21
	endif
20 	continue

21	fug=deltap*sumx*exp(ye(m))/ge(m)
	if (sumg.gt.0.) then
	if (fug/sumg.lt.1.e-4) goto 16
	endif
14	sumg=sumg+fug

15	continue
	
16	continue
	gcsth(l)=.75*sumg/exp(yg(l))
        enddo 
        
1000    format (1x,5(1pd12.4,1x))
	return
	end 
*************************************************************************

      subroutine bgel(x,ge,yg,ye,ne,ng,deltap,deltax,p)
! we construct the probability function P(ga,gb)--see equation (5.17) in Bloumenthal & Gould (1970)
      parameter (npmax=400,ntotal=800)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),ye(ne),p(npmax,npmax)
      dimension xa(2*ng),yga(2*ng)
      common/tfl/xnorm

! we divide the photons into 2*ng bins
	do 21 m=1,2*ng-1
21	xa(m)=x(1)*exp(deltap*(m-1)) 

	do 22 m=1,2*ng-1,2 
	n=(m+1)/2
22      yga(m)=yg(n)

	do 23 m=2,2*ng-2,2
	n=m/2 
23      yga(m)=.5*(yg(n)+yg(n+1))


	do 10 k=1,ne
	ga=ge(k)

	do 11 l=1,k
	gb=ge(l)

	if (k.eq.l) p(k,l)=0.
 
	xsc=ga-gb
	xmn=4./gb/8. !May 2017 -- original xmn=4./gb
	sumx=0.	 				
	
	do 12 m=1,2*ng-1
	if (xa(m).lt.xmn) then
	kmin=m
	goto 12
	endif

	if (ge(k).le.1.*xa(m)) goto110
	
	if (xa(m).ge.xsc) goto110

	e0mn=xa(m)/4./(ge(k)-xsc)/ge(k)
	e0mn=xsc/4./(ge(k)-xsc)/ge(k)
	if (xa(m).lt.e0mn) goto12
	
	Esc=xsc/ge(k)
	Gebg=4.*ge(k)*xa(m)

        qs = Esc/(Gebg*(1.-Esc))
	
        sa=2.*qs*log(qs)+(1.+2.*qs)*(1.-qs)
        sb=(gebg*qs)**2.*(1.-qs)/(2.*(1+Gebg*qs))
	Fg=sa+sb
	
	if (Fg.lt.0.) goto 110

	fux=deltap*exp(yga(m))*xnorm*Fg
	sumx=sumx+fux
!       stop integration when the relative change in sumx is less than 1e-4      
!!      CAUTION: this choice makes the code faster, but it may give wrong results when deep in KN. 
!!      In that case, comment out the following if statement
	if (sumx.gt.0.) then   
	if (fux/sumx.lt.1.e-4) goto 110
	endif
12 	continue
110	p(k,l)=.75*sumx/ge(k)**2.
11	continue

10	continue
	
1000    format (1x,5(1pd12.4,1x))
	return
	end 
*************************************************************************

        subroutine ggabsal(deltax,x,yg,ng,ygi)
        implicit real*8(a-h,o-z)
        dimension x(ng),yg(ng),ygi(ng)
!       uses the approximate expression from Coppi & Blandford (1990)
!       see also equation (55) in Mastichiadis & Kirk (1990)
	xmin=x(1)
	xmax=x(ng)

	do n=1,ng
	xinv=1./x(n)
	sum=0.
            do m=1,ng
	if (x(m).gt.xinv) then
	w=x(m)*x(n)
	r=0.652*(w**2.-1.)*log(w)/w**3.
	sum=sum+deltax*r*exp(yg(m))*x(m)
	end if
            enddo
	ygi(n)=sum
        enddo 
        return
	end
*************************************************************************

      subroutine ggcreal(deltax,x,ge,yg,ne,ng,yggcre)
      implicit real*8(a-h,o-z)
      dimension x(ng),ge(ne),yg(ng),yga(500),yggcre(ng)

	do 1 n=1,ne
	sum=0.
	if (ge(n)**2.lt.2.) then
        yggcre(n)=0.
	goto 1
	endif
	
	xa=ge(n)+sqrt(ge(n)**2-2.)
	xa=2.*ge(n)
	xb=1./xa

	if (xa.gt.x(ng).or.xa.lt.x(1)) then
	  yga(n)=-99.
	  goto 12
	end if

	if (xb.gt.x(ng).or.xb.lt.x(1)) then
	  yga(n)=-99.
	  goto 12
	end if

	do ka=1,ng
	 if (x(ka).gt.xa) then
	   ma=ka-1
	   goto 30
	 endif
        enddo

30        yga(n)=yg(ma+1)+(yg(ma)-yg(ma+1))*log(xa/x(ma+1))
     $    /log(x(ma)/x(ma+1))

	do m=1,ng
	if (x(m).gt.xb) then
	w=x(m)*xa
	r=0.652*(w**2.-1.)*log(w)/w**3.
	sum=sum+deltax*r*exp(yg(m))*x(m)
	end if
        enddo

12      yggcre(n)=exp(yga(n))*sum
!  	write (6,1000) ge(n),yggcre(n)

1       continue
1000    format (2x,4(1pd12.4,2x))
	end
*************************************************************************
        
        subroutine pgeeray(x,gp,ge,yg,yp,np,ne,ng,deltax,pgray) 
        parameter (npmax=400,ntotal=800)
	implicit real*8(a-h,o-z)
	dimension x(ng),gp(np),yg(ng),yp(np),pgray(ne),ge(ne),sray(40)
        dimension fpbh(npmax,npmax,npmax)
        common/bet/fpbh
        common/nonlin/factor
	common/sgm/sray

	do 10 i=1,ne
        sum2=0.
            do 11 j=1,np
	if (gp(j).lt.1.e3) goto 11
        sum1=0.
                do 12 k=1,ng 
       if (x(k).gt.1.) goto 14
       if (fpbh(j,i,k).eq.0.) goto 12
! here we perform a simple interpolation to get the value
! of Ray's angle average cross section
        y0=x(k)*gp(j)
        yl=log10(y0)
        do 33 m=1,40
        if (.1*m.lt.yl.and..1*(m+1).ge.yl) then
        sla=log10(sray(m))
        slb=log10(sray(m+1))
        xma=.1*m
        xmb=.1*(m+1)
        ysl=(slb-sla)/(xmb-xma)*(yl-xma)+sla
        csect=10.**ysl
        goto 34
        else
        csect=.0
        endif
33      continue
34      continue

	fa=deltax*x(k)*fpbh(j,i,k)*exp(yg(k))*csect
     $    *factor/ge(i)
     
        sum1=sum1+fa
12      continue
14      continue
        
        fb=deltax*gp(j)*exp(yp(j))*sum1/2.
        sum2=sum2+fb

11      continue
        pgray(i)=sum2
10      continue
1000    format (2x,5(1pd12.4,2x))
        end
*********************************************************************

	subroutine pgelos(x,gp,yg,np,ng,deltax,ypel)
        parameter (npmax=400,ntotal=800)
	implicit real*8(a-h,o-z)
	dimension x(ng),gp(np),yg(ng),yp(np),ypel(np),sray(40)
        dimension fpgloss(40)
	common/prl/fpgloss
	common/nonlin/factor
	common/sgm/sray

	do 30 i=1,np

	sumg=0.
	if (gp(i).lt.1.e2) goto 35
            do 31 j=1,ng
	if (x(j).gt.1.) goto 35 
	yu=log10(x(j)*gp(i))
	if (yu.lt..05.and.yu.gt.4.05) then
	goto 31
	else
                do 32 k=1,40
	yumn=.1*k-.05
	yumx=.1*(k+1)-.05
	if (yumn.le.yu.and.yu.lt.yumx) then
	yuint=yumn+.05
	xint=10.**yuint/gp(i)
	kint=int(10.02*yuint)
! here we perform a simple interpolation to get the value
! of Ray's angle average cross section
        y0=x(j)*gp(i)
        yl=log10(y0)
        do 33 m=1,40
        if (.1*m.lt.yl.and..1*(m+1).ge.yl) then
        sla=log10(sray(m))
        slb=log10(sray(m+1))
        xma=.1*m
        xmb=.1*(m+1)
        ysl=(slb-sla)/(xmb-xma)*(yl-xma)+sla
        csect=10.**ysl
        goto 34
        else
        csect=.0
        endif
33      continue
34      continue


	fg=deltax*x(j)*exp(yg(j))*csect*
     $    fpgloss(kint)*1836.1*factor
	sumg=sumg+fg
	endif
32	continue
	endif
31	continue
35	ypel(i)=sumg
30	continue
1000    format (2x,5(1pd12.4,2x)) 
	end	

************************************************* 

      subroutine PGPI_READER(Tpg0,Rpg01,Rpg02,Rpg03,Rpg04,Rpg05,Rpg13,
     $Rpg14,Rpg15,Rpg16,Rpg17,Rpg18,Rpm02,Rpm15,Rpm18,Tpgm,RatioM,
     $Tpgpi,RatioPi,Tpgk,RatioKph,RatioKe,RatioKm,RatioKpi,RatioKnt,
     $Rpp04,Rpp17,Rpk01,Rpk03,Rpk05,Rpk07,Rpk08,Rpk16,Rpk18,
     $Rpg07,Rpg08,Rpg09,Rpg10)
!     Pion-photo-production on arbitrary soft photon spectrum.
!     target is approximated by 161 isotropic monoenergetic photon
!     fields having energy 1.d-19*10.D0**((IE-1)/10.D0) GeV, and number
!     density TARGET_FACTORS(IE) cm^-3.  This subroutine combines data
!     for photopion-production on these monoenergetic fields
!     and outputs data for the composite field.  Distribution is for each
!     non-decayed particle, and is given the corresponding SOPHIA particle codes 1-18:
!    1 -- photon
!    2 -- positron
!    3 -- electron
!    4 -- mu+
!    5 -- mu-
!    6 -- pi0
!    7 -- pi+
!    8 -- pi-
!    9 -- K+
!   10 -- K-
!   11 -- K0l
!   12 -- K0s
!   13 -- proton
!   14 -- neutron
!   15 -- nu_e
!   16 -- bar nu_e
!   17 -- nu_mu
!   18 -- bar nu_mu
!   all energies in GeV
!     A.G.F. Reimer 29/11/04

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (M=300,M2=90000,N=120)

!! DIST() is the energy distribution in GeV of all the produced particles
!! EMID() is the proton energy for the population M
!! TARGET_FACTORS() is the photon density for each energy level

      DIMENSION DIST(M+3),TOT_DIST(18,M),EMID(M)
      DIMENSION TARGET_FACTORS(161), Tpg0(3,101), Tpgm(M), RatioM(M)
      DIMENSION Tpgpi(M), RatioPi(M), Tpgk(M), RatioKph(M)   !pk
      DIMENSION RatioKe(M), RatioKpi(M), RatioKm(M), RatioKnt(M)   !pk
      DIMENSION Rpg01(3,51,150)
      DIMENSION Rpg02(3,51,150)
      DIMENSION Rpg03(3,51,150)
      DIMENSION Rpg04(3,51,150)
      DIMENSION Rpg05(3,51,150)
      DIMENSION Rpg07(3,51,150)
      DIMENSION Rpg08(3,51,150)
      DIMENSION Rpg09(3,51,150)
      DIMENSION Rpg10(3,51,150)
      DIMENSION Rpg13(3,51,150)
      DIMENSION Rpg14(3,51,150)
      DIMENSION Rpg15(3,51,150)
      DIMENSION Rpg16(3,51,150)
      DIMENSION Rpg17(3,51,150)
      DIMENSION Rpg18(3,51,150)
      DIMENSION Rpm02(130,150)
      DIMENSION Rpm15(130,150)
      DIMENSION Rpm18(130,150)
      DIMENSION Rpp04(130,150)
      DIMENSION Rpp17(130,150)
      DIMENSION Rpk01(130,150)
      DIMENSION Rpk03(130,150)
      DIMENSION Rpk05(130,150)
      DIMENSION Rpk07(130,150)
      DIMENSION Rpk08(130,150)
      DIMENSION Rpk16(130,150)
      DIMENSION Rpk18(130,150)
      common/iflg/imsyn,ipsyn,iksyn,ineutron
		
!       write(*,*) 'initializing pgpi_reader'

!      DO I=1,300
!         EMID(I)=1.D-3*10.D0**((I-0.5D0)/20.D0) ! GeV
!      ENDDO
        DIST(301)=0
        E_DEN_SOFT=0.0D0       ! GeV/cm^3
        RATE=0.D0          ! s^-1
	RatioM(n)=0
	RatioPi(n)=0
	RatioKph(n)=0
	RatioKe(n)=0
	RatioKm(n)=0
	RatioKpi(n)=0
	RatioKnt(n)=0

!         EMID(IE0)=1.D-3*10.D0**((IE0-0.5D0)/20.D0) ! GeV

         IE_SOFT=0
         E_SOFT=0.0D0
                             IF (iksyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_EP002p.dat',status='old')
                             ELSEIF (ipsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP002p.dat',status='old')
                             ELSEIF (imsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP002m.dat',status='old')
                             ELSE
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP002a.dat',status='old')
                             END IF
!		write(*,*) 'EP_002 opened'
 1       READ(10,*,END=4)ITSOPHIA_0,N_INJ,IE_0,IE_SOFT_101,
     &        E_0,E_SOFT,DEN_SOFT,T_INT

         IE_SOFT=IE_SOFT_101+60   ! DATA FILE HAD ONLY 101 TARGET ENERGIES
! Interaction time stored in memory
!		write(*,*) 'EP_002 first line read'
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
         IE_SOFT_D = (IE_SOFT_101+1)/2
         Tpg0(1,IE_SOFT_D) = T_INT
        END IF

         DO ITSOPH=1,18
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
		IS=MIN
		DO WHILE(IS.LE.MAX)
               READ(10,*) DIST(IS)
		IS=IS+1
		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
!               IF (ITSOPHIA.EQ.1) THEN
!               Rpg01(1,IE_SOFT_101,IA) = DIST(IS) + DIST(IS+1)
!               ENDIF
               If (ITSOPHIA.EQ.2) THEN
               Rpg02(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.3) THEN
               Rpg03(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.4) THEN
               Rpg04(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.5) THEN
               Rpg05(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.7) THEN
               Rpg07(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.8) THEN
               Rpg08(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.9) THEN
               Rpg09(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.10) THEN
               Rpg10(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.13) THEN
               Rpg13(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.14) THEN
               Rpg14(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.15) THEN
               Rpg15(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.16) THEN
               Rpg16(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.17) THEN
               Rpg17(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.18) THEN
               Rpg18(1,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
        END IF
               END DO
!		write(*,*) 'EP_002 done with protons'
                   IF (mod(MIN,4).EQ.0) THEN
                   MIN1=MIN
                   ELSEIF(mod(MIN-1,4).EQ.0) THEN
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   ELSEIF(mod(MIN-2,4).EQ.0) THEN
                   MIN1=MIN-2
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
                   ELSE
                   MIN1=MIN-3
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
		   DIST(MIN1+2)=0
                   END IF
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
               DO IS=MIN1,MAX,4
                  IB = IS/4
               IF (ITSOPHIA.EQ.1) THEN
               Rpg01(1,IE_SOFT_D,IB) = DIST(IS) + DIST(IS+1) +
     $         DIST(IS+2) + DIST(IS+3)
               ENDIF
               END DO
            ENDIF
        END IF
         ENDDO

         GOTO 1
 4       CONTINUE
          IE_SOFT=0
         E_SOFT=0.0D0

!		write(*,*) 'all fine for EP_002'
                             IF (iksyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_EP062p.dat',status='old')
                             ELSEIF (ipsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP062p.dat',status='old')
                             ELSEIF (imsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP062m.dat',status='old')
                             ELSE
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP062a.dat',status='old')
                             END IF
 2       READ(10,*,END=5)ITSOPHIA_0,N_INJ,IE_0,IE_SOFT_101,
     &        E_0,E_SOFT,DEN_SOFT,T_INT
         IE_SOFT=IE_SOFT_101+60   ! DATA FILE HAD ONLY 101 TARGET ENERGIES
! Interaction time stored in memory
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
         IE_SOFT_D = (IE_SOFT_101+1)/2
         Tpg0(2,IE_SOFT_D) = T_INT
!	write(*,*) Tpg0(2,51)
        END IF
!	write(*,*) Tpg0(2,51)
!	pause
         DO ITSOPH=1,18
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
		IS=MIN
		DO WHILE(IS.LE.MAX)
               READ(10,*) DIST(IS)
		IS=IS+1
		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
!               If (ITSOPHIA.EQ.1) THEN
!               Rpg01(2,IE_SOFT_101,IA) = DIST(IS) + DIST(IS+1)
!               ENDIF
               If (ITSOPHIA.EQ.2) THEN
               Rpg02(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.3) THEN
               Rpg03(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.4) THEN
               Rpg04(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.5) THEN
               Rpg05(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.7) THEN
               Rpg07(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.8) THEN
               Rpg08(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.9) THEN
               Rpg09(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.10) THEN
               Rpg10(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.13) THEN
               Rpg13(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.14) THEN
               Rpg14(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.15) THEN
               Rpg15(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.16) THEN
               Rpg16(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.17) THEN
               Rpg17(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.18) THEN
               Rpg18(2,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
        END IF
               END DO
                   IF (mod(MIN,4).EQ.0) THEN
                   MIN1=MIN
                   ELSEIF(mod(MIN-1,4).EQ.0) THEN
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   ELSEIF(mod(MIN-2,4).EQ.0) THEN
                   MIN1=MIN-2
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
                   ELSE
                   MIN1=MIN-3
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
		   DIST(MIN1+2)=0
                   END IF
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
               DO IS=MIN1,MAX,4
                  IB = IS/4
               IF (ITSOPHIA.EQ.1) THEN
               Rpg01(2,IE_SOFT_D,IB) = DIST(IS) + DIST(IS+1) +
     $         DIST(IS+2) + DIST(IS+3)
               ENDIF
               END DO
            ENDIF
        END IF
         ENDDO

         GOTO 2
 5       CONTINUE
          IE_SOFT=0
         E_SOFT=0.0D0
!		write(*,*) 'all fine for EP_062'

                             IF (iksyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_EP120p.dat',status='old')
                             ELSEIF (ipsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP120p.dat',status='old')
                             ELSEIF (imsyn.EQ.1) THEN
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP120m.dat',status='old')
                             ELSE
         OPEN (UNIT=10, FILE=
     &   './pgpi_decd_EP120a.dat',status='old')
                             END IF
 3       READ(10,*,END=6)ITSOPHIA_0,N_INJ,IE_0,IE_SOFT_101,
     &        E_0,E_SOFT,DEN_SOFT,T_INT
         IE_SOFT=IE_SOFT_101+60   ! DATA FILE HAD ONLY 101 TARGET ENERGIES
! Interaction time stored in memory
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
         IE_SOFT_D = (IE_SOFT_101+1)/2
         Tpg0(3,IE_SOFT_D) = T_INT
        END IF

         DO ITSOPH=1,18
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
		IS=MIN
		DO WHILE(IS.LE.MAX)
               READ(10,*) DIST(IS)
		IS=IS+1
		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
!               If (ITSOPHIA.EQ.1) THEN
!               Rpg01(3,IE_SOFT_101,IA) = DIST(IS) + DIST(IS+1)
!               ENDIF
               If (ITSOPHIA.EQ.2) THEN
               Rpg02(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.3) THEN
               Rpg03(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.4) THEN
               Rpg04(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.5) THEN
               Rpg05(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.7) THEN
               Rpg07(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.8) THEN
               Rpg08(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.9) THEN
               Rpg09(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.10) THEN
               Rpg10(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.13) THEN
               Rpg13(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.14) THEN
               Rpg14(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.15) THEN
               Rpg15(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.16) THEN
               Rpg16(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.17) THEN
               Rpg17(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPHIA.EQ.18) THEN
               Rpg18(3,IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
               ENDIF
        END IF
               END DO
                   IF (mod(MIN,4).EQ.0) THEN
                   MIN1=MIN
                   ELSEIF(mod(MIN-1,4).EQ.0) THEN
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   ELSEIF(mod(MIN-2,4).EQ.0) THEN
                   MIN1=MIN-2
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
                   ELSE
                   MIN1=MIN-3
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
		   DIST(MIN1+2)=0
                   END IF
        IF (mod(IE_SOFT_101+1,2).EQ.0) THEN
               DO IS=MIN1,MAX,4
                  IB = IS/4
               IF (ITSOPHIA.EQ.1) THEN
               Rpg01(3,IE_SOFT_D,IB) = DIST(IS) + DIST(IS+1) +
     $         DIST(IS+2) + DIST(IS+3)
               ENDIF
               END DO
            ENDIF
        END IF
         ENDDO

         GOTO 3
 6       CONTINUE

                        IF (iksyn.EQ.1.OR.ipsyn.EQ.1.OR.imsyn.EQ.1) THEN
!        Muon decay
          OPEN (UNIT=10, FILE=
     &   './dec_mup.dat',status='old')


 7       READ(10,*,END=8)ITSOPHIA_0,xx,IE_m,xx,
     &        xx,xx,xx,T_INT

         IE_SOFT=IE_m   ! muon file
! Interaction time stored in memory

        IF (mod(IE_m,2).EQ.0) THEN
         IE_SOFT_D = (IE_m)/2
         Tpgm(IE_SOFT_D) = T_INT
        END IF
	p02=0
	p15=0
	p18=0
         DO ITSOPH=1,3
! 1=positron, 2=nu_e, 3= bar nu_mu
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
                READ(10,*)(DIST(IS),IS=MIN,MAX)
!		IS=MIN
!		DO WHILE(IS.LE.MAX)
!               READ(10,*) DIST(IS)
!		IS=IS+1
!		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_m+1,2).EQ.0) THEN

               If (ITSOPH.EQ.1) THEN
               Rpm02(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p02 = p02 + Rpm02(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.2) THEN
               Rpm15(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p15 = p15 + Rpm15(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.3) THEN
               Rpm18(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p18 = p18 + Rpm18(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
        END IF
               END DO
!		write(*,*) 'EP_002 done with protons'

            END IF
         ENDDO
        IF (mod(IE_m+1,2).EQ.0) THEN
		RatioM(IE_SOFT_D)=p02/(p02+p15+p18)
!		write(*,*)IE_SOFT_D,p02,p15,p18,RatioM(IE_SOFT_D)
!		pause
	END IF
         GOTO 7
8        Continue
                             END IF

                             IF (ipsyn.EQ.1.OR.iksyn.EQ.1) THEN
C        Pion decay
          OPEN (UNIT=10, FILE=
     &   './dec_pip.dat',status='old')


17       READ(10,*,END=18)ITSOPHIA_0,xx,IE_pi,xx,
     &        xx,xx,xx,T_INT

         IE_SOFT=IE_pi   ! pion file
! Interaction time stored in memory
        IF (mod(IE_pi,2).EQ.0) THEN
         IE_SOFT_D = (IE_pi-1)/2
         Tpgpi(IE_SOFT_D) = T_INT
        END IF
	p04=0
	p17=0
         DO ITSOPH=1,2
! 1=mu+, 2=nu_m
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
                READ(10,*)(DIST(IS),IS=MIN,MAX)
!		IS=MIN
!		DO WHILE(IS.LE.MAX)
!               READ(10,*) DIST(IS)
!		IS=IS+1
!		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_pi+1,2).EQ.0) THEN
               If (ITSOPH.EQ.1) THEN
               Rpp04(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p04 = p04 + Rpp04(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ttpt1= ttpt1+ DIST(IS) + DIST(IS+1)
               ENDIF
               If (ITSOPH.EQ.2) THEN
               Rpp17(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p17 = p17 + Rpp17(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ttpt2= ttpt2+ DIST(IS) + DIST(IS+1)
               ENDIF
        END IF
               END DO


            END IF
         ENDDO
        IF (mod(IE_pi+1,2).EQ.0) THEN
                If (p04+p17.GT.0) then
		RatioPi(IE_SOFT_D)=p04/(p04+p17)
		Else
		RatioPi(IE_SOFT_D)= 0
		End if
	END IF
        ptot =  10**((IE_SOFT_D-10.25)/10)
!		write(*,*)IE_SOFT_D,p04,p17,p04+p17,ptot
!		pause
         GOTO 17
18       Continue
                             END IF

                             IF (iksyn.EQ.1) THEN
C        Pion decay
          OPEN (UNIT=10, FILE=
     &   './dec_kam.dat',status='old')


27       READ(10,*,END=28)ITSOPHIA_0,xx,IE_k,xx,
     &        xx,xx,xx,T_INT

         IE_SOFT=IE_k   ! kaon file
! Interaction time stored in memory

        IF (mod(IE_k,2).EQ.0) THEN
         IE_SOFT_D = (IE_k-7)/2
         Tpgk(IE_SOFT_D) = T_INT
        END IF
	p01=0
	p03=0
	p05=0
	p07=0
	p08=0
	p16=0
	p18=0
         DO ITSOPH=1,18
! 1=photon, ..., 18=bar nu_mu
            READ(10,*)ITSOPHIA,MIN,MAX
            IF((MAX.GE.MIN).AND.(MIN.GT.0))THEN
                READ(10,*)(DIST(IS),IS=MIN,MAX)
!		IS=MIN
!		DO WHILE(IS.LE.MAX)
!               READ(10,*) DIST(IS)
!		IS=IS+1
!		END DO
			DIST(MAX+1)=0
			DIST(MAX+2)=0
			DIST(MAX+3)=0
                   IF (mod(MIN,2).EQ.0) THEN
                   MIN1=MIN
                   ELSE
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   END IF
               DO IS=MIN1,MAX,2
                  IA = IS/2
        IF (mod(IE_k+1,2).EQ.0) THEN

               If (ITSOPH.EQ.3) THEN
               Rpk03(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p03 = p03 + Rpk03(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.5) THEN
               Rpk05(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p05 = p05 + Rpk05(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.7) THEN
               Rpk07(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p07 = p07 + Rpk07(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.8) THEN
               Rpk08(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p08 = p08 + Rpk08(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.16) THEN
               Rpk16(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p16 = p16 + Rpk16(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF
               If (ITSOPH.EQ.18) THEN
               Rpk18(IE_SOFT_D,IA) = DIST(IS) + DIST(IS+1)
		p18 = p18 + Rpk18(IE_SOFT_D,IA)*
     $           10**((IA-30.5)/10)
               ENDIF

        END IF
               END DO

                   IF (mod(MIN,4).EQ.0) THEN
                   MIN1=MIN
                   ELSEIF(mod(MIN-1,4).EQ.0) THEN
                   MIN1=MIN-1
		   DIST(MIN1)=0
                   ELSEIF(mod(MIN-2,4).EQ.0) THEN
                   MIN1=MIN-2
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
                   ELSE
                   MIN1=MIN-3
		   DIST(MIN1)=0
		   DIST(MIN1+1)=0
		   DIST(MIN1+2)=0
                   END IF
        IF (mod(IE_k+1,2).EQ.0) THEN
               DO IS=MIN1,MAX,4
                  IB = IS/4
               IF (ITSOPHIA.EQ.1) THEN
               Rpk01(IE_SOFT_D,IB) = DIST(IS) + DIST(IS+1) +
     $         DIST(IS+2) + DIST(IS+3)
		p01 = p01 + Rpk01(IE_SOFT_D,IB)*
     $           10**((IB-15.25)/5)
               ENDIF
               END DO
            ENDIF
!        END IF

            END IF
         ENDDO
        IF (mod(IE_k+1,2).EQ.0) THEN
                If (p01+p03+p05+p07+p08+p16+p18.GT.0) then
		RatioKph(IE_SOFT_D)=p01/(p01+p03+p05+p07+p08+p16+p18)
		RatioKe(IE_SOFT_D)=p03/(p01+p03+p05+p07+p08+p16+p18)
		RatioKm(IE_SOFT_D)=p05/(p01+p03+p05+p07+p08+p16+p18)
	RatioKpi(IE_SOFT_D)=(p07+p08)/(p01+p03+p05+p07+p08+p16+p18)
	RatioKnt(IE_SOFT_D)=(p16+p18)/(p01+p03+p05+p07+p08+p16+p18)
		Else
		RatioKph(IE_SOFT_D)=0
		RatioKe(IE_SOFT_D)=0
		RatioKm(IE_SOFT_D)=0
	RatioKpi(IE_SOFT_D)=0
	RatioKnt(IE_SOFT_D)=0
		End if
!		write(*,*)IE_SOFT_D,p02,p15,p18,RatioM(IE_SOFT_D)
!		pause
	END IF
         GOTO 27
28       Continue
                             END IF

10      FORMAT(I4,I4)
         IE_SOFT=0
         E_SOFT=0.0D0

!	do IE_SOFT_D=1,130
!	write(*,*) IE_SOFT_D, Tpgk(IE_SOFT_D)
!     $  ,RatioKe(IE_SOFT_D), RatioKm(IE_SOFT_D), RatioKpi(IE_SOFT_D)
!	end do
!	pause
!	write(*,*) Rpm02(129,106),Rpm02(130,106)
!	pause
      CLOSE (UNIT=10)
      END
      
************************************************************************
      FUNCTION FSYNCH(T)
! 
!  THE FUNCTION F(X)=X * INTEGRAL(X TO INFINITY) K_5/3(Z) DZ
!  WHICH IS USED IN THE CALCULATION OF SYNCHROTRON EMISSIVITY.
!  (SEE  'RADIO ASTROPHYSICS' BY A.G.PACHOLCZYK, 1970, (FREEMAN,
!  SAN FRANCISCO).
!                                 R.J.PROTHEROE, ADELAIDE 20 OCT 1988.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X1(35),Y1(35),C1(34,3)
      DIMENSION X2(126),Y2(126),C2(125,3),C21(125),C22(125),C23(125)
      EQUIVALENCE (C2(1,1),C21(1)),(C2(1,2),C22(1)),(C2(1,3),C23(1))
      F1(X,A1)=A1*(X/2.D0)**(1.D0/3.D0)
      F2(X,A2)=A2*EXP(0.5D0*DLOG(X)-X)
!
      DATA A1/2.7082359D0/,A2/1.25331414/
      DATA X1/
     *  1.00D-4, 2.00D-4, 5.00D-4, 1.00D-3, 2.00D-3, 5.00D-3, 1.00D-2,
     *  2.00D-2, 3.00D-2, 4.00D-2, 5.00D-2, 6.00D-2, 7.00D-2, 8.00D-2,
     *  9.00D-2, 1.00D-1, 1.10D-1, 1.20D-1, 1.30D-1, 1.40D-1, 1.50D-1,
     *  1.60D-1, 1.70D-1, 1.80D-1, 1.90D-1, 2.00D-1, 2.10D-1, 2.20D-1,
     *  2.30D-1, 2.40D-1, 2.50D-1, 2.60D-1, 2.70D-1, 2.80D-1, 2.90D-1/
      DATA X2/
     *  2.90D-1, 3.00D-1, 3.10D-1, 3.20D-1, 3.30D-1, 3.40D-1, 3.50D-1,
     *  3.60D-1, 3.70D-1, 3.80D-1, 3.90D-1, 4.00D-1, 4.10D-1, 4.20D-1,
     *  4.30D-1, 4.40D-1, 4.50D-1, 4.60D-1, 4.70D-1, 4.80D-1, 4.90D-1,
     *  5.00D-1, 5.20D-1, 5.40D-1, 5.60D-1, 5.80D-1, 6.00D-1, 6.20D-1,
     *  6.40D-1, 6.60D-1, 6.80D-1, 7.00D-1, 7.20D-1, 7.40D-1, 7.60D-1,
     *  7.80D-1, 8.00D-1, 8.20D-1, 8.40D-1, 8.60D-1, 8.80D-1, 9.00D-1,
     *  9.20D-1, 9.40D-1, 9.60D-1, 9.80D-1, 1.00D+0, 1.05D+0, 1.10D+0,
     *  1.15D+0, 1.20D+0, 1.25D+0, 1.30D+0, 1.35D+0, 1.40D+0, 1.45D+0,
     *  1.50D+0, 1.55D+0, 1.60D+0, 1.65D+0, 1.70D+0, 1.75D+0, 1.80D+0,
     *  1.85D+0, 1.90D+0, 1.95D+0, 2.00D+0, 2.10D+0, 2.20D+0, 2.30D+0,
     *  2.40D+0, 2.50D+0, 2.60D+0, 2.70D+0, 2.80D+0, 2.90D+0, 3.00D+0,
     *  3.10D+0, 3.20D+0, 3.30D+0, 3.40D+0, 3.50D+0, 3.60D+0, 3.70D+0,
     *  3.80D+0, 3.90D+0, 4.00D+0, 4.10D+0, 4.20D+0, 4.30D+0, 4.40D+0,
     *  4.50D+0, 4.60D+0, 4.70D+0, 4.80D+0, 4.90D+0, 5.00D+0, 5.25D+0,
     *  5.50D+0, 5.75D+0, 6.00D+0, 6.25D+0, 6.50D+0, 6.75D+0, 7.00D+0,
     *  7.25D+0, 7.50D+0, 7.75D+0, 8.00D+0, 8.25D+0, 8.50D+0, 8.75D+0,
     *  9.00D+0, 9.25D+0, 9.50D+0, 9.75D+0, 1.00D+1, 1.20D+1, 1.40D+1,
     *  1.60D+1, 1.80D+1, 2.00D+1, 2.50D+1, 3.00D+1, 4.00D+1, 5.00D+1/
      DATA Y1/
     *  1.00D+0, 9.94D-1, 9.96D-1, 9.91D-1, 9.97D-1, 9.74D-1, 9.61D-1,
     *  9.37D-1, 9.18D-1, 9.02D-1, 8.86D-1, 8.71D-1, 8.58D-1, 8.44D-1,
     *  8.32D-1, 8.20D-1, 8.09D-1, 7.97D-1, 7.86D-1, 7.76D-1, 7.65D-1,
     *  7.56D-1, 7.47D-1, 7.37D-1, 7.28D-1, 7.19D-1, 7.11D-1, 7.02D-1,
     *  6.93D-1, 6.85D-1, 6.77D-1, 6.68D-1, 6.60D-1, 6.53D-1, 6.45D-1/
      DATA Y2/
     *  1.82D+0, 1.80D+0, 1.79D+0, 1.78D+0, 1.77D+0, 1.76D+0, 1.75D+0,
     *  1.73D+0, 1.72D+0, 1.71D+0, 1.70D+0, 1.70D+0, 1.69D+0, 1.68D+0,
     *  1.67D+0, 1.66D+0, 1.66D+0, 1.65D+0, 1.64D+0, 1.63D+0, 1.63D+0,
     *  1.62D+0, 1.61D+0, 1.59D+0, 1.58D+0, 1.57D+0, 1.56D+0, 1.55D+0,
     *  1.54D+0, 1.53D+0, 1.52D+0, 1.51D+0, 1.50D+0, 1.50D+0, 1.49D+0,
     *  1.48D+0, 1.47D+0, 1.46D+0, 1.46D+0, 1.45D+0, 1.44D+0, 1.44D+0,
     *  1.43D+0, 1.43D+0, 1.42D+0, 1.42D+0, 1.42D+0, 1.40D+0, 1.39D+0,
     *  1.38D+0, 1.37D+0, 1.36D+0, 1.35D+0, 1.34D+0, 1.33D+0, 1.32D+0,
     *  1.30D+0, 1.30D+0, 1.29D+0, 1.29D+0, 1.29D+0, 1.28D+0, 1.27D+0,
     *  1.27D+0, 1.27D+0, 1.27D+0, 1.25D+0, 1.23D+0, 1.23D+0, 1.23D+0,
     *  1.23D+0, 1.23D+0, 1.23D+0, 1.21D+0, 1.20D+0, 1.20D+0, 1.20D+0,
     *  1.20D+0, 1.18D+0, 1.19D+0, 1.20D+0, 1.19D+0, 1.19D+0, 1.19D+0,
     *  1.19D+0, 1.19D+0, 1.18D+0, 1.19D+0, 1.19D+0, 1.19D+0, 1.17D+0,
     *  1.15D+0, 1.17D+0, 1.14D+0, 1.16D+0, 1.13D+0, 1.13D+0, 1.13D+0,
     *  1.12D+0, 1.12D+0, 1.11D+0, 1.10D+0, 1.09D+0, 1.09D+0, 1.09D+0,
     *  1.08D+0, 1.09D+0, 1.08D+0, 1.08D+0, 1.10D+0, 1.08D+0, 1.07D+0,
     *  1.07D+0, 1.07D+0, 1.07D+0, 1.07D+0, 1.07D+0, 1.06D+0, 1.05D+0,
     *  1.04D+0, 1.04D+0, 1.04D+0, 1.03D+0, 1.02D+0, 1.02D+0, 1.00D+0/
      DATA C1/
     * -8.88D-16,-6.35D+1, 2.34D+1,-1.78D+1, 1.22D+1,-1.26D+1, 9.44D-1,
     * -3.09D+0,-1.50D+0,-1.58D+0,-1.56D+0,-1.43D+0,-1.31D+0,-1.35D+0,
     * -1.22D+0,-1.11D+0,-1.15D+0,-1.16D+0,-1.03D+0,-1.05D+0,-1.01D+0,
     * -9.17D-1,-9.23D-1,-9.09D-1,-9.24D-1,-8.81D-1,-8.41D-1,-8.83D-1,
     * -8.58D-1,-8.08D-1,-8.34D-1,-8.27D-1,-7.84D-1,-7.58D-1,-1.05D+6,
     *  4.13D+5,-1.24D+5, 4.16D+4,-1.16D+4, 3.30D+3,-5.82D+2, 1.78D+2,
     * -1.87D+1, 1.06D+1,-8.35D+0, 2.19D+1,-1.07D+1, 7.46D+0, 5.25D+0,
     *  5.26D+0,-8.76D+0, 7.92D+0, 4.74D+0,-7.09D+0, 1.12D+1,-1.64D+0,
     *  1.03D+0, 3.88D-1,-1.96D+0, 6.25D+0,-2.23D+0,-1.91D+0, 4.37D+0,
     *  5.97D-1,-3.21D+0, 3.95D+0, 3.28D-1, 2.29D+0, 4.87D+9,-5.97D+8,
     *  1.10D+8,-1.77D+7, 1.65D+6,-2.59D+5, 2.53D+4,-6.56D+3, 9.76D+2,
     * -6.31D+2, 1.01D+3,-1.09D+3, 6.05D+2,-7.34D+1, 7.24D-2,-4.67D+2,
     *  5.56D+2,-1.06D+2,-3.94D+2, 6.10D+2,-4.28D+2, 8.90D+1,-2.14D+1,
     * -7.81D+1, 2.74D+2,-2.83D+2, 1.05D+1, 2.10D+2,-1.26D+2,-1.27D+2,
     *  2.39D+2,-1.21D+2, 6.55D+1,-9.82D+1/
      DATA C21/
     * -1.32D+0,-1.30D+0,-1.29D+0,-1.20D+0,-1.08D+0,-1.12D+0,-1.15D+0,
     * -1.06D+0,-9.53D-1,-1.05D+0,-9.44D-1,-8.96D-1,-8.40D-1,-7.85D-1,
     * -7.50D-1,-6.57D-1,-7.86D-1,-6.62D-1,-7.83D-1,-7.48D-1,-5.46D-1,
     * -6.18D-1,-7.86D-1,-6.19D-1,-5.08D-1,-5.38D-1,-5.33D-1,-5.57D-1,
     * -5.27D-1,-4.25D-1,-3.95D-1,-4.52D-1,-3.89D-1,-4.52D-1,-4.37D-1,
     * -3.27D-1,-3.92D-1,-4.53D-1,-3.69D-1,-2.79D-1,-3.64D-1,-3.74D-1,
     * -2.08D-1,-2.01D-1,-4.12D-2,-4.82D-3,-2.72D-1,-2.72D-1,-2.32D-1,
     * -2.85D-1,-1.50D-1,-2.74D-1,-1.50D-1,-1.47D-1,-2.50D-1,-2.82D-1,
     * -2.30D-1,-4.62D-2,-1.09D-1, 9.37D-3,-1.56D-1,-1.79D-1,-9.21D-2,
     *  1.84D-2,-2.42D-3,-2.14D-1,-2.50D-1,-1.20D-1, 7.88D-2,-2.46D-2,
     * -1.51D-2,-2.43D-2,-7.13D-2,-1.61D-1,-6.90D-2, 4.10D-2,-1.40D-3,
     * -1.43D-1,-5.99D-2, 1.24D-1, 1.27D-2,-2.62D-2, 1.83D-2,-1.29D-2,
     * -7.50D-3,-1.06D-1,-5.70D-3, 1.06D-1,-1.84D-2,-5.20D-2,-2.98D-1,
     *  3.12D-2,-7.91D-2,-3.64D-2, 3.25D-2,-2.25D-1, 8.50D-2,-4.70D-2,
     * -1.03D-2,-2.57D-2,-4.80D-2,-3.25D-2,-7.97D-3,-1.54D-2,-1.28D-2,
     * -8.27D-3, 1.09D-2,-3.91D-2, 4.01D-2, 5.86D-3,-6.76D-2, 1.17D-2,
     * -1.24D-2,-1.12D-2, 1.96D-3,-2.00D-3,-1.31D-2,-2.62D-3,-3.72D-3,
     * -2.18D-3,-1.30D-3,-5.91D-4,-1.91D-3,-4.06D-4,-1.72D-3/
      DATA C22/
     *  4.06D+0,-2.76D+0, 4.47D+0, 3.79D+0, 8.80D+0,-1.26D+1, 8.78D+0,
     *  8.23D-1, 9.69D+0,-1.92D+1, 2.97D+1,-2.50D+1, 3.06D+1,-2.51D+1,
     *  2.87D+1,-1.95D+1, 6.59D+0, 5.84D+0,-1.79D+1, 2.15D+1,-1.36D+0,
     * -5.82D+0,-2.59D+0, 1.10D+1,-5.42D+0, 3.89D+0,-3.64D+0, 2.45D+0,
     * -9.30D-1, 6.01D+0,-4.49D+0, 1.61D+0, 1.57D+0,-4.71D+0, 5.46D+0,
     *  3.99D-2,-3.30D+0, 2.36D-1, 3.95D+0, 5.60D-1,-4.79D+0, 4.30D+0,
     *  4.01D+0,-3.66D+0, 1.16D+1,-9.82D+0,-3.55D+0, 3.56D+0,-2.77D+0,
     *  1.69D+0, 1.01D+0,-3.49D+0, 5.96D+0,-5.91D+0, 3.86D+0,-4.50D+0,
     *  5.55D+0,-1.87D+0, 6.19D-1, 1.75D+0,-5.05D+0, 4.60D+0,-2.86D+0,
     *  5.07D+0,-5.49D+0, 1.27D+0,-2.00D+0, 3.31D+0,-1.32D+0, 2.86D-1,
     * -1.90D-1, 9.75D-2,-5.67D-1,-3.30D-1, 1.25D+0,-1.49D-1,-2.75D-1,
     * -1.14D+0, 1.97D+0,-1.24D-1,-9.91D-1, 6.02D-1,-1.57D-1,-1.56D-1,
     *  2.10D-1,-1.19D+0, 2.20D+0,-1.08D+0,-1.63D-1,-1.73D-1,-2.29D+0,
     *  5.58D+0,-6.68D+0, 7.11D+0,-6.42D+0, 3.84D+0,-7.40D-1, 2.12D-1,
     * -6.59D-2, 4.31D-3,-9.35D-2, 1.56D-1,-5.74D-2, 2.76D-2,-1.70D-2,
     *  3.51D-2, 4.17D-2,-2.42D-1, 5.58D-1,-6.95D-1, 4.01D-1,-8.43D-2,
     * -1.21D-2, 1.71D-2, 3.54D-2,-5.12D-2, 6.95D-3,-1.72D-3, 1.17D-3,
     * -3.97D-4, 8.37D-4,-4.84D-4, 2.20D-4, 8.02D-5,-2.12D-4/
      DATA C23/
     * -2.27D+2, 2.41D+2,-2.29D+1, 1.67D+2,-7.14D+2, 7.13D+2,-2.65D+2,
     *  2.96D+2,-9.64D+2, 1.63D+3,-1.82D+3, 1.85D+3,-1.86D+3, 1.80D+3,
     * -1.61D+3, 8.69D+2,-2.48D+1,-7.93D+2, 1.31D+3,-7.61D+2,-1.49D+2,
     *  5.39D+1, 2.26D+2,-2.73D+2, 1.55D+2,-1.26D+2, 1.02D+2,-5.64D+1,
     *  1.16D+2,-1.75D+2, 1.02D+2,-7.04D-1,-1.05D+2, 1.69D+2,-9.03D+1,
     * -5.57D+1, 5.90D+1, 6.20D+1,-5.66D+1,-8.92D+1, 1.51D+2,-4.88D+0,
     * -1.28D+2, 2.55D+2,-3.58D+2, 1.04D+2, 4.74D+1,-4.22D+1, 2.97D+1,
     * -4.55D+0,-3.00D+1, 6.30D+1,-7.91D+1, 6.51D+1,-5.57D+1, 6.70D+1,
     * -4.95D+1, 1.66D+1, 7.52D+0,-4.53D+1, 6.43D+1,-4.97D+1, 5.29D+1,
     * -7.04D+1, 4.51D+1,-2.18D+1, 1.77D+1,-1.54D+1, 5.35D+0,-1.59D+0,
     *  9.59D-1,-2.21D+0, 7.89D-1, 5.27D+0,-4.66D+0,-4.20D-1,-2.88D+0,
     *  1.03D+1,-6.97D+0,-2.89D+0, 5.31D+0,-2.53D+0, 1.66D-3, 1.22D+0,
     * -4.68D+0, 1.13D+1,-1.09D+1, 3.05D+0,-3.17D-2,-7.04D+0, 2.62D+1,
     * -4.08D+1, 4.60D+1,-4.51D+1, 3.42D+1,-1.53D+1, 1.27D+0,-3.71D-1,
     *  9.37D-2,-1.30D-1, 3.32D-1,-2.84D-1, 1.13D-1,-5.94D-2, 6.95D-2,
     *  8.79D-3,-3.78D-1, 1.07D+0,-1.67D+0, 1.46D+0,-6.47D-1, 9.63D-2,
     *  3.90D-2, 2.43D-2,-1.16D-1, 7.76D-2,-1.45D-3, 4.82D-4,-2.61D-4,
     *  2.06D-4,-2.20D-4, 4.69D-5,-9.34D-6,-9.74D-6, 1.99D-5/
C
      IF(T.LE.0.D0)GOTO 103
      IF(T.LT.X1(1))GOTO 101
      IF(T.GE.X2(126))GOTO 102
      IF(T.GE.X1(35))GOTO 20
C
      I1=1
      IF(T.GE.X1(9))I1=9
      IF(T.GE.X1(18))I1=18
      IF(T.GE.X1(27))I1=27
      DO 10 K=I1,(I1+8)
      I=K
      IF(X1(I+1).GT.T)GOTO 11
   10 CONTINUE
   11 D=T-X1(I)
      S=((C1(I,3)*D+C1(I,2))*D+C1(I,1))*D+Y1(I)
      FSYNCH=F1(T,A1)*S
      RETURN
C
   20 I1=1
      IF(T.GE.X2(10))I1=10
      IF(T.GE.X2(20))I1=20
      IF(T.GE.X2(30))I1=30
      IF(T.GE.X2(40))I1=40
      IF(T.GE.X2(50))I1=50
      IF(T.GE.X2(60))I1=60
      IF(T.GE.X2(70))I1=70
      IF(T.GE.X2(80))I1=80
      IF(T.GE.X2(90))I1=90
      IF(T.GE.X2(100))I1=100
      IF(T.GE.X2(110))I1=110
      IF(T.GE.X2(120))I1=120
      DO 22 K=I1,(I1+9)
      I=K
      IF(X2(I+1).GT.T)GOTO 21
   22 CONTINUE
   21 D=T-X2(I)
      S=((C2(I,3)*D+C2(I,2))*D+C2(I,1))*D+Y2(I)
      FSYNCH=F2(T,A2)*S
      RETURN
  101 FSYNCH=F1(T,A1)
      RETURN
  102 FSYNCH=F2(T,A2)
      RETURN
  103 FSYNCH=0.D0
      RETURN
      END



*************************************************************************
      
