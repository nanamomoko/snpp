#+
# Name:
#	snpp
# PURPOSE:
#	calculate the S/N per pixel for CSST and simulate a noisy spectrum for any given template.
# CALLING SEQUENCE:
#   snpp,limitmag, repeatnum=10,obstime=300,targetmag=18,/skyperpixel,$
#         galtpl=,wavearr=wavearr,mockgal=mockgal,galflux=galflux
#   plot, wavearr, galflux  ; the input galaxy template
#   plot, wavearr, mockgal  ; the output spectrum with noise
#     
# INPUTS:
# OPTIONAL IUTPUTS:
#   darkcurrent    dark current, in e/s/pix, (defult: 0.0017)
#   deltal         the delta lambda per pixel, in unit of nm (defult: 0.1755555 nm)
#   fovp           diameter of fiber (or spaxel) in arcsec (defult: 0.2 arcsec)
#   filtera         the filter you chosed to estimate the S/N (defult: bessell_V)
#   galtpl         the filename of star-forming galaxy template you want to use. 
#                  They are in the ../obs/SFgal_tpl/ folder (default: SFgal_texp_FeH0_tau5_Ew10.fits)
#   lambdac        the noise at this wavelength wanted (defult: 550 nm)
#   npixel_width   the width of the spectrum on the CCD (defult: 3.0)
#   obstime        in seconds, single integration time (defult: 300s)
#   outfile        the output file name (defult: '../results/noise.dat' )
#   qinput         the throughput correct factor (defult: 1.0)
#   readnoise      read noise, in e/pix. (defult: 4.0)
#   redshift       the redshift of the target spectrum. (defult: 0.0)
#   repeatnum      repeat number (defult: 1.0)
#   skyperpixel    a second way of estimating the Sky, if know the sky photon number per pixel
#   skyv           V band sky brightness in Johnson V mag/arcsec^2 unit (defult: 22.5 mag/arcsec^2)
#  slitwidth      suit to the slit case. the length assumed to be 0.15 arcsec
#   snlimit        S/N limit (defult: 1.0)
#   specsample     pixels per spectral resolution element (defult: 2)
#   targetmag      the surface brightness of the target you want to calculate the S/N (defult: 22 .5 mag/arcsec^2)
#   teld           diameter of the telescope, in cm unit. (defult: d=200 cm)
# OUTPUTS:
#   limitmag       the Vband surface brightness needed to achieve the S/N limit (defult: 1.0)
# OPTIONAL OUTPUTS:
#   limitemi       the medien of Flambda*dlambda*sampling value of Ha line
#   limitemif      the limit detection of Ha flux 
#   snmean         the median S/N of the whole input target spectrum (mag_v=targetmag)
#   wavearr        the wave array (nm)
#   galflux        the input galaxy flux  (1e-13 erg/s/cm2/nm)
#   mockgal        the mocked galaxy flux  with noise (1e-13 erg/s/cm2/nm)
#
# v5: 15 August 2018      writen by Lei Hao, rivised by Jun Yin
# v7: 10 Sep 2019  by Jun Yin
#     1) remove the function im_filtermag, so do not need the Kcorrect package anymore.
#     2) 
#python 
# v7: 22 Sep 2019 by Mengting Ju
#-

#####################################################################################

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib
import pandas as pd
from scipy import interpolate
from sympy import *
import os

####################################################################################

def integral(x,y):
    
    nn=len(x)
    
    dx=x[1:]-x[:-1]
    yy=0.5*(y[1:]+y[:-1])
  
    return np.sum(dx*yy)

####################################################################################

class snpp(object):
    def __init__(self, limitmag=1.0, lambdac=550, deltal=0.1755555, qinput=1.0, fovp=0.2,
                 slitwidth=None,obstime=300, skyv=22.5, targetmag=None, repeatnum=1.0, 
                 outfile=False, spectype=None, teld=200, snmean=False, specsample=2, 
                 snlimit=1.0, readnoise=4.0, skyperpixel=None, npixel_width=3.0, 
                 limitemif=None,  darkcurrent=0.017, redshift=0.0, galtpl=False,
                 wavearr=None, galflux=False, mockgal=False, filtera=False):
        '''
        ; mydevice=!D.name
        ;   !p.font  = 0
        ;   !p.thick = 3
        ;   !x.thick = 3
        ;   !y.thick = 3
        ;   !p.charsize  = 1.0
        ;   !p.charthick = 8
        ; set_plot,'ps'  
        ;  device,file = '../graph/test.ps',/color,$
        ;   ysize=10.0,xsize=30.0,/iso, times,xoffset=0,yoffset=0
        ; loadct,39
        '''
        
        #Do extensive checking of possible input errors
        #
        self.limitmag=limitmag 
        self.lambdac=lambdac 
        self.deltal=deltal 
        self.qinput=qinput 
        self.fovp=fovp
        self.slitwidth=slitwidth
        self.obstime=obstime
        self.skyv=skyv
        self.targetmag=targetmag
        self.repeatnum=repeatnum
        self.teld=teld
        self.specsample=specsample
        self.snlimit=snlimit
        self.readnoise=readnoise
        self.npixel_width=npixel_width
        self.darkcurrent=darkcurrent
        self.redshift=redshift
        
        ###########################################################################
        #some basic unchanged parameters
        d=200.        # diameter of the telescope, in cm unit

        if self.teld:
            d=teld
        print('d:', d)
        
        obscure=0.0  #effective central obscuration, no unit
        telarea=3.14159/4.0*d*d*(1.0-obscure)  #effective area of the telescope, cm^2
        darkc=0.017   #dark current, in e/s/pix
        
        if self.darkcurrent:
            darkc=darkcurrent
        print('darkc:', darkc)
        
        rn=4.     #read noise, in e/pix

        if self.readnoise:
            rn=readnoise
        print('rn:', rn)
        
        planckh=6.626    # 10^{-27} erg*s
        cc=3.0   # speed of light, 10^{17} nm/s
        ####################################################################
        
        #load the filters
        if  filtera:
            
            filtersel=filtera
            
        else:
            filtersel='bessell_V.par'   #'../sdss_g0.par'
            
        filterpath='../obs/filters/'
        filterfile=filterpath+filtersel
        print(filterfile)
        # ;fluxfilter: max=1, min=0, no particular unit
        
        ia=0
        with open(filterfile,'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    ia=ia+1
                    continue

        band=pd.read_csv(filterfile,sep='\s+',header=None,skiprows=ia)
        wavefilter=np.array(band[0])
        fluxfilter=np.array(band[1])
        wavefilter=wavefilter/10.0  # in nm
        vmin=wavefilter[0]
        nw=len(wavefilter)
        vmax=wavefilter[nw-1]

        # find the central wavelength, effective wavelength, and FWHM of the given filter

        filtermid=(vmax-vmin)*0.5  #nm, central wavelength
        dwave=wavefilter[1:]-wavefilter[:-1]
        filtereff=np.nansum(dwave*wavefilter[1:]*fluxfilter[1:])/np.nansum(dwave*fluxfilter[1:]) #nm, effective wavelength
        rmax=np.max(fluxfilter)
        nnn=np.where(fluxfilter > 0.5*rmax)[0]
        FWHMmin=wavefilter[nnn[0]]
        FWHMmax=wavefilter[nnn[-1]]
        filterwid=FWHMmax-FWHMmin  #nm, FWHM
        
        ######################################################################
        
        # define wavelength array,
        #cover the range of 350nm to 1050nm, depend on the spectral resolution wanted. 
        
        #specr0=2000  ; no unit
        #if(keyword_set(specr)) then specr0=specr
        sampling=2.0    #pixels per spectral resolution element ?1D or 2D/linear or area?

        if self.specsample:
            sampling=specsample
        print('sampling:', sampling)
        
        #delta_lambda=500.0/specr0/sampling  ; in the same unit as lambda0
        delta_lambda=0.1755555
        
        if self.deltal:
            delta_lambda=deltal # has to be in unit of nm
        print('delta_lambda:', delta_lambda)
        
        narray=int((1000.0-350.0)/delta_lambda)  
        #figure out the wavelength array length, from 350nm to 1000nm, spacing at delta_lambda
        wavearr=350.0+delta_lambda*pl.frange(narray-1)
        # select out the array of V band filter
        ii=np.logical_and(wavearr >= vmin, wavearr <= vmax)
        wavetmp2=wavearr[ii]
        x=np.interp(wavetmp2,wavefilter,fluxfilter)
        integratef4=x*wavetmp2
        integconst=integral(wavetmp2, integratef4) # int(lambda*Rlambda*dlambda)
        
        #####################################################################
        
        # some less basic parameters, may change, but not often
        #qsys=0.10	; throughput of the whole system, should be a function of lambda
        
        throughput=pd.read_csv('../obs/IFU_throughput.dat',sep='\s+',header=None,skiprows=1)
        lambdaq=np.array(throughput[8])
        qtot=np.array(throughput[9]) #; throughput of the whole system,
        
        if not self.qinput:
            qinput=1.0
        print('qinput:', qinput)   

        qe=0.8 
        #;assuming the total throughput cannot reach the theory value, 0.3 is the upper limit. 
        qtot[qtot>=0.3]=0.3 
        q=qtot*qinput #*qe ;qtot of CSST already includes the CCD efficiency 
        fovsp=0.2  # diameter of fiber (or spaxel) in arcsec ?
        #fov2=3.14159/4.0*(0.2)^2      ; fiber area in (arcsec)^2
        fov2=(fovsp)**2     #*3.14159/4.0      ; fiber (or spaxel) area in (arcsec)^2
        if self.fovp:
            fov2=(fovp)**2 #*3.14159/4.0
        #   for slit (point source)
        if self.slitwidth:
            fov2=1 
        print('fov2:', fov2)
        
        slitunit=0.074  # arcsec. the length of slit which conresponds to a pixel length on IFU CCD 
        
        ##############################################################################
        
        # SKY
        #define V band sky brightness
        
        iskyv0=22.5  # in Johnson V mag/arcsec^2 unit
        if skyv:
            iskyv0=skyv
        print('iskyv0:', iskyv0)
        
        lambdav=filtereff   #in nm

        #sky brightness corresponding to this sky magnitude
        iskyv0_jy=3631.0*10**(-iskyv0/2.5+3.0)  # sky flux in V in mJy/arcsec^2 unit
        iskyv0_nm=iskyv0_jy*3.0/(lambdav/100.0)**2 #sky flux in V in 10^(-13)erg/s/cm^2/nm (/arcsec^2 ?)

        #readin the ground sky spectrum 
        skybg_50=pd.read_csv('../obs/skybg_50_10.dat',sep='\s+',header=None,skiprows=14)
        wavesky=np.array(skybg_50[0])
        fluxsky1=np.array(skybg_50[1])
        fluxsky2=fluxsky1/wavesky*1.98 #change the sky flux unit to 10^(-13)erg/s/cm^2/nm/arcsec^2

        #This fluxsky is in unit of phot/s/nm/arcsec^2/m^2, to convert it to F_lambda/arcsec^2, 
        #need to do fluxsky(phot/s/nm/arcsec^2/m^2)*h(6.625*10^{-27}erg.s)*nu(1/s)*10{-4}(m^2/cm^2)
        #=fluxsky*c(3.0*10^{17}nm/s)/lambda(nm)*6.6*10{-31} erg/s/cm^2/nm/arcsec^2
        #=fluxsky/lambda*1.98*10^{-13}erg/s/cm^2/nm/arcsec^2 

        #find out the normalization of the sky,
        #normalization=iskyv0_nm*(integrate(bandpass*lambda*dlambda)/integrate(bandpass*lambda*F_sky_lambda*dlambda))
        ii=np.logical_and(wavesky >= vmin, wavesky <= vmax)
        wavetmp=wavesky[ii]
        fluxtmp=fluxsky1[ii]

        x=np.interp(wavetmp,wavefilter,fluxfilter)
        vfluxtmp=x*fluxtmp*1.98  #bandpass*lambda*F_sky_lambda(fluxsky2)=bandpass*fluxsky*1.98, x10^(-13)
        skyintegrate=integral(wavetmp, vfluxtmp)
        skynorm=iskyv0_nm*integconst/skyintegrate 
        fluxsky3=np.interp(wavearr,wavesky,fluxsky2)
        
        fluxsky=fluxsky3*skynorm   
        # get the sky spectrum in wavearr grid, the unit should now be the same as fluxvega: 10^(-13) erg/s/nm/cm^2  (/arcsec^2 ?)

        fluxskypp=fluxsky
        print('fluxskypp:',fluxskypp)
        #a second way of estimating the Sky, if know the sky electron number per pixel

        if skyperpixel is not None:
                #since the numbers are given by the main survey, our detected Sky electron will be less, so scale a rough factor of 0.9
                scaletemp=0.9
                ii=np.logical_and(wavearr >= 255, wavearr <= 400)
                counta=len(np.where(ii==1)[0])
                fluxskypp[ii]=0.028/counta  
                ii=np.logical_and(wavearr >= 400, wavearr <= 600)
                countb=len(np.where(ii==1)[0])
                fluxskypp[ii]=0.229/countb  
                ii=np.logical_and(wavearr >= 600, wavearr <= 900)
                countc=len(np.where(ii==1)[0])
                fluxskypp[ii]=0.301/countc  
                ii=np.where(wavearr >= 900)[0]
                countd=len(ii)
                fluxskypp[ii]=0.301/countd 
                fluxskypp=fluxskypp/0.074^2*fov2*scaletemp
            
        print('fluxskypp:', fluxskypp)
        
        # define basic target brightness, parameters constantly change
        itarget=22.5  # in Johnson V mag/arcsec^2 unit
        if self.targetmag:
            itarget=targetmag
        print('itarget:', itarget)
        
        itarget_jy=3631.0*10**(-itarget/2.5+3.0)  # target flux in V in mJy/arcsec^2 unit
        itarget_nm=itarget_jy*3.0/(lambdav/100.0)**2 #target flux in V in 10^(-13)erg/s/cm^2/nm (/arcsec^2 ?)

        #readin the galaxy spectrum
        '''
            ; readcol, '../obs/allgalaxy.dat', wavegal, eflux, s0f, saf, sbf, scf, /silent
            ; wavegal=wavegal*1000.0
            ; spectype0=4	;default use Sb galaxy template spectrum
            ; if(keyword_set(spectype)) then spectype0=spectype  ;unless specified
            ; if(spectype0 eq 1) then galflux1=eflux
            ; if(spectype0 eq 2) then galflux1=s0f
            ; if(spectype0 eq 3) then galflux1=saf
            ; if(spectype0 eq 4) then galflux1=sbf
            ; if(spectype0 eq 5) then galflux1=scf
        '''

        
        if galtpl :
            galtpl=galtpl
        else:
            galtpl='SFgal_texp_FeH0_tau5_Ew10.fits'
            
        tplfile='../obs/SFgal_tpl/'+galtpl
        print('tplfile:',tplfile)   
        
        sfgal=fits.open(tplfile)
        wavegal=sfgal[1].data['wave']/10. #change A to nm
        galflux2=sfgal[1].data['flux']
        galflux1=np.interp(wavearr,wavegal,galflux2)

        #;normalize the galaxy spectrum to the V band magnitude specified.
        ii=np.logical_and(wavegal >= vmin, wavegal <= vmax)
        wavetmp=wavegal[ii]
        fluxtmp=galflux2[ii]
        x=np.interp(wavetmp,wavefilter,fluxfilter)
        vfluxtmp=x*wavetmp*fluxtmp #bandpass*lambda*F_gal_lambda
        galintegrate=integral(wavetmp,vfluxtmp)
        galnorm=itarget_nm*integconst/galintegrate
        galflux=galnorm*galflux1   # the unit should now be in 10^(-13)erg/s/nm/cm^2 (/arcsec^2 ?)
        
        ##########################################################################
        
        #define observation information, parameters constantly change
        obst=300.0  # in seconds, single integration time
        if self.obstime:
            obst=obstime
        print('obst:',obst)
        repn=1.0   # repeating time
        if self.repeatnum:
            repn=repeatnum
        print('repn:',repn)
        npixw=3.0
        if self.npixel_width:
            npixw=npixel_width
        print('npixw:',npixw)
        #   sky of slit area (slitwidth*npixw*slitlength) will go into the CCD
        if self.slitwidth:
            fluxskypp=fluxskypp*slitwidth*npixw*slitunit  
        print('fluxskypp:',fluxskypp)
        
        expf2=np.zeros(narray)
        expfemi=np.zeros(narray)
        snarray=np.zeros(narray)
        mockgal=np.zeros(narray)
        tmp=np.zeros(narray)
        lista=np.zeros(narray*10).reshape(narray,10)
        
        for i in range(narray):
            lambda0=wavearr[i]
            qlambda=np.interp(lambda0,lambdaq,q)
            hv=planckh*cc/lambda0 #;10^{-10}erg
            delta_hz=cc*delta_lambda/lambda0/lambda0*sampling #;10^17 1/s

            #now that many fluxes are in 10^(-13)erg/s/nm/cm^2, to convert it to Jy, need to multiple: 
            #lambda0^2/c(in nm)=lambda0^2(nm)/(3.*10^(17))*10^(-13)erg/s/Hz/cm^2
            #=lambda^2(nm)*3.33*10^(-31)erg/s/Hz/cm^2=lambda^2(nm)*3.33*10^(-8)Jy

            #find out sky value at lambda0    
            #calculate n_sky/pixel
            isky=fluxsky[i]*lambda0**2*0.0333*fov2   #in uJy/spaxel unit
            iskyall=isky*telarea/1000.0   #in 10-26 erg/s/Hz /spaxel
            fsky=qlambda*iskyall*delta_hz   #10^{-9} erg/s /spaxel
            nsky=fsky/hv*10.0   #in unit of #e/s /spaxel

            '''
            if(keyword_set(skyperpixel)) then begin
                nsky=fluxskypp[i]*sampling  ; #e/s in npixw*sampling pixels 
            endif
            ;print, "Sky electron counts", nsky, nsky0, fluxskypp[i]
            '''

            #calculate n_source/pixel
            isource=galflux[i]*lambda0**2*0.0333*fov2   #in uJy/spaxel unit
            isall=isource*telarea/1000.0   #in 10-26 erg/s/Hz /spaxel
            fs=qlambda*isall*delta_hz   #10^{-9} erg/s /spaxel
            ns=fs/hv*10.0   #in unit of #e/s /spaxel
            #print, "Source electron counts", ns

            darkn=(darkc*repn*obst*npixw*sampling)
            rnn2=rn**2*(repn*npixw*sampling)
            sourcenn=(ns*repn*obst)
            skynn=(nsky*repn*obst)
            tmp[i]=skynn

            #nn1=sqrt(2.0*rnn^2+2.0*darkn^2+sourcenn^2+2.0*skynn^2)
            #nn1=sqrt(rnn^2+darkn^2+sourcenn^2+skynn^2)
            nn1=np.sqrt(rnn2+darkn+skynn+sourcenn)  #total noise
            sn1=repn*ns*obst/nn1  #S/N
            snarray[i]=sn1
            #nn=sqrt(2.0*rnn^2+2.0*darkn^2+2.0*skynn^2)
            nn=np.sqrt(rnn2+darkn+skynn)  #system noise
            #print, "total noise, system noise, sn, readnoise, dark, source, sky", nn1, nn, sn1, rnn, darkn, sourcenn, skynn 

            #set the detection limit 
            detlimit=1.0
            #if(keyword_set(snlimit)) then detlimit=snlimit
            #N_{source}/sqrt(N_{source}+nn^2)=detlimit, ==m
            #N_{source}^2-m^2*N_{source}-m^2*nn^2=0 (solve the equation)
            #N_{source}=(m^2+sqrt(m^4+4m^2*nn^2))/2.0
            nntmp=detlimit**2+np.sqrt(detlimit**4+4.0*detlimit**2*nn**2)
            nntmp=nntmp/2.0

            #calculate detection limit in uJy and mag

            fnn=nntmp
            f1=fnn/obst/repn    #in e/s
            f2=f1*hv       #in 10^{-10} erg/s
            f3=f2/delta_lambda #in 10^{-10} erg/s/nm
            f1=f3/telarea #in 10^{-10} erg/s/nm/cm^2
            f2=f1/qlambda    #in 10^{-10} erg/s/nm/cm^2
            expf2[i]=f2/fov2*100000.   # in 10^{-15} erg/s/nm/cm^2/arcsec^2
            #expfemi[i]=expf2[i]*delta_lambda*sampling   # in 10^{-15} erg/s/cm^2/arcsec^2, the above multiplied by the spectral resolution
            #print, "detection limit is", f2,"microJy/arcsec^2"
            #print, "detection limit is", magf2, "AB mag / arcsec^2"
            mockgal[i]=galflux[i]+galflux[i]/snarray[i]*np.random.randn(1,1)[0][0]  #in 10^{-13} erg/s/nm/cm^2

            lista[i,:]=[lambda0, sn1, galflux[i], nn1,\
                        np.sqrt(sourcenn), nn, np.sqrt(rnn2),np.sqrt(darkn), \
                        np.sqrt(skynn), mockgal[i]]
            
        '''
        ;mockgal=galflux+galflux/snarray*randomn(seed,narray)
        ; plot,wavearr,galflux,xrange=[350,1000],xs=1,xtitle='lambda (nm)',ytitle='galaxy flux'
        ; oplot,wavearr,mockgal,color=250
        ; oplot,wavearr,galflux,thick=3
        ; label=strcompress(string(obst))+'s *'+strcompress(string(repn))+' times'
        ; xyouts,800,max(galflux),label,charsize=2
        ; plot,wavearr,expf2*0.01,xrange=[350,1000],xs=1,xtitle='lambda (nm)',ytitle='expected galaxy flux'
        ; oplot,wavearr,galflux,color=250

        ; plot,wavearr,fluxsky,xrange=[350,1000],xs=1,xtitle='lambda (nm)',ytitle='sky flux'
        ; plot,wavearr,tmp,xrange=[350,1000],xs=1,xtitle='lambda (nm)',ytitle='#photon'


        ; ii=where(wavearr ge vmin and wavearr le vmax)
        ; wavetmp=wavearr(ii)
        ; fluxtmp=expf2(ii)/100000.  ;10^{-10} erg/s/cm^2/arcsec^2
        ; x=interpol(fluxfilter, wavefilter, wavetmp)
        ; vfluxtmp=x*wavetmp*fluxtmp  ;bandpass*lambda*F_gal_lambda
        ; gexpintegrate=integral(wavetmp, vfluxtmp)
        ; magf2=-2.5*(alog10((gexpintegrate*(lambdav/100.0)^2)/(integconst*3631.0*3.0)))
        ; print,'magf2=',magf2
        ; plot,wavetmp,fluxtmp,xrange=[350,1000],xs=1,xtitle='lambda',ytitle='expected galaxy flux'
        '''
        
        ii=np.logical_and(wavearr >= FWHMmin , wavearr <= FWHMmax)
        wavetmp=wavearr[ii]
        if len(snarray[ii]) % 2 ==0:
            snm=sorted(list(snarray[ii]))[int(0.5*len(snarray[ii]))]
        else:
            snm=np.median(snarray[ii])    #the median SN of FWHM range to acchieve the sn limit
        
        im=np.where(snarray[ii] == snm)[0]
        
        fact=np.reshape(expf2[ii][im]*0.01/galflux[ii][im],1)
        fact=fact[0]
        limitmag=-2.5*np.log10(fact)+itarget  

        #print,limitmag
        #oplot,wavetmp, fluxtmp, color=250
        snmean=np.median(snarray)
        
        z=0.0
        if self.redshift:
            z=redshift
        print('z:',z)
        
        waveha=656.3*(1.0+z)
        ii=np.logical_and(wavearr >= (waveha-0.5) , wavearr < (waveha+0.5))  #1nm ,10A
        nii=len(np.where(ii==1)[0])
        ii_1=np.logical_and(wavearr >= (waveha-10) , wavearr <= (waveha-5))
        ii_2=np.logical_and(wavearr <= (waveha+10) , wavearr >= (waveha+5))
        icont=np.logical_or(ii_1==1,ii_2==1)

        contrms=np.sqrt(np.sum(mockgal[icont]**2)/len(mockgal[icont]))
        h=3*contrms*np.sqrt(nii) # hight>3 con
        w=1.  # width of 10A
        limitemif=np.sqrt(2*np.pi)*h*w  #h=3*cont, w=10A 
        
        ############################################################################################
        
        # write file
        namedat=np.array(['lambda','S/N','tar_flux','tot_noise','sc_noise', \
                          'sys_noise', 'readnoise','dark_noise', 'sky_noise', 'mockgal'])
        unit=np.array(['nm', ' ','1e-13 erg/s/cm2/nm',\
                       '#e','#e','#e','#e','#e','#e', '1e-13 erg/s/cm2/nm'])
        
        hdr=fits.Header()
        for i in range(len(namedat)):
            hdr[str(i)]=unit[i]
        hun1=fits.PrimaryHDU(header=hdr)
        hun2=fits.BinTableHDU.from_columns([fits.Column(name=namedat[i],array=np.array(lista[:,i]),format='1E') for i in range(len(namedat))])
        hdulist = fits.HDUList([hun1,hun2])
        
        if(os.path.exists('../results/noise_'+str(filtersel)+'_'+str(galtpl)+'.fits'))==1:
            os.remove('../results/noise_'+str(filtersel)+'_'+str(galtpl)+'.fits')
            
        hdulist.writeto('../results/noise_'+str(filtersel)+'_'+str(galtpl)+'.fits')
        
        #####################################################################
        
        print('The END!')
        
#######################################################################################################

snpp()
snpp(filtera='sdss_i0.par')
snpp(filtera='sdss_g0.par',galtpl='SFgal_texp_FeH0_tau1_Ewd.fits')

