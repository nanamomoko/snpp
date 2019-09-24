;+
; Name:
;	snpp
; PURPOSE:
;	calculate the S/N per pixel for CSST and simulate a noisy spectrum for any given template.
; CALLING SEQUENCE:
;   snpp,limitmag, repeatnum=10,obstime=300,targetmag=18,/skyperpixel,$
;         galtpl=,wavearr=wavearr,mockgal=mockgal,galflux=galflux
;   plot, wavearr, galflux  ; the input galaxy template
;   plot, wavearr, mockgal  ; the output spectrum with noise
;     
; INPUTS:
; OPTIONAL IUTPUTS:
;   darkcurrent    dark current, in e/s/pix, (defult: 0.0017)
;   deltal         the delta lambda per pixel, in unit of nm (defult: 0.1755555 nm)
;   fovp           diameter of fiber (or spaxel) in arcsec (defult: 0.2 arcsec)
;   filter         the filter you chosed to estimate the S/N (defult: SDSS g)
;   galtpl         the filename of star-forming galaxy template you want to use. 
;                  They are in the ../obs/SFgal_tpl/ folder (default: SFgal_texp_FeH0_tau5_Ew10.fits)
;   lambdac        the noise at this wavelength wanted (defult: 550 nm)
;   npixel_width   the width of the spectrum on the CCD (defult: 3.0)
;   obstime        in seconds, single integration time (defult: 300s)
;   outfile        the output file name (defult: '../results/noise.dat' )
;   qinput         the throughput correct factor (defult: 1.0)
;   readnoise      read noise, in e/pix. (defult: 4.0)
;   redshift       the redshift of the target spectrum. (defult: 0.0)
;   repeatnum      repeat number (defult: 1.0)
;   skyperpixel    a second way of estimating the Sky, if know the sky photon number per pixel
;   skyv           V band sky brightness in Johnson V mag/arcsec^2 unit (defult: 22.5 mag/arcsec^2)
;   slitwidth      suit to the slit case. the length assumed to be 0.15 arcsec
;   snlimit        S/N limit (defult: 1.0)
;   specsample     pixels per spectral resolution element (defult: 2)
;   targetmag      the surface brightness of the target you want to calculate the S/N (defult: 22.5 mag/arcsec^2)
;   teld           diameter of the telescope, in cm unit. (defult: d=200 cm)
; OUTPUTS:
;   limitmag       the Vband surface brightness needed to achieve the S/N limit (defult: 1.0)
; OPTIONAL OUTPUTS:
;   limitemi       the medien of Flambda*dlambda*sampling value of Ha line
;   limitemif      the limit detection of Ha flux 
;   snmean         the median S/N of the whole input target spectrum (mag_v=targetmag)
;   wavearr        the wave array (nm)
;   galflux        the input galaxy flux  (1e-13 erg/s/cm2/nm)
;   mockgal        the mocked galaxy flux  with noise (1e-13 erg/s/cm2/nm)
;
; v5: 15 August 2018      writen by Lei Hao, rivised by Jun Yin
; v7: 10 Sep 2019  by Jun Yin
;     1) remove the function im_filtermag, so do not need the Kcorrect package anymore.
;     2) 
;-

PRO snpp,limitmag, lambdac=lambdac, deltal=deltal, qinput=qinput, fovp=fovp, $
	slitwidth=slitwidth,obstime=obstime, skyv=skyv, targetmag=targetmag, repeatnum=repeatnum, $
	outfile=outfile, spectype=spectype, teld=teld, snmean=snmean, specsample=specsample, $
	snlimit=snlimit, readnoise=readnoise, skyperpixel=skyperpixel, npixel_width=npixel_width, $
	limitemif=limitemif, darkcurrent=darkcurrent,redshift=redshift,galtpl=galtpl,$
	wavearr=wavearr,galflux=galflux,mockgal=mockgal,filter=filter

; mydevice=!D.name
;   !p.font  = 0
;   !p.thick = 3
;   !x.thick = 3
;   !y.thick = 3
;   !p.charsize  = 1.0
;   !p.charthick = 8
; set_plot,'ps'  
;  device,file = '../graph/test.ps',/color,$
;   ysize=10.0,xsize=30.0,/iso,/times,xoffset=0,yoffset=0
; loadct,39

    outfile0='../results/noise.dat'
    if(keyword_set(outfile)) then outfile0=outfile
    openw, uu, outfile0, /get_lun
    printf, uu, F='(A1,A11,11A12)','#','lambda','S/N','tar_flux','tot_noise','sc_noise', $
       'sys_noise', 'readnoise','dark_noise', 'sky_noise', 'mockgal'
    printf, uu, F='(A1,A11,A12,A18,A6,5A12,A24,A18)','#','nm', ' ','1e-13 erg/s/cm2/nm',$
       '#e','#e','#e','#e','#e','#e', '1e-13 erg/s/cm2/nm'

;*******************************
;some basic unchanged parameters
    d=200.        ; diameter of the telescope, in cm unit
    if(keyword_set(teld)) then d=teld
    obscure=0.0  ;effective central obscuration, no unit
    telarea=3.14159/4.0*d*d*(1.0-obscure)  ; effective area of the telescope, cm^2
    darkc=0.017   ;dark current, in e/s/pix
    if(keyword_set(darkcurrent)) then darkc=darkcurrent
    rn=4     ;read noise, in e/pix
    if(keyword_set(readnoise)) then rn=readnoise
    planckh=6.626    ; 10^{-27} erg*s
    cc=3.0   ; speed of light, 10^{17} nm/s

;*******************************
;load the filters    
;    filtersel='./bessell/bessell_V.dat' ;'./sdss/sdss_g0.res'
    filtersel='./bessell_V.par' ;'./sdss_g0.par'
    if(keyword_set(filter)) then filtersel=filter
    filterpath='../obs/filters/'
    filterfile=filterpath+filtersel
    print,filterfile
    readcol, filterfile, wavefilter, fluxfilter;, /silent ;fluxfilter: max=1, min=0, no particular unit
    wavefilter=wavefilter/10.0  ; in nm
    vmin=wavefilter[0]
    vmax=wavefilter[n_elements(wavefilter)-1]
  ; find the central wavelength, effective wavelength, and FWHM of the given filter
      filtermid=(vmax-vmin)*0.5  ;nm, central wavelength
      dwave=wavefilter-shift(wavefilter,+1)
      nw=n_elements(wavefilter)
      filtereff=total(dwave[1:nw-1]*wavefilter[1:nw-1]*fluxfilter[1:nw-1])$
      /total(dwave[1:nw-1]*fluxfilter[1:nw-1]) ;nm, effective wavelength
      rmax=max(fluxfilter)
      nnn=where(fluxfilter gt 0.5*rmax)
      FWHMmin=wavefilter[nnn[0]]
      FWHMmax=wavefilter[nnn[n_elements(nnn)-1]]
      filterwid=FWHMmax-FWHMmin  ;nm, FWHM


;*******************************
;define wavelength array, cover the range of 350nm to 1050nm, depend on the spectral resolution wanted. 

    ;specr0=2000  ; no unit
    ;if(keyword_set(specr)) then specr0=specr
    sampling=2.0    ;pixels per spectral resolution element ?1D or 2D/linear or area?
    if(keyword_set(specsample)) then sampling=specsample
    ;delta_lambda=500.0/specr0/sampling  ; in the same unit as lambda0
    delta_lambda=0.1755555
    if(keyword_set(deltal)) then delta_lambda=deltal ; has to be in unit of nm
    narray=long((1000.0-350.0)/delta_lambda)  ;figure out the wavelength array length, from 350nm to 1000nm, spacing at delta_lambda
    wavearr=350.0+delta_lambda*findgen(narray)

    ; select out the array of V band filter
    ii=where(wavearr ge vmin and wavearr le vmax)
    wavetmp2=wavearr(ii)  
    x=interpol(fluxfilter, wavefilter,wavetmp2)
    integratef4=x*wavetmp2
    integconst=integral(wavetmp2, integratef4) ; int(lambda*Rlambda*dlambda)

;*******************************
;some less basic parameters, may change, but not often
;qsys=0.10	; throughput of the whole system, should be a function of lambda
    readcol, '../obs/IFU_throughput.dat', a,a,a,a,a,a,a,a,lambdaq,qtot,/silent ; throughput of the whole system,
    if(not keyword_set(qinput)) then qinput=1.0
    qe=0.8 
    qtot=qtot<0.3 ;assuming the total throughput cannot reach the theory value, 0.3 is the upper limit. 
    q=qtot*qinput;*qe ;qtot of CSST already includes the CCD efficiency 
    fovsp=0.2  ; diameter of fiber (or spaxel) in arcsec ?
    ;fov2=3.14159/4.0*(0.2)^2      ; fiber area in (arcsec)^2
    fov2=(fovsp)^2;*3.14159/4.0      ; fiber (or spaxel) area in (arcsec)^2
    if(keyword_set(fovp)) then fov2=(fovp)^2;*3.14159/4.0
;   for slit (point source)
    if(keyword_set(slitwidth)) then fov2=1 
    slitunit=0.074  ; arcsec. the length of slit which conresponds to a pixel length on IFU CCD 
;**************************************
;SKY
;define V band sky brightness 
    iskyv0=22.5  ; ; in Johnson V mag/arcsec^2 unit
    if(keyword_set(skyv)) then iskyv0=skyv
    lambdav=filtereff   ;in nm

    ;sky brightness corresponding to this sky magnitude
    iskyv0_jy=3631.0*10^(-iskyv0/2.5+3.0)  ; sky flux in V in mJy/arcsec^2 unit
    iskyv0_nm=iskyv0_jy*3.0/(lambdav/100.0)^2 ;sky flux in V in 10^(-13)erg/s/cm^2/nm (/arcsec^2 ?)

;readin the ground sky spectrum 
    readcol, '../obs/skybg_50_10.dat', wavesky, fluxsky1, /silent
    fluxsky2=fluxsky1/wavesky*1.98	;change the sky flux unit to 10^(-13)erg/s/cm^2/nm/arcsec^2
    ;This fluxsky is in unit of phot/s/nm/arcsec^2/m^2, to convert it to F_lambda/arcsec^2, 
    ;need to do fluxsky(phot/s/nm/arcsec^2/m^2)*h(6.625*10^{-27}erg.s)*nu(1/s)*10{-4}(m^2/cm^2)
    ;=fluxsky*c(3.0*10^{17}nm/s)/lambda(nm)*6.6*10{-31} erg/s/cm^2/nm/arcsec^2
    ;=fluxsky/lambda*1.98*10^{-13}erg/s/cm^2/nm/arcsec^2 

    ;find out the normalization of the sky,
    ;normalization=iskyv0_nm*(integrate(bandpass*lambda*dlambda)/integrate(bandpass*lambda*F_sky_lambda*dlambda))
    ii=where(wavesky ge vmin and wavesky le vmax)
    wavetmp=wavesky(ii)
    fluxtmp=fluxsky1(ii)
    x=interpol(fluxfilter, wavefilter, wavetmp)
    vfluxtmp=x*fluxtmp*1.98  ;bandpass*lambda*F_sky_lambda(fluxsky2)=bandpass*fluxsky*1.98, x10^(-13)
    skyintegrate=integral(wavetmp, vfluxtmp)
    skynorm=iskyv0_nm*integconst/skyintegrate 
    
    fluxsky3=interpol(fluxsky2, wavesky, wavearr)
    fluxsky=fluxsky3*skynorm   ; get the sky spectrum in wavearr grid, the unit should now be the same as fluxvega: 10^(-13) erg/s/nm/cm^2  (/arcsec^2 ?)
    
    fluxskypp=fluxsky

;a second way of estimating the Sky, if know the sky electron number per pixel
    if(keyword_set(skyperpixel)) then begin
    ;since the numbers are given by the main survey, our detected Sky electron will be less, so scale a rough factor of 0.9
        scaletemp=0.9
        ii=where(wavearr ge 255 and wavearr le 400, counta)
        fluxskypp[ii]=0.028/counta  
        ii=where(wavearr ge 400 and wavearr le 600, countb)
        fluxskypp[ii]=0.229/countb  
        ii=where(wavearr ge 600 and wavearr le 900, countc)
        fluxskypp[ii]=0.301/countc  
        ii=where(wavearr ge 900, countd)
        fluxskypp[ii]=0.301/countd 
        fluxskypp=fluxskypp/0.074^2*fov2*scaletemp
    end 
;********************************************
;define basic target brightness, parameters constantly change
    itarget=22 ; in Johnson V mag/arcsec^2 unit
    if(keyword_set(targetmag)) then itarget=targetmag
    itarget_jy=3631.0*10^(-itarget/2.5+3.0)  ; target flux in V in mJy/arcsec^2 unit
    itarget_nm=itarget_jy*3.0/(lambdav/100.0)^2 ;target flux in V in 10^(-13)erg/s/cm^2/nm (/arcsec^2 ?)

;readin the galaxy spectrum
    ; readcol, '../obs/allgalaxy.dat', wavegal, eflux, s0f, saf, sbf, scf, /silent
    ; wavegal=wavegal*1000.0
    ; spectype0=4	;default use Sb galaxy template spectrum
    ; if(keyword_set(spectype)) then spectype0=spectype  ;unless specified
    ; if(spectype0 eq 1) then galflux1=eflux
    ; if(spectype0 eq 2) then galflux1=s0f
    ; if(spectype0 eq 3) then galflux1=saf
    ; if(spectype0 eq 4) then galflux1=sbf
    ; if(spectype0 eq 5) then galflux1=scf

    tplfile='../obs/SFgal_tpl/SFgal_texp_FeH0_tau5_Ew10.fits'
    if(keyword_set(galtpl)) then tplfile='../obs/SFgal_tpl/'+galtpl
    sfgal=mrdfits(tplfile,1,hd)
    wavegal=sfgal.wave/10.  ;change A to nm
    galflux2=sfgal.flux   
    galflux1=interpol(galflux2,wavegal, wavearr)  
    
    ;normalize the galaxy spectrum to the V band magnitude specified. 
    ii=where(wavegal ge vmin and wavegal le vmax)
    wavetmp=wavegal(ii)
    fluxtmp=galflux2(ii)
    x=interpol(fluxfilter, wavefilter, wavetmp)
    vfluxtmp=x*wavetmp*fluxtmp  ;bandpass*lambda*F_gal_lambda
    galintegrate=integral(wavetmp, vfluxtmp)
    galnorm=itarget_nm*integconst/galintegrate
    galflux=galnorm*galflux1   ; the unit should now be in 10^(-13)erg/s/nm/cm^2 (/arcsec^2 ?)

;*********************************
;define observation information, parameters constantly change
    obst=300.0  ; in seconds, single integration time
    if(keyword_set(obstime)) then obst=obstime
    repn=1.0   ; repeating time
    if(keyword_set(repeatnum)) then repn=repeatnum
    npixw=3.0
    if(keyword_set(npixel_width)) then npixw=npixel_width
;   sky of slit area (slitwidth*npixw*slitlength) will go into the CCD
    if(keyword_set(slitwidth)) then fluxskypp=fluxskypp*slitwidth*npixw*slitunit  
    
    expf2=fltarr(narray)
    expfemi=fltarr(narray)
    snarray=fltarr(narray)
    mockgal=fltarr(narray)
    tmp=fltarr(narray)

    for i=0L, long(narray-1) do begin
        lambda0=wavearr[i]
        qlambda=interpol(q,lambdaq,lambda0)
        hv=planckh*cc/lambda0    ;10^{-10}erg
        delta_hz=cc*delta_lambda/lambda0/lambda0 * sampling ;10^17 1/s
    
        ;now that many fluxes are in 10^(-13)erg/s/nm/cm^2, to convert it to Jy, need to multiple: 
        ;lambda0^2/c(in nm)=lambda0^2(nm)/(3.*10^(17))*10^(-13)erg/s/Hz/cm^2
        ;=lambda^2(nm)*3.33*10^(-31)erg/s/Hz/cm^2=lambda^2(nm)*3.33*10^(-8)Jy
    
        ;find out sky value at lambda0    
        ;calculate n_sky/pixel
    	isky=fluxsky[i]*lambda0^2*0.0333*fov2   ;in uJy/spaxel unit
    	iskyall=isky*telarea/1000.0   ;in 10-26 erg/s/Hz /spaxel
    	fsky=qlambda*iskyall*delta_hz   ;10^{-9} erg/s /spaxel
    	nsky=fsky/hv*10.0   ;in unit of #e/s /spaxel
            if(keyword_set(skyperpixel)) then begin
                nsky=fluxskypp[i]*sampling  ; #e/s in npixw*sampling pixels 
            endif
    	;print, "Sky electron counts", nsky, nsky0, fluxskypp[i]

	    ;calculate n_source/pixel
	    isource=galflux[i]*lambda0^2*0.0333*fov2   ;in uJy/spaxel unit
	    isall=isource*telarea/1000.0   ;in 10-26 erg/s/Hz /spaxel
	    fs=qlambda*isall*delta_hz   ;10^{-9} erg/s /spaxel
	    ns=fs/hv*10.0   ;in unit of #e/s /spaxel
	    ;print, "Source electron counts", ns
    
	    darkn=(darkc*repn*obst*npixw*sampling)
	    rnn2=rn^2*(repn*npixw*sampling)
	    sourcenn=(ns*repn*obst)
	    skynn=(nsky*repn*obst)
	    tmp[i]=skynn

	    ;;nn1=sqrt(2.0*rnn^2+2.0*darkn^2+sourcenn^2+2.0*skynn^2)
	    ;nn1=sqrt(rnn^2+darkn^2+sourcenn^2+skynn^2)
	    nn1=sqrt(rnn2+darkn+skynn+sourcenn)  ;total noise
	    sn1=repn*ns*obst/nn1  ;S/N
	    snarray[i]=sn1
	    ;nn=sqrt(2.0*rnn^2+2.0*darkn^2+2.0*skynn^2)
	    nn=sqrt(rnn2+darkn+skynn)  ;system noise
	    ;print, "total noise, system noise, sn, readnoise, dark, source, sky", nn1, nn, sn1, rnn, darkn, sourcenn, skynn 
    
	    ;set the detection limit 
	    detlimit=1.0
	    if(keyword_set(snlimit)) then detlimit=snlimit
	    ;N_{source}/sqrt(N_{source}+nn^2)=detlimit, ==m
	    ;N_{source}^2-m^2*N_{source}-m^2*nn^2=0 (solve the equation)
	    ;N_{source}=(m^2+sqrt(m^4+4m^2*nn^2))/2.0
	    nntmp=detlimit^2+sqrt(detlimit^4+4.0*detlimit^2*nn^2)
	    nntmp=nntmp/2.0

    	;calculate detection limit in uJy and mag
    
        fnn=nntmp
        f1=fnn/obst/repn    ;in e/s
        f2=f1*hv       ;in 10^{-10} erg/s
        f3=f2/delta_lambda ;in 10^{-10} erg/s/nm
        f1=f3/telarea ;in 10^{-10} erg/s/nm/cm^2
        f2=f1/qlambda    ;in 10^{-10} erg/s/nm/cm^2
        expf2[i]=f2/fov2*100000.   ; in 10^{-15} erg/s/nm/cm^2/arcsec^2
        ;expfemi[i]=expf2[i]*delta_lambda*sampling   ; in 10^{-15} erg/s/cm^2/arcsec^2, the above multiplied by the spectral resolution
        ;print, "detection limit is", f2,"microJy/arcsec^2"
        ;print, "detection limit is", magf2, "AB mag / arcsec^2"
        mockgal[i]=galflux[i]+galflux[i]/snarray[i]*randomn(seed,1)  ;in 10^{-13} erg/s/nm/cm^2

        printf, uu, F='(F12.3,2E12.3,7F12.5,2E12.3)',lambda0,sn1,galflux[i],nn1,sqrt(sourcenn), nn, $ 
             sqrt(rnn2),sqrt(darkn), sqrt(skynn), mockgal[i]
    endfor


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
    

    ii=where(wavearr ge FWHMmin and wavearr le FWHMmax)
    wavetmp=wavearr(ii)
    snm=median(snarray[ii])    ;the median SN of FWHM range to acchieve the sn limit
    im=where(snarray[ii] eq snm)
    fact=reform(expf2[ii[im]]*0.01/galflux[[ii[im]]])
    fact=fact[0]
    limitmag=-2.5*alog10(fact)+itarget  

    ;print,limitmag
    ;oplot,wavetmp, fluxtmp, color=250
    snmean=median(snarray)
    
    z=0.0
    if(keyword_set(redshift)) then z=redshift
    waveha=656.3*(1.0+z)
    ii=where(wavearr ge (waveha-0.5) and wavearr le (waveha+0.5),nii)  ;1nm ,10A
    icont=where((wavearr ge (waveha-10) and wavearr le (waveha-5)) or (wavearr le (waveha+10) and wavearr ge (waveha+5)))
    contrms=sqrt(total(mockgal(icont)^2)/n_elements(mockgal(icont)))
    h=3*contrms*sqrt(nii) ; hight>3 con
    w=1.  ; width of 10A
    limitemif=sqrt(2*!PI)*h*w  ;h=3*cont, w=10A 
    ;print,limtmag,limitemif
    
free_lun, uu

 ; device,/close
 ; set_plot,mydevice
print,'The END!'
 
end
