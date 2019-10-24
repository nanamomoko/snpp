
# V1.03alpha
# V1.03beta
#v1.04beta: 17 Oct 2019 by Mengting Ju

##################################################################
from __future__ import print_function
from astropy.io import fits
import numpy as np
from snpp import snpp
import inwf
import inmg

#################################################################

def snpp_example():
    
    ######################################
    
    
    #select model and magnitude
    targetmag=19.
    galtpl='SFgal_texp_FeH0_tau5_Ew10.fits'
    deltal=1.755555
    filtera='sdss_g0.par'
    sampling=2.0
    
    result=inmg.input_mag_model(targetmag,galtpl,deltal, filtera,sampling)
    wavearr=result[0]   #A
    galflux=result[1]   #10^-12 erg/s/A/cm2
    integconst=result[2]
    delta_lambda=np.array([deltal]*len(wavearr))    
    '''
    ####################################
    
    #select put in wave and flux    
    filee=fits.open('MockGal-M21Z0.01-FOV6R0.1-W350n1000n.fits')
    fluxx=np.sum(np.sum(filee[1].data[:20,:20],1),0)  #erg/s/A/cm2
    wavee=filee[2].data    #A
    redshift=0.01
    filtera='sdss_r0.par'
    sampling=2.0
    
    result=inwf.input_wave_flux(wavee,fluxx,redshift,filtera)
    wavearr=result[0]  #A
    galflux=result[1]  #10^-12 erg/s/A/cm2
    integconst=result[2]
    delta_lambda=result[3]
    
    '''
    #####################################
    #filename='test.fits'
    filename='mg_weak_30001_24_1_19.fits'
    ss=snpp(wavearr=wavearr,galflux=galflux,integconst=integconst,
            deltal=delta_lambda,filename=filename,
            specsample=sampling,filtera=filtera,
            readnoise=4.0,fovp=0.2,npixel_width=2.0,
            obstime=300,repeatnum=1,skyv=22.5,qinput=1.0, 
            slitwidth=None,skyperpixel=True,
            snlimit=1.0,targetmaglimit=19.)

#########################################################

if __name__=='__main__':
    snpp_example()
    
