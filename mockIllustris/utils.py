from platform import python_version

print(python_version())

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
#import sedpy
from sedpy import observate
from astropy.cosmology import WMAP9, FlatLambdaCDM, Planck13, Planck15, Planck18, z_at_value
import astropy.units as u
import astropy.constants as con
import os
import scipy
from scipy import ndimage
from astropy.convolution import convolve
#import webbpsf
from astropy.wcs import WCS
from astropy import units as u
from astropy.nddata import Cutout2D
import math
import pandas as pd


os.environ['WEBBPSF_PATH']='/Users/rebeccanevin/Documents/CfA_Code/mockIllustris/webbpsf-data'


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx
# Plot all of the different viewpoints

def plot_all_angles(file):
    # file is output'+str(snap)+'_'+str(id)+'
    fig = plt.figure(figsize = (15,10))
    ax = fig.add_subplot(231)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img1_total.fits')#41130
    ax.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    ax1 = fig.add_subplot(232)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img2_total.fits')#41130
    ax1.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    ax2 = fig.add_subplot(233)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img3_total.fits')#41130
    ax2.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    ax3 = fig.add_subplot(234)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img4_total.fits')#41130
    ax3.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    ax4 = fig.add_subplot(235)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img5_total.fits')#41130
    ax4.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    ax5 = fig.add_subplot(236)
    im_plot = pyfits.open('SKIRT_files/'+file+'/dusty_img6_total.fits')#41130
    ax5.imshow(np.sum(im_plot[0].data, axis=0), norm=matplotlib.colors.LogNorm())

    plt.show()

def combine(wav_rest, im, z, redshift_observation, cen_x, cen_y, cut_arc, filter_id, filter_name,  brighten, bg_limiting_m,size_detector_pix_arc, l_eff, verbose):
    
    
    
    
    # Convert to Jy from MJy/sr
    Jy = im[0].data*1e6
    size = np.shape(Jy)[1]
    size_pix_arc_obs = im[0].header['CDELT1']
    kpc_arcmin_obs=Planck15.kpc_proper_per_arcmin(redshift_observation)#insert the redshift to get the kpc/arcmin scaling                      
    
    if verbose:
        plt.clf()
        plt.imshow(abs(np.sum(Jy, axis=0)),  cmap = 'magma', norm=matplotlib.colors.LogNorm())#vmin=1e2,vmax=1e6))#, 
                   #norm = matplotlib.colors.LogNorm(vmax = np.max(img_cut[idx_lambda,:,:])/2))
        #plt.scatter(size/2, size/2, color='red')
        plt.colorbar( label = 'MJy/sr')
        locs = [0,size/2,size]
        labels = [round((size/2 - x)*size_pix_arc_obs,1) for x in locs]
        labels_x = [-round((size/2 - x)*size_pix_arc_obs*(kpc_arcmin_obs.value/60),1) for x in locs]

        plt.xticks(locs,labels_x)
        plt.yticks(locs,labels)
        plt.title('Observed Broadband Image')
        plt.xlabel('kpc')
        plt.ylabel('Arcsec')
        plt.show()

        print('size of pixels at 10pMpc', size_pix_arc_obs, 'arcsec', size_pix_arc_obs*(kpc_arcmin_obs.value/60), 'kpc')

    # Figure out 
    kpc_arcsec_at_z = Planck15.kpc_proper_per_arcmin(z).value/60 #kpc/arcsec
    print('kpc/arcsec at z', kpc_arcsec_at_z)
    # Okay so we know how many kpc it is across, 
    # lets figure out how many arcsec that is at this redshift
    
    kpc_across = size*size_pix_arc_obs*(kpc_arcmin_obs.value/60)
    arc_across = kpc_across / kpc_arcsec_at_z
    size_pix_arc = arc_across / size
    
    print('size pic in arc at z', size_pix_arc)
    
    
    # Now, time to start moving the observation to higher redshift by first accounting 
    # for the physical size of the pixels at the 10 p Mpc distance
    '''
    STOP
    
    
    
    
    
    
    
    # I think I actually should get rid of steradians RIGHT AWAY
    Flux_Jy = (im[0].data*1e6) #Jy # also used to be (pixelscale_pc/1e6)**2) [to get rid of SR]
    
    
    
    
    pixelscale_kpc = size_pix_arc * (kpc_arcmin.value/60)
    print('size of a pixel in kpc', pixelscale_kpc)
    print('size of a pixel in arc', size_pix_arc)
    
    
    size_pix = int(cut_arc/size_pix_arc)
    
    
    if size_pix > np.shape(Flux_Jy)[1]/2: # this would happen if the image isn't big enough
        size_pix = int(np.shape(Flux_Jy)[1]/2)
        cut_arc = size_pix*size_pix_arc
    # so this is cutting it up to an arcsec size
    Jy_cut = Flux_Jy[:,cen_x-size_pix:cen_x+size_pix,cen_y-size_pix:cen_y+size_pix]
    

    idx_lambda = find_nearest(wav_rest, l_eff)[1]
    
    
    
    

    
    
    '''
    
    # Specific intensity (I_nu) is W/m^2/sr/Hz
    # is invariant to distance
    
    # Code from Ben
    from astropy import units as u
    from astropy.constants import c
    
    # Convert back to a surface brightness
    B_nu_rest = Jy#/size_pix_arc**2

    pixel_size = 1 * u.sr #size_pix_arc * u.arcsec # replace with output pixel size, or with 1 * u.sr if you want to keep output units as sr^{-1}

    # Add units
    B_nu_rest *= u.Jy / u.sr#(u.arcsec*u.arcsec)#**2
    wav_rest *= u.micron
    
    if verbose:
        plt.clf()
        plt.plot(wav_rest, np.sum(B_nu_rest, axis=(1,2)))
        plt.title('Is this the rest wavelength?')
        plt.axvline(x = 656.6*10**(-3))
        #plt.xlim([0.5,1.5])
        plt.show()

        print(np.shape(B_nu_rest), np.shape(wav_rest))
        print('this would be the factor you need to brighten by to get intrinsic surface brightness')
        print((1+redshift_observation)**4)
    

    
    
    
    
    '''cut_lambda = 3*(1+redshift)#3
    cut_idx_wav = find_nearest(wave_obs.value, cut_lambda)[1] # This gives you the index where you need to cut
    
    cut_wav = wave_obs[1:cut_idx_wav+1]
    cut_img_cut = I_nu_obs[:cut_idx_wav,:,:]
    
    cut_lambda = 3#*(1+redshift)#3
    cut_idx_rest = find_nearest(wav_rest[1:].value, cut_lambda)[1]
    cut_wav_rest = wav_rest[1:cut_idx_rest+1]
    cut_img_rest = B_nu_rest[:cut_idx_rest,:,:]
    
    # While you're at it just plot the f*ing SED
    plt.clf()
    plt.plot(cut_wav, np.sum(cut_img_cut, axis=(1,2)), label='Observed')
    #plt.plot(cut_wav_rest, np.sum(cut_img_rest, axis=(1,2)), label='Rest frame, pre dimmed')
    plt.title('SED + Filters')
    plt.xlabel('Wavelength, Microns')
    
    # go and grab the filter to overplot plz
    ylim = 2*np.max(np.sum(cut_img_cut.value, axis=(1,2)))
    plot_filter(filter_id, filter_name, filtercolor, ylim)
    ##index_filter = find_nearest(cut_wav, l_eff)[1]
    #plt.axvline(x=l_eff)
    plt.legend()
    plt.ylabel('Flux [MJy/sr]')
    plt.ylim([0, ylim])
    #plt.xlim([3000,10000])
    plt.xlim([0.7, 5])
    plt.show()'''
    
    #if verbose:
        #kpc_arcmin_obs=cosmo.kpc_proper_per_arcmin(redshift_observation)
        #plot_SED(redshift, wave_obs, I_nu_obs, size*size_pix_arc)
        
    
   
    
    
    '''
    plt.clf()
    plt.plot(wave_obs, np.sum(I_nu_obs, axis=(1,2))/size**2)
    plt.title('Average Spectrum for the FOV')
    plt.show()

    
    
    # First, do it all at once, summed up
    f_nu_obs_in_one_pixel = np.sum(I_nu_obs, axis=(1,2))  #* (pixel_size**2)

    

    # unit conversions, switch to f_lambda, and strip the astropy units again
    f_lambda_obs_per_area = f_nu_obs_in_one_pixel.to(u.erg/u.s/u.Hz /u.cm**2 / u.sr ).value * c.to(u.AA/u.s).value / wave_obs.to(u.AA).value**2
    #/u.arcsec**2

    # filter projection
    mags_per_area = observate.getSED(wave_obs.to(u.AA).value, f_lambda_obs_per_area, filterlist=bandlist)
    plt.clf()
    plt.plot(wave_obs, f_nu_obs_in_one_pixel)
    plt.title('Summed Spectrum for the FOV mag = '+str(mags_per_area))
    plt.show()
    '''
    
    # Okay split it up
    # Do the redshifting
    a = (1+z)
    wave_obs = wav_rest * a
    # wave_obs = wav_obs[1:]
    # Source for cosmological surface brightness dimming: 
    # https://iopscience.iop.org/article/10.1088/0004-637X/796/2/102/meta
    # B_nu is the surface brightness and I_nu is the specific intensity 
    # (which has been corrected for cosmo dimming)
    I_nu_obs = B_nu_rest / (a/(1+redshift_observation))**3# 
    
    
    '''
    # Compare this to Jacob's method - he uses 1/D_L**2
    #Flux = fnu*(1/D_L)**2*extinction*(1+z)* (pixel*pixel in pc)/(1e6)**2 * 1e6 #units are Jy
    luminosity_distance = Planck15.luminosity_distance(z).value
    print('luminosity distance', luminosity_distance, 'Mpc')
    size_pix_pc = 60.081283552314665
    # Multiply surface brightness by omega (unit solid angle) = A/D_A**2 to get flux
    # D_L = (1+z)^2 D_A
    dim_factor_Jacob = a*(10e6/luminosity_distance*1e6)**2 * (size_pix_pc**2)/(10e6)**2
    print('redshift dimming factor Jacob', dim_factor_Jacob)
    
    print('what is 1+z', a)
    '''
    
    
    if verbose:
        plt.clf()
        plt.plot(wave_obs, np.sum(I_nu_obs, axis=(1,2)))
        plt.title('Observed wavelength')
        plt.axvline(x = 656.6*10**(-3))

        plt.show()

        plt.clf()
        plt.imshow(np.sum(I_nu_obs, axis=0), norm=matplotlib.colors.LogNorm())
        plt.colorbar()
        plt.title('Dimmed SB (1+z)^4')
        plt.show()
    
    
    bandlist = observate.load_filters(filter_name)
    


    
    img_mag = np.zeros((np.shape(I_nu_obs)[1], np.shape(I_nu_obs)[2]))
    

    
    
    for i in range(np.shape(I_nu_obs)[1]):
        for j in range(np.shape(I_nu_obs)[2]):
            # the output of this is the value of the pixel in AB mags/arcsec^2
            
            # convert to flux density per pixel (or per sr if pixel_size = 1 * u.sr)
            f_nu_obs_in_one_pixel = I_nu_obs[:,i,j]  #* (pixel_size**2)


            # unit conversions, switch to f_lambda, and strip the astropy units again
            f_lambda_obs_per_area = f_nu_obs_in_one_pixel.to(u.erg/u.s/u.Hz /u.cm**2 / u.sr  ).value * c.to(u.AA/u.s).value / wave_obs.to(u.AA).value**2
            #/u.arcsec**2
            
            # filter projection
            img_mag[i,j] = observate.getSED(wave_obs.to(u.AA).value, f_lambda_obs_per_area, filterlist=bandlist)


            #img_mag[i,j] = mags_per_area
            
    
    # Redshift dimming applied now
    
    Jy_per_sr = 3631 * 10**(-0.4 * img_mag) #/ (a/(1+redshift_observation))**4

   

    # We are assuming the PSF has the same pixelscale as the JWST image for that filter
    # (Well actually it gives you a couple of extensions, one of them has the same sampling (ext 1))
    '''
    nc = webbpsf.NIRCam()
    nc.filter=filter_id
    psf = nc.calc_psf(fov_arcsec=cut_arc) # was 2 arcsec
    '''
    # From this site, a good approx of the PSF is a 0.16 FWHM (@ max)
    # https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-predicted-performance/nircam-point-spread-functions
    
    # Only works for wide filters
    dict_filter_fwhm = {'070':0.03,'090':0.034,'115':0.04,'150':0.05,'200':0.066,'277':0.091,'356':0.115,'444':0.145}
    lambda_filter = str.split(str.split(filter_id,'F')[1],'W')[0]
    fwhm = dict_filter_fwhm[lambda_filter] # in arcsec
        
    
    size_pix = size_detector_pix_arc #arcsec
    
    # Let's make the thing about 10 times as large as the FWHM
    # size = 10*0.16 = 1.6" across
    # how many pixels is this?
    factor = 10
    size_across_arc = factor*fwhm
    size_across_pix = int(size_across_arc / size_pix)
    
    x = np.linspace(0,size_across_pix - 1,size_across_pix)
    y = np.linspace(0,size_across_pix - 1,size_across_pix)
    
    def twod_gauss(x, y, mu_x, mu_y, sigma):
        # How the hell is this not a single number?
        value = np.exp(- ( (x - mu_x)**2 + (y - mu_y)**2 ) / (2 * sigma**2))
        if value < 0.001:
            value ==0
        return value
        
    
    twod = np.zeros((size_across_pix, size_across_pix))
    for i in range(len(x)):
        for j in range(len(y)):
            twod[i, j] = twod_gauss(x[i], y[j], size_across_pix/2, size_across_pix/2, (fwhm/2.54)/size_pix)
    
    if verbose:
        plt.clf()
        plt.imshow(twod)
        plt.colorbar()
        plt.show()




        print('size before rebinning', size_pix_arc)
        print('this will be size after rebinning', size_detector_pix_arc)
    factor = size_pix_arc/size_detector_pix_arc
    # you have to use the difference between the current size of pixels (in the sim image)
    # and the size of pixels on the JWST detector
    rebin = scipy.ndimage.zoom(Jy_per_sr, factor, order=0)

    
    if verbose:
        plt.clf()
        plt.imshow(rebin, cmap='magma_r')#, vmin=min_mag, vmax=max_mag)
        plt.colorbar()
        plt.title('Rebinned Jy/arcsec2')
        size = np.shape(rebin)[1]
        locs = [0,size/2,size]
        labels = [-round((size/2 - x)*size_detector_pix_arc,1) for x in locs]
        labels_y = [-round((size/2 - x)*size_detector_pix_arc*kpc_arcsec_at_z,1) for x in locs]
        plt.xticks(locs,labels)
        plt.yticks(locs,labels_y)
        plt.xlabel('Arcsec')
        plt.ylabel('kpc')
        plt.show()

    '''if verbose:
        plt.clf()
        plt.imshow(rebin, cmap='magma_r')#, vmin=min_mag, vmax=max_mag)
        plt.colorbar()
        plt.title('Rebinned Jy/arcsec2')
        size = np.shape(rebin)[1]
        locs = [0,size/2,size]
        labels = [round((size/2 - x)*size_detector_pix_arc,1) for x in locs]
        plt.xticks(locs,labels)
        plt.yticks(locs,labels)
        plt.show()'''
        
        
    
    try:
        result = (convolve(rebin, twod))# psf[1].data
    except:
        # this is an error that occurs when the psf doesn't have odd dimensions
        
        
        result = convolve(rebin, twod[1:,1:])
        
        
    
    if verbose:
        plt.clf()
        plt.imshow(result, cmap='magma_r')
        plt.colorbar()
        plt.title('Convolved and Rebinned')
        plt.show()

    # 1 sr = 4.25 x 10^10 arcsec2
    nJy = ((result*10**9)/(4.25*10**10))*size_detector_pix_arc**2
    
    

    # the last step is to add background
    nJy_bg = background_resids(bg_limiting_m, size_detector_pix_arc, l_eff, rebin, verbose)

    
    # Alternately, just add in the JADES background + bg galaxies using Guitarra
    #nJy_bg = Guitarra_bg(cut_arc, filter_id)
   
    
    AB = convert_to_AB(nJy)
    AB_bg = convert_to_AB(nJy_bg)
    if verbose:
        plt.clf()
        plt.imshow(nJy_bg)
        plt.colorbar()
        plt.title('Just background nJy '+str(filter_name))
        plt.show()



        plt.clf()
        plt.imshow(nJy)
        plt.colorbar()
        plt.title('Just image nJy '+str(filter_name))
        plt.show()

        '''plt.clf()
        plt.imshow(AB)
        plt.colorbar()
        plt.title('Just image AB '+str(filter_name))
        plt.show()'''

        plt.clf()
        plt.hist(nJy_bg.flatten(), bins=100)
        plt.xlabel('')
        plt.show()
    
    plt.clf()
    plt.imshow(abs(nJy+nJy_bg), norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.title('nJy image '+str(filter_name))
    plt.show()
    
    '''plt.clf()
    plt.imshow(abs(nJy+nJy_bg), norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.title('nJy image log '+str(filter_name))
    plt.show()'''
    
    
    
    
    # output is image with added bg, then size of the pixels in arc, then size of the pixels in kpc 
    return nJy+nJy_bg, size_detector_pix_arc, size_detector_pix_arc*kpc_arcsec_at_z



def plot_filters(id_list, filtercolor, ymax_data):
    
    
    for j in range(len(id_list)):
        filtertxt = np.loadtxt('greg_JWST_filters/nircam_'+str(id_list[j]))

        # use the max of the data

        plt.fill_between(filtertxt[:,0]/1e4, 0, filtertxt[:,1]*(ymax_data/4), color=filtercolor[j], alpha=0.5)
        plt.plot(filtertxt[:,0]/1e4, filtertxt[:,1]*(ymax_data/4), color=filtercolor[j])
def plot_SED(redshift, wave_obs, I_nu_obs,arcsec_side ):
    cut_lambda = 3*(1+redshift)#3
    cut_idx_wav = find_nearest(wave_obs.value, cut_lambda)[1] # This gives you the index where you need to cut
    
    cut_wav = wave_obs[1:cut_idx_wav+1]
    cut_img_cut = I_nu_obs[:cut_idx_wav,:,:]
    
    # While you're at it just plot the f*ing SED
    plt.clf()
    fig = plt.figure(figsize = (8,5))
    
    #plt.plot(cut_wav_rest, np.sum(cut_img_rest, axis=(1,2)), label='Rest frame, pre dimmed')
    #plt.title('SED + Filters')
    plt.xlabel(r'Observed Wavelength ($\mu$m)')
    
    # go and grab the filter to overplot plz
    ylim = 2*np.max(np.sum(cut_img_cut.value, axis=(1,2)))
    
    filter_list = ['F070W','F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W']
    filtercolor = ['#0D69AB', '#6E99C9','#9FC3E9','#6874AC',
                   '#342B75','#6B327B','#923978','#DC9095','#E4ADC8']
    plot_filters(filter_list, filtercolor, ylim)
    ##index_filter = find_nearest(cut_wav, l_eff)[1]
    #plt.axvline(x=l_eff)
    #plt.legend()
    plt.ylabel('Integrated Surface Brightness [MJy/sr]')
    plt.step(cut_wav, np.sum(cut_img_cut, axis=(1,2)),  color='black')
    plt.ylim([0.01e10, ylim*0.7])
    #plt.annotate('Star-Forming Galaxy at $z = $'+str(redshift), xy=(0.02,0.9), xycoords='axes fraction', color='black')
    #plt.xlim([3000,10000])
    plt.xlim([0.7, 5])
    plt.show()
    
    
    
    plt.clf()
    fig = plt.figure(figsize = (8,5))
    #plot_filters(filter_list, filtercolor, ylim)
    plt.ylabel('Integrated Flux [nJy]')
    flux = ((((np.sum(cut_img_cut, axis=(1,2)))/4.25e10)*arcsec_side**2)/1e6)# in Jy
    
    ylim = 2*np.max(flux)
    plot_filters(filter_list, filtercolor, ylim)
    plt.step(cut_wav, 1e9*flux,  color='black')
    #plt.ylim([0.01e10, ylim*0.7])
    #plt.annotate('Star-Forming Galaxy at $z = $'+str(redshift), xy=(0.02,0.9), xycoords='axes fraction', color='black')
    #plt.xlim([3000,10000])
    plt.xlim([0.7, 5])
    plt.show()
    
    
    # Now convert to AB mags:
    plt.clf()
    fig = plt.figure(figsize = (8,5))
    #plot_filters(filter_list, filtercolor, ylim)
    plt.ylabel('Integrated Flux [AB mag]')
    
    AB = -2.5 * np.log10(np.array(flux.value) / 3631)
    
    ylim = 2*np.max(AB)
    plot_filters(filter_list, filtercolor, ylim)
    plt.step(cut_wav, AB,  color='black')
    plt.ylim([30,20])
    # compare to Christina's paper: https://iopscience.iop.org/article/10.3847/1538-4365/aabcbb/pdf
    
    #plt.annotate('Star-Forming Galaxy at $z = $'+str(redshift), xy=(0.02,0.9), xycoords='axes fraction', color='black')
    #plt.xlim([3000,10000])
    plt.xlim([0.7, 5])
    plt.show()
    
    #plt.savefig('mock_panels/SED.png', dpi=500)
    
    
def background_resids(bg_limiting_magnitude,size_detector_pix_arc, l_eff, rebin, verbose):
    # This section takes the limiting magnitude and converts to a surface brightness 
    # in units of 
    
    # Super useful ref from Paul Green: https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf
    
    bg_limiting_nJy = (3631*10**(bg_limiting_magnitude/(-2.5)))*1e9

    print('bg limit Jy', bg_limiting_nJy)

    

    bg_sigma = 0.2*bg_limiting_nJy
    
    
    

    background = np.zeros(np.shape(rebin))
    mu, sigma = 0, bg_sigma # mean and standard deviation

    for i in range(np.shape(background)[0]):
        for j in range(np.shape(background)[1]):
            s = np.random.normal(mu, sigma, 1)
            background[i,j] = s

    
    '''if verbose:
        plt.clf()
        plt.imshow(background, cmap='magma_r')
        plt.colorbar(label='erg/s/cm2/A')
        plt.title('Background')
        plt.show()

        plt.clf()
        plt.imshow(rebin+background, cmap='magma_r')
        plt.colorbar(label='erg/s/cm2/A')
        plt.title('Background + image')
        plt.show()
        
        
        # Now convert back into Jy and then AB mags
        # convert from AB mags to Janskies
        # Janskies = 10**(img_mag/(-2.5))*3631# the AB mag zeropt is defined to be 3631 Jy
    
        # Now convert back to erg/s/cm2/A
        # per_A = 2.99792458E-05 * Janskies / l_eff**2# l_eff is the effective wavelength

        Jy = l_eff**2*(rebin+background)/2.99792458E-05
        
        xcen = 45
        ycen = 35
        plt.clf()
        plt.imshow(Jy, cmap = 'magma_r')
        plt.colorbar(label='Jy/arcsec**2')
        plt.scatter(xcen, ycen, color='black', marker='x')
        plt.show()
        
        plt.clf()
        plt.imshow(Jy[xcen-3:xcen+4,ycen-3:ycen+4], cmap = 'magma_r')
        plt.colorbar(label='Jy/arcsec**2')
        #plt.scatter(45,35, color='black', marker='x')
        plt.show()
        
        ap_sum_Jy = np.sum(Jy[xcen-3:xcen+4,ycen-3:ycen+4])
        ap_avg_Jy = np.mean(Jy[xcen-3:xcen+4,ycen-3:ycen+4])
        print('Sum Jy of aperture', ap_sum_Jy)
        print('Avg Jy of aperture', ap_avg_Jy)
        
        
        img_mag = -2.5*np.log10(abs(Jy/3631))
        bg_mag = -2.5*np.log10(abs(l_eff**2*(background)/2.99792458E-05)/3631)
        
        print('mag of this aperture', -2.5*np.log10(ap_sum_Jy/3631))
        print('mag of this aperture', -2.5*np.log10(ap_avg_Jy/3631))
        
        plt.clf()
        plt.imshow(bg_mag, cmap = 'magma_r')
        plt.colorbar(label='mag/arcsec**2')
        plt.title('Background in mag/arcsec**2')
        plt.show()
        
        plt.clf()
        plt.imshow(img_mag, cmap = 'magma_r')
        plt.colorbar(label='mag/arcsec**2')
        plt.show()
        
        # Now calculate just the average bg mag - 
        Jy_bg_mean = l_eff**2*(bg_limiting_flambda)/2.99792458E-05
        print('mag of bg sum',-2.5*np.log10(Jy_bg_mean*int(np.sqrt(0.2/size_detector_pix_arc**2))/3631))
        print('mag of bg mean',-2.5*np.log10(Jy_bg_mean/3631))'''
        
        
    return background
def convert_to_nJy(convolved, pix_arc, wavelength):
    #input is a surface brightness units are [erg/s/cm2/A/arcsec2]
    
    erg_s_cm2_A = convolved*pix_arc**2 # convert from a surface brightness to a specific flux?
    
    # Then from here ( https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf )
    erg_s_cm2_Hz = erg_s_cm2_A*wavelength**2*3.34e-19
    Jy = erg_s_cm2_Hz/1e-23
    return Jy*1e9

def convert_to_AB(nJy):
    # Convert to AB mags :)
    #bg_limiting_nJy = (3631*10**(bg_limiting_magnitude/(-2.5)))*1e9
    
    return -2.5*np.log10((nJy/1e9)/3631)

