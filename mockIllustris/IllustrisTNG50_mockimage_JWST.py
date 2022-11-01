# Nov 1 2022 - updating to be a .py file
### Goal: Create a mock image from a TNG50 SKIRT datacube
#  This notebook will walk through the steps to go from the SKIRT output cube into a reasonable mock of the user's choice of filter (from JWST or HST). 

#  The user can also set the redshift of the image and the size of the image cutout.

#  Things this notebook does:
#  1. Imports SED from SKIRT and describes the file (from Xuejian Shen).
#  2. Applies cosmological dimming and IGM absorption (from Xuejian Shen).
#  3. Uses SEDPY to create an image of a given filter.
#  4. Applies an appropriate PSF.
#  5. Rebins to the pixelscale of the instrument.
#  6. Introduces appropriate background residual.

#  Things I haven't done yet: 
#  1. Centering software to center on the galaxy of choice.
#  2. Determine (if we need) and implement background galaxies.
#  3. Potentially create an error image.

#  Things that need to be double checked:
#  1. The cosmological dimming.
#  2. Adding the background residuals in the appropriate units.
import utils
import os
import astropy
from astropy.cosmology import WMAP9, FlatLambdaCDM, Planck13, Planck15, Planck18, z_at_value
from astropy import units as u


os.environ['WEBBPSF_PATH']='/Users/rebeccanevin/Documents/CfA_Code/mockIllustris/webbpsf-data'


verbose=True

# Make a dictionary for redshifts and snapshots
# So you don't need to look this up all the time
redshift_snap = {84:0.2,50:1,33:2}
size_snap = {84:50.11341574807246798, 50:30.04064159676324562}
snap = 50
redshift = redshift_snap[snap]

# The SKIRT camera is placed at a distance of 10 pMpc, so first use this to calculate the physical
# size of the pixels



# Okay but whats the difference between:
# angular diameter distance (d_A = x/theta, where x is the physical size and theta is angular size)
# 29.357047613319295
# comoving distance (29.290825794667786)
# luminosity distance (Full FOV in kpc 29.22505026845403)

# For small distances, z = v/c = d/D_=0.6774
# Hubble distance D_H = 9.26e25 1/h m, so in our case 1.367 e 26 m
# so in our case z = 10 Mpc / 1.367e26 m 
# z = 0.002257

redshift_observation = z_at_value(Planck15.angular_diameter_distance,10 * u.Mpc, zmax=1.5)#Planck15.angular_diameter_distance
#print('redshift according to comoving', redshift_observation)

#redshift_observation = 0.002257
#print(redshift_observation)

kpc_arcmin_camera = Planck15.kpc_proper_per_arcmin(redshift_observation)#insert the redshift to get the kpc/arcmin scaling                      

print(kpc_arcmin_camera)


id = 41130
merger = True

im = astropy.io.fits.open('SKIRT_files/output'+str(snap)+'_'+str(id)+'/dusty_img3_total.fits')#41130
hdr = im[0].header

FOV_arc = hdr['CDELT1']*np.shape(im[0].data)[1]
FOV_kpc = size_snap[snap]
print('Full FOV in arcsec', FOV_arc)
print('Full FOV in kpc', FOV_arc*(kpc_arcmin_camera.value/60))
print('FOV according to SKIRT', FOV_kpc)

#theta = d/D
redshift_based_on_size = (FOV_kpc/np.shape(im[0].data)[1])/(10*100)
print('redshift based on size', redshift_based_on_size)

print('size of a single pixel in kpc', (hdr['CDELT1'])*(kpc_arcmin_camera.value/60))

#print(repr(hdr))
wav_array = im[1].data
wav = []
for j in range(len(wav_array)):
    wav.append(wav_array[j][0])
wav = np.array(wav)

wav_rest = wav / (1+redshift_observation)

if merger:
    table = pd.read_table('Tables/mergers_all_nodup_'+str(snap)+'_L35n2160TNG.txt')
    info = table[table['Subfind_ID']==id]
else:
    table = pd.read_table('Tables/nonmergers_snap_'+str(snap)+'_L35n2160TNG.txt')
    info = table[table['Subfind_ID']==id]
    
plt.clf()
plt.imshow(im[0].data[10,:,:], norm=matplotlib.colors.LogNorm())
plt.annotate(str(info), xy=(0,50), xycoords='data', color='white', size=8)
plt.title('Plot at a single wavelength')
plt.show()

# Run multiple different galaxies through and plot them

verbose=False

brighten = 1


cut_arc = 100#12.6

# I'm just goint to center on the snapshot for now
cenx =int(np.shape(im[0].data)[1]/2)
ceny = int(np.shape(im[0].data)[2]/2)

print('redshift', redshift)


cut_arc = 200#







# Run through multiple galaxies

holder_images = np.zeros((len(os.listdir('SKIRT_files/')), np.shape(im[0].data)[2], np.shape(im[0].data)[2]))
info_list = []

img_f444w_list = []
img_f356w_list = []
img_f277w_list = []
img_f200w_list = []
img_f115w_list = []
img_f090w_list = []

pix_f444w_arc_list = []
pix_f444w_kpc_list = []

pix_f356w_arc_list = []
pix_f356w_kpc_list = []

pix_f277w_arc_list = []
pix_f277w_kpc_list = []

pix_f200w_arc_list = []
pix_f200w_kpc_list = []

pix_f115w_arc_list = []
pix_f115w_kpc_list = []

pix_f090w_arc_list = []
pix_f090w_kpc_list = []

redshift_list = []

counter = 0
for file in os.listdir('SKIRT_files/'):
    print(file)
    # Okay figure out the id and the snap and the redshift
    out_1 = file.split('_')[0]
    snap = int(float(out_1.split('output')[1]))
    id = float(file.split('_')[1])
    redshift = redshift_snap[snap]
    print('snap', snap, 'id', id, 'redshift', redshift)
    
    # Get info about this galaxy
    if merger:
        table = pd.read_table('Tables/mergers_all_nodup_'+str(snap)+'_L35n2160TNG.txt')
        info = table[table['Subfind_ID']==id]
    else:
        table = pd.read_table('Tables/nonmergers_snap_'+str(snap)+'_L35n2160TNG.txt')
        info = table[table['Subfind_ID']==id]
    print('info', info)
    info_list.append(info)
    try:
        plot_all_angles(file)
    except:# FileNotFoundError:
        print('doesnt have the .fits files yet')
        continue
    
    im = pyfits.open('SKIRT_files/'+file+'/dusty_img1_total.fits')#41130
    hdr = im[0].header

    FOV_arc = hdr['CDELT1']*np.shape(im[0].data)[1]
    

    #print(repr(hdr))
    wav_array = im[1].data
    wav = []
    for j in range(len(wav_array)):
        wav.append(wav_array[j][0])
    wav = np.array(wav)

    wav_rest = wav / (1+redshift_observation)
    
    
    filter_id = 'F444W'
    size_detector_pix_arc = 0.031
    l_eff = 0.814#e4
    # need to convert from microns to A
    bg_limiting_sb_flambda = 7.57e-18 # this is flambda/arcsec^2
    bg_limiting_mag = 29.8#29.8# from the CEERs white paper 
    filter_name = ['acs_wfc_f814w']
    
    mock_f444w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                     brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f444w_list.append(mock_f444w)
    pix_f444w_arc_list.append(pix_arc)
    pix_f444w_kpc_list.append(pix_kpc)
    
    
    
    filter_id = 'F444W'
    size_detector_pix_arc = 0.063
    l_eff = 4.44#e4
    # need to convert from microns to A
    bg_limiting_sb_flambda = 7.57e-18 # this is flambda/arcsec^2
    bg_limiting_mag = 29.8#29.8# from the CEERs white paper 
    filter_name = ['jwst_f444w']
    
    mock_f444w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                     brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f444w_list.append(mock_f444w)
    pix_f444w_arc_list.append(pix_arc)
    pix_f444w_kpc_list.append(pix_kpc)
    
    
    
    filter_id = 'F356W'
    size_detector_pix_arc = 0.063
    l_eff = 3.56#e4
    bg_limiting_mag = 30.1# from the CEERs white paper 
    filter_name = ['jwst_f356w']
    mock_f356w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                         brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f356w_list.append(mock_f356w)
    pix_f356w_arc_list.append(pix_arc)
    pix_f356w_kpc_list.append(pix_kpc)
    
    filter_id = 'F277W'
    size_detector_pix_arc = 0.063
    l_eff = 2.77#e4
    bg_limiting_mag = 30.2# from the CEERs white paper 
    filter_name = ['jwst_f277w']
    mock_f277w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                         brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    
    img_f277w_list.append(mock_f277w)
    pix_f277w_arc_list.append(pix_arc)
    pix_f277w_kpc_list.append(pix_kpc)
    
    filter_id = 'F200W'
    size_detector_pix_arc = 0.031#0.031
    l_eff = 2.00#e4
    bg_limiting_mag = 30.6#29.7#30.6
    filter_name = ['jwst_f200w']
    mock_f200w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                         brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f200w_list.append(mock_f200w)
    pix_f200w_arc_list.append(pix_arc)
    pix_f200w_kpc_list.append(pix_kpc)
    
    filter_id = 'F115W'
    l_eff = 1.15#e4#70e4
    size_detector_pix_arc = 0.031#0.063#0.031
    bg_limiting_mag = 30.5#30.5#29.8#28.7 Unclear if this is /arcsec^2 or not
    filter_name = ['jwst_f115w']
    mock_f115w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                         brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f115w_list.append(mock_f115w)
    pix_f115w_arc_list.append(pix_arc)
    pix_f115w_kpc_list.append(pix_kpc)
    
    filter_id = 'F090W'
    size_detector_pix_arc = 0.031#0.031
    l_eff = 0.9#e4
    bg_limiting_mag = 30.2
    filter_name = ['jwst_f090w']
    mock_f090w, pix_arc, pix_kpc = combine(wav_rest, im, redshift, redshift_observation, cenx, ceny, cut_arc, filter_id, filter_name, 
                         brighten, bg_limiting_mag,size_detector_pix_arc, l_eff, verbose)
    img_f090w_list.append(mock_f090w)
    pix_f090w_arc_list.append(pix_arc)
    pix_f090w_kpc_list.append(pix_kpc)
    
    redshift_list.append(redshift)
    
    counter+=1
    

STOP









'''filter_id = 'F090W'
size_detector_pix_arc = 0.031#0.031
l_eff = 0.9#e4
bg_limiting_mag = 30.2
filter_name = ['jwst_f090w']
mock_f090w = combine(wav, im5, redshift, cenx, ceny, cut_arc, filter_id, filter_name, 
                     brighten, bg_limiting_mag,
              size_detector_pix_arc, l_eff, verbose)'''







 
