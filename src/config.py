import numpy as np

# Path definitions

mockdir = '/share/des/sv/BCC-SVA1-WL-v3.0/'
golddir = '/share/des/disc2/y1/gold_v101/'
pzdir = '/share/des/sv/photoz/DES_PDF_Stacker/'
pztestdir = '/home/troxel/cosmosis/cosmosis-des-library/photoztests/y1/'
cosmosisrootdir = '/home/troxel/cosmosis/cosmosis-des-library/wl/y1prep/'
wcsfile = '/share/des/disc2/y1/y1a1_image_wcs_info.txt'
spointsfile = 'y1a1_special_field_points.fits'
pointingfile = '/home/troxel/catcutDES/y1a1_telradec2.txt'
y1sysmapdir = '/share/des/disc2/y1/sysmaps/'
svsysmapdir = '/share/des/sv/systematics_maps/'
redmagicdir = '/share/des/sv/redmagicv6.3.3/'

# Dictionaries

cfg = {
  
  'lbins':10,
  'sbins':2,
  'slop':0.1,
  'tbins':8,
  'cbins':5,
  'sep':np.array([1.,400.]),
  'num_patch':126,
  'num_reg':150,
  'bs':False,
  'wt':False,
  'pzrw':False

}

# Change lookup names at own risk - will need to modify any direct references throught the code, too.

ng_col_lookup = {
  
  'coadd':'coadd_objects_id',
  'ra':'ra',
  'dec':'dec',
  'e1':'exp_e_1',
  'e2':'exp_e_2',
  'psf1':'psfrec_e_1',
  'psf2':'psfrec_e_2',
  'psffwhm':'psfrec_t',
  's1':'exp_e_sens_1',
  's2':'exp_e_sens_2',
  'cov11':'exp_e_cov_1_1',
  'cov12':'exp_e_cov_1_2',
  'cov22':'exp_e_cov_2_2',
  'radius':'exp_t',
  'tsnr':'exp_t_s2n',
  'snr':'exp_s2n_w'
  
}

i3_col_lookup = {
  
  'coadd':'coadd_objects_id',
  'ra':'ra',
  'dec':'dec',
  'e1':'e1',
  'e2':'e2',
  'psf1':'mean_hsm_psf_e1_sky',
  'psf2':'mean_hsm_psf_e2_sky',
  'hsmpsf1':'mean_hsm_psf_e1_sky',
  'hsmpsf2':'mean_hsm_psf_e2_sky',
  'psffwhm':'mean_psf_fwhm',
  'm':'nbc_m',
  'c1':'nbc_c1',
  'c2':'nbc_c2',
  'w':'w',
  'rgp':'mean_rgpp_rp',
  'snr':'snr',
  'ra_off':'ra_as',
  'dec_off':'dec_as',
  'rad':'radius',
  'bamp':'bulge_a',
  'damp':'disc_a',
  'bflux':'bulge_flux',
  'dflux':'disc_flux',
  'ratflux':'flux_ratio',
  'resmin':'min_residuals',
  'resmax':'max_residuals',
  'modmin':'model_min',
  'modmax':'model_max',
  'like':'likelihood',
  'evals':'levmar_like_evals',
  'iter':'levmar_iterations',
  'flux':'mean_flux',
  'cov11':'covmat_2_2',
  'cov22':'covmat_3_3',
  'cov12':'covmat_2_3',
  'nexp':'n_exposure',
  'stamp':'stamp_size',
  'flagi':'sex_flag_i',
  'flagr':'sex_flag_r',
  'info':'info_flag',
  'error':'error_flag',
  'chi2pix':'chi2_pixel',
  'nrej':'nreject',
  'tile':'tilename',
  'rsnr':'round_snr',
  'maskfrac':'mean_mask_fraction',
  'fluxfrac':'mean_unmasked_flux_frac',
  'modelmu':'mean_model_edge_mu',
  'modelsig':'mean_model_edge_sigma',
  'mu':'mean_edge_mu',
  'sig':'mean_edge_sigma',
  'xoff':'x',
  'yoff':'y',
  'row':'orig_row',
  'col':'orig_col',
  'psf1_exp':'psf_e1_sky',
  'psf2_exp':'psf_e2_sky',
  'psf1b_exp':'psf_e1',
  'psf2b_exp':'psf_e2',
  'hsmpsf1_exp':'hsm_psf_e1_sky',
  'hsmpsf2_exp':'hsm_psf_e2_sky',
  'expnum':'expnum',
  'ccd':'ccd',
  'g':'mag_auto_g',
  'r':'mag_auto_r',
  'i':'mag_auto_i',
  'z':'mag_auto_z',
  'pz':'desdm_pz',
  'gold_flag':'gold_flag',
  'gold_mask':'gold_mask'

}


gold_col_lookup = {
  

  'coadd':True,
  'ra':True,
  'dec':True,
  'sva1_gold_flags':True,
  'sva1_gold_flags':True,
  'sva1_spte_flags':True,
  'sva1_gold_mag_flags':True,
  'photoz_bin':True,
  'ngmix_flags':True,
  'im3shape_flags':True,
  'mag_auto_r':True,
  'mag_auto_g':True,
  'mag_auto_i':True,
  'mag_auto_z':True,

}

gal_col_lookup = {

  'coadd':'COADD_OBJECTS_ID',
  'ra':'RA',
  'dec':'DEC',
  'zp':'ZREDMAGIC',
  'l':'ZLUM',
  'error':'REDMAGICFLAG'

}


log_val = {
  
  'ra':False,
  'dec':False,
  'e1':False,
  'e2':False,
  'psf1':False,
  'psf2':False,
  'psffwhm':False,
  'm':False,
  'c1':False,
  'c2':False,
  'w':False,
  'rgp':False,
  'snr':True,
  'ra_off':False,
  'dec_off':False,
  'rad':False,
  'bamp':False,
  'damp':False,
  'bflux':False,
  'dflux':False,
  'ratflux':False,
  'resmin':False,
  'resmax':False,
  'modmin':False,
  'modmax':False,
  'like':True,
  'nlike':True,
  'evals':True,
  'iter':True,
  'flux':True,
  'cov11':False,
  'cov22':False,
  'cov12':False,
  'nexp':False,
  'stamp':False,
  'flagi':False,
  'flagr':False,
  'info':False,
  'error':False,
  'chi2pix':False,
  'nrej':False,
  'tile':False,
  'rsnr':True,
  'maskfrac':False,
  'fluxfrac':False,
  'modelmu':False,
  'modelsig':False,
  'mu':False,
  'sig':False,
  'pos':False,
  'psfpos':False,
  'invfluxfrac':True,
  'airmass':False,
  'exptime':False,
  'fwhm':False,
  'maglimit':False,
  'skybrite':False,
  'skysigma':False

}

lbl = {
  
  'coadd':'Coadd Objects ID',
  'ra':'RA',
  'dec':'Dec',
  'e1':r'$e_1$',
  'e2':r'$e_2$',
  'psf1':r'PSF $e_1$',
  'psf2':r'PSF $e_2$',
  'psffwhm':'PSF FWHM',
  'm':'m',
  'c1':r'$c_1$',
  'c2':r'$c_2$',
  'w':'Weight',
  'rgp':r'$R_{gp}/R_{p}$',
  'snr':'Signal-to-Noise',
  'ra_off':'RA offset',
  'dec_off':'Dec offset',
  'rad':'Radius',
  'bamp':'Bulge Amp.',
  'damp':'Disc Amp.',
  'bflux':'Bulge Flux',
  'dflux':'Disc Flux',
  'ratflux':'Flux Ratio',
  'resmin':'Min Model Res.',
  'resmax':'Max Model Res.',
  'modmin':'Model Min.',
  'modmax':'Model Max.',
  'like':'Likelihood',
  'nlike':'Neg. Likelihood',
  'evals':'Levmar Evals',
  'iter':'Iterations',
  'flux':'Flux',
  'cov11':'Cov 11',
  'cov22':'Cov 22',
  'cov12':'Cov 12',
  'nexp':'Num. Exposures',
  'stamp':'Stamp Size',
  'flagi':'Flags i',
  'flagr':'Flags r',
  'info':'Info Flag',
  'error':'Error Flag',
  'chi2pix':r'$\chi^{2}$ pix$^{-1}$',
  'nrej':'Num. Reject Exp.',
  'tile':'Tile Name',
  'rsnr':'Round SNR',
  'maskfrac':'Mask Frac.',
  'fluxfrac':'UnMask Flux Frac.',
  'modelmu':r'Mdoel Edge $\mu$',
  'modelsig':r'Model Edge $\sigma$',
  'mu':r'$Edge \mu$',
  'sig':r'$Edge \sigma$',
  'pos':'Pos. Angle',
  'psfpos':'PSF Pos. Angle',
  'invfluxfrac':'Inv. Flux Frac.',
  'airmass':'Air Mass',
  'exptime':'Exp. Time',
  'fwhm':'FWHM',
  'maglimit':'Mag. Limit',
  'skybrite':'Sky Brightness',
  'skysigma':r'Sky $\sigma$',
  'dpsf':r'PSF $e_1-$ PSF $e_2$',
  'g':r'$g$ magauto',
  'r':r'$r$ magauto',
  'i':r'$i$ magauto',
  'z':r'$z$ magauto',
  'pz':r'$z$'

}

map_name_y1 = {
  
  'airmass':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',
  'exptime':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz',
  'fwhm':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',
  'maglimit':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit2__.fits.gz',
  'skybrite':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',
  'skysigma':y1sysmapdir+'Y1A1_SPT_and_S82_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz'
  
}

map_name_sv = {

  'airmass':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_AIRMASS_coaddweights_mean.fits.gz',
  'ebv':svsysmapdir+'Planck_EBV_2048r_Q.fits',
  'exptime':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz',
  'fwhm':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_FWHM_coaddweights_mean.fits.gz',
  'maglimit':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_maglimit__.fits.gz',
  'skybrite':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz',
  'skysigma':svsysmapdir+'SVA1_IMAGE_SRC_band_r_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz'
  
}

