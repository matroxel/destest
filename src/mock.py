import numpy as np
import healpy as hp
import fitsio as fio
import catalog
import config

'''
Modified code by Oliver to turn flask catalogs by Lucas 
(the ones called kgg-s*-f2z*_c*.fits) into source catalogs 
at Y1 position (with optional pixel-to-pixel weighting).

'''

neff_mcal = {
  1 : 1.736599,
  2 : 1.464735,
  3 : 1.62532,
  4 : 0.811633
}

neff_i3 = {
  1 : 1.21746,
  2 : 0.98196,
  3 : 0.983865,
  4 : 0.36885
}

sig_mcal = {
  1 : 0.394925,
  2 : 0.394925,
  3 : 0.394925,
  4 : 0.394925
}

sig_i3 = {
  1 : 0.3867957,
  2 : 0.3867957,
  3 : 0.3867957,
  4 : 0.3867957
}

euler_angle_1 = {
  1 : 0.0,
  2 : 180.0,
  3 : 0.0,
  4 : 90.0,
  5 : 180.0,
  6 : 270.0,
  7 : 0.0,
  8 : 180.0
}

euler_angle_2 = {
  1 : 0.0,
  2 : 0.0,
  3 : 60.0,
  4 : 60.0,
  5 : 60.0,
  6 : 60.0,
  7 : 120.0,
  8 : 120.0
}

class methods(object):

  @staticmethod
  def rotate_mock_rescale_nsigma( zbin,                 # tomographic bin index - 0 to n
                                  rlsn,                 # mock cutout (realisation)
                                  seed,                 # mock seed 
                                  wfile,                # input weight file (optional)
                                  out_file='tmp.fits',  # output file name + path
                                  neff_orig=neff_mcal,  # dictionary for neff
                                  sig_orig=sig_mcal,    # dictionary for sigma_e
                                  neff_ratio=1.0,       # ratio of original to new neff (default half density)
                                  nside=4096):          # nside of maps (default 4096)

    # out = mock.methods.rotate_mock_rescale_nsigma(3, 1, 1, wfile='text/pzrw_metacalibration_snr_0.fits.gz')

    npix     = hp.nside2npix(nside)
    neff_pix = 1. / (hp.nside2pixarea(nside, degrees=True)*3600.)
    neff_new = neff_orig[zbin]*neff_ratio

    mapfile = '/global/cscratch1/sd/seccolf/y1_patch/seed'+str(seed)+'/kgg-s'+str(seed)+'-f2z'+str(zbin)+'_c'+str(rlsn)+'.fits'
    fmap = fio.FITS(mapfile)[-1].read(columns=['PIXEL','Q_STOKES','U_STOKES'])
    theta, phi         = hp.pix2ang(nside,fmap['PIXEL']) #theta and phi of the footprint pixels
    pix_rotator        = hp.Rotator(deg=False, rot=[euler_angle_1[int(rlsn)]*np.pi/180., euler_angle_2[int(rlsn)]*np.pi/180.])
    theta_rot, phi_rot = pix_rotator(theta,phi)
    rot_pix            = hp.ang2pix(nside,theta_rot,phi_rot)

    w        = fio.FITS(wfile)[-1].read()
    s1,s2    = catalog.CatalogMethods.sort2(rot_pix,w['pix']) # cut to intersection of data / mock catalogs
    rot_pix  = rot_pix[s1]

    wpix     = w['pix'][s2]
    w1       = w['weight'][s2]
    w2       = w['weightsq'][s2]
    w        = None

    xsorted  = np.argsort(wpix) # translate wpix to rot_pix
    mask     = xsorted[np.searchsorted(wpix[xsorted], rot_pix)]

    map_ra   = phi[s1]/np.pi*180.0
    map_dec  = (np.pi/2.0 - theta[s1])/np.pi*180.0
    map_g1   = fmap['Q_STOKES'][s1]
    map_g2   = fmap['U_STOKES'][s1]
    map_w    = w1[mask]
    map_sige = np.sqrt((sig_orig[zbin]**2/2.)*(neff_new/neff_pix)*(w2/w1))[mask]
    fmap     = None

    n        = np.random.poisson(neff_new/neff_pix,size=len(map_ra))
    print n,n.min(),n.max()
    print map_sige,map_sige.min(),map_sige.max()

    out = np.zeros(len(map_ra),dtype=[('ra','f4')]+[('dec','f4')]+[('e1','f4')]+[('e2','f4')]+[('w','f4')])      
    out['ra']  = map_ra
    out['dec'] = map_dec
    out['e1']  = map_g1 + np.random.randn(len(map_ra))*map_sige/np.sqrt(n)
    out['e2']  = map_g2 + np.random.randn(len(map_ra))*map_sige/np.sqrt(n)
    out['w']   = map_w*n
    # fio.write(out_file,out,clobber=True)

    return out # original pixel positions (not DES pixel positions)

  @staticmethod
  def save_weights(cat,val,w,bins):

    for i in xrange(cat.sbins):
      if cat.cat=='mcal':
        mask = bins[i]
      else:
        mask = bins==i
      pix  = catalog.CatalogMethods.radec_to_hpix(cat.ra[mask],cat.dec[mask],nside=4096,nest=False)
      upix = np.unique(pix)
      w0   = w[mask]
      w1   = np.bincount(pix,weights=w0)
      w2   = np.bincount(pix,weights=w0*w0)
      mask = np.where(w1!=0)[0]
      w1   = w1[mask]
      w2   = w2[mask]

      out  = np.empty(len(upix),dtype=[('pix',int)]+[('weight','f4')]+[('weightsq','f4')])
      out['pix']      = upix
      out['weight']   = w1
      out['weightsq'] = w2
      fio.write('text/pzrw_'+cat.name+'_'+val+'_'+str(i)+'.fits.gz',out,clobber=True)

    return

