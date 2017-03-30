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
  1 : 1.30157,
  2 : 1.41888,
  3 : 1.45387,
  4 : 0.77419
}

neff_i3 = {
  1 : 1.30157,
  2 : 1.41888,
  3 : 1.45387,
  4 : 0.77419
}

sig_mcal = {
  1 : 0.39,
  2 : 0.39,
  3 : 0.39,
  4 : 0.39
}

sig_i3 = {
  1 : 0.39,
  2 : 0.39,
  3 : 0.39,
  4 : 0.39
}

euler_angle_1 = {
  0 : 0.0,
  1 : 180.0,
  2 : 0.0,
  3 : 90.0,
  4 : 180.0,
  5 : 270.0,
  6 : 0.0,
  7 : 180.0
}

euler_angle_2 = {
  0 : 0.0,
  1 : 0.0,
  2 : 60.0,
  3 : 60.0,
  4 : 60.0,
  5 : 60.0,
  6 : 120.0,
  7 : 120.0
}

class methods(object):

  @staticmethod
  def rotate_mock_rescale_nsigma( zbin,                 # tomographic bin index - 0 to n
                                  rlsn,                 # mock cutout (realisation)
                                  mapfile,              # input shear map file name + path
                                  out_file,             # output file name + path
                                  wfile=None,           # input weight file (optional)
                                  neff_orig=neff_mcal,  # dictionary for neff
                                  sig_orig=sig_mcal,    # dictionary for sigma_e
                                  neff_ratio=0.5,       # ratio of original to new neff (default half density)
                                  nside=4096):          # nside of maps (default 4096)

    neff_new = neff_orig[zbin]*neff_ratio
    npix     = hp.nside2npix(nside)
    map_g1   = np.zeros(npix)
    map_g2   = np.zeros(npix)
    map_w    = np.zeros(npix)
    map_sige = np.zeros(npix)
    map_g1.fill(-100.)
    map_g2.fill(-100.)
    map_w.fill(-100.)
    map_sige.fill(-100.)

    fmap = fio.FITS(mapfile)[-1].read(columns=['PIXEL','Q_STOKES','U_STOKES'])
    pix  = fmap['PIXEL']

    if wfile is not None:
      w     = fio.FITS(wfile)[-1].read()
      wpix  = w['pix']
      s1,s2 = catalog.CatalogMethods.sort2(pix,wpix) # cut to intersection of data / mock catalogs
      fmap  = fmap[s1]
      pix   = pix[s1]
      wpix  = wpix[s2]
      w     = w['weight'][s2]
      w2    = w['weightsq'][s2]
    else:
      wpix  = pix
      w     = np.ones(len(pix))
      w2    = np.ones(len(pix))

    theta, phi         = hp.pix2ang(nside,pix) #theta and phi of the footprint pixels
    pix_rotator        = hp.Rotator(deg=False, rot=[euler_angle_1[int(rlsn)]*np.pi/180., euler_angle_2[int(rlsn)]*np.pi/180.])
    theta_rot, phi_rot = pix_rotator(theta,phi)
    rot_pix            = hp.ang2pix(nside,theta_rot,phi_rot)

    neff_pix          = 1. / (hp.nside2pixarea(nside, degrees=True)*3600.)
    map_g1[rot_pix]   = fmap['Q_STOKES']
    map_g2[rot_pix]   = fmap['U_STOKES']
    map_w[rot_pix]    = w
    map_sige[rot_pix] = np.sqrt((sig_orig[zbin]**2/2.)*(neff_new/neff_pix)*(w2/w))

    n       = np.random.poisson(neff_new/neff_pix,size=len(rot_pix))
    gal_pix = np.repeat(rot_pix,n)

    out = np.zeros(len(gal_pix),dtype=[('ra','f4')]+[('dec','f4')]+[('e1','f4')]+[('e2','f4')]+[('w','f4')])      
    out['dec'], out['ra'] = catalog.CatalogMethods.hpix_to_radec(th, ph, nside=nside, nest=False)
    out['e1'] = map_g1[gal_pix] + np.random.randn(len(gal_pix))*map_sige[gal_pix]
    out['e2'] = map_g2[gal_pix] + np.random.randn(len(gal_pix))*map_sige[gal_pix]
    out['w']  = map_w[gal_pix]
    fio.write(out_file,out,clobber=True)

    return 

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

