import numpy as np
import healpy as hp
import fitsio as fio
import sys_split
import catalog
import config
import treecorr
import time
import pickle
import glob
import twopoint as tp
from scipy.optimize import curve_fit

'''
Modified code by Oliver to turn flask catalogs by Lucas 
(the ones called kgg-s*-f2z*_c*.fits) into source catalogs 
at Y1 position (with optional pixel-to-pixel weighting).

'''

def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)

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

    if wfile is None:

      map_ra   = phi/np.pi*180.0
      map_dec  = (np.pi/2.0 - theta)/np.pi*180.0
      map_g1   = fmap['Q_STOKES']
      map_g2   = fmap['U_STOKES']
      map_w    = np.ones(len(map_ra))
      map_sige = np.sqrt((sig_orig[zbin]**2/2.)*(neff_new/neff_pix))

    else:

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

    out = np.zeros(len(map_ra),dtype=[('ra','f4')]+[('dec','f4')]+[('e1','f4')]+[('e2','f4')]+[('w','f4')])      
    out['ra']  = map_ra
    out['dec'] = map_dec
    out['e1']  = map_g1 + np.random.randn(len(map_ra))*map_sige/np.sqrt(n)
    out['e2']  = map_g2 + np.random.randn(len(map_ra))*map_sige/np.sqrt(n)
    out['w']   = map_w*n
    out        = out[n!=0]
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

class run(object):

  @staticmethod
  def loop_2pt(zbin,catname,val):

    t0=time.time()

    cnt=0
    for j in range(100):
      for i in range(8):
        for k in range(3):
          print j,i,k,time.time()-t0
          if k==0:
            wfile = None
          else:
            wfile='text/pzrw_'+catname+'_'+val+'_'+str(k-1)+'.fits.gz'

          out = methods.rotate_mock_rescale_nsigma(zbin, i+1, j+1, wfile=wfile)
          cat = treecorr.Catalog(g1=out['e1'], g2=out['e2'], w=out['w'], ra=out['ra'], dec=out['dec'], ra_units='deg', dec_units='deg')
          gg  = treecorr.GGCorrelation(nbins=20, min_sep=2.5, max_sep=250., sep_units='arcmin', bin_slop=0.2, verbose=0)
          gg.process(cat)

          d = {
            'theta' : np.exp(gg.meanlogr),
            'xip' : gg.xip,
            'xim' : gg.xim,
            'err' : np.sqrt(gg.varxi)
          }

          save_obj(d,'text/flask_GG_'+catname+'_'+val+'_'+str(cnt)+'_'+str(k)+'.cpickle')
        cnt+=1

    return

  @staticmethod
  def get_data_cov(zbin):

    if config.cov.get('path') is None:
      return None
    else:
      try:
        cov=tp.TwoPointFile.from_fits(config.cov.get('path')).covmat_info
      except:
        return None

      ind0 = {1:0,2:80,3:140,4:180}
      ind1 = {1:20,2:100,3:160,4:199}

      xip = cov.covmat[ind0[zbin]:ind1[zbin],ind0[zbin]:ind1[zbin]]
      xim = cov.covmat[ind0[zbin]:ind1[zbin],ind0[zbin]:ind1[zbin]]

      print xip
      print xim

      return xip,xim

  @staticmethod
  def amp_fit(x,y,cov):

    def func(x,a):
      return a*x

    print np.shape(cov),len(x),len(y)

    params=curve_fit(func,x,y,p0=(0.1),sigma=np.diagonal(cov),bounds=(-2,2))

    return params[0]

  @staticmethod
  def get_amp_cov(zbin,catname,val):

    a=[]
    b=[]
    c=[]
    covp,covm = run.get_data_cov(zbin)
    for i in range(800):
      try:
        d0 = load_obj('text/flask_GG_'+str(i)+'_0.cpickle')
        d1 = load_obj('text/flask_GG_'+str(i)+'_1.cpickle')
        d2 = load_obj('text/flask_GG_'+str(i)+'_2.cpickle')
      except IOError:
        continue

      a.append( run.amp_fit(d0['xip'],d1['xip'],covp) )
      b.append( run.amp_fit(d0['xip'],d2['xip'],covp) )
      c.append( run.amp_fit(d1['xip'],d2['xip'],covp) )

    a=np.array(a)
    b=np.array(b)
    c=np.array(c)
    acov=np.sum((a-np.mean(a))*(a-np.mean(a)))*(len(a)-1.)/len(a)
    bcov=np.sum((b-np.mean(b))*(b-np.mean(b)))*(len(b)-1.)/len(b)
    ccov=np.sum((c-np.mean(c))*(c-np.mean(c)))*(len(c)-1.)/len(c)

    print 'a',np.mean(a),acov
    print 'b',np.mean(b),bcov
    print 'c',np.mean(c),ccov

    return
