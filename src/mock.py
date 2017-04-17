import numpy as np
import healpy as hp
import fitsio as fio
import sys_split
import catalog
import config
try:
  import treecorr
except:
  print "No treecorr"
  treecorr=None
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

zbounds = {
  1 : [0.2, 0.43],
  2 : [0.43,0.63],
  3 : [0.63,0.9],
  4 : [0.9,1.3]
}

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
                                  neff_ratio=0.5,       # ratio of original to new neff (default half density)
                                  nside=4096):          # nside of maps (default 4096)

    # out = mock.methods.rotate_mock_rescale_nsigma(3, 1, 1, wfile='text/pzrw_metacalibration_snr_0.fits.gz')

    npix     = hp.nside2npix(nside)
    neff_pix = 1. / (hp.nside2pixarea(nside, degrees=True)*3600.)
    neff_new = neff_orig[zbin]*neff_ratio

    mapfile = '/fs/scratch/cond0083/flask/kgg-s'+str(seed)+'-f2z'+str(zbin)+'_c'+str(rlsn)+'.fits'
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
      map_neff = neff_new/neff_pix#*map_w/np.mean(map_w)
      map_sige = np.sqrt((sig_orig[zbin]**2/2.)*(map_neff))

    else:

      w        = fio.FITS(wfile)[-1].read()
      ind      = np.argsort(w['pix'])
      w        = w[ind]
      ind      = np.argsort(rot_pix)
      rot_pix  = rot_pix[ind]
      fmap     = fmap[ind]
      theta    = theta[ind]
      phi      = phi[ind]

      w        = w[np.in1d(w['pix'],rot_pix,assume_unique=False)]
      mask     = np.in1d(rot_pix,w['pix'],assume_unique=False)
      fmap     = fmap[mask]
      rot_pix  = rot_pix[mask]
      theta    = theta[mask]
      phi      = phi[mask]

      diff     = np.diff(rot_pix)
      diff     = np.where(diff!=0)[0]+1
      diff     = np.append([0],diff)
      diff     = np.append(diff,[None])

      w1       = np.zeros(len(rot_pix))
      # w2       = np.zeros(len(rot_pix))
      for i in xrange(len(diff)-1):
        w1[diff[i]:diff[i+1]] = w['weight'][i]
        # w2[diff[i]:diff[i+1]] = w['weightsq'][i]
      w        = None

      map_ra   = phi/np.pi*180.0
      map_dec  = (np.pi/2.0 - theta)/np.pi*180.0
      map_g1   = fmap['Q_STOKES']
      map_g2   = fmap['U_STOKES']
      map_w    = w1
      map_neff = neff_new/neff_pix#*map_w/np.mean(map_w)
      # map_sige = np.sqrt((sig_orig[zbin]**2/2.)*(map_neff)*(w2/w1))
      map_sige = np.sqrt((sig_orig[zbin]**2/2.)*(map_neff))*np.ones(len(map_ra))

    fmap     = None
    n        = np.random.poisson(map_neff,size=len(map_ra))

    out = np.zeros(np.sum(n),dtype=[('ra','f4')]+[('dec','f4')]+[('e1','f4')]+[('e2','f4')]+[('w','f4')])
    out['ra']  = np.repeat(map_ra,n)
    out['dec'] = np.repeat(map_dec,n)
    out['e1']  = np.repeat(map_g1,n)+np.random.randn(len(out))*np.repeat(map_sige,n)
    out['e2']  = np.repeat(map_g2,n)+np.random.randn(len(out))*np.repeat(map_sige,n)
    out['w']   = np.repeat(map_w,n)
    # out['w']   = np.ones(len(out))
    # fio.write(out_file,out,clobber=True)

    return out # original pixel positions (not DES pixel positions)

  @staticmethod
  def save_weights(cat,val,zbin,w,bins,mask0):

    for i in xrange(cat.sbins):
      if cat.cat=='mcal':
        mask = bins[i][np.in1d(bins[i],np.where(w!=0)[0],assume_unique=False)]
      else:
        mask = (bins==i)&mask0&(w!=0)
      pix  = catalog.CatalogMethods.radec_to_hpix(cat.ra[mask],cat.dec[mask],nside=4096,nest=False)
      upix = np.unique(pix)
      w0   = w[mask]
      w1   = np.bincount(pix,weights=w0)
      w2   = np.bincount(pix,weights=w0*w0)
      mask = np.where(w1!=0)[0]
      # upix = upix[mask]
      w1   = w1[mask]
      w2   = w2[mask]

      out  = np.empty(len(upix),dtype=[('pix',int)]+[('weight','f4')]+[('weightsq','f4')])
      out['pix']      = upix
      out['weight']   = w1
      out['weightsq'] = w2
      fio.write('text/pzrw_'+cat.name+'_'+val+'_'+str(zbin+1)+'_'+str(i)+'.fits.gz',out,clobber=True)

    return

class run(object):

  @staticmethod
  def loop_2pt_data(cat,col):

    if col=='snr':
      no2pt=None
    else:
      no2pt=1
    for i in range(4):
      if cat.cat == 'mcal':
        catalog.CatalogMethods.add_cut_sheared(cat,'pz',cmin=zbounds[i+1][0],cmax=zbounds[i+1][1],remove=False)
        xip,xim,gt,split,edge=sys_split.split_methods.split_gals_2pt_along(cat,None,col,blind=False,plot=False,no2pt=no2pt,zbin=i)
        catalog.CatalogMethods.add_cut_sheared(cat,'pz',cmin=zbounds[i+1][0],cmax=zbounds[i+1][1],remove=True)
      else:
        mask = (cat.pz>zbounds[i+1][0])&(cat.pz<zbounds[i+1][1])
        xip,xim,gt,split,edge=sys_split.split_methods.split_gals_2pt_along(cat,None,col,mask=mask,blind=False,plot=True,no2pt=no2pt,zbin=i)

    return

  @staticmethod
  def loop_2pt(catname,val):

    t0=time.time()

    cnt=0
    for j in range(36):
      for i in range(8):
        for zbin in range(4):
          for k in range(3):
            print j,i,k,time.time()-t0
            if k==0:
              wfile = None
              if val!='snr':
                continue
            else:
              wfile='text/pzrw_'+catname+'_'+val+'_'+str(zbin+1)+'_'+str(k-1)+'.fits.gz'

            out = methods.rotate_mock_rescale_nsigma(zbin+1, i+1, j+1, wfile=wfile)
            cat = treecorr.Catalog(g1=out['e1'], g2=out['e2'], w=out['w'], ra=out['ra'], dec=out['dec'], ra_units='deg', dec_units='deg')
            gg  = treecorr.GGCorrelation(nbins=20, min_sep=2.5, max_sep=250., sep_units='arcmin', bin_slop=0.2, verbose=0)
            gg.process(cat)

            d = {
              'theta' : np.exp(gg.meanlogr),
              'xip' : gg.xip,
              'xim' : gg.xim,
              'err' : np.sqrt(gg.varxi)
            }

            save_obj(d,'text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(cnt)+'_'+str(k)+'.cpickle')
        cnt+=1

    return

  @staticmethod
  def loop_2pt_noweight(catname,ii):

    t0=time.time()

    cnt=0
    for j in range(36):
      for i in range(8):
        for zbin in range(4):
          for k in range(3):
            if cnt>225:
              continue
            if (cnt+1)%5!=ii:
              continue
            print j,i,k,time.time()-t0
            out = methods.rotate_mock_rescale_nsigma(zbin+1, i+1, j+1, wfile=None)
            cat = treecorr.Catalog(g1=out['e1'], g2=out['e2'], w=out['w'], ra=out['ra'], dec=out['dec'], ra_units='deg', dec_units='deg')
            gg  = treecorr.GGCorrelation(nbins=20, min_sep=2.5, max_sep=250., sep_units='arcmin', bin_slop=0.2, verbose=0)
            gg.process(cat)

            d = {
              'theta' : np.exp(gg.meanlogr),
              'xip' : gg.xip,
              'xim' : gg.xim,
              'err' : np.sqrt(gg.varxi)
            }

            save_obj(d,'text/flask_GG_'+catname+'_noweight_'+str(zbin)+'_'+str(cnt)+'.cpickle')
        cnt+=1

    return

  @staticmethod
  def get_data_cov(zbin,full=False):

    if config.cov.get('path') is None:
      return None
    else:
      try:
        cov=tp.TwoPointFile.from_fits(config.cov.get('path')).covmat_info
      except:
        return None

      ind0 = {0:0,1:80,2:140,3:180}
      ind1 = {0:20,1:100,2:160,3:200}

      if full:
        xip = np.zeros((80,80))
        xim = np.zeros((80,80))
        for i in range(4):
          xip[i*20:(i+1)*20,i*20:(i+1)*20] = cov.covmat[cov.starts[0]+ind0[zbin]:cov.starts[0]+ind1[zbin],cov.starts[0]+ind0[zbin]:cov.starts[0]+ind1[zbin]]
          xim[i*20:(i+1)*20,i*20:(i+1)*20] = cov.covmat[cov.starts[1]+ind0[zbin]:cov.starts[1]+ind1[zbin],cov.starts[1]+ind0[zbin]:cov.starts[1]+ind1[zbin]]
        return xip,xim

      xip = cov.covmat[cov.starts[0]+ind0[zbin]:cov.starts[0]+ind1[zbin],cov.starts[0]+ind0[zbin]:cov.starts[0]+ind1[zbin]]
      xim = cov.covmat[cov.starts[1]+ind0[zbin]:cov.starts[1]+ind1[zbin],cov.starts[1]+ind0[zbin]:cov.starts[1]+ind1[zbin]]

      return xip,xim

  @staticmethod
  def amp_fit(xip,dxip,cov):

    chi2st=999999
    amp=0.
    for a in xrange(-200,200):
      x=dxip-a*xip/100.
      chi2=np.dot(x,np.dot(np.linalg.inv(cov),x))
      if chi2<chi2st:
        chi2st=chi2
        amp=a/100.

    return amp

  @staticmethod
  def get_chi2(xip,cov):

    chi2=np.dot(xip,np.dot(np.linalg.inv(cov),xip))/(len(xip)-1)
    print xip,np.diagonal(cov)

    return chi2

  @staticmethod
  def get_amp_cov(zbin,catname,val,xi,full=False):

    a=[]
    b=[]
    c=[]
    ddd=[]
    covp,covm = run.get_data_cov(zbin,full)
    for i in range(800):
      dd0 = []
      dd1 = []
      dd2 = []
      if full:
        try:
          for zbin in range(4):
            d0 = load_obj('text/flask_GG_'+catname+'_'+'snr'+'_'+str(zbin)+'_'+str(i)+'_0.cpickle')
            d1 = load_obj('text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(i)+'_1.cpickle')
            d2 = load_obj('text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(i)+'_2.cpickle')
            dd0 = np.append(dd0,d0[xi])
            dd1 = np.append(dd1,d1[xi])
            dd2 = np.append(dd2,d2[xi])
          dd0 = np.array(dd0)
          dd1 = np.array(dd1)
          dd2 = np.array(dd2)
        except IOError:
          continue

        if len(dd0)<80:
          print i

      else:
        try:
          d0 = load_obj('text/flask_GG_'+catname+'_'+'snr'+'_'+str(zbin)+'_'+str(i)+'_0.cpickle')
          d1 = load_obj('text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(i)+'_1.cpickle')
          d2 = load_obj('text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(i)+'_2.cpickle')
        except IOError:
          continue

        dd0 = d0[xi]
        dd1 = d1[xi]
        dd2 = d2[xi]

      if xi == 'xip':
        cov0 = covp
      else:
        cov0 = covm
      a.append( run.amp_fit(dd0,dd2-dd0,cov0) )
      b.append( run.amp_fit(dd0,dd0-dd1,cov0) )
      c.append( run.amp_fit(dd0,dd2-dd1,cov0) )

      ddd.append(dd2-dd1)

    print np.array(dd2)-np.array(dd1)
    a=np.array(a)
    b=np.array(b)
    c=np.array(c)
    ddd=np.array(ddd)
    acov=np.mean((a-np.mean(a))*(a-np.mean(a)))*(len(a)-1)/(len(a)-1-1)
    bcov=np.mean((b-np.mean(b))*(b-np.mean(b)))*(len(b)-1)/(len(b)-1-1)
    ccov=np.mean((c-np.mean(c))*(c-np.mean(c)))*(len(c)-1)/(len(c)-1-1)
    cov = np.zeros((len(ddd[0,:]),len(ddd[0,:])))
    for i in range(len(cov[0,:])):
      for j in range(len(cov[0,:])):
        cov[i,j]=np.mean((ddd[:,i]-np.mean(ddd[:,i]))*(ddd[:,j]-np.mean(ddd[:,j])))*(len(ddd)-1)/(len(ddd)-len(ddd[0,:])-1)

    # print 'a',np.mean(a),acov
    # print 'b',np.mean(b),bcov
    # print 'c',np.mean(c),ccov

    return acov, bcov, ccov, cov

  @staticmethod
  def cat_2pt_results(catname,full=False):

    if catname == 'im3shape':
      vals = ['snr','psf1','psf2','rgp','ebv','skybrite','fwhm','airmass','maglim','colour']      
    else:
      vals = ['snr','psf1','psf2','size','ebv','skybrite','fwhm','airmass','maglim','colour']

    print catname
    for xi in ['xip','xim']:
      for val in vals:
        # if val=='snr':
        #   continue
        if full:
          dd0 = []
          dd1 = []
          dd2 = []
        for zbin in range(4):
          covp,covm = run.get_data_cov(zbin,full)
          d0 = load_obj('text/data_GG_'+catname+'_'+str(zbin)+'.cpickle')
          d1 = load_obj('text/data_GG_'+catname+'_'+val+'_'+str(zbin)+'_1.cpickle')
          d2 = load_obj('text/data_GG_'+catname+'_'+val+'_'+str(zbin)+'_2.cpickle')
          if full:
            dd0 = np.append(dd0,d0[xi])
            dd1 = np.append(dd1,d1[xi])
            dd2 = np.append(dd2,d2[xi])
          else:
            dd0 = d0[xi]
            dd1 = d1[xi]
            dd2 = d2[xi]
          if xi == 'xip':
            cov0=covp
          else:
            cov0=covm
          if full&(zbin<3):
            continue
          if full:
            dd0 = np.array(dd0)
            dd1 = np.array(dd1)
            dd2 = np.array(dd2)
          a = run.amp_fit(dd0,dd2-dd0,cov0)
          b = run.amp_fit(dd0,dd0-dd1,cov0)
          c = run.amp_fit(dd0,dd2-dd1,cov0)
          acov,bcov,ccov,cov = run.get_amp_cov(zbin,catname,val,xi,full)
          chi2 = run.get_chi2(dd2-dd1,cov)
          print xi, val, zbin, 'a = '+str(np.around(a,2))+' +- '+str(np.around(np.sqrt(acov),2))
          print xi, val, zbin, 'b = '+str(np.around(b,2))+' +- '+str(np.around(np.sqrt(bcov),2))
          print xi, val, zbin, 'c = '+str(np.around(c,2))+' +- '+str(np.around(np.sqrt(ccov),2))
          print xi, val, zbin, 'red. chi2 = ',str(np.around(chi2,2))

    return

# a = np.zeros((168,20))
# for i in range(168):
#   d0 = load_obj('text/flask_GG_'+catname+'_'+val+'_'+str(zbin)+'_'+str(i)+'_1.cpickle')
#   a[i,:] = d0['xip']

# print np.sqrt(np.sum((a-np.mean(a,axis=0))*(a-np.mean(a,axis=0)),axis=0)*((168)-1.)/(168)*(168-20-1)/(168-1))



# (np.sum(out['w']**2)/(np.sum(out['w']))**2)**-1/(len(out)*0.73766043036137707)

# (np.sum(out['e1']**2*out['w']**2)/np.sum(out['w']**2))

# for zbin in range(4):
#   a = np.zeros((224,20))
#   for i in range(224):
#     d = load_obj('text/flask_GG_metacalibration_noweight_'+str(zbin)+'_'+str(i)+'.cpickle')
#     a[i,:] = d['xip']
#   print zbin, np.sqrt(np.sum((a-np.mean(a,axis=0))*(a-np.mean(a,axis=0)),axis=0)/len(a))*(len(a)-20-1)/(len(a)-1)

# import src.config as config
# config.cov['path']='../des-mpp/y1_mcal/2pt_fits/2pt_NG.fits'
# import src.mock as mock
# mock.run.cat_2pt_results('im3shape',full=True)

