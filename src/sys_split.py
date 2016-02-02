import numpy as np
import fitsio as fio
import healpy as hp

import catalog
import config
import fig
import txt
import lin
import corr


class split(object):

  @staticmethod
  def cat_splits_lin(vals,cat,mask=None):
    """
    Loop over array names in vals of CatalogStore cat with optional masking. Calls split_methods.split_gals_lin_along().
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    txt.write_methods.heading('Linear Splits',cat,label='linear_splits',create=True)

    for x in vals:

      txt.write_methods.heading(x,cat,label='linear_splits',create=False)

      arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=split_methods.split_gals_lin_along(cat,x,mask=mask,log=config.log_val.get(x,None),plot=True)
      
      txt.write_methods.write_append(x+'  '+str(arr1)+'  '+str(arr1err),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('e  '+str(e1)+'  '+str(e2),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('e err  '+str(e1err)+'  '+str(e2err),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('slope  '+str(m1)+'  '+str(m2),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('slope err  '+str(m1err)+'  '+str(m2err),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('intercept  '+str(b1)+'  '+str(b2),cat,label='linear_splits',create=False)
      txt.write_methods.write_append('intercept err  '+str(b1err)+'  '+str(b2err),cat,label='linear_splits',create=False)

    return

  @staticmethod
  def cat_splits_2pt(vals,cat,cat2=None,mask=None):
    """
    Loop over array names in vals of CatalogStore cat with optional masking. Calls split_methods.split_gals_2pt_along().
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    if cat2 is None:
      cat2=cat        

    txt.write_methods.heading('2pt Splits',cat,label='2pt_splits',create=True)

    for x in vals:

      txt.write_methods.heading(x,cat,label='2pt_splits',create=False)

      xi,gt,split=split_methods.split_gals_2pt_along(cat,cat2,x,mask=mask,mean=True,jkon=False,mock=False,log=False,plot=True)

      txt.write_methods.write_append(x+'  '+str(np.min(getattr(cat,x)[mask]))+'  '+str(split)+'  '+str(np.max(getattr(cat,x)[mask])),cat,label='2pt_splits',create=False)

      txt.write_methods.write_append('theta  '+str(xi[0]),cat,label='2pt_splits',create=False)

      txt.write_methods.write_append('xip  '+str(xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xip  '+str((xi[1][0]-xi[2][0])/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xip  '+str((xi[3][0]-xi[2][0])/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xip  '+str((xi[3][0]-xi[1][0])/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xip err  '+str(xi[7][0]/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xip err  '+str(xi[9][0]/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xip err  '+str(xi[8][0]/xi[2][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('amp xip  '+str(xi[10][0])+'  '+str(xi[11][0])+'  '+str(xi[12][0]),cat,label='2pt_splits',create=False)

      txt.write_methods.write_append('xim  '+str(xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xim  '+str((xi[1][1]-xi[2][1])/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xim  '+str((xi[3][1]-xi[2][1])/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xim  '+str((xi[3][1]-xi[1][1])/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xim err  '+str(xi[7][1]/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xim err  '+str(xi[9][1]/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xim err  '+str(xi[8][1]/xi[2][1]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('amp xim  '+str(xi[10][1])+'  '+str(xi[11][1])+'  '+str(xi[12][1]),cat,label='2pt_splits',create=False)

      txt.write_methods.write_append('theta  '+str(gt[0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('gt  '+str(gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('lower delta gt  '+str((gt[1][0]-gt[2][0])/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('upper delta gt  '+str((gt[3][0]-gt[2][0])/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('delta gt  '+str((gt[3][0]-gt[1][0])/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('lower delta gt err  '+str(gt[7][0]/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('upper delta gt err  '+str(gt[9][0]/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('delta gt err  '+str(gt[8][0]/gt[2][0]),cat,label='2pt_tan_splits',create=False)
      txt.write_methods.write_append('amp gt  '+str(gt[10][0])+'  '+str(gt[11][0])+'  '+str(gt[12][0]),cat,label='2pt_tan_splits',create=False)

    return


class split_methods(object):

  @staticmethod
  def load_maps(cat,maps=None):
    """
    Load list of 'systematic' maps (see Boris et al) and match to positions in cat. If maps=None, load all in dictionary. See config.py.
    """

    if not hasattr(maps, '__len__'):
      if cat.release=='y1':
        maps=np.array(list(config.map_name_y1.keys()))
      elif cat.release=='sv':
        maps=np.array(list(config.map_name_sv.keys()))
    print maps
    for i,x in enumerate(maps):
      print i,x
      if x=='ebv':
        setattr(cat,x,split_methods.get_maps(cat.ra,cat.dec,x,release=cat.release,nside=2048,map=True))
      else:
        setattr(cat,x,split_methods.get_maps(cat.ra,cat.dec,x,release=cat.release))

    return  

  @staticmethod
  def get_maps(ra,dec,x,release='y1',nside=4096,map=False,nested=False):
    """
    Match 'systematic' map to catalog positions.
    """

    if release=='y1':
      xx=config.map_name_y1.get(x,None)
    elif release=='sv':
      xx=config.map_name_sv.get(x,None)

    if map:
      tmp = hp.read_map(xx)
      array = tmp[hp.ang2pix(nside, np.pi/2.-np.radians(dec),np.radians(ra), nest=nested)]
    else:
      fits=fio.FITS(xx)
      tmp=fits[-1].read()
      map = np.zeros(12*nside**2)
      map[tmp['PIXEL']] = tmp['SIGNAL']
      array = map[hp.ang2pix(nside, np.pi/2.-np.radians(dec),np.radians(ra), nest=nested)]

    return array

  @staticmethod
  def split_gals_lin_along(cat,val,mask=None,jkon=True,mock=False,log=False,label='',plot=False,fit=True):
    """
    Split catalog cat into cat.lbins equal (weighted) parts along val. Plots mean shear in bins and outputs slope and error information.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    array=getattr(cat,val)

    if log:
      array=np.log10(array)

    arr1,arr1err,e1,e1err,e2,e2err=lin.linear_methods.bin_means(array,cat,mask=mask,mock=mock,log=log)

    if fit:
      m1,m2,b1,b2,m1err,m2err,b1err,b2err=lin.fitting.sys_lin_fit(array,cat,mask=mask,log=False,noe=False)
    else:
      m1,m2,b1,b2,m1err,m2err,b1err,b2err=0.,0.,0.,0.,0.,0.,0.,0.

    # if (jkon)&(cat.use_jk==1):
    #   me1err,me2err,slp1err,slp2err,b1err,b2err=jackknife_methods.lin_err0(array,cat,label,mask0=mask,parallel=parallel)
    # elif (jkon)&(cat.use_jk==2):
    #   me1err,me2err,slp1err,slp2err,b1err,b2err=BCC_Methods.jk_iter_lin(array,cat,label,parallel=parallel)

    if plot:
      fig.plot_methods.plot_lin_split(arr1,e1,e2,e1err,e2err,m1,m2,b1,b2,cat,val,log=log,label=label)

    return arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err

  @staticmethod
  def split_gals_2pt_along(cat,cat2,val,mask=None,jkon=False,mock=False,log=False,plot=False):
    """
    Calculates xi and tangential shear for halves of catalog split along val. Optionally reweight each half by redshift distribution (cat.pzrw).
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    array=getattr(cat,val)

    if log:
      array=np.log10(array)

    s=config.lbl.get(val,None)
    if config.log_val.get(val,None):
      s='log '+s

    bins,w,edge=split_methods.get_mask_wnz(cat,array,cat.pz,mask=mask,label=s,plot=plot)
    return

    for j in xrange(2):
      if j==0:
        theta,out0,err0,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask)
        theta,outa,erra,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask&m,wa=w1)
        theta,outb,errb,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask&(~m),wa=w2)
      else:
        theta,out0,err0,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',maskb=mask,ran=False)
        theta,outa,erra,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',maskb=mask&m,wb=w1,ran=False)
        theta,outb,errb,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',maskb=mask&(~m),wb=w2,ran=False) 

      # if (jkon)&(cat.use_jk==1):
      #   me1err,me2err,slp1err,slp2err,b1err,b2err=jackknife_methods.lin_err0(array,cat,label,mask0=mask,parallel=parallel)
      # elif (jkon)&(cat.use_jk==2):
      #   me1err,me2err,slp1err,slp2err,b1err,b2err=BCC_Methods.jk_iter_lin(array,cat,label,parallel=parallel)

      a0=np.zeros((2))
      a1=np.zeros((2))
      a2=np.zeros((2))
      if ~jkon:
        derra=erra
        derr0=err0
        derrb=errb
        for i in xrange(2):
          if (j==1)&(i==1):
            continue
          a0[i]=split_methods.amp_shift(out0[i],outb[i]-outa[i],np.sqrt(erra[i]*errb[i]))
          a1[i]=split_methods.amp_shift(out0[i],outa[i]-out0[i],erra[i])
          a2[i]=split_methods.amp_shift(out0[i],outb[i]-out0[i],errb[i])

      if j==0:
        xi=(theta,outa,out0,outb,erra,err0,errb,derra,derr0,derrb,a1,a0,a2)
      else:
        gt=(theta,outa,out0,outb,erra,err0,errb,derra,derr0,derrb,a1,a0,a2)

    if plot:
      fig.plot_methods.plot_2pt_split(xi,gt,cat,val,split,log)

    return xi,gt,split

  @staticmethod
  def get_mask_wnz(cat,array,nz,mask=None,label='',plot=False):
    """
    Calculate splitting and redshift reweighting. 
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    w=np.ones(len(cat.coadd))
    edge=lin.linear_methods.find_bin_edges(array[mask],cat.sbins,w=cat.w[mask])
    bins=np.digitize(array[mask],edge)-1

    if cat.pzrw:
      w=split_methods.pz_weight(nz[mask],bins)
    else:
      w=np.ones(np.sum([mask]))

    if plot:
      fig.plot_methods.plot_pzrw(cat,nz,bins,w,label,edge)

    return bins,w,edge

  @staticmethod
  def pz_weight(nz,bins,binnum=200,pdf=False):
    """
    Reweight portions of galaxy population to match redshift distributions to that of the whole population.
    """

    w=np.ones(len(nz))
    if pdf:
      print 'transfer pdf support'
      return
    else:
      h0,b0=np.histogram(nz,bins=binnum)
      for j in xrange(np.max(bins)+1):
        binmask=(bins==j)
        h,b=np.histogram(nz[binmask],bins=b0)
        w[binmask]=1.*h0/h

      print np.sum(np.isnan(w)),np.sum(np.isinf(w))

    return w

  @staticmethod
  def amp_shift(xip,dxip,dxiperr):
    """
    Calculate amplitude shift for correlation splits.
    """

    chi2st=999999
    amp=0.
    for a in xrange(-200,200):
      chi2=np.sum((dxip-a*xip/100.)**2./dxiperr**2.)
      if chi2<chi2st:
        chi2st=chi2
        amp=a/100.

    return amp
