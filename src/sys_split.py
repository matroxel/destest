import numpy as np
import fitsio as fio
import healpy as hp
import multiprocessing
import os
import numpy.random as ran

import catalog
import config
import fig
import txt
import lin
import corr


class split(object):

  @staticmethod
  def cat_splits_lin_e(cat,cols=None,mask=None,p=None):
    """
    Loop over array names in cols of CatalogStore cat with optional masking. Calls split_methods.split_gals_lin_along().
    """

    print '!!!!!!!!!!!!!!!!!!!! check if plotting at mean of log(val) or log(mean(val))...'

    if p is not None:
      jobs=[]
      p=multiprocessing.Pool(processes=config.cfg.get('proc',32),maxtasksperchild=config.cfg.get('task',None))

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask,p=p)
    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    txt.write_methods.heading('Linear Splits',cat,label='linear_splits',create=True)

    for val in cols:

      array=getattr(cat,val)
      name=fig.plot_methods.get_filename_str(cat)

      if p is not None:
        job=p.apply_async(split_gals_lin_along_base,[[cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name],{'log':config.log_val.get(val,False),'plot':True})
        jobs.append(job)
      else:
        tmp,tmp,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=split_gals_lin_along_base(cat,val,array,mask,name,log=config.log_val.get(val,False),plot=True)
        if arr1 is None:
          continue
        txt.write_methods.heading(val,cat,label='linear_splits',create=False)
        # txt.write_methods.write_append(x+'  '+str(arr1)+'  '+str(arr1err),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e  '+str(e1)+'  '+str(e2),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e err  '+str(e1err)+'  '+str(e2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope  '+str(m1)+'  '+str(m2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope err  '+str(m1err)+'  '+str(m2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept  '+str(b1)+'  '+str(b2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept err  '+str(b1err)+'  '+str(b2err),cat,label='linear_splits',create=False)

    if p is not None:
      for job in jobs:
        val,tmp,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=job.get()
        txt.write_methods.heading(val,cat,label='linear_splits',create=False)
        # txt.write_methods.write_append(x+'  '+str(arr1)+'  '+str(arr1err),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e  '+str(e1)+'  '+str(e2),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e err  '+str(e1err)+'  '+str(e2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope  '+str(m1)+'  '+str(m2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope err  '+str(m1err)+'  '+str(m2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept  '+str(b1)+'  '+str(b2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept err  '+str(b1err)+'  '+str(b2err),cat,label='linear_splits',create=False)

        p.close()
        p.join()
      
    return

  @staticmethod
  def cat_splits_lin_full(cat,cols=None,mask=None,p=None):
    """
    Loop over array names in cols of CatalogStore cat with optional masking. Calls split_methods.split_gals_lin_along().
    """

    if p is not None:
      jobs=[]
      p=multiprocessing.Pool(processes=config.cfg.get('proc',32),maxtasksperchild=config.cfg.get('task',None))

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    txt.write_methods.heading('Linear Splits Full',cat,label='linear_splits_full',create=True)

    for val in cols:
      for val2 in cols:
        if (val==val2)|(val in ['e1','e2'])|(val2 in ['e1','e2']):
          continue

        array=getattr(cat,val)
        array2=getattr(cat,val2)
        name=fig.plot_methods.get_filename_str(cat)

        if p is not None:
          job=p.apply_async(split_gals_lin_along_base,[[cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name],{'log':config.log_val.get(val,False),'log2':config.log_val.get(val2,False),'plot':True,'e':False,'val2':val2,'array2':array2})
          jobs.append(job)
        else:
          tmp,tmp,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=split_gals_lin_along_base([cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name,log=config.log_val.get(val,False),log2=config.log_val.get(val2,False),plot=True,e=False,val2=val2,array2=array2)
          txt.write_methods.heading(val+'  '+val2,cat,label='linear_splits',create=False)
          # txt.write_methods.write_append(x+'  '+str(arr1)+'  '+str(arr1err),cat,label='linear_splits',create=False)
          # txt.write_methods.write_append('e  '+str(e1)+'  '+str(e2),cat,label='linear_splits',create=False)
          # txt.write_methods.write_append('e err  '+str(e1err)+'  '+str(e2err),cat,label='linear_splits',create=False)
          txt.write_methods.write_append('slope  '+str(m1)+'  '+str(m2),cat,label='linear_splits',create=False)
          txt.write_methods.write_append('slope err  '+str(m1err)+'  '+str(m2err),cat,label='linear_splits',create=False)
          txt.write_methods.write_append('intercept  '+str(b1)+'  '+str(b2),cat,label='linear_splits',create=False)
          txt.write_methods.write_append('intercept err  '+str(b1err)+'  '+str(b2err),cat,label='linear_splits',create=False)

    if p is not None:
      for job in jobs:
        val,val2,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=job.get()
        txt.write_methods.heading(val+'  '+val2,cat,label='linear_splits',create=False)
        # txt.write_methods.write_append(x+'  '+str(arr1)+'  '+str(arr1err),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e  '+str(e1)+'  '+str(e2),cat,label='linear_splits',create=False)
        # txt.write_methods.write_append('e err  '+str(e1err)+'  '+str(e2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope  '+str(m1)+'  '+str(m2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('slope err  '+str(m1err)+'  '+str(m2err),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept  '+str(b1)+'  '+str(b2),cat,label='linear_splits',create=False)
        txt.write_methods.write_append('intercept err  '+str(b1err)+'  '+str(b2err),cat,label='linear_splits',create=False)

        p.close()
        p.join()
       
    return

  @staticmethod
  def cat_splits_2pt(cat,cat2=None,cols=None,mask=None):
    """
    Loop over array names in cols of CatalogStore cat with optional masking. Calls split_methods.split_gals_2pt_along().
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    txt.write_methods.heading('2pt Splits',cat,label='2pt_splits',create=True)

    for x in cols:

      txt.write_methods.heading(x,cat,label='2pt_splits',create=False)

      xip,xim,gt,split,edge=split_methods.split_gals_2pt_along(cat,cat2,x,mask=mask,jkon=False,mock=False,log=False,plot=True)

      txt.write_methods.write_append(x+'  '+str(np.min(getattr(cat,x)[mask]))+'  '+str(split)+'  '+str(np.max(getattr(cat,x)[mask])),cat,label='2pt_splits',create=False)

      txt.write_methods.write_append('theta  '+str(xip[0]),cat,label='2pt_splits',create=False)

      # txt.write_methods.write_append('xip  '+str(xi[2]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xip  '+str((xip[1][1]-xip[1][0])/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xip  '+str((xip[1][2]-xip[1][0])/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xip  '+str((xip[1][2]-xip[1][1])/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xip err  '+str(xip[2][1]/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xip err  '+str(xip[2][2]/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xip err  '+str(xip[2][0]/xip[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('amp xip  '+str(xip[4][0])+'  '+str(xip[4][1])+'  '+str(xip[4][2]),cat,label='2pt_splits',create=False)

      # txt.write_methods.write_append('xim  '+str(xi[2]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xim  '+str((xim[1][1]-xim[1][0])/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xim  '+str((xim[1][2]-xim[1][0])/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xim  '+str((xim[1][2]-xim[1][1])/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('lower delta xim err  '+str(xim[2][1]/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('upper delta xim err  '+str(xim[2][2]/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('delta xim err  '+str(xim[2][0]/xim[1][0]),cat,label='2pt_splits',create=False)
      txt.write_methods.write_append('amp xim  '+str(xim[4][0])+'  '+str(xim[4][1])+'  '+str(xim[4][2]),cat,label='2pt_splits',create=False)

      if cat2 is not None:
        # txt.write_methods.write_append('gt  '+str(xi[2]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('lower delta gt  '+str((gt[1][1]-gt[1][0])/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('upper delta gt  '+str((gt[1][2]-gt[1][0])/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('delta gt  '+str((gt[1][2]-gt[1][1])/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('lower delta gt err  '+str(gt[2][1]/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('upper delta gt err  '+str(gt[2][2]/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('delta gt err  '+str(gt[2][0]/gt[1][0]),cat,label='2pt_splits',create=False)
        txt.write_methods.write_append('amp gt  '+str(gt[4][0])+'  '+str(gt[4][1])+'  '+str(gt[4][2]),cat,label='2pt_splits',create=False)

    return

class split_methods(object):

  @staticmethod
  def split_gals_lin_along(cat,val,mask=None,mock=False,log=False,log2=False,label='',plot=False,fit=True,e=True,val2=None,trend=True):
    """
    Split catalog cat into cat.lbins equal (weighted) parts along val. Plots mean shear in bins and outputs slope and error information.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    array=getattr(cat,val)
    if val2 is not None:
      array2=getattr(cat,val2)

    name=fig.plot_methods.get_filename_str(cat)

    arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=split_gals_lin_along_base([cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name,mock=mock,log=log,log2=log2,label=label,plot=plot,fit=fit,e=e,val2=val2,array2=array2,trend=trend)
    return arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err

  @staticmethod
  def split_gals_2pt_along(cat,cat2,val,mask=None,jkon=False,mock=False,log=False,plot=False,blind=True):
    """
    Calculates xi and tangential shear for halves of catalog split along val. Optionally reweight each half by redshift distribution (cat.pzrw).
    """

    def get_a_st(cat, theta, out, err):

      a=[]
      derr=[]
      for i in xrange(cat.sbins+1):
        if ~jkon:
          derr.append(err[i])
      a.append(split_methods.amp_shift(out[0],out[2]-out[1],np.sqrt(err[1]*err[2])))
      a.append(split_methods.amp_shift(out[0],out[1]-out[0],err[1]))
      a.append(split_methods.amp_shift(out[0],out[2]-out[0],err[2]))

      if blind:
        r0=ran.rand()
        for i in xrange(cat.sbins+1):
          out[i]*=(r0*0.1+1.)
          err[i]*=(r0*0.1+1.)
          derr[i]*=(r0*0.1+1.)

      return (theta,out,err,derr,a)

    if cat.cat!='mcal':
      mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    else:
      mask=np.ones(len(cat.coadd)).astype(bool)

    array=getattr(cat,val)

    # if log:
    #   array=np.log10(array)

    s=config.lbl.get(val,val)
    if config.log_val.get(val,False):
      s='log '+s

    bins,w,edge=split_methods.get_mask_wnz(cat,array,val,cat.zp,mask=mask,label=val,plot=plot)
    print 'edge',edge,w

    print '00000000000000', -1
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG')
    xip=[out[0]]
    xiperr=[err[0]]
    xim=[out[1]]
    ximerr=[err[1]]
    print '---------'
    for i in xrange(cat.sbins):
      if cat.cat!='mcal':
        print 'why'
        theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask&(bins==i),wa=w)
        xip.append(out[0])
        xiperr.append(err[0])
        xim.append(out[1])
        ximerr.append(err[1])
      else:
        catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=False)
        print 'before',w
        print '00000000000000', i
        theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat,corr='GG',wa=w)
        catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=True)
        xip.append(out[0])
        xiperr.append(err[0])
        xim.append(out[1])
        ximerr.append(err[1])

    xip=get_a_st(cat, theta, xip, xiperr)
    xim=get_a_st(cat, theta, xim, ximerr)

    if cat2 is not None:
      print 'test in cat2?'
      theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',ran=False)    
      gt=[out[0]]
      gterr=[err[0]]
      for i in xrange(cat.sbins):
        if cat.cat!='mcal':
          theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',maskb=mask&(bins==0),wb=w,ran=False)
          gt.append(out[0])
          gterr.append(err[0])
        else:
          catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=False)
          theta,out,err,chi2=corr.xi_2pt.xi_2pt(cat2,catb=cat,corr='NG',wb=w,ran=False)
          gt.append(out[0])
          gterr.append(err[0])
          catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=True)

      gt=get_a_st(cat, theta, gt, gterr)
    else:
      gt=None

    # if (jkon)&(cat.use_jk==1):
    #   me1err,me2err,slp1err,slp2err,b1err,b2err=jackknife_methods.lin_err0(array,cat,label,mask0=mask,parallel=parallel)
    # elif (jkon)&(cat.use_jk==2):
    #   me1err,me2err,slp1err,slp2err,b1err,b2err=BCC_Methods.jk_iter_lin(array,cat,label,parallel=parallel)

    if plot:
      fig.plot_methods.plot_2pt_split(xip,xim,gt,cat,val,edge[1],log)

    return xip,xim,gt,split,edge

  @staticmethod
  def load_maps(cat,maps=None):
    """
    Load list of 'systematic' maps (see Boris et al) and match to positions in cat. If maps=None, load all in dictionary. See config.py.
    """

    if maps is None:
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
  def get_mask_wnz(cat,array,val,nz,mask=None,label='',plot=False):
    """
    Calculate splitting and redshift reweighting. 
    """

    if cat.cat!='mcal':
      mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    else:
      mask0=catalog.CatalogMethods.get_cuts_mask(cat,full=True)
      mask=mask0[0]

    w=np.ones(len(cat.coadd))
    if cat.wt:
      edge=lin.linear_methods.find_bin_edges(array[mask],cat.sbins,w=cat.w[mask])
    else:
      edge=lin.linear_methods.find_bin_edges(array[mask],cat.sbins)

    if cat.cat!='mcal':
      bins=np.digitize(array[mask],edge)-1
    else:
      bins=[]
      for i in range(cat.sbins):
        catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=False)
        bins.append(catalog.CatalogMethods.get_cuts_mask(cat,full=True))
        catalog.CatalogMethods.add_cut_sheared(cat,val,cmin=edge[i],cmax=edge[i+1],remove=True)

      print 'wnz',edge,bins,np.sum(bins[0]),np.sum(bins[1])

    if cat.pzrw:
      if cat.cat=='mcal':
        w=[]
        for i in range(5):
          w.append(split_methods.pz_weight(cat,nz,mask0[i],[bins[0][i],bins[1][i]]))
      else:
        w=split_methods.pz_weight(cat,nz,mask0,bins)
    else:
      if cat.cat!='mcal':
        w=np.ones(np.sum([mask]))
      else:
        w=[]
        for i in range(5):
          w.append(np.ones(len(nz)))

    print 'before fig',len(nz),len(w),len(bins),len(bins[0]),len(bins[1])

    if plot:
      if cat.cat=='mcal':
        fig.plot_methods.plot_pzrw(cat,nz,mask,[bins[0][0],bins[1][0]],w[0],label,edge)
      else:
        fig.plot_methods.plot_pzrw(cat,nz,mask,bins,w,label,edge)

    return bins,w,edge

  @staticmethod
  def pz_weight(cat,nz,mask,bins,binnum=100,pdf=False):
    """
    Reweight portions of galaxy population to match redshift distributions to that of the whole population.
    """

    if cat.cat!='mcal':
      print 'am i stupid?'
      nz=nz[mask]

    if pdf:
      print 'transfer pdf support'
      return
    else:
      h0,b0=np.histogram(nz,bins=binnum)
      for j in range(cat.sbins):
        if cat.cat!='mcal':
          binmask=(bins==j)
        else:
          binmask=bins[j]
          print 'binmask',j,binmask
        h,b=np.histogram(nz[binmask],bins=b0)
        w=np.ones(len(nz))
        print 'w0',len(w)
        for k in range(binnum):
          binmask2=(nz>b[k])&(nz<=b[k+1])
          if h[k]<0.01*h0[k]:
            w[binmask&binmask2]=0.
          else:
            w[binmask&binmask2]=0.5*h0[k]/h[k]

        print 'max/min/mean weight', k,np.max(w),np.min(w),np.mean(w[binmask])

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

def split_gals_lin_along_base(cat,val,array,mask,name,mock=False,log=False,log2=False,label='',plot=False,fit=True,e=True,val2=None,array2=None,trend=True):

  if cat.cat=='mcal':
    if not cat.tablesheared.get(val,False):
      print 'Not registered as sheared value:  ', val

  if val2 is not None:
    if log2:
      array2=np.log10(array2)

  if log:
    array=np.log10(array)

  print 'not passing w to bin_means - fix'

  if e:
    arr1,arr1err,e1,e1err,e2,e2err=lin.linear_methods.bin_means(array,val,cat,w=None,mask=mask,mock=mock,log=log)
  else:
    arr1,arr1err,e1,e1err,e2,e2err=lin.linear_methods.bin_means(array,val,cat,w=None,mask=mask,mock=mock,log=log,noe=True,y=array2)

  if arr1 is None:
    return None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None

  if fit:
    m1,b1,m1err,b1err=lin.fitting.lin_fit(arr1,e1,e1err)
    if e:
      m2,b2,m2err,b2err=lin.fitting.lin_fit(arr1,e2,e2err)
    else:
      m2,b2,m2err,b2err=0.,0.,0.,0.
  else:
    m1,m2,b1,b2,m1err,m2err,b1err,b2err=0.,0.,0.,0.,0.,0.,0.,0.

  if plot:
    if e:
      fig.plot_methods.plot_lin_split(arr1,e1,e2,e1err,e2err,m1,m2,b1,b2,name,val,log=log,label=label,e=True,val2=None,trend=trend)
    else:
      fig.plot_methods.plot_lin_split(arr1,e1,e2,e1err,e2err,m1,m2,b1,b2,name,val,log=log,label=label,e=False,val2=val2,trend=trend)

  return val,val2,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err
