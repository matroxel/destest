import numpy as np
from scipy.optimize import curve_fit

import catalog
import config
import fig
import txt

class UseError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class linear_methods(object):
  """
  Linear statistics and support functions for modifying shear.
  """

  @staticmethod
  def get_lin_e_w_ms(cat,xi=False,mock=False,mask=None,w1=None):
    """
    For a general CatalogStore object cat, return properly corrected ellipticity and weight, m/s values. Used in many functions.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if mock:

      return cat.e1[mask], cat.e2[mask], cat.w[mask], np.ones(np.sum(mask))

    elif cat.cat=='gal':

      if cat.wt:
        return None, None, cat.w[mask], None
      else:
        return None, None, np.ones(np.sum(mask)), None

    elif cat.cat==None:

      print 'None catalog type. Assuming no e, bias/sensitivity corections.'
      if cat.wt:
        return np.ones(np.sum(mask)), np.ones(np.sum(mask)), cat.w[mask], np.ones(np.sum(mask))
      else:
        return np.ones(np.sum(mask)), np.ones(np.sum(mask)), np.ones(np.sum(mask)), np.ones(np.sum(mask))

    elif cat.cat not in ['i3','ng']:

      print 'Unknown catalog type. Assuming no bias/sensitivity corections.'
      if cat.wt:
        return cat.e1[mask], cat.e2[mask], cat.w[mask], np.ones(np.sum(mask))
      else:
        return cat.e1[mask], cat.e2[mask], np.ones(np.sum(mask)), np.ones(np.sum(mask))

    elif cat.cat=='i3':

      if cat.bs:
        e1=cat.e1-cat.c1
        e2=cat.e2-cat.c2
        ms=cat.m+1.
      else:
        e1=cat.e1
        e2=cat.e2
        ms=np.ones(len(cat.coadd))     

    elif cat.cat=='ng':

      e1=cat.e1
      e2=cat.e2
      if cat.bs:
        if xi:
          ms=np.sqrt(cat.s1*cat.s2)
        else:
          ms=(cat.s1+cat.s2)/2.
      else:  
        ms=np.ones(len(cat.coadd))

    if cat.wt:
      w=cat.w
    else:
      w=np.ones(len(cat.coadd))

    if w1 is not None:
      w=np.sqrt(w*w1)

    return e1[mask],e2[mask],w[mask],ms[mask]

  @staticmethod
  def calc_mean_stdev_rms_e(cat,mask=None,mock=False,full=True):
    """
    For a general CatalogStore object cat, return mean, std dev, rms of ellipticities.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    e1,e2,w,ms=linear_methods.get_lin_e_w_ms(cat,mock=mock,mask=mask)
    wms=np.sum(w*ms)
    ww=np.sum(w**2)
    mean1=np.sum(w*e1)/wms
    mean2=np.sum(w*e2)/wms
    if full:
      std1=np.sqrt(np.sum(w*(e1-mean1)**2)/wms)
      std2=np.sqrt(np.sum(w*(e2-mean2)**2)/wms)
      rms1=np.sqrt(np.sum((w*e1)**2)/ww)
      rms2=np.sqrt(np.sum((w*e2)**2)/ww)

      return mean1,mean2,std1,std2,rms1,rms2
    else:
      return mean1,mean2

    return

  @staticmethod
  def calc_mean_stdev_rms(cat,x,mask=None,mock=False,full=True):
    """
    For a general CatalogStore object cat, return mean, std dev, rms of array x in cat.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if isinstance(x,str):
      x1=getattr(cat,x)
    else:
      x1=x

    w=linear_methods.get_lin_e_w_ms(cat,mock=mock,mask=mask)[2]

    mean=np.sum(w*x1[mask])/np.sum(w)
    if full:
      std=np.sqrt(np.sum(w*(x1[mask]-mean)**2)/np.sum(w))
      rms=np.sqrt(np.sum((w*x1[mask])**2)/np.sum(w**2))
      return mean,std,rms
    else:
      return mean

    return


  @staticmethod    
  def find_bin_edges(x,nbins,w=None):
    """
    For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins.
    """

    if w is None:
      xs=np.sort(x)
      r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
      return xs[r.astype(int)]

    fail=False
    ww=np.sum(w)/nbins
    i=np.argsort(x)
    k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
    k=k.astype(int)
    r=np.zeros((nbins+1))
    ist=0
    for j in xrange(1,nbins):
      print k[j],r[j-1]
      if k[j]<r[j-1]:
        print 'Random weight approx. failed - attempting brute force approach'
        fail=True
        break
      w0=np.sum(w[i[ist:k[j]]])
      if w0<=ww:
        for l in xrange(k[j],len(x)):
          w0+=w[i[l]]
          if w0>ww:
            r[j]=x[i[l]]
            ist=l
            break
      else:
        for l in xrange(k[j],0,-1):
          w0-=w[i[l]]
          if w0<ww:
            r[j]=x[i[l]]
            ist=l
            break

    if fail:

      ist=np.zeros((nbins+1))
      ist[0]=0
      for j in xrange(1,nbins):
        wsum=0.
        for k in xrange(ist[j-1].astype(int),len(x)):
          wsum+=w[i[k]]
          if wsum>ww:
            r[j]=x[i[k-1]]
            ist[j]=k
            break

    r[0]=x[i[0]]
    r[-1]=x[i[-1]]

    return r

  @staticmethod
  def binned_mean_e(bin,cat,mask=None,mock=False):

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    y_mean1=[]
    y_err1=[]
    y_mean2=[]
    y_err2=[]

    for i in xrange(cat.lbins):
      mask0=(bin==i)&mask
      mean1,mean2,std1,std2,rms1,rms2=linear_methods.calc_mean_stdev_rms_e(cat,mask0,mock=mock)
      y_mean1=np.append(y_mean1,mean1)
      y_err1=np.append(y_err1,std1/np.sqrt(np.sum(mask0)))
      y_mean2=np.append(y_mean2,mean2)
      y_err2=np.append(y_err2,std2/np.sqrt(np.sum(mask0)))

    return y_mean1,y_err1,y_mean2,y_err2

  @staticmethod
  def binned_mean_x(bin,x,cat,mask=None,mock=False):

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    x_mean=[]
    x_err=[]

    for i in xrange(cat.lbins):
      mask0=(bin==i)&mask
      mean,std,rms=linear_methods.calc_mean_stdev_rms(cat,x,mask0,mock=mock)
      x_mean=np.append(x_mean,mean)
      x_err=np.append(x_err,std/np.sqrt(np.sum(mask0)))

    return x_mean,x_err


  @staticmethod
  def bin_means(x,cat,mask=None,mock=False,log=False,noe=False,y=None):
    """
    For array x in CatalogStore object cat, calculate the means of shear in equal bins of x. Returns errors in both x and y directions.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cat.wt:
      edge=linear_methods.find_bin_edges(x[mask],cat.lbins,cat.w[mask])
    else:
      edge=linear_methods.find_bin_edges(x[mask],cat.lbins)

    xbin=np.digitize(x,edge)-1

    x_mean,x_err=linear_methods.binned_mean_x(xbin,x,cat,mask,mock=mock)
    if noe:
      if y is None:
        return x_mean,x_err
      else:
        e1_mean,e1_err=linear_methods.binned_mean_x(xbin,y,cat,mask,mock=mock)
        e2_mean=np.zeros(len(e1_mean))
        e2_err=np.zeros(len(e1_mean))
    else:
      e1_mean,e1_err,e2_mean,e2_err=linear_methods.binned_mean_e(xbin,cat,mask,mock=mock)

    return x_mean,x_err,e1_mean,e1_err,e2_mean,e2_err

class fitting(object):

  @staticmethod
  def lin_fit(x,y,sig):
    """
    Find linear fit with errors for two arrays.
    """

    def func(x,m,b):
      return m*x+b

    params=curve_fit(func,x,y,p0=(0.,0.),sigma=sig)

    m,b=params[0]
    merr,berr=np.sqrt(np.diagonal(params[1]))

    return m,b,merr,berr

  @staticmethod
  def sys_lin_fit(x,cat,mask=None,log=False,noe=False,mock=False):
    """
    Find linear fit corresponding to mean shears in bins of x.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cat.wt:
      edge=linear_methods.find_bin_edges(x[mask],cat.lbins,cat.w[mask])
    else:
      edge=linear_methods.find_bin_edges(x[mask],cat.lbins)
    xbin=np.digitize(x,edge)-1

    x_mean,x_err=linear_methods.binned_mean_x(xbin,x,cat,mask,mock=mock)
    if cat.wt:
      w_mean,w_err=linear_methods.binned_mean_x(xbin,cat.w,cat,mask,mock=mock)
    else:
      w_mean,w_err=linear_methods.binned_mean_x(xbin,np.ones(len(cat.coadd)),cat,mask,mock=mock)
    e1_mean,e1_err,e2_mean,e2_err=linear_methods.binned_mean_e(xbin,cat,mask,mock=mock)

    m1,b1,m1err,b1err=fitting.lin_fit(x_mean,e1_mean,1./w_mean)
    m2,b2,m2err,b2err=fitting.lin_fit(x_mean,e2_mean,1./w_mean)

    return m1,m2,b1,b2,m1err,m2err,b1err,b2err


class hist(object):

  @staticmethod
  def hist_tests(vals,cat,mask=None):
    """
    Loop over array vals, containing stored catalog column variables in CatalogStore object cat. Optionally mask the elements used.

    Produces plots of 1D histograms for each element in vals.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    for x in vals:

      print 'hist',x

      x1=getattr(cat,x)[mask]
      if config.log_val.get(x,False):
        x1=np.log10(x1)

      fig.plot_methods.plot_hist(x1,name=cat.name,label=x)

    return

  @staticmethod
  def hist_comp_tests(vals,cat,cat2,mask=None,mask2=None,bins=50):
    """
    Loop over array vals, containing stored catalog column variables in CatalogStore object cat. Optionally mask the elements used.

    Produces plots of 1D histograms for each element in vals.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    mask2=catalog.CatalogMethods.check_mask(cat2.coadd,mask2)

    for x in vals:

      print 'hist',x

      x1=getattr(cat,x)[mask]
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      x2=getattr(cat2,x)[mask2]
      if config.log_val.get(x,False):
        x2=np.log10(x2)

      fig.plot_methods.plot_comp_hist(x1,x2,name=cat.name,name2=cat2.name,bins=bins,label=x)

    return

  @staticmethod
  def hist_2D_tests(valsx,valsy,cat,cat2=None,mask=None,mask2=None,match_col=None):
    """
    Loop over array vals(x|y), containing stored catalog column variables in CatalogStore object cat (and cat2). Optionally mask the elements used. 

    Produces plots of 2D histograms for each cross pair in valsx and valsy. If both cat and cat2 provided, optionally provide an array name in cat with which to match the two catalogs (default value of None indicates to use the coadd ids). This is useful for matching between data releases.
    """


    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cat2!=None:
      mask2=catalog.CatalogMethods.check_mask(cat2.coadd,mask2)

      if match_col is None:
        coadd=cat.coadd
      else:
        coadd=getattr(cat,match_col)
        mask=mask&(coadd>0)

      yname=cat2.name

      m1,s1,m2,s2=catalog.CatalogMethods.sort(coadd[mask],cat2.coadd[mask2])
    
    else:
      yname=cat.name

    for ix,x in enumerate(valsx):
      for iy,y in enumerate(valsy):   

        if cat2 is None:
          if (ix>=iy):
            continue
        else:
          if ix!=iy:
            continue

        print 'hist 2D',x,y

        x1=getattr(cat,x)[mask]
        if config.log_val.get(x,False):
          x1=np.log10(x1)

        if cat2 is None:
          y1=getattr(cat,y)[mask]
        else:
          y1=getattr(cat2,y)[mask2]

          x1=x1[m1]
          y1=(y1[m2])[s2]

        if config.log_val.get(y,False):
          y1=np.log10(y1)

        fig.plot_methods.plot_2D_hist(x1,y1,xname=cat.name,yname=yname,xlabel=x,ylabel=y)

    return

  @staticmethod
  def tile_tests(vals,cat):
    """
    Loop over array vals, containing stored catalog column variables in CatalogStore object cat that have had means computed tile-by-tile with summary_stats.tile_stats(). This expects the output file from that function to exist.

    Produces plots of 1D histograms for each mean value tile-by-tile.
    """

    s='S12,f8'+',f8'*len(vals)*3
    tmp=np.genfromtxt(txt.write_methods.get_file(cat,label='tile_stats'),names=True,dtype=s)

    for i,x in enumerate(vals):

      print 'tile hist',x

      x1=tmp[x]
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='mean')

      x1=tmp[x+'_std']
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='std')

      x1=tmp[x+'_rms']
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='rms')

    return

  @staticmethod
  def tile_tests_2D(valsx,valsy,cat,cat2=None):
    """
    Loop over array vals(x|y), containing stored catalog column variables in CatalogStore object cat (and cat2) that have had means computed tile-by-tile with summary_stats.tile_stats(). This expects the output file from that function to exist.

    Produces plots of 2D histograms for each cross pair in valsx and valsy tile-by-tile.
    """

    s='S12,f8'+',f8'*len(valsx)*3
    tmp=np.genfromtxt(txt.write_methods.get_file(cat,label='tile_stats'),names=True,dtype=s)
    if cat2!=None:
      s='S12,f8'+',f8'*len(valsy)*3
      tmp2=np.genfromtxt(txt.write_methods.get_file(cat2,label='tile_stats'),names=True,dtype=s)
      yname=cat2.name
    else:
      tmp2=tmp
      yname=cat.name

    for ix,x in enumerate(valsx):
      for iy,y in enumerate(valsy):

        if cat2 is None:
          if (ix>=iy):
            continue
        else:
          if ix!=iy:
            continue

        print 'tile hist 2D',x,y

        x1=tmp[x]
        if config.log_val.get(x,False):
          x1=np.log10(x1)
        y1=tmp2[y]
        if config.log_val.get(y,False):
          y1=np.log10(y1)

        fig.plot_methods.plot_2D_hist(x1,y1,bins=20,xname=cat.name,yname=yname,xlabel=x,ylabel=y,xtile='mean',ytile='mean')

    return

class footprint(object):

  @staticmethod
  def hexbin_tests(vals,cat,mask=None):
    """
    Produces a set of hexbin plots (mean value in cells across the survey footprint) for each value listed in vals and stored in CatalogeStore object cat. Optionally mask the arrays listed in vals.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    for x in vals:

      print 'hexbin',x
      x1=getattr(cat,x)[mask]

      fig.plot_methods.plot_hexbin(x1,cat,mask=mask,name=cat.name,label=x)

    return

  @staticmethod
  def tile_tests(vals,cat,mask=None):
    """
    A version of hexbin_tests that maps mean value tile-by-tile instead of by hexbin cell.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    s='S12,f8'+',f8'*len(vals)*3
    tmp=np.genfromtxt(txt.write_methods.get_file(cat,label='tile_stats'),names=True,dtype=s)

    for j,x in enumerate(vals):

      x1=np.zeros(len(cat.coadd))
      for i in xrange(len(tmp)):
        mask0=(cat.tile==tmp['tile'][i])
        x1[mask0]=np.ones(np.sum(mask0))*tmp[x][i]

      fig.plot_methods.plot_hexbin(x1[mask],cat,mask=mask,name=cat.name,label=x,tile='mean')

      x1=np.zeros(len(cat.coadd))
      for i in xrange(len(tmp)):
        mask0=(cat.tile==tmp['tile'][i])
        x1[mask0]=np.ones(np.sum(mask0))*tmp[x+'_std'][i]

      fig.plot_methods.plot_hexbin(x1[mask],cat,mask=mask,name=cat.name,label=x,tile='std')

      x1=np.zeros(len(cat.coadd))
      for i in xrange(len(tmp)):
        mask0=(cat.tile==tmp['tile'][i])
        x1[mask0]=np.ones(np.sum(mask0))*tmp[x+'_rms'][i]

      fig.plot_methods.plot_hexbin(x1[mask],cat,mask=mask,name=cat.name,label=x,tile='rms')

    return

  @staticmethod
  def footprint_tests(cat,vals,mask=None,bins=100,label='',cap=None):
    """
    If vals==[], produces a galaxy density plot over the survey fottprint for the catalog with mask. If vals contains a list of flag columns, it maps the density of objects that fail each flag value.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if vals==[]:
      fig.plot_methods.plot_footprint(cat,mask=mask,label=label,bins=bins,cap=cap)

    else:
      for val in vals:
        flag=getattr(cat,val)
        for i in xrange(summary_stats.n_bits_array(cat,val)):
          if np.sum((flag & 2**i) != 0)>1000:
            print 'footprint flag',val,i
            if (val=='info')&(i<21):
              continue
            fig.plot_methods.plot_footprint(cat,mask=((flag & 2**i) != 0)&mask,label=getattr(config,val+'_name').get(i),bins=bins,cap=cap)

    return


class summary_stats(object):

  @staticmethod
  def n_bits_array(cat,arr,mask=None):
    """
    Convenience function to return maximum flag bit.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    arr1=getattr(cat,arr)

    return len(bin(np.max(arr1[mask])))-2

  @staticmethod
  def i3_flags_dist(cat,mask=None):
    """
    Produce a summary of error properties in the catalog. Identifies nan or inf values in arrays and lists the distribution of objects that fail info/error flags. Needs to be updated to accept the name of the flag columns to generalise from im3shape.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)


    # txt.write_methods.heading('error flags',cat,label='flags_dist',create=True)
    # for i in xrange(summary_stats.n_bits_array(cat,'error')):
    #   txt.write_methods.write_append(str(i)+'  '+str(np.sum((cat.error & 2**i) != 0))+'  '+str(np.sum(cat.error == 2**i)),cat,label='flags_dist',create=False)

    # txt.write_methods.heading('info flags',cat,label='flags_dist',create=False)
    # for i in xrange(summary_stats.n_bits_array(cat,'info')):
    #   txt.write_methods.write_append(str(i)+'  '+str(np.sum((cat.info[mask] & 2**i) != 0))+'  '+str(np.sum(cat.info[mask] == 2**i)),cat,label='flags_dist',create=False)

    txt.write_methods.heading('checking for bad values',cat,label='flags_dist',create=False)

    for x in dir(cat):
      obj = getattr(cat,x)
      if isinstance(obj,np.ndarray):
        if len(obj)==len(cat.coadd):
          if isinstance(obj[0], str)|('__' in x):
            continue
          line=x
          if np.sum(np.isnan(obj[mask]))>0:
            line+='  !!!NANS!!!'
          if np.sum(np.isinf(obj[mask]))>0:
            line+='  !!!INF!!!'
          txt.write_methods.write_append(line,cat,label='flags_dist',create=False)

    return

  @staticmethod
  def e1_psf_stats(cat,mask=None):
    """
    Produce a summary of basic shear and psf properties in the catalog. Writes out the mean, std dev, and rms of e1, e2, psf e1, psf e2 and psf fwhm, along with number of galaxies for some mask. Needs to be updated to be more flexible about what information should be incldued in the summary (which columns).
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    txt.write_methods.heading('Summary',cat,label='summary',create=True)
    txt.write_methods.heading('num gals  '+str(np.sum(mask)),cat,label='summary',create=False)
    txt.write_methods.heading('e1, e2',cat,label='summary',create=False)

    mean1,mean2,std1,std2,rms1,rms2=linear_methods.calc_mean_stdev_rms_e(cat,mask)

    txt.write_methods.write_append('mean  '+str(mean1)+'  '+str(mean2),cat,label='summary',create=False)
    txt.write_methods.write_append('stdev  '+str(std1)+'  '+str(std2),cat,label='summary',create=False)
    txt.write_methods.write_append('rms  '+str(rms1)+'  '+str(rms2),cat,label='summary',create=False)

    txt.write_methods.heading('psf e1, psf e2, psf fwhm',cat,label='summary',create=False)

    mean1,std1,rms1=linear_methods.calc_mean_stdev_rms(cat,'psf1',mask)
    mean2,std2,rms2=linear_methods.calc_mean_stdev_rms(cat,'psf2',mask)
    mean3,std3,rms3=linear_methods.calc_mean_stdev_rms(cat,'psffwhm',mask)

    txt.write_methods.write_append('mean  '+str(mean1)+'  '+str(mean2)+'  '+str(mean3),cat,label='summary',create=False)
    txt.write_methods.write_append('stdev  '+str(std1)+'  '+str(std2)+'  '+str(std3),cat,label='summary',create=False)
    txt.write_methods.write_append('rms  '+str(rms1)+'  '+str(rms2)+'  '+str(rms3),cat,label='summary',create=False)

    return

  @staticmethod
  def tile_stats(cat,vals,mask=None):
    """
    Produce a summary of basic properties tile-by-tile in the catalog. Writes out the mean, std dev, and rms of the values listed in vals for some mask. This output file is used in the other tile functions in module lin.
    """
    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    line='#tile  numgal  '
    for val in vals:
      line+=val+'  '+val+'_std  '+val+'_rms  '

    txt.write_methods.write_append(line,cat,label='tile_stats',create=True)

    for i,x in enumerate(np.unique(cat.tile)):
      if np.sum(mask&(cat.tile==x))>0:
        line=x+'  '+str(np.sum(mask&(cat.tile==x)))+'  '

        if ('e1' in vals)|('e2' in vals):
          e1,e2,e1_std,e2_std,e1_rms,e2_rms=linear_methods.calc_mean_stdev_rms_e(cat,mask&(cat.tile==x))
          line+=str(np.around(e1,5))+'  '+str(np.around(e2,5))+'  '+str(np.around(e1_std,5))+'  '+str(np.around(e2_std,5))+'  '+str(np.around(e1_rms,5))+'  '+str(np.around(e2_rms,5))+'  '
        for val in vals:
          if val in ['e1','e2']:
            continue

          x1=linear_methods.calc_mean_stdev_rms(cat,val,mask&(cat.tile==x))
          line+=str(np.around(x1[0],5))+'  '+str(np.around(x1[1],5))+'  '+str(np.around(x1[2],5))+'  '

        txt.write_methods.write_append(line,cat,label='tile_stats',create=False)

    return

