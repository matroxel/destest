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
      # print k[j],r[j-1]
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

class hist(object):

  @staticmethod
  def hist_tests(cat,cat2=None,cols=None,mask=None):
    """
    Loop over array cols, containing stored catalog column variables in CatalogStore object cat. Optionally mask the elements used.

    Produces plots of 1D histograms for each element in cols.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)
    if cat2 is not None:
      mask2=catalog.CatalogMethods.check_mask(cat2.coadd,mask2)

    for val in cols:

      print 'hist',val

      x1=getattr(cat,val)[mask]
      if config.log_val.get(val,False):
        x1=np.log10(x1)

      if cat2 is None:
        if cat.wt:
          fig.plot_methods.plot_hist(x1,name=cat.name,label=val,w=cat.w[mask])
        else:
          fig.plot_methods.plot_hist(x1,name=cat.name,label=val)
      else:
        x2=getattr(cat2,val)[mask2]
        if config.log_val.get(val,False):
          x2=np.log10(x2)

        if cat.wt:
          fig.plot_methods.plot_comp_hist(x1,x2,name=cat.name,name2=cat2.name,bins=bins,label=x)
        else:
          fig.plot_methods.plot_comp_hist(x1,x2,name=cat.name,name2=cat2.name,bins=bins,label=x)

    return

  @staticmethod
  def hist_2D_tests(cat,cat2=None,colsx=None,colsy=None,mask=None,mask2=None,match_col=None):
    """
    Loop over array cols(x|y), containing stored catalog column variables in CatalogStore object cat (and cat2). Optionally mask the elements used. 

    Produces plots of 2D histograms for each cross pair in valsx and colsy. If both cat and cat2 provided, optionally provide an array name in cat with which to match the two catalogs (default value of None indicates to use the coadd ids). This is useful for matching between data releases.
    """


    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if colsx is None:
      colsx=catalog.CatalogMethods.get_cat_colnames(cat)


    if cat2!=None:
      mask2=catalog.CatalogMethods.check_mask(cat2.coadd,mask2)
      if colsy is None:
        colsy=catalog.CatalogMethods.get_cat_colnames(cat2)

      if match_col is None:
        coadd=cat.coadd
      else:
        coadd=getattr(cat,match_col)
        mask=mask&(coadd>0)

      yname=cat2.name

      s1,s2=catalog.CatalogMethods.sort2(coadd[mask],cat2.coadd[mask2])
    
    else:
      yname=cat.name
      colsy=colsx

    for ix,x in enumerate(colsx):
      for iy,y in enumerate(colsy):   

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

          x1=x1[s1]
          y1=y1[s2]

        if config.log_val.get(y,False):
          y1=np.log10(y1)

        fig.plot_methods.plot_2D_hist(x1,y1,xname=cat.name,yname=yname,xlabel=x,ylabel=y)

    return

  @staticmethod
  def tile_tests(cat):
    """
    Loop over cols that have had means computed tile-by-tile with summary_stats.tile_stats(). This expects the output from that function to exist.

    Produces plots of 1D histograms for each mean value tile-by-tile.
    """

    try:
      cat.tilecols
    except NameError:
      print 'you must first call lin.summary_stats.tile_stats(cat)'
      return

    for i,x in enumerate(cat.tilecols):

      print 'tile hist',x

      x1=cat.tilemean[:,i]
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='mean')

      x1=cat.tilestd[:,i]
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='std')

      x1=cat.tilerms[:,i]
      if config.log_val.get(x,False):
        x1=np.log10(x1)
      fig.plot_methods.plot_hist(x1,name=cat.name,label=x,bins=20,tile='rms')

    return

  @staticmethod
  def tile_tests_2D(cat,cat2=None):
    """
    Loop over cols that have had means computed tile-by-tile with summary_stats.tile_stats(). This expects the output from that function to exist.

    Produces plots of 2D histograms for each cross pair in valsx and valsy tile-by-tile.
    """

    try:
      cat.tilecols
    except NameError:
      print 'you must first call lin.summary_stats.tile_stats(cat)'
      return

    if cat2 is not None:
      try:
        cat2.tilecols
      except NameError:
        print 'you must first call lin.summary_stats.tile_stats(cat)'
        return

      if cat2.tilecols!=cat.tilecols:
        print 'tile col lists do not agree between cat and cat2'
        return
    else:
      cat2=cat

    for ix,x in enumerate(cat.tilecols):
      for iy,y in enumerate(cat2.tilecols):

        if cat2 is cat:
          if (ix>=iy):
            continue
        else:
          if ix!=iy:
            continue

        print 'tile hist 2D',x,y

        x1=cat.tilemean[:,ix]
        if config.log_val.get(x,False):
          x1=np.log10(x1)
        y1=cat2.tilemean[:,iy]
        if config.log_val.get(y,False):
          y1=np.log10(y1)

        fig.plot_methods.plot_2D_hist(x1,y1,bins=20,xname=cat.name,yname=cat2.name,xlabel=x,ylabel=y,xtile='mean',ytile='mean')

    return

class footprint(object):

  @staticmethod
  def hexbin_tests(cat,cols=None,mask=None):
    """
    Produces a set of hexbin plots (mean value in cells across the survey footprint) for each value listed in cols and stored in CatalogeStore object cat. Optionally mask the arrays listed in cols.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    for x in cols:

      print 'hexbin',x
      x1=getattr(cat,x)[mask]
      if config.log_val.get(x,False):
        x1=np.log10(x1)

      fig.plot_methods.plot_hexbin(x1,cat,mask=mask,name=cat.name,label=x)

    return

  @staticmethod
  def tile_tests(cat,mask=None):
    """
    A version of hexbin_tests that maps mean value tile-by-tile instead of by hexbin cell.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    try:
      cat.tilecols
    except NameError:
      print 'you must first call lin.summary_stats.tile_stats(cat)'
      return

    for j,x in enumerate(cat.tilecols):

      x1=np.zeros(np.sum(mask))
      for i in xrange(len(cat.tilelist)):
        mask0=(cat.tile[mask]==cat.tilelist[i])
        x0=cat.tilemean[i,j]
        if config.log_val.get(x,False):
          x0=np.log10(x0)
        x1[mask0]=np.ones(np.sum(mask0))*x0

      fig.plot_methods.plot_hexbin(x1,cat,mask=mask,name=cat.name,label=x,tile='mean')

      x1=np.zeros(np.sum(mask))
      for i in xrange(len(cat.tilelist)):
        mask0=(cat.tile[mask]==cat.tilelist[i])
        x0=cat.tilestd[i,j]
        if config.log_val.get(x,False):
          x0=np.log10(x0)
        x1[mask0]=np.ones(np.sum(mask0))*x0

      fig.plot_methods.plot_hexbin(x1,cat,mask=mask,name=cat.name,label=x,tile='std')

      x1=np.zeros(np.sum(mask))
      for i in xrange(len(cat.tilelist)):
        mask0=(cat.tile[mask]==cat.tilelist[i])
        x0=cat.tilerms[i,j]
        if config.log_val.get(x,False):
          x0=np.log10(x0)
        x1[mask0]=np.ones(np.sum(mask0))*x0

      fig.plot_methods.plot_hexbin(x1,cat,mask=mask,name=cat.name,label=x,tile='rms')

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
  def i3_flags_vals_check(cat,flags=['error','info'],mask=None):
    """
    Produce a summary of error properties in the catalog. Identifies nan or inf values in arrays and lists the distribution of objects that fail flags.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    for flag in flags:

      txt.write_methods.heading(flag+' flags',cat,label='flags_dist',create=True)
      for i in xrange(summary_stats.n_bits_array(cat,flag)):
        total=np.sum((cat.error & 2**i) != 0)
        unique=np.sum(cat.error == 2**i)
        txt.write_methods.write_append(str(i)+'  '+str(total)+'  '+str(unique)+'  '+str(np.around(1.*total/len(cat.coadd),5)),cat,label='flags_dist',create=False)

    txt.write_methods.heading('checking for bad values',cat,label='flags_dist',create=False)

    for x in dir(cat):
      obj = getattr(cat,x)
      if isinstance(obj,np.ndarray):
        if len(obj)==len(cat.coadd):
          if isinstance(obj[0], str)|('__' in x):
            continue
          line=x
          if np.isnan(obj[mask]).any():
            line+='  !!!NANS!!!'
          if np.isinf(obj[mask]).any():
            line+='  !!!INF!!!'
          txt.write_methods.write_append(line,cat,label='flags_dist',create=False)

    return

  @staticmethod
  def val_stats(cat,cols=None,mask=None):
    """
    Produce a summary of basic properties in the catalog. Writes out the mean, std dev, and rms of cols, along with number of galaxies for some optional mask.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    txt.write_methods.heading('Summary',cat,label='summary',create=True)
    txt.write_methods.heading('num gals  '+str(np.sum(mask)),cat,label='summary',create=False)
    txt.write_methods.heading('mean, std, rms, min, max',cat,label='summary',create=False)

    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    for val in cols:
      mean,std,rms=linear_methods.calc_mean_stdev_rms(cat,val,mask)
      txt.write_methods.write_append(val+'  '+str(mean)+'  '+str(std)+'  '+str(rms)+'  '+str(np.min(getattr(cat,val)))+'  '+str(np.max(getattr(cat,val))),cat,label='summary',create=False)

    return

  @staticmethod
  def tile_stats(cat,cols=None,mask=None):
    """
    Produce a summary of basic properties tile-by-tile in the catalog. Writes out the mean, std dev, and rms of the values listed in cols for some mask.
    """
    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if cols is None:
      cols=catalog.CatalogMethods.get_cat_colnames(cat)

    line='#tile  numgal  '
    for val in cols:
      line+=val+'  '+val+'_std  '+val+'_rms  '

    txt.write_methods.write_append(line,cat,label='tile_stats',create=True)

    cat.tilelist=np.unique(cat.tile)
    cat.tilecols=cols
    cat.tilenums=np.zeros(len(cat.tilelist))
    cat.tilemean=np.zeros((len(cat.tilelist),len(cols)))
    cat.tilestd=np.zeros((len(cat.tilelist),len(cols)))
    cat.tilerms=np.zeros((len(cat.tilelist),len(cols)))

    for i,x in enumerate(cat.tilelist):
      if np.sum(mask&(cat.tile==x))>0:
        line=x+'  '+str(np.sum(mask&(cat.tile==x)))+'  '
        cat.tilenums[i]=np.sum(mask&(cat.tile==x))

        for j,val in enumerate(cols):
          if val=='e1':
            e1,e2,e1_std,e2_std,e1_rms,e2_rms=linear_methods.calc_mean_stdev_rms_e(cat,mask&(cat.tile==x))
            line+=str(np.around(e1,5))+'  '+str(np.around(e1_std,5))+'  '+str(np.around(e1_rms,5))+'  '
            cat.tilemean[i,j]=e1
            cat.tilestd[i,j]=e1_std
            cat.tilerms[i,j]=e1_rms
          elif val=='e2':
            e1,e2,e1_std,e2_std,e1_rms,e2_rms=linear_methods.calc_mean_stdev_rms_e(cat,mask&(cat.tile==x))
            line+=str(np.around(e2,5))+'  '+str(np.around(e2_std,5))+'  '+str(np.around(e2_rms,5))+'  '
            cat.tilemean[i,j]=e2
            cat.tilestd[i,j]=e2_std
            cat.tilerms[i,j]=e2_rms
          else:
            x1=linear_methods.calc_mean_stdev_rms(cat,val,mask&(cat.tile==x))
            line+=str(np.around(x1[0],5))+'  '+str(np.around(x1[1],5))+'  '+str(np.around(x1[2],5))+'  '
            cat.tilemean[i,j]=x1[0]
            cat.tilestd[i,j]=x1[1]
            cat.tilerms[i,j]=x1[2]

        txt.write_methods.write_append(line,cat,label='tile_stats',create=False)

    return

