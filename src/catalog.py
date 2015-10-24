import numpy as np
import glob
import fitsio as fio
import pandas as pd
import config
import pickle

noval=999999

class CatalogStore(object):
  """
  A flexible class that reads and stores shear catalog information. An example is given for im3shape in the main testsuite.py.

  Use:

  :name:      str - A label for identifying the catalog in output. Multiple catalog classes can be defined, but each should have a unique name.
  :setup:     bool - Read in data and setup class structure? If false, only basic structure is created. Used for custom catalogs to feed into functions that expect a CatalogStore object.
  :cutfunc:   func - Pre-defined function for determining which objects in the catalog to store. See examples in CatalogMethods below.
  :cattype:   str - Currently one of 'im', 'ng', or 'gal' for im3shape, ngmix or a galaxy/lens catalog, respectively. This determines the dictionary in config.py used to match shortnames for quantities with column names in the catalog files.
  :cols:      [str] - Array of shortnames of columns to save from files. If None, saves all columns in dictionary.
  :catdir:    str - Directory to look for catalog files. Currently looks for files with *fits* extensions. Supports multiple files.
  :goldfile:  str - Common 'gold' file in combined, flatcat SV format for ngmix, im3shape. Not fully migrated from SV version of code. Will update when common file format for Y1 agreed upon.
  :catfile:   str - Point to a single catalog file. Also expects fits currently, but simple to add option for text files with column headers or h5 if useful.
  :ranfile:   str - Random file to be used with catalog type 'gal'. Currently rigid format, expecting numpy array of ra,dec.
  :jkbuild:   bool - Build jackknife regions from kmeans? Currently partly implemented - 150 regions. Can be more flexible if needed.
  :jkload:    bool - Load jackknife regions from file. Currently not fully migrated from SV code.
  :tiles:     [str] - Array of DES tile names to load, maching to multi-file fits format for catalog. Could be made more flexible to match strings in filenames generally.
  :release:   str - Release version, e.g. 'y1', 'sv'. Used for identifying release-specific files like systematics maps.

  """

  def __init__(self,name,setup=True,cutfunc=None,cattype=None,cols=None,catdir=None,goldfile=None,catfile=None,ranfile=None,jkbuild=False,jkload=False,tiles=None,release='y1'):

    if setup:

      if cattype=='i3':
        table=config.i3_col_lookup
      elif cattype=='ng':
        table=config.ng_col_lookup
      elif cattype=='gal':
        table=config.gal_col_lookup
      else:
        raise CatValError('No catalog type cattype specified.')

      if cutfunc is None:
        print 'Assuming no mask of catalog.'
        cutfunc=CatalogMethods.final_null_cuts()

      if cols is None:
        print 'Assuming all columns in dictionary.'
        cols=np.array(list(table.keys()))

      if goldfile!=None:
        print 'not ready to use flat catalog format in this version'
        if (i3file is None)|(ngfile is None):
          raise CatValError('Assumed flat catalog style and no im3shape or ngmix file specified.')

        cols1=[table.get(x,None) for x in cols]
        for i,x in enumerate(CatalogMethods.get_cat_cols_matched(catdir,goldfile,catfile,cols1,table,cuts,full=full,ext=ext)):
          setattr(self,cols[i],x)

      elif (catfile!=None)|(catdir!=None):
        if (catfile!=None)&(catdir!=None):
          raise CatValError('Both catfile and catdir specified.')
        if catdir is None:
          catdir=catfile
        else:
          catdir=catdir+'*fits*'

        cols1=[table.get(x,None) for x in cols]
        for i,x in enumerate(CatalogMethods.get_cat_cols(catdir,cols1,table,cutfunc,tiles)):
          setattr(self,cols[i],x)

      else:
        raise CatValError('Please specify the source of files: catfile, catdir, or goldfile/i3file/ngfile')

      if cattype in ['i3','ng']:
        if ('e1' in cols)&('e2' in cols):
          self.pos=0.5*np.arctan2(self.e2,self.e1)+np.pi/2.
          self.e=np.sqrt(self.e1**2.+self.e2**2.)
        if ('psf1' in cols)&('psf2' in cols):
          self.psfpos=0.5*np.arctan2(self.psf2,self.psf1)+np.pi/2.
          self.dpsf=self.psf1-self.psf2
        if ('psf1_exp' in cols)&('psf2_exp' in cols):
          self.psfpos=0.5*np.arctan2(self.psf2_exp,self.psf1_exp)+np.pi/2.
          self.dpsf=self.psf1_exp-self.psf2_exp
        if 'fluxfrac' in cols:
          self.invfluxfrac=1.001-self.fluxfrac
        if 'ccd' in cols:
          self.ccd-=1
        if 'like' in cols:
          self.nlike=-self.like

      self.name=name
      self.cat=cattype
      self.release=release
      #Number of bins in linear split functions (see e.g., sys_split module).
      self.lbins=config.cfg.get('lbins',None)
      #Number of bins to split signal in for systematics null tests.
      self.sbins=config.cfg.get('sbins',None)
      #Binslop for use in Treecorr.
      self.slop=config.cfg.get('slop',None)
      #Number of separation bins in Treecorr.
      self.tbins=config.cfg.get('tbins',None)
      #Number of cosebi bins (not migrated from SV yet)
      self.cbins=config.cfg.get('cbins',None)
      #Separation [min,max] for Treecorr.
      self.sep=config.cfg.get('sep',None)
      self.calc_err=False
      #Number of simulation patches for covariances.
      self.num_patch=config.cfg.get('num_patch',None)
      #Whether to use weighting and bias/sensitivity corrections in calculations.
      if cattype=='gal':
        self.bs=False
        self.wt=config.cfg.get('wt',None)
      else:
        self.bs=config.cfg.get('bs',None)
        self.wt=config.cfg.get('wt',None)
      #Whether to reweight n(z) when comparing subsets of data for null tests.
      self.pzrw=config.cfg.get('pzrw',None) 

      if (cattype=='gal')&(ranfile is not None):
        tmp=np.load(ranfile)
        self.ran_ra=tmp[:,0]
        self.ran_dec=tmp[:,1]
        mask=(self.ran_ra>60)&(self.ran_ra<95)&(self.ran_dec>-61)&(self.ran_dec<-42)
        self.ran_ra=self.ran_ra[mask]
        self.ran_dec=self.ran_dec[mask] 

      if jkbuild:
        X=np.vstack((self.ra,self.dec)).T
        km0 = km.kmeans_sample(X, config.cfg.get('num_reg',None), maxiter=100, tol=1.0e-5)
        self.regs=km0.labels
        #Number of jackknife regions
        self.num_reg=config.cfg.get('num_reg',None)
      else:
        self.num_reg=0
        self.regs=np.ones(len(self.coadd))


      # if jkload:
      #   self.num_reg,self.regs=jackknife_methods.load_jk_regs(self.ra,self.dec,64,jkdir)
      # else:
      #   self.num_reg=0
      #   self.regs=np.ones(len(self.coadd))


    else:

      self.name=name
      self.cat=cattype
      self.lbins=10
      self.sbins=2
      self.slop=0.1
      self.tbins=8
      self.cbins=5
      self.sep=np.array([1,400])
      self.calc_err=False
      self.num_patch=126
      self.bs=False
      self.wt=False
      self.pzrw=False 
      self.release=release

    return


class MockCatStore(object):
  """
  A flexible class that reads and stores mock catalog information. Used for covarianace calculations. Mock catalog module not yet migrated from SV code.
  """

  coadd=[]
  ra=[]
  dec=[]
  z=[]
  A00=[]
  A01=[]
  A10=[]
  A11=[]
  w=[]
  e1=[]
  e2=[]
  g1=[]
  g2=[]
  lbins=10
  cbins=5
  tbins=9
  sep=np.array([1,400])
  slop=0.15
  use_jk=False

  def __init__(self,setup=True,filenum=0,mocktype='sva1_gold'):

    if setup:
      for ifile,file in enumerate(glob.glob(config.mockdir+mocktype+'/*.fit')):
        if ifile==filenum:
          print 'simfile',ifile,file
          try:
            fits=fio.FITS(file)
            #tmparray=apt.Table.read(file)
          except IOError:
            print 'error loading fits file: ',file
            return
          tmparray=fits[-1].read()
          self.coadd=np.arange(len(tmparray['id']))
          self.ra=tmparray['ra']
          self.dec=tmparray['dec']
          self.z=tmparray['z']
          self.A00=tmparray['A00']
          self.A01=tmparray['A01']
          self.A10=tmparray['A10']
          self.A11=tmparray['A11']
          self.w=[]
          self.e1=[]
          self.e2=[]
          self.prop=[]
          self.g1=tmparray['g1']
          self.g2=tmparray['g2']
          break


class PZStore(object):
  """
  A flexible class that reads and stores photo-z catalog information. An example is given in the main testsuite.py.

  Use:

  :name:      str - A label for identifying the catalog in output. Multiple catalog classes can be defined, but each should have a unique name.
  :setup:     bool - Read in data and setup class structure? If false, only basic structure is created. Used for custom catalogs to feed into functions that expect a PZStore object.
  :pztype:    str - This identifies which photo-z code is stored in the class. In SV was used for dealing with multiple formats. Hopefully deprecated in Y1 analysis. For now identifies h5 table name.
  :filetype:  bool - Type of file to be read. dict - standard dict file for passing n(z) for spec validation, h5 - h5 file of pdfs, fits - fits file of pdfs. Non-dict file support not migrated from SV code yet - can use setup=False to store pdf information in class manually, though.
  :file:      str - File to be read in.
  """

  def __init__(self,name,setup=True,pztype='skynet',filetype=None,file=None):

    if setup:

      if file is None:
        raise CatValError('Need a source file.')

      if filetype=='dict':
        d=CatalogMethods.load_spec_test_file(file)
        self.tomo=len(d['phot'])
        self.boot=len(d['phot'][0])-1
        self.bins=len(d['phot'][0][0])
        self.bin=d['binning']
        self.binlow=d['binning']-(d['binning'][1]-d['binning'][0])/2.
        self.binhigh=d['binning']+(d['binning'][1]-d['binning'][0])/2.
        self.spec=np.zeros((self.tomo,self.bins))
        self.pz=np.zeros((self.tomo,self.bins))
        self.bootspec=np.zeros((self.tomo,self.boot,self.bins))
        self.bootpz=np.zeros((self.tomo,self.boot,self.bins))
        for i in xrange(self.tomo):
          self.spec[i,:]=d['spec'][i][0]
          self.pz[i,:]=d['phot'][i][0]
          for j in xrange(self.boot):
            self.bootspec[i,j,:]=d['spec'][i][j+1]
            self.bootpz[i,j,:]=d['phot'][i][j+1]
      elif filetype=='h5':
        store=pd.HDFStore(config.pzdir+file, mode='r')
        pz0=store[pztype]
        self.coadd=pz0.index.values
        zm0=pz0.z_mean.values
        self.z_mean_full=zm0.astype('float32')
        zm0=pz0.z_peak.values
        self.z_peak_full=zm0.astype('float32')
        pz0=pz0[['pdf_{}'.format(i) for i in xrange(200)]].values
        self.bins=200
        self.bin=(np.linspace(0.005, 1.8, 201)[1:] + np.linspace(0.005, 1.8, 201)[0:-1])/2.0
        self.pz_full=pz0.astype('float32')
        self.w=np.ones(len(self.z_mean_full))
        store.close()
      elif filetype=='fits':
        return

      self.pztype=pztype
      self.name=name

      #self.pz_full
      #self.z_mean_full
      #self.z_median_full
      #self.coadd

    else:

      self.pztype=pztype
      self.name=name

    return

class ColError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class CatValError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class CatalogMethods(object):

  @staticmethod
  def get_cat_cols(dir,cols,table,cuts,tiles=None,maxiter=999999,hdu=-1):

    import fitsio as fio
    import re

    noval=999999

    lenst=0
    for ifile,file in enumerate(glob.glob(dir)):
      if ifile>maxiter:
        break
      if hasattr(tiles, '__len__'):
        m=re.search('.*/(DES\d\d\d\d-\d\d\d\d).*',file)
        if m:
          if m.group(1) not in tiles:
            continue
      print 'file',ifile,file
      try:
        fits=fio.FITS(file)
      except IOError:
        print 'error loading fits file: ',file
        return

      tmparray=fits[hdu].read()
      lenst+=len(tmparray)
      print 'len',len(tmparray),lenst

      colex,colist=CatalogMethods.col_exists(cols,fits[hdu].get_colnames())
      if colex<1:
        raise ColError('columns '+colist+' do not exist in file: '+file)

      cutcols=[table.get(x,None) for x in cuts['col']]
      colex,colist=CatalogMethods.col_exists(cutcols,fits[hdu].get_colnames())
      if colex<1:
        raise ColError('cut columns '+colist+' do not exist in file: '+file)

      mask=np.array([])
      for icut,cut in enumerate(cuts): 
        mask=CatalogMethods.cuts_on_col(mask,tmparray,cutcols[icut],cut['min'],cut['eq'],cut['max'])

      tmparray=fits[hdu].read(columns=cols)

      if lenst==len(tmparray):
        array=tmparray[mask]
      else:
        array=np.append(array,tmparray[mask],axis=0)

      fits.close()

    return [(array[col]) for col in cols]



  # @staticmethod
  # def get_cat_cols_matched(dir,goldfile,catfile,cols,table,cuts,full=False,maxiter=999999,hdu=-1,ext='*.fits*'):

  #   print 'file',goldfile
  #   try:
  #     fits=fio.FITS(dir+goldfile)
  #   except IOError:
  #     print 'error loading fits file: ',goldfile
  #     return

  #   gold=fits[hdu].read()

  #   cols1=cols[[gold_col_lookup.get(x,False) for x in cols]]

  #   colex,colist=CatalogMethods.col_exists(cols1,fits[hdu].get_colnames())
  #   if colex<1:
  #     raise ColError('columns '+colist+' do not exist in file: '+file)
  #     print 'there are column name(s) not in file: ',file
  #     return

  #   print cuts

  #   colex,colist=CatalogMethods.col_exists(cuts['col'],fits[hdu].get_colnames())
  #   if colex<1:
  #     print 'there are cut column name(s) not in file: ',file
  #     raise ColError('cut columns '+colist+' do not exist in file: '+file)
  #     return

  #   mask=np.array([])
  #   for icut,cut in enumerate(cuts): 
  #     mask=CatalogMethods.cuts_on_col(mask,gold,cut['col'],cut['min'],cut['eq'],cut['max'])

  #   gold=fits[hdu].read(columns=cols1)
  #   gold=gold[mask]

  #   print 'file',catfile
  #   try:
  #     fits=fio.FITS(dir+catfile)
  #   except IOError:
  #     print 'error loading fits file: ',catfile
  #     return

  #   cat=fits[hdu].read()

  #   cols1=~cols1

  #   colex,colist=CatalogMethods.col_exists(cols1,fits[hdu].get_colnames())
  #   if colex<1:
  #     raise ColError('columns '+colist+' do not exist in file: '+file)
  #     print 'there are column name(s) not in file: ',file
  #     return

  #   colex,colist=CatalogMethods.col_exists(cuts['col'],fits[hdu].get_colnames())
  #   if colex<1:
  #     print 'there are cut column name(s) not in file: ',file
  #     raise ColError('cut columns '+colist+' do not exist in file: '+file)
  #     return

  #   mask=np.array([])
  #   for icut,cut in enumerate(cuts): 
  #     mask=CatalogMethods.cuts_on_col(mask,cat,cut['col'],cut['min'],cut['eq'],cut['max'])

  #   cat=fits[hdu].read(columns=cols1)
  #   cat=cat[mask]

  #   array=np.column_stack((gold,cat))

  #   if full:
  #     return np.copy(array)
  #   else:
  #     return [(array[col]) for col in cols]


  @staticmethod
  def col_exists(cols,colnames):

    colist=''
    exists=np.in1d(cols,colnames)
    for i,val in enumerate(exists):
      if val==0:
        colist+=' '+cols[i]

    return np.sum(exists)/len(cols),colist

  @staticmethod
  def cuts_on_col(mask,array,col,valmin,valeq,valmax):
    noval=999999

    if mask.size==0:
      mask=np.ones((len(array[col])), dtype=bool)

    if (valmin==noval) & (valmax==noval):
      if valeq==noval:
        print 'warning, no range or equality set in cut on column '+col
      else:
        mask=mask & (array[col]==valeq)
    elif (valmin!=noval) & (valmax!=noval):
      if valeq!=noval:
        print 'cannot have both equality and range cut on column '+col
      else:
        mask=mask & (valmin<array[col]) & (array[col]<valmax)
    elif (valmin!=noval):
      mask=mask & (valmin<array[col])
    else:
      mask=mask & (array[col]<valmax)
    return mask

  @staticmethod
  def col_view(array, cols):

    dtype2 = np.dtype({name:array.dtype.fields[name] for name in cols})
    return np.ndarray(array.shape, dtype2, array, 0, array.strides)

  @staticmethod
  def add_cut(cuts,col,min,eq,max):
    
    if cuts.size==0:
      cuts=np.zeros((1), dtype=[('col',np.str,20),('min',np.float64),('eq',np.float64),('max',np.float64)])
      cuts[0]['col']=col
      cuts[0]['min']=min
      cuts[0]['eq']=eq
      cuts[0]['max']=max
    else:
      cuts0=np.zeros((1), dtype=[('col',np.str,20),('min',np.float64),('eq',np.float64),('max',np.float64)])
      cuts0[0]['col']=col
      cuts0[0]['min']=min
      cuts0[0]['eq']=eq
      cuts0[0]['max']=max
      cuts=np.append(cuts,cuts0,axis=0)

    return cuts

  @staticmethod
  def load_spec_test_file(file):

    f=open(file, 'r')
    d=pickle.load(f)
    f.close()

    return d

  @staticmethod
  def sort(a1,a2):
    """
    Sorts and matches two arrays of unique object ids (in DES this is coadd_objects_id).

    len(a1[mask1])==len(a2[mask2])
    (a1[mask1])[sort1]==a2[mask2]
    a1[mask1]==(a2[mask2])[sort2]

    """

    mask1=np.in1d(a1,a2,assume_unique=True)
    mask2=np.in1d(a2,a1,assume_unique=True)
    sort1=np.argsort(a1[mask1])[np.argsort(np.argsort(a2[mask2]))]
    sort2=np.argsort(a2[mask2])[np.argsort(np.argsort(a1[mask1]))]

    return mask1,sort1,mask2,sort2

  @staticmethod
  def get_new_nbcw(cat,file,w=True,prune=False):
    """
    Convenience function to match and add new nbc values from a fits file to a CatalogeStore object. Takes a CatalogStore object to modify.
    """

    fits=fio.FITS(file)
    tmp=fits[-1].read()

    m1,s1,m2,s2=CatalogMethods.sort(cat.coadd,tmp['coadd_objects_id'])
    if prune:
      CatalogMethods.match_cat(cat,m1)
    elif np.sum(m1)<len(cat.coadd):
      print 'missing ids in file'
      return

    tmp=tmp[m2]
    tmp=tmp[s2]

    cat.c1=tmp['c1']
    cat.c2=tmp['c2']
    cat.m=tmp['m']
    if w:
      cat.w=tmp['w']

    return

  @staticmethod
  def final_null_cuts():
    """
    Masking functions for use in CatalogStore initialisation. 

    Use:

    Each entry of CatalogMethods.add_cut(array,col,a,b,c) adds to array a structured definition of the mask to apply for a given column in the catalog, col. a,b,c are limiting values. If be is set, value in column must be equal to b. Otherwise it must be greater than a and/or less than c.
    """

    cuts=CatalogMethods.add_cut(np.array([]),'coadd',0,noval,noval)

    return cuts

  @staticmethod
  def redmagic():
    """
    Masking functions for use in CatalogStore initialisation. 

    Use:

    Each entry of CatalogMethods.add_cut(array,col,a,b,c) adds to array a structured definition of the mask to apply for a given column in the catalog, col. a,b,c are limiting values. If be is set, value in column must be equal to b. Otherwise it must be greater than a and/or less than c.
    """

    cuts=CatalogMethods.add_cut(np.array([]),'coadd',0,noval,noval)

    return cuts

  @staticmethod
  def default_ngmix009_cuts():
    """
    Masking functions for use in CatalogStore initialisation. 

    Use:

    Each entry of CatalogMethods.add_cut(array,col,a,b,c) adds to array a structured definition of the mask to apply for a given column in the catalog, col. a,b,c are limiting values. If be is set, value in column must be equal to b. Otherwise it must be greater than a and/or less than c.
    """

    cuts=CatalogMethods.add_cut(np.array([]),'exp_flags',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'flags_i',noval,noval,4)
    cuts=CatalogMethods.add_cut(cuts,'exp_arate',0.4,noval,0.6)
    cuts=CatalogMethods.add_cut(cuts,'exp_s2n_w',10,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'exp_t_s2n',3,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',60.,noval,95.)
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)

    return cuts

  @staticmethod
  def default_im3shapev8_cuts():
    """
    Masking functions for use in CatalogStore initialisation. 

    Use:

    Each entry of CatalogMethods.add_cut(array,col,a,b,c) adds to array a structured definition of the mask to apply for a given column in the catalog, col. a,b,c are limiting values. If be is set, value in column must be equal to b. Otherwise it must be greater than a and/or less than c.
    """

    cuts=CatalogMethods.add_cut(np.array([]),'error_flag',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'info_flag',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',60.,noval,95.)
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)
    # Add gold match cut...

    return cuts

  @staticmethod
  def default_rm_cuts():
    """
    Masking functions for use in CatalogStore initialisation. 

    Use:

    Each entry of CatalogMethods.add_cut(array,col,a,b,c) adds to array a structured definition of the mask to apply for a given column in the catalog, col. a,b,c are limiting values. If be is set, value in column must be equal to b. Otherwise it must be greater than a and/or less than c.
    """

    cuts=CatalogMethods.add_cut(np.array([]),'error',noval,0,noval)
    cuts=CatalogMethods.add_cut(cuts,'l',0,noval,noval)
    cuts=CatalogMethods.add_cut(cuts,'ra',60.,noval,95.)
    cuts=CatalogMethods.add_cut(cuts,'dec',-61,noval,-42.)
    # Add gold match cut...

    return cuts

  @staticmethod
  def ngmix_weight_calc(cat,sn=0.16):
    """
    Calculates ngmix weight columns from covariance. Takes a CatalogStore object to modify.
    """

    w=1./(2.*sn**2.+cat.cov11+cat.cov22+2.*cat.cov12)
    print w[w<0]

    return w


  @staticmethod
  def match_cat(cat,mask):
    """
    Takes a CatalogStore object to modify. Masks all arrays of length the array of unique object ids - cat.coadd.
    """

    for x in dir(cat):
      obj = getattr(cat,x)
      if isinstance(obj,np.ndarray):
        if len(obj)==len(cat.coadd):
          setattr(cat,x,obj[mask])

    return

  @staticmethod
  def check_mask(array,mask):
    """
    Convenience function to return true array for mask if None supplied in other functions.
    """

    if mask is None:
      return np.ones(len(array)).astype(bool)
    else:
      return mask


  @staticmethod
  def info_flag(cat):
    """
    Takes properly constructed im3shape CatalogStore object and adds info_flag values.
    """

    import healpy as hp
    gdmask=hp.read_map(config.goldir+'y1a1_gold_1.0.1_wide_footprint_4096.fit')
    badmask=hp.read_map(config.goldir+'y1a1_gold_1.0.1_wide_badmask_4096.fit')

    pix=hp.ang2pix(4096, np.pi/2.-np.radians(i3.dec),np.radians(i3.ra), nest=False)
    i3.gold_mask=(gdmask[pix] >=1)
    i3.gold_flag=badmask[pix]

    info_cuts =[
        'i3.gold_mask==False',
        'i3.gold_flag>0',
        'i3.modest!=1',
        'i3.maskfrac>.75', #0.75
        'i3.evals>10000',
        'i3.flagr==1',
        'i3.flagr==2',
        'i3.fluxfrac<.75', #.7
        'i3.snr<10.', #10
        'i3.snr>10000.', #10000
        'i3.rgp<1.1', #1.1
        'i3.rgp>3.5', #3.5
        'i3.rad>5', #10
        'i3.rad<.1', 
        'np.sqrt(i3.ra_off**2.+i3.dec_off**2.)>1', #1
        'i3.chi2pix<.5', #.8
        'i3.chi2pix>1.5', #2
        'i3.resmin<-0.2',#-2
        'i3.resmax>0.2',#2
        'i3.psffwhm>7',
        'i3.psffwhm<0',
        'i3.error!=0'
        # model image?
    ]

    i3.info = np.zeros(len(i3.coadd), dtype=np.int64)
    for i,cut in enumerate(info_cuts):
      mask=eval(cut).astype(int)
      print i,cut,np.sum(mask)
      j=1<<i
      flags=mask*j
      i3.info|=flags


    return

  @staticmethod
  def nbc_sv(cat):
    """
    Constructs nbc values for im3shape based on SV greatDES simulation results. Takes CatalogStore object.
    """

    def basis_m(snr,rgp):
      snr1=snr/100.
      rgp1=(rgp-1.)/10.

      func=np.zeros((len(snr1),18))
      func[:,0]=1/snr1**2*1/rgp1**2
      func[:,1]=1/snr1**3*1/rgp1**3
      func[:,2]=1/snr1**3*1/rgp1**2
      func[:,3]=1/snr1**2*1/rgp1**3
      func[:,4]=1/snr1**4*1/rgp1**3
      func[:,5]=1/snr1**4*1/rgp1**4
      func[:,6]=1/snr1**3*1/rgp1**4
      func[:,7]=1/snr1**2.5*1/rgp1**2.5
      func[:,8]=1/snr1**2.5*1/rgp1**3
      func[:,9]=1/snr1**3*1/rgp1**2.5
      func[:,10]=1/snr1**1.5*1/rgp1**1.5
      func[:,11]=1/snr1**2.5*1/rgp1**1.5
      func[:,12]=1/snr1**1.5*1/rgp1**2.5
      func[:,13]=1/snr1**1.5*1/rgp1**2
      func[:,14]=1/snr1**2*1/rgp1**1.5
      func[:,15]=1/snr1**1.25*1/rgp1**1.75
      func[:,16]=1/snr1**1.75*1/rgp1**1.25
      func[:,17]=1/snr1**4*1/rgp1**4
      return func

    def basis_a(snr,rgp):
      snr1=snr/100.
      rgp1=(rgp-1.)/10.

      func=np.zeros((len(snr1),18))
      func[:,0] =1/snr1**2*1/rgp1**2
      func[:,1] =1/snr1**3*1/rgp1**3
      func[:,2] =1/snr1**3*1/rgp1**2
      func[:,3] =1/snr1**2*1/rgp1**3
      func[:,4] =1/snr1**4*1/rgp1**3
      func[:,5] =1/snr1**4*1/rgp1**4
      func[:,6] =1/snr1**3*1/rgp1**4
      func[:,7] =1/snr1**2.5*1/rgp1**2.5
      func[:,8] =1/snr1**2.5*1/rgp1**3
      func[:,9] =1/snr1**3* 1/rgp1**2.5
      func[:,10]=1/snr1**1.5*1/rgp1**1.5
      func[:,11]=1/snr1**2.5*1/rgp1**1.5
      func[:,12]=1/snr1**1.5*1/rgp1**2.5
      func[:,13]=1/snr1**3*1/rgp1**5
      func[:,14]=1/snr1**5*1/rgp1**3
      func[:,15]=1/snr1**5*1/rgp1**5
      func[:,16]=1/snr1**5*1/rgp1**4
      func[:,17]=1/snr1**4*1/rgp1**5
      return func

    wm=np.array([-1.05e-03,1.47e-06,8.10e-05,6.73e-06,0.00e+00,0.00e+00,0.00e+00,6.56e-05,-6.16e-06,-2.09e-05,-7.63e-03,-1.37e-03,-1.08e-04,1.63e-03,5.37e-03,1.63e-04,2.28e-03,-2.73e-11])
    wa=np.array([-1.67283612e-04, 1.09715332e-06, 5.95801408e-05, 6.39015150e-07, 2.97121531e-08, -3.60228146e-10, 4.73608639e-09, 4.05443791e-05, -3.52379986e-06, -1.95627195e-05, 8.97549111e-04, -3.23420375e-04, -1.91942923e-06, -6.57971727e-12, -1.41012000e-09, 1.61504257e-15, 2.36381064e-11, -1.76498862e-12])

    b=basis_m(cat.snr,cat.rgp)
    cat.m = np.dot(b, wm)
    b=basis_a(cat.snr,cat.rgp)
    a = np.dot(b, wa)
    cat.c1=a*cat.psf1
    cat.c2=a*cat.psf2

    return

  @staticmethod
  def write_output(i3,indir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/spte_sv_v1_partial/bord/main_cats/',outdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/tmp/'):
    """
    Matches extra columns to CatalogStore object (i3) and rewrites info to fits files it was built from.
    """

    import glob
    for ifile,file in enumerate(glob.glob(indir+'*')):
      print 'file',ifile,file
      fits=fio.FITS(outdir+file[82:],'rw')
      tmp=fits[-1].read()
      m1,s1,m2,s2=CatalogMethods.sort(i3.coadd,tmp['coadd_objects_id'])

      fits[-1].insert_column('mag_auto_g', (i3.g[m1])[s1])
      fits[-1].insert_column('mag_auto_r', (i3.r[m1])[s1])
      fits[-1].insert_column('mag_auto_i', (i3.i[m1])[s1])
      fits[-1].insert_column('mag_auto_z', (i3.z[m1])[s1])
      fits[-1].insert_column('desdm_pz', (i3.pz[m1])[s1])
      fits[-1].insert_column('modest', (i3.modest[m1])[s1])
      fits[-1].insert_column('gold_mask', (i3.gold_mask[m1])[s1].astype(int))
      fits[-1].insert_column('gold_flag', (i3.gold_flag[m1])[s1].astype(int))
      tmp=fits[-1].read()
      tmp['info_flag']=(i3.info[m1])[s1]
      fits.write(tmp)
      fits.close()

    return

  @staticmethod
  def write_output2(i3,indir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/spte_sv_v1_partial/bord/main_cats/',outdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/tmp/'):
    """
    Matches extra columns to CatalogStore object (i3) and rewrites info to fits files it was built from.
    """

    import glob
    for ifile,file in enumerate(glob.glob(indir+'*')):
      print 'file',ifile,file
      fits=fio.FITS(outdir+file[82:],'rw')
      tmp=fits[-1].read()
      m1,s1,m2,s2=CatalogMethods.sort(i3.coadd,tmp['coadd_objects_id'])

      fits[-1].insert_column('m', (i3.m[m1])[s1])
      fits[-1].insert_column('c1', (i3.c2[m1])[s1])
      fits[-1].insert_column('c2', (i3.c1[m1])[s1])
      fits.close()
      
    return
