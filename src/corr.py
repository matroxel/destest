import numpy as np
import treecorr
import os.path

import catalog
import fig
import txt
import lin

class UseError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class xi_2pt(object):

  @staticmethod
  def xi_2pt(cata,catb=None,k=None,ga=None,gb=None,corr='GG',maska=None,maskb=None,wa=None,wb=None,ran=True,mock=False,erron=True,jkmask=None,label0='',plot=False,conj=False):
    """
    This is a flexible convenience wrapper for interaction with treecorr to work on CatalogStore objects. Some basic examples are given in corr_tests() of the main testsuite.py. g1, g2 correctly by c1, c2 if ellipticities and cat.bs is true. Correction by sensitivity, 1+m applied if cat.bs=True. Weighting applied if cat.wt is true. Other config properties for treecorr stored in CatalogStore object. See catalog.py or config.py. Not all correlation types fully integrated or tested. For example, only one kappa value is currently possible. Will be updated in future as useful.

    Use:

    :cata, catb:    CatalogStore - Must supply both cata, catb (can be same reference) if NG or NK correlation. Otherwise catb is optional.
    :k:             str - Array name in cata, catb to use for kappa correlation. 
    :ga, gb:        str - Array names for g1, g2 treecorr inputs. If None assume e1, e2.
    :corr:          str - Type of correlation for treecorr.
    :maska, maskb:  [bool] - Masking array to apply to input catalogs.
    :wa, wb:        [float] - Additional weights to apply after cat.w is used. Combined as e.g., w=sqrt(cat.w*wa).
    :ran:           bool - Use randoms in correlation calculation. If True, assumes cat.ran_ra, cat.ran_dec exist.
    :mock:          bool - If mock catalog from sims. Used when calculating covariances from sims, not currently migrated from SV code.
    :erron:         bool - Calculate jackknife or sim cov errors. If False, uses treecorr error outputs. Not currently migrated from SV code. When implemented requires cat.calc_err in ('jk', 'mock').
    :jkmask:        [bool] - For jk, mock cov calculation loop over regions/sims.
    :label0:        str - Additional (optional) label string used in some outputs.
    :plot:          bool - Plot output?

    Output (len cat.tbins):

    :theta:         [float] - Treecorr np.exp(meanlogr)
    :out:           ([float]x4) - Output of signal e.g., (xi+,xi-,xi+im,x-im). For correlations with only one xi output, (xi,0.,xi_im,0.).
    :err:           ([float]x4) - Same but for sqrt(var).
    :chi2:          ([float]x4) - Same but for chi^2 if using jk or sim covariance.
    :conj:          bool - Conjugate calculation.

    """

    maska=catalog.CatalogMethods.check_mask(cata.coadd,maska)
    jkmask=catalog.CatalogMethods.check_mask(cata.coadd,jkmask)

    maska0=maska&jkmask

    if wa is None:
      wa=np.ones(len(cata.coadd))

    e1,e2,w,ms=lin.linear_methods.get_lin_e_w_ms(cata,xi=True,mock=mock,mask=maska0,w1=wa)

    if catb is None:
      if corr not in ['GG','NN','KK']:
        raise UseError('Must supply both cata,catb for NG,NK correlations.')

    if ga is not None:
      e1=getattr(cata,ga+'1')[maska]
      e2=getattr(cata,ga+'2')[maska]
    else:
      ga='e'
    if catb is None:
      gb=ga
    if conj:
      e2=-e2

    if (corr=='GG')|((catb!=None)&(corr=='KG')):
      catxa=treecorr.Catalog(g1=e1, g2=e2, w=w, ra=cata.ra[maska0], dec=cata.dec[maska0], ra_units='deg', dec_units='deg')
      catma=treecorr.Catalog(k=ms, w=w, ra=cata.ra[maska0], dec=cata.dec[maska0], ra_units='deg', dec_units='deg')

    elif (corr=='NN')|((catb!=None)&(corr in ['NG','NK'])):
      catxa=treecorr.Catalog(w=w, ra=cata.ra[maska0], dec=cata.dec[maska0], ra_units='deg', dec_units='deg')
      if ran:
        catra=treecorr.Catalog(w=w, ra=cata.ran_ra[maska0], dec=cata.ran_dec[maska0], ra_units='deg', dec_units='deg')

    elif corr=='KK':
      if k is None:
        raise UseError('Must specify k for KK correlation.')
      if k not in dir(cata):
        raise UseError('Unknown k field specified.')
      catxa=treecorr.Catalog(k=getattr(cata, k)[maska0], w=w, ra=cata.ra[maska0], dec=cata.dec[maska0], ra_units='deg', dec_units='deg')

    if catb is not None:

      maskb=catalog.CatalogMethods.check_mask(catb.coadd,maskb)

      if wb is None:
        wb=np.ones(len(catb.coadd))

      e1,e2,w,ms=lin.linear_methods.get_lin_e_w_ms(catb,xi=True,mock=mock,mask=maskb,w1=wb)

      if gb is not None:
        e1=getattr(catb,gb+'1')[maskb]
        e2=getattr(catb,gb+'2')[maskb]
      else:
        gb='e'
      if conj:
        e2=-e2


      if corr in ['GG','NG','KG']:
        catxb=treecorr.Catalog(g1=e1, g2=e2, w=w, ra=catb.ra[maskb], dec=catb.dec[maskb], ra_units='deg', dec_units='deg')
        catmb=treecorr.Catalog(k=ms, w=w, ra=catb.ra[maskb], dec=catb.dec[maskb], ra_units='deg', dec_units='deg')
      elif corr=='NN':
        catxb=treecorr.Catalog(w=w, ra=catb.ra[maskb], dec=catb.dec[maskb], ra_units='deg', dec_units='deg')
        if ran:
          catrb=treecorr.Catalog(w=w, ra=catb.ran_ra[maskb], dec=catb.ran_dec[maskb], ra_units='deg', dec_units='deg')
      elif corr in ['KK','NK']:
        if k is None:
          raise UseError('Must specify k for KK correlation.')
        if k not in dir(catb):
          raise UseError('Unknown k field specified.')
        catxb=treecorr.Catalog(k=getattr(catb, k)[maskb], w=w, ra=catb.ra[maskb], dec=catb.dec[maskb], ra_units='deg', dec_units='deg')

    xim=None
    xip_im=None
    xim_im=None
    ximerr=None
    xiperr_im=None
    ximerr_im=None
    if corr=='GG':
      gg = treecorr.GGCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      kk = treecorr.KKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      if catb is None:
        gg.process(catxa)
        kk.process(catma)
      else:
        gg.process(catxa,catxb)
        kk.process(catma,catmb)

      xip = gg.xip/kk.xi
      xim = gg.xim/kk.xi
      xiperr = ximerr = np.sqrt(gg.varxi)
      xip_im = gg.xip_im/kk.xi
      xim_im = gg.xim_im/kk.xi
      theta = np.exp(gg.meanlogr)

    elif corr=='NN':
      nn = treecorr.NNCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      if ran:
        nr = treecorr.NNCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
        rr = treecorr.NNCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)

      if catb is None:
        nn.process(catxa)
        xip=nn.npairs
        xiperr=np.sqrt(nn.npairs)
        if ran:
          nr.process(catxa,catra)
          rr.process(catra)
        xip,xiperr=nn.calculateXi(rr,nr)
        xiperr=np.sqrt(xiperr)
      else:
        rn = treecorr.NNCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
        nn.process(catxa,catxb)
        xip=nn.npairs
        xiperr=np.sqrt(nn.npairs)
        if ran:
          nr.process(catxa,catrb)
          nr.process(catra,catxb)
          rr.process(catra,catrb)
        xip,xiperr=nn.calculateXi(rr,nr,rn)
        xiperr=np.sqrt(xiperr)
      theta=np.exp(nn.meanlogr)

    elif corr=='KK':

      kk = treecorr.KKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      if catb is None:
        kk.process(catxa)
      else:
        kk.process(catxa,catxb)
      xip=kk.xi
      xiperr=np.sqrt(kk.varxi)
      theta=np.exp(kk.meanlogr)

    elif corr=='KG':

      kg = treecorr.KGCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      kk = treecorr.KKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      kg.process(catxa,catxb)
      kk.process(catxa,catmb)
      xip=kg.xi/kk.xi
      xiperr=np.sqrt(kg.varxi)
      xip_im=kg.xi_im/kk.xi
      theta=np.exp(kg.meanlogr)

    elif corr=='NG':

      ng = treecorr.NGCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      nk = treecorr.NKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      ng.process(catxa,catxb)
      nk.process(catxa,catmb)
      xip=ng.xi/nk.xi
      xiperr=np.sqrt(ng.varxi)
      xip_im=ng.xi_im/nk.xi
      if ran:
        rg = treecorr.NGCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
        rk = treecorr.NKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
        rg.process(catra,catxb)
        rk.process(catra,catmb)
        xip,xip_im,xiperr=ng.calculateXi(rg)
        tmpa,tmp=nk.calculateXi(rk)
        if np.sum(tmpa)==0:
          tmpa=np.ones(len(xip))
        xip/=tmpa
        xiperr=np.sqrt(xiperr)
        xip_im/=tmpa
      theta=np.exp(ng.meanlogr)

    elif corr=='NK':

      nk = treecorr.NKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
      nk.process(catxa,catxb)
      xip=nk.xi
      xiperr=np.sqrt(nk.varxi)
      if ran:
        rk = treecorr.NKCorrelation(nbins=cata.tbins, min_sep=cata.sep[0], max_sep=cata.sep[1], sep_units='arcmin',bin_slop=cata.slop,verbose=0)
        rk.process(catra,catxb)
        xip,xiperr=nk.calculateXi(rk)
        xiperr=np.sqrt(xiperr)
      theta=np.exp(nk.meanlogr)

    out=[xip,xim,xip_im,xim_im]
    err=[xiperr,ximerr,xiperr,ximerr]
    chi2=[0.,0.,0.,0.]

    if erron:
      kwargs={'catb':catb,'k':k,'corr':corr,'maska':maska,'maskb':maskb,'wa':wa,'wb':wb,'ran':ran}
      if catb is None:
        if corr in ['KK','NK','KG']:
          label='xi_2pt_'+cata.name+'_'+k+'_'+corr+'_'+label0
        else:
          label='xi_2pt_'+cata.name+'_'+corr+'_'+label0
      else:
        if corr in ['KK','NK','KG']:
          label='xi_2pt_'+cata.name+'-'+catb.name+'_'+k+'_'+corr+'_'+label0
        else:
          label='xi_2pt_'+cata.name+'-'+catb.name+'_'+corr+'_'+label0
      if cata.calc_err=='jk':
        err,chi2=jackknife_methods.jk(cata,xi_2pt.xi_2pt,[xip,xim,xip_im,xim_im],label,**kwargs)
      elif cata.calc_err=='mock':
        ggperr,ggmerr,chi2p,chi2m,ceerr,cberr,cechi2,cbchi2=BCC_Methods.jk_iter_xi(cat,ggp,ggm,ce,cb,mask,w,cosebi=cosebi,parallel=parallel)

    if plot:
      fig.plot_methods.fig_create_xi(cata,catb,corr,theta,out,err,k,ga,gb)

    return theta,out,err,chi2


  @staticmethod
  def psf_alpha(cat,mask=None):

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    theta,gpout,gperr,chi2=xi_2pt.xi_2pt(cat,corr='GG',ga='psf',maska=mask,plot=False)
    theta,ppout,pperr,chi2=xi_2pt.xi_2pt(cat,cat,corr='GG',gb='psf',maska=mask,plot=False)
    theta,gpout2,gperr2,chi2=xi_2pt.xi_2pt(cat,corr='GG',ga='psf',maska=mask,plot=False,conj=True)
    theta,ppout2,pperr2,chi2=xi_2pt.xi_2pt(cat,cat,corr='GG',gb='psf',maska=mask,plot=False,conj=True)

    e1,e2=lin.linear_methods.calc_mean_stdev_rms_e(cat,mask=mask,full=False)
    psf1=lin.linear_methods.calc_mean_stdev_rms(cat,'psf1',mask=mask,full=False)
    psf2=lin.linear_methods.calc_mean_stdev_rms(cat,'psf2',mask=mask,full=False)
    e=e1+1j*e2
    psf=psf1+1j*psf2

    a=ppout[0]+1j*ppout[2]-np.abs(e)**2
    b=ppout2[0]+1j*ppout2[2]-e**2
    c=gpout[0]+1j*gpout[2]-np.conj(e)*psf
    d=gpout2[0]+1j*gpout2[2]-e*psf

    alphap=(a*np.conj(c)-np.conj(b)*d)/(np.conj(a)*a-np.conj(b)*b)
    alpham=(np.conj(a)*d-b*np.conj(c))/(np.conj(a)*a-np.conj(b)*b)

    alpha0=xi_2pt.calc_alpha(gpout[0],ppout[0],e1,e2,psf1,psf2)

    txt.write_methods.heading('Xi PSF alpha calculation results',cat,label='xi_alpha',create=True)
    txt.write_methods.write_append('theta  '+str(theta),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('gp+  '+str(gpout[0]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('gp-  '+str(gpout[1]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('pp+  '+str(ppout[0]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('pp-  '+str(ppout[1]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('gp err  '+str(gperr[0]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('pp err  '+str(pperr[0]),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('alpha+ (real)  '+str(alphap.real),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('alpha+ (imag)  '+str(alphap.imag),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('alpha- (real)  '+str(alpham.real),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('alpha- (imag)  '+str(alpham.imag),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('alpha0  '+str(alpha0),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('e1  '+str(e1),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('e2  '+str(e2),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('psf1  '+str(psf1),cat,label='xi_alpha',create=False)
    txt.write_methods.write_append('psf2  '+str(psf2),cat,label='xi_alpha',create=False)

    fig.plot_methods.fig_create_xi_alpha(cat,theta,gpout,ppout,gperr,pperr,alphap,alpham,alpha0)

    return 

  @staticmethod
  def calc_alpha(gp,pp,e1,e2,psf1,psf2):

    return (gp-e1*psf1-e2*psf2)/(pp-psf1**2-psf2**2)

  @staticmethod
  def create_shear_rm_cat(cat,cat2):
    rm10s=catalog.CatalogStore('y1_rm_shear',setup=False)

    m1,s1,m2,s2=catalog.CatalogMethods.sort(cat.coadd,cat2.coadd)
    rm10s.e1=(cat2.e1[m2])[s2]
    rm10s.e2=(cat2.e2[m2])[s2]
    rm10s.ra=(cat.ra[m1])[s1]
    rm10s.dec=(cat.dec[m1])[s1]
    rm10s.zp=(cat.zp[m1])[s1]
    rm10s.coadd=(cat.coadd[m1])[s1]

    return rm10s

  @staticmethod
  def ia_estimatorb(cat,cat2,dlos=100.,rbins=5,rmin=.1,rmax=200.,logr=False,lum=0.,comm=None,rank=None,size=None,output=False):

    import numpy.random as rand
    from mpi4py import MPI

    def chi(z,omegam=0.27,c=299792.458,H=70,n=1000):
      
      a=1.718
      b=0.315

      return z/np.sqrt(1.+a*omegam*z+b*np.sqrt(omegam)*z*z)*c/H

    def rote(cat,cat2,i,j):

      x=np.sin(cat.dec[i])*np.cos(cat2.dec[j])-np.sin(cat2.dec[j])*np.cos(cat.dec[i])*np.cos(cat.ra[i]-cat2.ra[j])
      y=np.cos(cat.dec[i])*np.sin(cat.ra[i]-cat2.ra[j])
      tdphi=2.*np.arctan2(y,x)
      
      if cat.bs:
        e1=np.dot(cat2.e1[j]-cat2.c1[j],np.cos(tdphi))+np.dot(cat2.e2[j]-cat2.c2[j],np.sin(tdphi))
        e2=np.dot(cat2.e1[j]-cat2.c1[j],np.sin(tdphi))-np.dot(cat2.e2[j]-cat2.c2[j],np.cos(tdphi))
      else:
        e1=np.dot(cat2.e1[j],np.cos(tdphi))+np.dot(cat2.e2[j],np.sin(tdphi))
        e2=np.dot(cat2.e1[j],np.sin(tdphi))-np.dot(cat2.e2[j],np.cos(tdphi))

      return e1,e2

    def pairs0(cat,cat2,dlos,r):
      print 'inside pairs'

      from sklearn.neighbors import KDTree
      # from mpi4py import MPI
      import time
      # import multiprocessing as mp

      # comm = MPI.COMM_WORLD
      # rank = comm.Get_rank()
      # size = comm.Get_size()

      jpairs=[]
      jbins=[]
      tree=KDTree(np.vstack((cat.ra,cat.dec)).T, leaf_size=2)
      i=tree.query_radius(np.vstack((cat2.ra,cat2.dec)).T, r=.05)
      t1=time.time()
      print 'done tree'
      for x,ind in enumerate(i):
        if x%10000==0:
          print 'pairs',x,time.time()-t1
        # zmask=(np.abs(cat2.r[x]-cat.r[ind])<=dlos)
        sep=physsep(cat2,cat,x,ind)
        bins=np.digitize(sep,r)-1
        mask=(bins!=-1)&(bins!=len(r)-1)&(cat.r[ind]-cat2.r[x]>dlos)
        jpairs.append(ind[mask])
        jbins.append(bins[mask])

        if x%100==0:
          print 'pairs',x,len(ind),np.sum(mask),time.time()-t1

      return jpairs,jbins


    def pairs(cat,cat2,dlos,r,mpibins,rank,lum=0.,output=True):
      #print 'inside pairs'

      from sklearn.neighbors import KDTree
      # from mpi4py import MPI
      import time
      # import multiprocessing as mp

      # comm = MPI.COMM_WORLD
      # rank = comm.Get_rank()
      # size = comm.Get_size()

      jpairs=[]
      jbins=[]
      x=cat.r*np.cos(cat.dec)*np.cos(cat.ra)
      y=cat.r*np.cos(cat.dec)*np.sin(cat.ra)
      z=cat.r*np.sin(cat.dec)
      tree=KDTree(np.vstack((x,y,z)).T, leaf_size=2)
      x=cat2.r*np.cos(cat2.dec)*np.cos(cat2.ra)
      y=cat2.r*np.cos(cat2.dec)*np.sin(cat2.ra)
      z=cat2.r*np.sin(cat2.dec)
      mask=np.arange(mpibins[rank],mpibins[rank+1]).astype(int)
      i=tree.query_radius(np.vstack((x[mask],y[mask],z[mask])).T, r=np.sqrt(r[-1]**2+dlos**2))
      t1=time.time()
      #print 'done tree'
      for x,ind in enumerate(i):
        if output:
          if x%10000==0:
            print 'pairs',x,time.time()-t1
        if len(ind)==0:
          jpairs.append([])
          jbins.append([])
          continue
        sep=physsep(cat2,cat,x,ind)
        bins=np.digitize(sep,r)-1
        mask=(bins!=-1)&(bins!=len(r)-1)&(np.abs(cat2.r[x]-cat.r[ind])<=dlos)&(cat.lum[ind]>=lum)
        jpairs.append(ind[mask])
        jbins.append(bins[mask])

        # if x%100==0:
        #   print 'pairs',x,len(ind),np.sum(mask),time.time()-t1

      return jpairs,jbins

    def angsep(cat,cat2,i,j):

      import math
      
      sep=(np.cos(cat.dec[i])*np.cos(cat.ra[i])-np.cos(cat2.dec[j])*np.cos(cat2.ra[j]))**2+(np.cos(cat.dec[i])*np.sin(cat.ra[i])-np.cos(cat2.dec[j])*np.sin(cat2.ra[j]))**2+(np.sin(cat.dec[i])-np.sin(cat2.dec[j]))**2
      
      return [2.*math.asin(np.sqrt(sep[x])/2.) for x in range(len(j))]

    def physsep(cat,cat2,i,j): 
      
      return np.dot(np.mean(cat.r[i]+cat2.r[j]),angsep(cat,cat2,i,j))

    def angsep2(cat,cat2,i,j):

      import math
      
      sep=(np.cos(cat.dec[i])*np.cos(cat.ra[i])-np.cos(cat2.dec[j])*np.cos(cat2.ra[j]))**2+(np.cos(cat.dec[i])*np.sin(cat.ra[i])-np.cos(cat2.dec[j])*np.sin(cat2.ra[j]))**2+(np.sin(cat.dec[i])-np.sin(cat2.dec[j]))**2
      
      return 2.*math.asin(np.sqrt(sep)/2.)

    def physsep2(cat,cat2,i,j): 
      
      return (cat.r[i]+cat2.r[j])/2.*angsep2(cat,cat2,i,j)

    if logr:
      r=np.logspace(np.log(rmin),np.log(rmax),rbins+1,base=np.exp(1))
    else:
      r=np.linspace(rmin,rmax,rbins+1)

    # m1,s1,m2,s2=catalog.CatalogMethods.sort(cat.coadd,cat2.coadd)
    # cat.e1=-99.*np.ones(len(cat.coadd))
    # cat.e2=-99.*np.ones(len(cat.coadd))
    # cat.e1[m1]=(cat2.e1[m2])[s2]
    # cat.e2[m1]=(cat2.e2[m2])[s2]
    cat.r=chi(cat.zp)
    catr=catalog.CatalogStore('y1_rm10_ran',setup=False,cattype='gal')
    catr.coadd=np.arange(len(cat.ran_ra))
    catr.ra=cat.ran_ra
    catr.dec=cat.ran_dec
    catr.r=np.random.choice(cat.r,size=len(catr.ra))
    catr.regs=cat.ran_regs
    catr.lum=10.*np.ones(len(cat.ran_dec))

    # cat.ra=cat.ra/180.*np.pi
    # cat.dec=cat.dec/180.*np.pi
    # catr.ra=catr.ra/180.*np.pi
    # catr.dec=catr.dec/180.*np.pi

    d0=len(cat.coadd)
    d1=np.sum(cat.e1!=0.)
    r0=len(catr.coadd)

    print d0,d1,r0

    #print 'catalog build done'

    mm=np.zeros((cat.num_regs+1,rbins))
    dm=np.zeros((cat.num_regs+1,rbins))
    rm=np.zeros((cat.num_regs+1,rbins))
    dd=np.zeros((cat.num_regs+1,rbins))
    dr=np.zeros((cat.num_regs+1,rbins))
    dr0=np.zeros((cat.num_regs+1,rbins))
    rr=np.zeros((cat.num_regs+1,rbins))
    ee=np.zeros((cat.num_regs+1,rbins))
    xx=np.zeros((cat.num_regs+1,rbins))
    de=np.zeros((cat.num_regs+1,rbins))
    dx=np.zeros((cat.num_regs+1,rbins))
    re=np.zeros((cat.num_regs+1,rbins))
    rx=np.zeros((cat.num_regs+1,rbins))

    if size is not None:
      mpibins=np.floor(np.linspace(0,len(cat.coadd),size+1))
    else:
      mpibins=np.array([0,len(cat.coadd)])
      rank=0
    j,bin=pairs(cat,cat,dlos,r,mpibins,rank,lum=lum,output=output)
    #print 'pairs dd done'
    for i in xrange(int(mpibins[rank+1]-mpibins[rank])):
      if cat.wt:
        w=cat.w[i]
      else:
        w=1.
      if output:
        if i%10000==0:
          print 'pairs dd',i
      if len(j[i])>0:
        for x,ind in enumerate(j[i]):
          if cat.wt:
            w2=cat.w[ind]
          else:
            w2=1.
          if cat.e1[ind]!=0:
            dd[0,bin[i][x]]+=1
            dd[cat.regs[i]+1,bin[i][x]]+=1
            e1,e2=rote(cat,cat,i,ind)
            de[0,bin[i][x]]+=e1*w
            if cat.bs:
              dm[0,bin[i][x]]+=(1.+cat.m[i])*w
            dx[0,bin[i][x]]+=e2*w
            de[cat.regs[i]+1,bin[i][x]]+=e1*w
            dx[cat.regs[i]+1,bin[i][x]]+=e2*w
            if cat.e1[i]!=0:
              e1a,e2a=rote(cat,cat,ind,i)
              ee[0,bin[i][x]]+=e1*e1a*w*w2
              if cat.bs:
                mm[0,bin[i][x]]+=(1.+cat.m[i])*(1.+cat.m[ind])*w*w2
              xx[0,bin[i][x]]+=e1*e2a*w*w2
              ee[cat.regs[i]+1,bin[i][x]]+=e1*e1a*w*w2
              xx[cat.regs[i]+1,bin[i][x]]+=e1*e2a*w*w2

    if size is not None:
      mpibins=np.floor(np.linspace(0,len(catr.coadd),size+1))
    else:
      mpibins=np.array([0,len(catr.coadd)])
      rank=0
    j,bin=pairs(cat,catr,dlos,r,mpibins,rank,lum=lum,output=output)
    #print 'pairs dr done',len(j)
    for i in xrange(int(mpibins[rank+1]-mpibins[rank])):
      if cat.wt:
        w=cat.w[ind]
      else:
        w=1.
      if output:
        if i%10000==0:
          print 'pairs dr',i
      if len(j[i])>0:
        for x,ind in enumerate(j[i]):
          dr0[0,bin[i][x]]+=1
          dr0[catr.regs[i]+1,bin[i][x]]+=1
          if cat.e1[ind]!=0:
            dr[0,bin[i][x]]+=1
            dr[catr.regs[i]+1,bin[i][x]]+=1
            e1,e2=rote(catr,cat,i,ind)
            re[0,bin[i][x]]+=e1*w
            rx[0,bin[i][x]]+=e2*w
            re[catr.regs[i]+1,bin[i][x]]+=e1*w
            if cat.bs:
              rm[0,bin[i][x]]+=(1.+cat.m[ind])*w
            rx[catr.regs[i]+1,bin[i][x]]+=e2*w

    if size is not None:
      mpibins=np.floor(np.linspace(0,len(catr.coadd),size+1))
    else:
      mpibins=np.array([0,len(catr.coadd)])
      rank=0
    j,bin=pairs(catr,catr,dlos,r,mpibins,rank,lum=lum,output=output)
    #print 'pairs rr done'    
    for i in xrange(int(mpibins[rank+1]-mpibins[rank])):
      if output:
        if i%10000==0:
          print 'pairs rr',i
      if len(j[i])>0:
        for x,ind in enumerate(j[i]):
          rr[0,bin[i][x]]+=1
          rr[catr.regs[i]+1,bin[i][x]]+=1

    if size is not None:
      comm.Barrier()

      if rank==0:
        mma = np.zeros_like(mm)
        dma = np.zeros_like(dm)
        rma = np.zeros_like(rm)
        dda = np.zeros_like(dd)
        dra = np.zeros_like(dr)
        dr0a = np.zeros_like(dr0)
        rra = np.zeros_like(rr)
        eea = np.zeros_like(ee)
        xxa = np.zeros_like(xx)
        dea = np.zeros_like(de)
        rea = np.zeros_like(re)
        dxa = np.zeros_like(dx)
        rxa = np.zeros_like(rx)
      else:
        mma = None
        dma = None
        rma = None
        dda = None
        dra = None
        dr0a = None
        rra = None
        eea = None
        xxa = None
        dea = None
        dxa = None
        rea = None
        rxa = None

      comm.Reduce([mm, MPI.DOUBLE],[mma, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([dm, MPI.DOUBLE],[dma, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([rm, MPI.DOUBLE],[rma, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([dd, MPI.DOUBLE],[dda, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([dr, MPI.DOUBLE],[dra, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([dr0, MPI.DOUBLE],[dr0a, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([rr, MPI.DOUBLE],[rra, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([ee, MPI.DOUBLE],[eea, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([xx, MPI.DOUBLE],[xxa, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([de, MPI.DOUBLE],[dea, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([dx, MPI.DOUBLE],[dxa, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([re, MPI.DOUBLE],[rea, MPI.DOUBLE],op=MPI.SUM,root=0)
      comm.Reduce([rx, MPI.DOUBLE],[rxa, MPI.DOUBLE],op=MPI.SUM,root=0)

      if rank==0:
        mm = mma
        dm = dma
        rm = rma
        dd = dda
        dr = dra
        dr0 = dr0a
        rr = rra
        ee = eea
        xx = xxa
        de = dea
        dx = dxa
        re = rea
        rx = rxa        

    if rank==0:

      gp=np.zeros_like(de)
      gx=np.zeros_like(dx)
      gee=np.zeros_like(ee)
      gxx=np.zeros_like(xx)
      gpm=np.zeros_like(dm)
      gmm=np.zeros_like(mm)

      gp[0]=(de[0]/d1/d0-re[0]/d1/r0)/rr[0]*r0*r0
      gx[0]=(dx[0]/d1/d0-rx[0]/d1/r0)/rr[0]*r0*r0
      gee[0]=ee[0]/d1/d1/rr[0]*r0*r0
      gxx[0]=xx[0]/d1/d1/rr[0]*r0*r0
      if cat.bs:
        gpm[0]=(dm[0]/d1/d0-rm[0]/d1/r0)/rr[0]*r0*r0
        gmm[0]=mm[0]/d1/d1/rr[0]*r0*r0
      for i in xrange(cat.num_regs):
        gp[i+1]=((de[0]-de[i+1])/d1/d0-(re[0]-re[i+1])/d1/r0)/(rr[0]-rr[i+1])*r0*r0
        gx[i+1]=((dx[0]-dx[i+1])/d1/d0-(rx[0]-rx[i+1])/d1/r0)/(rr[0]-rr[i+1])*r0*r0
        gee[i+1]=(ee[0]-ee[i+1])/d1/d1/(rr[0]-rr[i+1])*r0*r0
        gxx[i+1]=(xx[0]-xx[i+1])/d1/d1/(rr[0]-rr[i+1])*r0*r0
        if cat.bs:
          gpm[i+1]=((dm[0]-dm[i+1])/d1/d0-(rm[0]-rm[i+1])/d1/r0)/(rr[0]-rr[i+1])*r0*r0
          gmm[i+1]=(mm[0]-ee[i+1])/d1/d1/(rr[0]-rr[i+1])*r0*r0

      if cat.bs:
        gp/=gpm
        gx/=gpm
        gee/=gmm
        gxx/=gmm

      covgp=np.zeros((rbins,rbins))
      covgx=np.zeros((rbins,rbins))
      covee=np.zeros((rbins,rbins))
      covxx=np.zeros((rbins,rbins))

      for i in xrange(rbins):
        for j in xrange(rbins):
          covgp[i,j]=np.sum((gp[1:,i]-np.mean(gp[1:,i]))*(gp[1:,j]-np.mean(gp[1:,j])))*(cat.num_regs-1.)/cat.num_regs
          covgx[i,j]=np.sum((gx[1:,i]-np.mean(gx[1:,i]))*(gx[1:,j]-np.mean(gx[1:,j])))*(cat.num_regs-1.)/cat.num_regs
          covee[i,j]=np.sum((gee[1:,i]-np.mean(gee[1:,i]))*(gee[1:,j]-np.mean(gee[1:,j])))*(cat.num_regs-1.)/cat.num_regs
          covxx[i,j]=np.sum((gxx[1:,i]-np.mean(gxx[1:,i]))*(gxx[1:,j]-np.mean(gxx[1:,j])))*(cat.num_regs-1.)/cat.num_regs

    # dd/=d0**2
    # ee/=d1**2
    # xx/=d1**2
    # de/=d1*d0
    # dx/=d1*d0
    
    # dr/=d1*r0
    # dr/=d0*r0
    # re/=d1*r0
    # rx/=d1*r0

    # rr/=r0*r0

    # cat.ra=cat.ra*180./np.pi
    # cat.dec=cat.dec*180./np.pi
    # catr.ra=catr.ra*180./np.pi
    # catr.dec=catr.dec*180./np.pi

    if rank==0:
      np.save('r.npy',r)
      np.save('gp.npy',gp)
      np.save('gx.npy',gx)
      np.save('ee.npy',ee)
      np.save('xx.npy',xx)
      np.save('covgp.npy',covgp)
      np.save('covgx.npy',covgx)
      np.save('covee.npy',covee)
      np.save('covxx.npy',covxx)
      return r,gp[0],gx[0],gee[0],gxx[0],np.sqrt(np.diag(covgp)),np.sqrt(np.diag(covgx)),np.sqrt(np.diag(covee)),np.sqrt(np.diag(covxx)),[d0,d1,r0,dd,dr,dr0,rr,ee,xx,de,dx,re,rx]
    else:
      return None,None,None,None,None,None,None,None,None,None
