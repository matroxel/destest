import numpy as np
import treecorr
import os.path

import catalog
import fig
import lin

class UseError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class xi_2pt(object):

  @staticmethod
  def xi_2pt(cata,catb=None,k=None,ga=None,gb=None,corr='GG',maska=None,maskb=None,wa=None,wb=None,ran=True,mock=False,erron=True,jkmask=None,label0='',plot=False):

    maska=catalog.CatalogMethods.check_mask(cata.coadd,maska)
    jkmask=catalog.CatalogMethods.check_mask(cata.coadd,jkmask)

    maska0=maska&jkmask

    if wa is None:
      wa=np.ones(len(cata.coadd))

    e1,e2,w,ms=lin.linear_methods.get_lin_e_w_ms(cata,xi=True,mock=mock,mask=maska0,w1=wa)

    if catb is None:
      if corr not in ['GG','NN','KK']:
        raise UseError('Must supply both cata,catb for NG,NK correlations.')

    if ga!=None:
      e1=getattr(cata,ga+'1')[maska]
      e2=getattr(cata,ga+'2')[maska]
    else:
      ga='e'
    if catb is None:
      gb=ga

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

    if catb!=None:

      maskb=catalog.CatalogMethods.check_mask(catb.coadd,maskb)

      if wb is None:
        wb=np.ones(len(catb.coadd))

      e1,e2,w,ms=lin.linear_methods.get_lin_e_w_ms(catb,xi=True,mock=mock,mask=maskb,w1=wb)

      if gb!=None:
        e1=getattr(cata,gb+'1')[maskb]
        e2=getattr(cata,gb+'2')[maskb]
      else:
        gb='e'

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
      fig.plot_methods.fig_create_xi(cata,corr,theta,out,err,k,ga,gb)

    return theta,out,err,chi2

