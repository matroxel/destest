import numpy as np
import os

import catalog
import lin
import fig
import txt
import cosmosis


class pz_methods(object):

  @staticmethod
  def build_nofz_bins(pz0,pzref=None,pzlow=0.0,pzhigh=2.5,cat=None,bins=3,split='mean',nztype='pdf',pzmask=None,catmask=None,spec=False,boot=None):
    """
    Build n(z) for a PZStore object that contains full pdf information. Optionally match to a catalog. Masking for both pz's and catalog optional. split determines which point estimate is used for binning. bins is number of tomographic bins. pzlow, pzhigh are bounds of reliable photo-zs. Stores result in PZStore object.
    """

    pzmask=catalog.CatalogMethods.check_mask(pz0.coadd,pzmask)
    if pzref is None:
      pzref=pz0

    if cat is not None:
      catmask=catalog.CatalogMethods.check_mask(cat.coadd,catmask)
      m1,s1,m2,s2=catalog.CatalogMethods.sort(pz0.coadd,cat.coadd)
      s2=s2[catmask]
      s2=np.arange(len(cat.coadd))
      catmask=np.ones(len(cat.coadd)).astype(bool)
      pzmask=pzmask&m1
      catmask=catmask&m2
      if pz0.wt:
        w=pz0.w*(cat.w[catmask])[s2]
        pointw=pzref.w*(cat.w[catmask])[s2]
      else:
        w=(cat.w[catmask])[s2]
        pointw=(cat.w[catmask])[s2]
    else:
      if pz0.wt:
        w=pz0.w
        pointw=pzref.w
      else:
        w=np.ones(len(pz0.w))
        pointw=np.ones(len(pzref.w))

    m1,s1,m2,s2=catalog.CatalogMethods.sort(pz0.coadd,pzref.coadd)
    if split=='mean':
      pointz=pzref.z_mean_full[m2][s2]
    elif split=='peak':
      pointz=pzref.z_peak_full[m2][s2]
    elif split=='spec':
      pointz=pzref.spec_full[m2][s2]
    else:
      print 'need split type'
      return

    mask0=(pointz>=pzlow)&(pointz<=pzhigh)

    if nztype=='mean':
      pzdist=pz0.z_mean_full
      point=True
    elif nztype=='peak':
      pzdist=pz0.z_peak_full
      point=True
    elif nztype=='mc':
      pzdist=pz0.z_mc_full
      point=True
    else:
      if pz0.pdftype=='sample':
        pzdist=pz0.pz_full
        point=True
      else:
        point=False
        pzdist=None

    pz0.tomo=bins+1

    nofz,specnofz=pz_methods.nofz_bins(pz0,pointz,pzdist,mask0,w,pointw,bins=bins,pzmask=pzmask,spec=spec,point=point)

    pz0.pz=nofz
    if spec:
      pz0.spec=specnofz

    if boot is not None:
      pz0.boot=boot
      pz0.bootspec=np.zeros((pz0.tomo,pz0.boot,pz0.bins))
      pz0.bootpz=np.zeros((pz0.tomo,pz0.boot,pz0.bins))
      for i in xrange(boot):
        bootmask=np.zeros(len(mask0),dtype=int)
        bootmask[np.random.choice(np.where(mask0)[0],np.sum(mask0),replace=True)]=1

        nofz,specnofz=pz_methods.nofz_bins(pz0,pointz,pzdist,bootmask.astype(bool),w,bins=bins,pzmask=pzmask,spec=spec,point=point)

        pz0.bootpz[:,i,:]=nofz
        if spec:
          pz0.bootspec[:,i,:]=specnofz

    return

  @staticmethod
  def nofz_bins(pz0,pointz,pzdist,mask0,w,pointw,bins=3,pzmask=None,spec=False,point=False):

    edge=lin.linear_methods.find_bin_edges(pointz[mask0],bins,pointw[mask0])
    print edge
    xbins=np.digitize(pointz,edge)-1

    nofz=np.zeros((bins+1,pz0.bins))

    if point:
      nofz[0,:],b=np.histogram(pzdist[pzmask&mask0],bins=np.append(pz0.binlow,pz0.binhigh[-1]),weights=w[pzmask&mask0])
      nofz[0,:]/=np.sum(nofz[0,:])*(pz0.bin[1]-pz0.bin[0])

      for i in xrange(bins):
        mask=(xbins==i)
        nofz[i+1,:],b=np.histogram(pzdist[pzmask&mask0&mask],bins=np.append(pz0.binlow,pz0.binhigh[-1]),weights=w[pzmask&mask&mask0])
        nofz[i+1,:]/=np.sum(nofz[i+1,:])*(pz0.bin[1]-pz0.bin[0])

    else:
      print np.sum(pzmask&mask0),len(pzmask),len(mask0)
      nofz[0,:]=np.sum((pz0.pz_full[pzmask&mask0].T*w[pzmask&mask0]).T,axis=0)
      nofz[0,:]/=np.sum(nofz[0,:])*(pz0.bin[1]-pz0.bin[0])

      for i in xrange(bins):
        mask=(xbins==i)
        nofz[i+1,:]=np.sum((pz0.pz_full[pzmask&mask0&mask].T*w[pzmask&mask&mask0]).T,axis=0)
        nofz[i+1,:]/=np.sum(nofz[i+1,:])*(pz0.bin[1]-pz0.bin[0])

    specnofz=np.zeros((bins+1,pz0.bins))

    if spec:
      from weighted_kde import gaussian_kde
      tmp=gaussian_kde(pz0.spec_full[pzmask&mask0],weights=w[pzmask&mask0],bw_method='scott')
      specnofz[0,:]=tmp(pz0.bin)
      specnofz[0,:]/=np.sum(specnofz[0,:])*(pz0.bin[1]-pz0.bin[0])

      for i in xrange(bins):
        mask=(xbins==i)
        tmp=gaussian_kde(pz0.spec_full[pzmask&mask0&mask],weights=w[pzmask&mask0&mask],bw_method='scott')
        specnofz[i+1,:]=tmp(pz0.bin)
        specnofz[i+1,:]/=np.sum(specnofz[i+1,:])*(pz0.bin[1]-pz0.bin[0])

    return nofz,specnofz

  @staticmethod
  def pdf_mvsk(x,pdf0,dm,dv,ds,dk,p6=False):
    from numpy.polynomial.hermite import hermval
    from scipy.stats import norm
    from astropy.convolution import convolve

    v,m=norm.fit(pdf0)

    if dv>=1.:
      if dv==1:
        dv+=1e-4
      pdf2=norm.pdf(x,loc=np.median(x)-dm,scale=np.sqrt((dv*v)**2-v**2))
      pdf=convolve(pdf0,pdf2,normalize_kernel=True)
    else:
      v2=np.sqrt(1./(1./(v*dv)**2-1./v**2))
      m2=v2**2*((m+dm)/(v*dv)**2-m/v**2)
      pdf2=norm.pdf(x,loc=m2,scale=v2)
      norm=1./np.sqrt(2.*np.pi*v**2*v2**2/(v*dv)**2)*np.exp(-(v*dv)**2*(m-m2)**2/v2**2/v**2/2.)
      pdf=pdf0*pdf2/norm
    if p6:
      return hermval(x-m,[1.,0.,0.,ds/6.,dk/24.,0.,ds**2./72.])*pdf
    else:
      return hermval(x-m,[1.,0.,0.,ds/6.,dk/24.])*pdf


class pz_spec_validation(object):

  @staticmethod
  def calc_bootstrap(pz,test,dir,tomobins,notomo=False):
    """
    Calculate bootstrap for correlation functions in spec validation tests.
    """

    ratio=np.zeros((tomobins,pz.boot))
    var=np.zeros((tomobins))

    if notomo:
      data0=np.loadtxt(dir+test+'/out/sim_data_notomo_spec_'+pz.name+'/data.txt')[2:]
    else:
      data0=np.loadtxt(dir+test+'/out/sim_data_spec_'+pz.name+'/data.txt')[:,2:]

    for i in xrange(pz.boot):
      if notomo:
        data=np.loadtxt(dir+test+'/out/sim_data_notomo_'+pz.name+'_'+str(i)+'/data.txt')[2:]
      else:
        data=np.loadtxt(dir+test+'/out/sim_data_'+pz.name+'_'+str(i)+'/data.txt')[:,2:]

      for bin in range(tomobins):
        if notomo:
          if bin>0:
            continue
          ratio[bin,i]=np.mean((data-data0)/data0)
        else:
          ratio[bin,i]=np.mean((data[bin,:]-data0[bin,:])/data0[bin,:])

    for bin in range(tomobins):
      if notomo&(bin>0):
        continue
      var[bin]=np.sum((ratio[bin,:]-np.mean(ratio[bin,:]))*(ratio[bin,:]-np.mean(ratio[bin,:])))/(pz.boot-1.)

    if notomo:
      return np.sqrt(var[0])

    print var

    return np.sqrt(var)

  @staticmethod
  def calc_bootstrap_sig8(pz,test,dir,param,notomo=False):
    """
    Calculate bootstrap for cosmological parameters in spec validation tests.
    """

    from astropy.table import Table

    mean=[]
    for i in xrange(50):
      if notomo:
        mean0=Table.read(dir+test+'/out/sim_data_notomo_'+pz.name+'_'+str(i)+'/means.txt', format='ascii')
      else:
        mean0=Table.read(dir+test+'/out/sim_data_'+pz.name+'_'+str(i)+'/means.txt', format='ascii')

      for row in mean0:
        if param in row['parameter']:
          mean=np.append(mean,row['mean'])
    print mean

    var=np.sum((mean-np.mean(mean))*(mean-np.mean(mean)))/(pz.boot-1.)

    print var

    return np.sqrt(var)

  @staticmethod
  def bpz_prior(tpf='/home/troxel/bpz/bpz/build_prior_info2.txt',mf='/home/troxel/bpz/bpz/build_prior_info.txt',specf='/home/troxel/spec_cat_0.fits.gz'):

    import fitsio as fio
    from scipy import ndimage
    
    z=np.arange(0.0100,3.0100,0.0100)
    m0=np.arange(10,40,.25)
    spec=fio.FITS(specf)[-1].read()['z_spec']
    wt=fio.FITS(specf)[-1].read()['weights']
    tp=np.loadtxt(tpf)
    tp=tp.reshape((len(spec), 64))
    mf=np.loadtxt(mf)
    coadd=mf[:,0]
    m=mf[:,1]
    m1,s1,m2,s2=catalog.CatalogMethods.sort(coadd,fio.FITS(specf)[-1].read()['coadd_objects_id'])
    spec=spec[s2]
    wt=wt[s2]

    priors1=np.zeros((len(z),len(m0),64))
    priors2=np.zeros((len(z),len(m0),64))
    for i in range(len(z)):
      for mi in range(len(m0)):
        mask=(spec>z[i]-(z[1]-z[0])/2.)&(spec<=z[i]+(z[1]-z[0])/2.)&(m[i]>m0[mi]-(m0[1]-m0[0])/2.)&(m[i]<=m0[mi]+(m0[1]-m0[0])/2.)
        if np.sum(mask)==0:
          priors1[i,mi,:]=.000001*np.ones((64))
          priors1[i,mi,:]/=np.sum(priors1[i,:])
          priors2[i,mi,:]=.000001*np.ones((64))
          priors2[i,mi,:]/=np.sum(priors2[i,:])
        else:
          priors1[i,mi,:]=np.mean(tp[mask],axis=0)
          priors1[i,mi,:]/=np.sum(priors1[i,:])
          priors2[i,mi,:]=np.average(tp[mask],axis=0,weights=wt[mask])+.0001*np.ones((64))
          priors2[i,mi,:]/=np.sum(priors2[i,:])

    np.save('bpz_prior1.npy',priors1)
    np.save('bpz_prior2.npy',priors2)
    gauss1=ndimage.gaussian_filter(priors1, (10, 10,10))
    gauss2=ndimage.gaussian_filter(priors2, (10, 10,10))
    np.save('bpz_prior1_smooth.npy',gauss1)
    np.save('bpz_prior2_smooth.npy',gauss2)

    return priors1, priors2

  @staticmethod
  def zbin_tests(pz0):

    pz_methods.build_nofz_bins(pz0,pzlow=0.0,pzhigh=5.,cat=None,bins=3,split='mean',pzmask=None,catmask=None,spec=False,point=False)
    pz0.spec=np.copy(pz0.pz)
    cosmosis.make.nofz(pz0,'test_zbins_3bin_251',zmax=5.3)
    cosmosis.make.nofz(pz0,'test_zbins_3bin_200',zmax=5.01,zbins=200)
    cosmosis.make.nofz(pz0,'test_zbins_3bin_150',zmax=5.01,zbins=150)
    cosmosis.make.nofz(pz0,'test_zbins_3bin_100',zmax=5.01,zbins=100)
    cosmosis.make.nofz(pz0,'test_zbins_3bin_50',zmax=5.01,zbins=50)

    pz_methods.build_nofz_bins(pz0,pzlow=0.0,pzhigh=5.,cat=None,bins=6,split='mean',pzmask=None,catmask=None,spec=False,point=False)
    pz0.spec=np.copy(pz0.pz)
    cosmosis.make.nofz(pz0,'test_zbins_6bin_251',zmax=5.3)
    cosmosis.make.nofz(pz0,'test_zbins_6bin_200',zmax=5.01,zbins=200)
    cosmosis.make.nofz(pz0,'test_zbins_6bin_150',zmax=5.01,zbins=150)
    cosmosis.make.nofz(pz0,'test_zbins_6bin_100',zmax=5.01,zbins=100)
    cosmosis.make.nofz(pz0,'test_zbins_6bin_50',zmax=5.01,zbins=50)

    return


  @staticmethod
  def zmax_tests(pz0):

    pz_methods.build_nofz_bins(pz0,pzlow=0.0,pzhigh=5.,cat=None,bins=3,split='mean',pzmask=None,catmask=None,spec=False,point=False)
    pz0.spec=np.copy(pz0.pz)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_50',zmax=5.3)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_30',zmax=3.)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_25',zmax=2.5)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_20',zmax=2.)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_175',zmax=1.75)
    cosmosis.make.nofz(pz0,'test_zmax_3bin_15',zmax=1.5)

    pz_methods.build_nofz_bins(pz0,pzlow=0.0,pzhigh=5.,cat=None,bins=6,split='mean',pzmask=None,catmask=None,spec=False,point=False)
    pz0.spec=np.copy(pz0.pz)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_50',zmax=5.3)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_30',zmax=3.)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_25',zmax=2.5)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_20',zmax=2.)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_175',zmax=1.75)
    cosmosis.make.nofz(pz0,'test_zmax_6bin_15',zmax=1.5)

    return

  @staticmethod
  def spec_comp(pz0,ref=0):

    for j,split in enumerate(['mean','peak','spec']):
      for i,x in enumerate(['pdf','mc','mean','peak']):
        for j,n in enumerate([3,6]):
          for ipz,pz in enumerate(pz0):
            print x,n,pz
            pz.wt=True
            pz_methods.build_nofz_bins(pz,pzref=pz0[ref],pzlow=0.2,pzhigh=1.3,cat=None,bins=n,split=split,nztype=x,spec=True)
          fig.plot_methods.plot_nofz_comp_pz(pz0,label='spec_'+split+'_'+x,spec=True)
          for ipz,pz in enumerate(pz0):
            pz.wt=False
            pz_methods.build_nofz_bins(pz,pzref=pz0[ref],pzlow=0.2,pzhigh=1.3,cat=None,bins=n,split=split,nztype=x,spec=True)
          fig.plot_methods.plot_nofz_comp_pz(pz0,label='spec_'+split+'_'+x,spec=True)

    return

  @staticmethod
  def sig_crit_inv(z_l,z_s,c=299792.458,G=4.302e-3,omegam=0.27,H0=70,const=6.01389688186074e-19):

    import cosmology

    c0=cosmology.Cosmo(H0=H0,omega_m=omegam)

    zs,zl=np.meshgrid(z_s,z_l)

    invsigcrit=c0.Da(0.,zl)*c0.Da(zl,zs)/c0.Da(0.,zs)*4.*np.pi*G/c/c/1e6
    invsigcrit[np.isnan(invsigcrit)]=0.
    invsigcrit[np.isinf(invsigcrit)]=0.
    invsigcrit[invsigcrit<0]=0.

    return np.reshape(invsigcrit,(len(z_l),len(z_s)))

  @staticmethod
  def sig_crit_spec(pz0,bins=3,pzlow=0.,pzhigh=2.5,load=False,bootstrap=None,point=None,lensbins=[.1,1.,20]):

    import matplotlib.gridspec as gridspec

    for ipz,pz in enumerate(pz0):

      pzmask=(pz.spec_full>pzlow)&(pz.spec_full<pzhigh)

      if load:
        if pz.wt:
          edge=lin.linear_methods.find_bin_edges(pz.spec_full[pzmask],bins,pz.w[pzmask])
        else:
          edge=lin.linear_methods.find_bin_edges(pz.spec_full[pzmask],bins)
      else:

        z_l=np.linspace(lensbins[0],lensbins[1],lensbins[2])

        if point is not None:
          sigcritpz0=pz_spec_validation.sig_crit_inv(z_l,getattr(pz,'z_'+point+'_full')).T
        else:
          pz_full=pz.pz_full/np.sum(pz.pz_full,axis=1, keepdims=True)
          sigcritpz0=np.dot(pz_full,pz_spec_validation.sig_crit_inv(z_l,pz.bin).T)

        sigcritspec0=pz_spec_validation.sig_crit_inv(z_l,pz.spec_full).T

        for ibin in xrange(bins+1):
          if pz.wt:
            edge=lin.linear_methods.find_bin_edges(pz.spec_full[pzmask],bins,pz.w[pzmask])
          else:
            edge=lin.linear_methods.find_bin_edges(pz.spec_full[pzmask],bins)
          xbins=np.digitize(pz.spec_full,edge)-1
          if ibin==0:
            mask=np.ones(len(xbins)).astype(bool)
          else:
            mask=(xbins==ibin)&pzmask

          if pz.wt:
            if np.sum(pz.w[mask])==0:

              print pz.name,ibin,' ---- Weights sum to zero???'

              data=np.zeros(len(z_l))

            else:

              sigcritpz=np.average(sigcritpz0[mask,:],weights=pz.w[mask],axis=0)

              sigcritspec=np.average(sigcritspec0[mask,:],weights=pz.w[mask],axis=0)
          else:
            sigcritpz=np.average(sigcritpz0[mask,:],axis=0)

            sigcritspec=np.average(sigcritspec0[mask,:],axis=0)

          data=(sigcritpz-sigcritspec)/sigcritspec

          if bootstrap is not None:
            bootspec=np.zeros((bootstrap,len(z_l)))
            bootpz=np.zeros((bootstrap,len(z_l)))
            for boot in xrange(bootstrap):
              if ibin==0:
                bootmask=np.random.choice(np.sum(mask),np.sum(mask),replace=True)
              else:
                bootmask=np.random.choice(np.where(mask)[0],np.sum(mask),replace=True)

              if pz.wt:
                if np.sum(pz.w[bootmask])==0:

                  print pz.name,ibin,' ---- Weights sum to zero???'

                  bootpz[boot,:]=np.ones(len(z_l))
                  bootspec[boot,:]=np.ones(len(z_l))

                else:
                  bootpz[boot,:]=np.average(sigcritpz0[bootmask,:],weights=pz.w[bootmask],axis=0)
                  bootspec[boot,:]=np.average(sigcritspec0[bootmask,:],weights=pz.w[bootmask],axis=0)
              else:
                bootpz[boot,:]=np.average(sigcritpz0[bootmask,:],axis=0)
                bootspec[boot,:]=np.average(sigcritspec0[bootmask,:],axis=0)

            bootdata=(bootpz-bootspec)/bootspec

            var=np.sqrt(np.mean((bootdata-np.mean(bootdata,axis=0))*(bootdata-np.mean(bootdata,axis=0)),axis=0)*bootstrap/(bootstrap-1.))
          else:
            var=np.zeros(len(data))

          if pz.wt:
            np.save('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_dat_weighted.npy',data)
            np.save('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_var_weighted.npy',var)
          else:
            np.save('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_dat.npy',data)
            np.save('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_var.npy',var)
          np.save('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_z.npy',z_l)

    fig.plot_methods.sig_crit_spec(pz0,bins,bootstrap,point,lensbins,edge)

    for ipz,pz in enumerate(pz0):

      if not load:

        z_l=(np.linspace(0.,2.02,101)[:-1]+np.linspace(0.,2.02,101)[1:])/2.
        z_s=np.linspace(0.,2.02,101)

        if point is not None:
          sigcritpz0=pz_spec_validation.sig_crit_inv(z_l,getattr(pz,'z_'+point+'_full')).T
        else:
          pz_full=pz.pz_full/np.sum(pz.pz_full,axis=1, keepdims=True)
          sigcritpz0=np.dot(pz_full,pz_spec_validation.sig_crit_inv(z_l,pz.bin).T)

        sigcritspec0=pz_spec_validation.sig_crit_inv(z_l,pz.spec_full).T

        data=np.zeros((100,100))

        for i in range(100):
          mask=(pz.spec_full>z_s[i])&(pz.spec_full<=z_s[i+1])

          if pz.wt:
            if np.sum(pz.w[mask])==0:

              print pz.name,ibin,' ---- Weights sum to zero???'

            else:

              sigcritpz=np.average(sigcritpz0[mask,:],weights=pz.w[mask],axis=0)

              sigcritspec=np.average(sigcritspec0[mask,:],weights=pz.w[mask],axis=0)
          else:
            sigcritpz=np.average(sigcritpz0[mask,:],axis=0)

            sigcritspec=np.average(sigcritspec0[mask,:],axis=0)

          data[:,i]=(sigcritpz-sigcritspec)/sigcritspec

        data[np.isinf(data)]=np.nan
        data[data>1]=np.nan
        np.save('text/spec_'+pz.name+'_invsigcrit_zl.npy',z_l)
        if pz.wt:
          np.save('text/spec_'+pz.name+'_invsigcrit_dat_weighted.npy',data)
        else:
          np.save('text/spec_'+pz.name+'_invsigcrit_dat.npy',data)

    fig.plot_methods.sig_crit_spec2(pz0,point)

    return