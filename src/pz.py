import numpy as np
import os

import catalog
import lin
import fig
import txt


class pz_methods(object):

  @staticmethod
  def build_nofz_bins(pz0,pzlow,pzhigh,cat=None,bins=3,split='mean',pzmask=None,catmask=None):

    pzmask=catalog.CatalogMethods.check_mask(pz0.coadd,pzmask)

    if cat is not None:
      catmask=catalog.CatalogMethods.check_mask(cat.coadd,catmask)
      m1,s1,m2,s2=catalog.CatalogMethods.sort(pz0.coadd,cat.coadd)
      s2=s2[catmask]
      s2=np.arange(len(cat.coadd))
      catmask=np.ones(len(cat.coadd)).astype(bool)
      pzmask=pzmask&m1
      catmask=catmask&m2
      w=pz0.w[pzmask]*(cat.w[catmask])[s2]
    else:
      w=pz0.w[pzmask]

    if split=='mean':
      pointz=pz0.z_mean_full[pzmask]
    elif split=='median':
      pointz=pz0.z_median_full[pzmask]
    else:
      print 'need split type'
      return

    edge=lin.linear_methods.find_bin_edges(pointz,bins,w)
    xbins=np.digitize(pointz,edge)-1

    nofz=np.zeros((bins+1,pz0.bins))
    nofz[0,:]=np.sum(pz0.pz_full[pzmask],axis=0)
    nofz[0,:]/=np.sum(nofz[0,:])

    for i in xrange(bins):
      mask=(xbins==i)&(pointz>=pzlow)&(pointz<=pzhigh)
      nofz[i+1,:]=np.sum(pz0.pz_full[pzmask&mask],axis=0)
      nofz[i+1,:]/=np.sum(nofz[i+1,:])

    pz0.pz=nofz
    pz0.tomo=bins+1

    return

class pz_spec_validation(object):

  @staticmethod
  def plot_nofz(pz0,test):

    import matplotlib
    matplotlib.use ('agg')
    import matplotlib.pyplot as plt
    plt.style.use('/home/troxel/SVA1/SVA1StyleSheet.mplstyle')

    col=['k','r','g','b','r','g','b']
    for i in xrange(len(pz0.pz)-1):
      plt.plot(pz0.bin,pz0.pz[i+1,:],color=col[i+1],linestyle='-',linewidth=1.,drawstyle='steps-mid',label='')
      plt.plot(pz0.bin,pz0.spec[i+1,:],color=col[i+1],linestyle=':',linewidth=3.,drawstyle='steps-mid',label='')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$n(z)$')
    plt.savefig('pz_nofz_'+test+'.png')
    plt.close()

    return

  @staticmethod
  def calc_bootstrap(test,dir,tomobins,notomo=False):

    ratio=np.zeros((tomobins,50))
    var=np.zeros((tomobins))

    if notomo:
      data0=np.loadtxt(dir+test+'/out/sim_data_notomo_spec_skynet/data.txt')[2:]
    else:
      data0=np.loadtxt(dir+test+'/out/sim_data_spec_skynet/data.txt')[:,2:]

    for i in xrange(50):
      if notomo:
        data=np.loadtxt(dir+test+'/out/sim_data_notomo_skynet_'+str(i)+'/data.txt')[2:]
      else:
        data=np.loadtxt(dir+test+'/out/sim_data_skynet_'+str(i)+'/data.txt')[:,2:]

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
      var[bin]=np.mean((ratio[bin,:]-np.mean(ratio[bin,:]))*(ratio[bin,:]-np.mean(ratio[bin,:])))*50./49.

    if notomo:
      return np.sqrt(var[0])

    print var

    return np.sqrt(var)

  @staticmethod
  def calc_bootstrap_sig8(test,dir,param,notomo=False):

    from astropy.table import Table

    mean=[]
    for i in xrange(50):
      if notomo:
        mean0=Table.read(dir+test+'/out/sim_data_notomo_skynet_'+str(i)+'/means.txt', format='ascii')
      else:
        mean0=Table.read(dir+test+'/out/sim_data_skynet_'+str(i)+'/means.txt', format='ascii')

      for row in mean0:
        if param in row['parameter']:
          mean=np.append(mean,row['mean'])
    print mean

    var=np.mean((mean-np.mean(mean))*(mean-np.mean(mean)))*49./50.

    print var

    return np.sqrt(var)
