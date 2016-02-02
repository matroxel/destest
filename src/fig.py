
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
plt.style.use('/home/troxel/SVA1/SVA1StyleSheet.mplstyle')
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import healpy as hp

import catalog
import config
import pz

class plot_methods(object):
  """
  Plotting routines used by various modules.
  """

  @staticmethod
  def plot_hist(x1,bins=500,name='',label='',tile=''):

    plt.figure()
    plt.hist(x1,bins=bins)
    plt.ylabel(r'$n$')
    s=config.lbl.get(label,None)
    if config.log_val.get(label,None):
      s='log '+s
    plt.xlabel(s)
    plt.minorticks_on()
    if tile!='':
      name='tile_'+tile+'_'+name
    plt.savefig('plots/hist/hist_'+name+'_'+label+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_2D_hist(x1,y1,bins=500,xname='',yname='',xlabel='',ylabel='',xtile='',ytile=''):

    plt.figure()
    plt.hist2d(x1,y1,bins=bins)
    s=config.lbl.get(xlabel,None)
    if config.log_val.get(xlabel,None):
      s='log '+s
    plt.xlabel(s)
    s=config.lbl.get(ylabel,None)
    if config.log_val.get(ylabel,None):
      s='log '+s
    plt.ylabel(s)
    plt.minorticks_on()
    if xtile!='':
      xname='tile_'+xtile+'_'+xname
    plt.savefig('plots/hist/hist_2D_'+xname+'_'+yname+'_'+xlabel+'_'+ylabel+'.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.hist2d(x1,y1,bins=bins,norm=LogNorm())
    s=config.lbl.get(xlabel,None)
    if config.log_val.get(xlabel,None):
      s='log '+s
    plt.xlabel(s)
    s=config.lbl.get(ylabel,None)
    if config.log_val.get(ylabel,None):
      s='log '+s
    plt.ylabel(s)   
    plt.minorticks_on()
    plt.savefig('plots/hist/hist_2D_'+xname+'_'+yname+'_'+xlabel+'_'+ylabel+'_log.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_hexbin(x1,cat,mask=None,bins=20,name='',label='',tile=''):

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)
    s82=mask&(cat.dec>-10)
    spta=mask&(~s82)&(cat.ra<0)
    sptc=mask&(~s82)&(cat.ra>50)
    sptb=mask&(~s82)&(~spta)&(~sptc)

    plot_methods.plot_hexbin_base(x1,cat,mask=mask&s82,label=label,bins=bins,part='s82',name=name,tile=tile)
    plot_methods.plot_hexbin_base(x1,cat,mask=mask&spta,label=label,bins=bins,part='spta',name=name,tile=tile)
    plot_methods.plot_hexbin_base(x1,cat,mask=mask&sptb,label=label,bins=bins,part='sptb',name=name,tile=tile)
    plot_methods.plot_hexbin_base(x1,cat,mask=mask&sptc,label=label,bins=bins,part='sptc',name=name,tile=tile)

    return

  @staticmethod
  def plot_hexbin_base(x1,cat,mask=None,bins=20,name='',label='',part='',tile=''):

    ra1=np.max(cat.ra[mask])
    ra0=np.min(cat.ra[mask])
    dec1=np.max(cat.dec[mask])
    dec0=np.min(cat.dec[mask])

    plt.figure()

    plt.hexbin(cat.ra[mask],cat.dec[mask],x1,gridsize=(int((ra1-ra0)*bins),int((dec1-dec0)*bins)), cmap=plt.cm.afmhot,linewidth=0)
    cb = plt.colorbar(orientation='horizontal')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.xlim((ra0-0.1*(ra1-ra0),ra1+0.1*(ra1-ra0)))
    plt.ylim((dec0-0.1*(dec1-dec0),dec1+0.1*(dec1-dec0)))
    plt.minorticks_on()
    s=config.lbl.get(label,None)
    if config.log_val.get(label,None):
      s='log '+s
    if tile!='':
      name='tile_'+tile+'_'+name
    cb.set_label(s)
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/hexbin_'+name+'_'+label+'_'+part+'.png', dpi=500,bbox_inches='tight')
    plt.close()

    # plt.figure()
    # plt.hexbin(cat.ra[mask],cat.dec[mask],x1,gridsize=bins,bins='log', cmap=plt.cm.afmhot,linewidth=0)
    # cb = plt.colorbar()
    # plt.xlabel('RA')
    # plt.ylabel('Dec')
    # plt.minorticks_on()
    # s=config.lbl.get(label,None)
    # if config.log_val.get(label,None):
    #   s='log '+s
    # cb.set_label(s)
    # plt.gca().set_aspect('equal', 'box')
    # plt.savefig('plots/footprint/hexbin_'+name+'_'+label+'_'+part+'_log.png', dpi=500,bbox_inches='tight')
    # plt.close()

    return


  @staticmethod
  def plot_field_footprint(cat,mask=None,label='',bins=100):

    dra=np.max(cat.ra)-np.min(cat.ra)
    ddec=np.max(cat.dec)-np.min(cat.dec)

    plt.figure()
    a=plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=((np.min(cat.ra[mask])-.1*dra,np.max(cat.ra[mask])+.1*dra),(np.min(cat.dec[mask])-.1*ddec,np.max(cat.dec[mask])+.1*ddec)),normed=True,cmax=1.2e-5,cmin=0.5e-5, cmap=plt.cm.afmhot,interpolation='nearest')
    cb = plt.colorbar(orientation='horizontal')
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/field_'+cat.name+'_'+label+'.png', dpi=500, bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=((np.min(cat.ra[mask])-.1*dra,np.max(cat.ra[mask])+.1*dra),(np.min(cat.dec[mask])-.1*ddec,np.max(cat.dec[mask])+.1*ddec)),normed=True,norm=LogNorm(), cmap=plt.cm.afmhot)
    #cb = plt.colorbar()
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/field_'+cat.name+'_'+label+'_log.png', dpi=500, bbox_inches='tight')
    plt.close()    

    return

  @staticmethod
  def plot_footprint(cat,mask=None,label='',bins=100,cap=None):

    s82=mask&(cat.dec>-10)
    spta=mask&(~s82)&(cat.ra<0)
    sptc=mask&(~s82)&(cat.ra>50)
    sptb=mask&(~s82)&(~spta)&(~sptc)

    plot_methods.plot_footprint_base(cat,mask=mask&s82,label=label,bins=bins,part='s82',cap=cap)
    plot_methods.plot_footprint_base(cat,mask=mask&spta,label=label,bins=bins,part='spta',cap=cap)
    plot_methods.plot_footprint_base(cat,mask=mask&sptb,label=label,bins=bins,part='sptb',cap=cap)
    plot_methods.plot_footprint_base(cat,mask=mask&sptc,label=label,bins=bins,part='sptc',cap=cap)

    return

  @staticmethod
  def plot_footprint_base(cat,mask=None,label='',bins=100,part='',cap=None):

    # if not hasattr(cat, 'gdmask'):
    #   cat.gdmask=hp.reorder(hp.read_map(config.golddir+'y1a1_gold_1.0.1_wide_footprint_4096.fit'),inp='ring',out='nested')
    #   cat.badmask=hp.reorder(hp.read_map(config.golddir+'y1a1_gold_1.0.1_wide_badmask_4096.fit'),inp='ring',out='nested')
    # if not hasattr(cat,'pix'):
    #   cat.pix=radec_to_hpix(ra,dec,nside=4096,nest=True)

    # mask0=(cat.gdmask>=1)&(cat.badmask==0)

    # hpmap=np.ones((12*4096**2))*hp.UNSEEN
    # hpmap[(cat.gdmask>=1)]=0
    # cnt=np.zeros(12*4096**2)
    # cnt[:np.max(cat.pix[mask])+1]=np.bincount(cat.pix[mask])
    # hpmap[mask0]=cnt[mask0]
    # hp.cartview(hpmap,latra=config.dec_lim.get(part),lonra=config.ra_lim.get(part),nest=True,xsize=10000,title=label)
    # plt.savefig('plots/footprint/footprint_'+cat.name+'_'+label+'_'+part+'.png', dpi=1000,bbox_inches='tight')
    # plt.close()

    dra=config.ra_lim.get(part)[1]-config.ra_lim.get(part)[0]
    ddec=config.dec_lim.get(part)[1]-config.dec_lim.get(part)[0]

    plt.figure()
    tmp=config.ra_lim.get(part),config.dec_lim.get(part)
    a=plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=tmp,normed=True, cmap=plt.cm.afmhot,cmax=cap)

    cb = plt.colorbar(orientation='horizontal')
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/footprint_'+cat.name+'_'+label+'_'+part+'.png', dpi=500, bbox_inches='tight')
    plt.close()

    # plt.figure()
    # plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=((np.min(cat.ra[mask])-.1*dra,np.max(cat.ra[mask])+.1*dra),(np.min(cat.dec[mask])-.1*ddec,np.max(cat.dec[mask])+.1*ddec)),normed=True,norm=LogNorm(), cmap=plt.cm.afmhot)
    # #cb = plt.colorbar()
    # plt.gca().set_aspect('equal', 'box')
    # plt.savefig('plots/footprint/footprint_'+cat.name+'_'+label+'_log.png', dpi=500, bbox_inches='tight')
    # plt.close()    

    return

  @staticmethod
  def plot_whisker(x,y,e1,e2,name='',label='',scale=.01,key='',chip=False):

    plt.figure()
    Q = plt.quiver(x,y,e1,e2,units='width',pivot='middle',headwidth=0,width=.0005)
    if chip:
      plt.quiverkey(Q,0.2,0.125,scale,str(scale)+' '+key,labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
      plt.xlim((-250,4250))
      plt.ylim((-200,2100))
    else:
      plt.quiverkey(Q,0.2,0.2,scale,str(scale)+' '+key,labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
    plt.savefig('plots/footprint/whisker_'+name+'_'+label+'.png', dpi=500,bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_lin_split(x,e1,e2,e1err,e2err,m1,m2,b1,b2,cat,val,log=False,label=''):

    plt.figure()
    plt.errorbar(x,e1,yerr=e1err,marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
    plt.errorbar(x+(x[1]-x[0])/5.,e2,yerr=e2err,marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
    plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})
    plt.errorbar(x,m1*x+b1,marker='',linestyle='-',color='r')
    plt.errorbar(x,m2*x+b2,marker='',linestyle='-',color='b')
    plt.axhline(.004,color='k')
    plt.axhline(-.004,color='k')
    plt.ylabel(r'$\langle e \rangle$')
    if config.log_val.get(val,None):
      plt.xlabel('log '+config.lbl.get(val,None))
    else:
      plt.xlabel(config.lbl.get(val,None))
    y1=np.min(np.minimum(e1,e2))   
    y2=np.max(np.maximum(e1,e2))
    plt.ylim((np.min([y1-(y2-y1)/10.,-.005]),np.max([y2+(y2-y1)/10.,.005])))
    plt.minorticks_on()
    plt.savefig('plots/split/lin_split_'+cat.name+'_'+val+'_bs-'+str(cat.bs)+label+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_2pt_split_sub(cat,val,split,n,yl,xi,i,log):

      plt.figure(0)
      ax=plt.subplot(3,3,n)
      ax.fill_between([1,500],-1,1,facecolor='gray',alpha=0.25)
      ax.fill_between([1,500],-2,2,facecolor='gray',alpha=0.2)
      plt.errorbar(xi[0],np.zeros((len(xi[0]))),marker='',linestyle='-',color='k',alpha=.8)
      plt.errorbar(xi[0],xi[10][i]*np.ones((len(xi[0]))),marker='',linestyle='-',color='b')
      plt.errorbar(xi[0],(xi[1][i]-xi[2][i])/xi[2][i],yerr=xi[7][i]/xi[2][i],marker='v',linestyle='',color='g')
      plt.errorbar(xi[0],xi[12][i]*np.ones((len(xi[0]))),marker='',linestyle='-',color='g')
      plt.errorbar(xi[0]*1.2,(xi[3][i]-xi[2][i])/xi[2][i],yerr=xi[9][i]/xi[2][i],marker='^',linestyle='',color='b')
      plt.errorbar(xi[0],xi[11][i]*np.ones((len(xi[0]))),marker='',linestyle='-',color='r')
      plt.errorbar(xi[0]*1.1,(xi[3][i]-xi[1][i])/xi[2][i],yerr=xi[8][i]/xi[2][i],marker='o',linestyle='',color='r')
      plt.xlabel(r'$\theta$ (arcmin)')
      plt.xscale('log')
      plt.ylim(-3,3)
      plt.xlim(1,500)
      ax.set_xticklabels([])
      plt.ylabel(r'$\Delta '+yl+r'/\sigma_{\Delta '+yl+r'}$')

      ax=plt.subplot(3,3,3+n)
      plt.errorbar(xi[0],xi[0]*xi[2][i],yerr=xi[0]*xi[5][i],marker='o',linestyle='',color='r',label='All (upper-lower in top)')
      s=config.lbl.get(val,None)
      if log:
        s='log '+s
      plt.errorbar(xi[0]*1.1,xi[0]*xi[1][i],yerr=xi[0]*xi[4][i],marker='v',linestyle='',color='g',label=s+r'$<$'+str(np.around(split,2)))
      plt.errorbar(xi[0]*1.2,xi[0]*xi[3][i],yerr=xi[0]*xi[6][i],marker='^',linestyle='',color='b',label=s+r'$>$'+str(np.around(split,2)))
      plt.xlabel(r'$\theta$ (arcmin)')
      plt.xscale('log')
      if n==1:
        leg=plt.legend(loc='upper left',ncol=1, frameon=False,prop={'size':12},framealpha=0.2)
      ax.set_xticklabels([])
      plt.ylabel(r'$\theta\times'+yl+r'$')
      plt.xlim(1,500)

      ax=plt.subplot(3,3,6+n)
      plt.errorbar(xi[0],xi[2][i],yerr=xi[5][i],marker='o',linestyle='',color='r')
      plt.errorbar(xi[0]*1.1,xi[1][i],yerr=xi[4][i],marker='v',linestyle='',color='g')
      plt.errorbar(xi[0]*1.2,xi[3][i],yerr=xi[6][i],marker='^',linestyle='',color='b')
      plt.xlabel(r'$\theta$ (arcmin)')
      plt.ylabel(r'$'+yl+r'$')
      plt.xscale('log')
      plt.yscale('log')
      plt.xlim(1,500)

      return

  @staticmethod
  def plot_2pt_split(xi,gt,cat,val,split,log):

    plt.figure(0,figsize=(15,10))

    for i in range(3):

      if i==0:
        plot_methods.plot_2pt_split_sub(cat,val,split,i+1,r'\xi_{+}',xi,0,log)
      elif i==1:
        plot_methods.plot_2pt_split_sub(cat,val,split,i+1,r'\xi_{-}',xi,1,log)
      elif i==2:
        plot_methods.plot_2pt_split_sub(cat,val,split,i+1,r'\gamma_{t}',gt,0,log) 

    plt.minorticks_on()
    plt.subplots_adjust(hspace=0,wspace=.4)
    plt.savefig('plots/split/2pt_split_'+cat.name+'_'+val+'.png', bbox_inches='tight')
    plt.close(0)

    return

  @staticmethod
  def plot_field_corr(cat,theta,out,err,label):

    plt.figure()
    plt.errorbar(theta,theta*out[0],yerr=theta*err[0],marker='o',linestyle='',color='r',label=r'$\gamma_{t}$')
    plt.errorbar(theta,theta*out[2],yerr=theta*err[2],marker='o',linestyle='',color='b',label=r'$\gamma_{x}$')
    if 'chip' not in label:
      plt.axvline(x=5.25*60, linewidth=1, color='k')
    elif 'corner' in label:
      plt.axvline(x=0.75*60, linewidth=1, color='k')
      plt.axvline(x=0.15*60, linewidth=1, color='k')
      plt.axvline(x=0.765*60, linewidth=1, color='k')
    elif 'centre' in label:
      plt.axvline(x=0.75*60/2., linewidth=1, color='k')
      plt.axvline(x=0.15*60/2., linewidth=1, color='k')
      plt.axvline(x=0.765*60/2., linewidth=1, color='k')
    plt.ylabel(r'$\theta \times\gamma$')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylim((-.005,.005))
    plt.xscale('log')
    plt.minorticks_on()
    plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
    plt.savefig('plots/xi/field_'+label+'_'+cat.name+'.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.errorbar(theta[out[0]>0],out[0][out[0]>0],yerr=err[0][out[0]>0],marker='o',linestyle='',color='r',label=r'$\gamma_{t}$')
    plt.errorbar(theta[out[2]>0],out[2][out[2]>0],yerr=err[2][out[2]>0],marker='o',linestyle='',color='b',label=r'$\gamma_{x}$')
    if np.sum(out[0]<0):
      plt.errorbar(theta[out[0]<0],-out[0][out[0]<0],yerr=err[0][out[0]<0],marker='x',linestyle='',color='r',label='')
    if np.sum(out[2]<0):
      plt.errorbar(theta[out[2]<0],-out[2][out[2]<0],yerr=err[2][out[2]<0],marker='x',linestyle='',color='b',label='')
    if 'chip' not in label:
      plt.axvline(x=5.25*60, linewidth=1, color='k')
    elif 'corner' in label:
      plt.axvline(x=0.75*60, linewidth=1, color='k')
      plt.axvline(x=0.15*60, linewidth=1, color='k')
      plt.axvline(x=0.765*60, linewidth=1, color='k')
    elif 'center' in label:
      plt.axvline(x=0.75*60/2., linewidth=1, color='k')
      plt.axvline(x=0.15*60/2., linewidth=1, color='k')
      plt.axvline(x=0.765*60/2., linewidth=1, color='k')
    plt.ylabel(r'$\gamma$')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.yscale('log')
    plt.xscale('log')
    plt.minorticks_on()
    plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
    plt.savefig('plots/xi/field_'+label+'_'+cat.name+'_log.png', bbox_inches='tight')
    plt.close()

    return


  @staticmethod
  def plot_field_corr2(cat,theta,out,err,out2,err2,label):

    plt.figure()
    plt.errorbar(theta,theta*out[0],yerr=theta*err[0],marker='o',linestyle='',color='r',label=r'$e_1$')
    plt.errorbar(theta,theta*out2[0],yerr=theta*err2[0],marker='o',linestyle='',color='b',label=r'$e_2$')
    if 'chip' not in label:
      plt.axvline(x=5.25*60, linewidth=1, color='k')
    elif 'corner' in label:
      plt.axvline(x=0.75*60, linewidth=1, color='k')
      plt.axvline(x=0.15*60, linewidth=1, color='k')
      plt.axvline(x=0.765*60, linewidth=1, color='k')
    elif 'centre' in label:
      plt.axvline(x=0.75*60/2., linewidth=1, color='k')
      plt.axvline(x=0.15*60/2., linewidth=1, color='k')
      plt.axvline(x=0.765*60/2., linewidth=1, color='k')
    plt.ylabel(r'$\langle e \rangle$')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.ylim((-.005,.005))
    plt.xscale('log')
    plt.minorticks_on()
    plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
    plt.savefig('plots/xi/field_'+label+'_'+cat.name+'_mean_e.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def fig_create_xi(cat,catb,corr,theta,out,err,k,ga,gb):

    plt.figure()
    if (corr=='GG')|(corr=='NG'):
      plt.errorbar(theta,theta*out[0],yerr=theta*err[0],marker='o',linestyle='',color='r',label=r'$\xi_{+}$')
      if corr=='GG':
        plt.errorbar(theta,theta*out[1],yerr=theta*err[1],marker='o',linestyle='',color='b',label=r'$\xi_{-}$')
    plt.ylabel(r'$\theta \xi$')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.xscale('log')

    plt.minorticks_on()
    plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
    if catb is None:
      plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_bs-'+str(cat.bs)+'.png', bbox_inches='tight')
    else:
      plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_'+catb.name+'_bs-'+str(cat.bs)+'.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    if (corr=='GG')|(corr=='NG'):
      plt.errorbar(theta[out[0]>0],out[0][out[0]>0],yerr=err[0][out[0]>0],marker='o',linestyle='',color='r',label=r'$\xi_{+}$')
      if corr=='GG':
        plt.errorbar(theta[out[1]>0],out[1][out[1]>0],yerr=err[1][out[1]>0],marker='o',linestyle='',color='b',label=r'$\xi_{-}$')
      if np.sum(out[0]<0):
        plt.errorbar(theta[out[0]<0],-out[0][out[0]<0],yerr=err[0][out[0]<0],marker='x',linestyle='',color='r',label='')
      if corr=='GG':
        if np.sum(out[1]<0):
          plt.errorbar(theta[out[1]<0],-out[1][out[1]<0],yerr=err[1][out[1]<0],marker='x',linestyle='',color='b',label='')
    plt.ylabel(r'$\theta \xi$')
    plt.xlabel(r'$\theta$ (arcmin)')
    plt.yscale('log')
    plt.xscale('log')

    plt.minorticks_on()
    plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
    if catb is None:
      plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_bs-'+str(cat.bs)+'_log.png', bbox_inches='tight')
    else:
      plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_'+catb.name+'_bs-'+str(cat.bs)+'_log.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def fig_create_xi_alpha(cat,theta,gpout,ppout,gperr,pperr,alphap,alpham,alpha0):

    ax=plt.subplot(2,1,1)
    plt.errorbar(theta,gpout[0],yerr=gperr[0],marker='.',linestyle='',color='r',label=r'$\xi^{gp}$')
    plt.errorbar(theta,ppout[0],yerr=pperr[0],marker='.',linestyle='',color='b',label=r'$\xi^{pp}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((.1,500))
    # plt.ylim(5e-7,6e-4)
    ax.set_xticklabels([])
    ax.minorticks_on()
    plt.legend(loc='upper right',ncol=2, frameon=True,prop={'size':12})
    plt.ylabel(r'$|\xi_{+}|$')

    ax=plt.subplot(2,1,2)
    plt.errorbar(theta,alphap.real,marker='',linestyle='-',color='b',label=r'$\alpha_{p}$')
    plt.errorbar(theta,alphap.imag,marker='',linestyle=':',color='b',label=r'$\alpha_{p}$')
    plt.errorbar(theta,alpham.real,marker='',linestyle='-',color='r',label=r'$\alpha_{m}$')
    plt.errorbar(theta,alpham.imag,marker='',linestyle=':',color='r',label=r'$\alpha_{m}$')
    plt.errorbar(theta,alpha0,marker='',linestyle='-',color='k',label=r'$\alpha_{0}$')
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim((.1,500))
    plt.ylim(-.5,.5)
    plt.ylabel(r'$\alpha$')
    plt.xlabel(r'$\theta$ (arcmin)')    
    plt.legend(loc='upper left',ncol=2, frameon=True,prop={'size':12})
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('plots/xi/xi_alpha_'+cat.name+'_bs-'+str(cat.bs)+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_pz_sig8(test,pz0,label='',boot=False,ylim=0.3,noparam=True):

    from astropy.table import Table
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    colors=['k','r','g','b','c','y']
    col=['r','g','b','c','y','b']

    fig, ax = plt.subplots()
    for i,pz1 in enumerate(pz0):

      theta=np.loadtxt(config.pztestdir+test+'/out/sim_data_notomo_spec_'+pz1.name+'/ell.txt')

      cov=np.loadtxt(config.pztestdir+test+'/out/sim_data_notomo_spec_'+pz1.name+'/covmat.txt')
      cov=np.sqrt(cov[4:])/np.sqrt(len(theta))
      xi00=np.loadtxt(config.pztestdir+test+'/out/sim_data_notomo_spec_'+pz1.name+'/data.txt')[2:]
      sig_notomo=np.mean(cov/xi00)

      xi00=np.loadtxt(config.pztestdir+test+'/out/sim_data_spec_'+pz1.name+'/data.txt')[:,2:]
      tomobins=len(xi00)
      cov=np.loadtxt(config.pztestdir+test+'/out/sim_data_spec_'+pz1.name+'/covmat.txt')
      tmp=cov[(cov[:,2]==cov[:,0])&(cov[:,3]==cov[:,1])]
      cov=np.zeros((tomobins,len(theta)))
      for k in range(tomobins):
        cov[k,:]=np.sqrt(tmp[k,4:])/np.sqrt(len(theta))
      sig_tomo=np.mean(cov/xi00,axis=1)

      ratio=np.zeros(tomobins+1)
      sig=np.zeros(tomobins+1)
      data=np.loadtxt(config.pztestdir+test+'/out/sim_data_notomo_'+pz1.name+'/data.txt')[2:]
      data0=np.loadtxt(config.pztestdir+test+'/out/sim_data_notomo_spec_'+pz1.name+'/data.txt')[2:]
      ratio[0]=np.mean((data[:]-data0[:])/data0[:])
      if boot:
        sig[0]=pz.pz_spec_validation.calc_bootstrap(pz1,test,config.pztestdir,tomobins,notomo=True)
      else:
        sig[0]=0. 

      data=np.loadtxt(config.pztestdir+test+'/out/sim_data_'+pz1.name+'/data.txt')[:,2:]
      data0=np.loadtxt(config.pztestdir+test+'/out/sim_data_spec_'+pz1.name+'/data.txt')[:,2:]
      for bin in range(tomobins):
        ratio[bin+1]=np.mean((data[bin,:]-data0[bin,:])/data0[bin,:])
      if boot:
        sig[1:]=pz.pz_spec_validation.calc_bootstrap(pz1,test,config.pztestdir,tomobins,notomo=False)
      else:
        sig[1:]=0.

      plt.errorbar(np.arange(tomobins+1),ratio,yerr=sig,marker='o',linestyle='',color=col[i],label=pz1.name)

    plt.fill_between(0-.4+.8*np.arange(100)/100.,-sig_notomo*np.ones(100),sig_notomo*np.ones(100),interpolate=True,color='k',alpha=0.2)
    for i in range(tomobins):
      plt.fill_between(i+1-.4+.8*np.arange(100)/100.,-sig_tomo[i]*np.ones(100),sig_tomo[i]*np.ones(100),interpolate=True,color='k',alpha=0.2)
    plt.plot(np.arange(tomobins+3)-1,np.zeros((tomobins+3)), marker='', linestyle='-',color='k',label='')
    ax.xaxis.set_major_locator(MultipleLocator(1.))
    ax.yaxis.set_major_locator(MultipleLocator(ylim/3.))
    plt.xticks(np.arange(tomobins+3)-1,np.append([' ','2D','11','21','22','31','32','33','41','42','43','44','51','52','53','54','55','61','62','63','64','65','66'][:tomobins+2],[' ']))
    plt.ylim((-ylim,ylim))
    plt.xlim((-1,tomobins+1))
    plt.ylabel(r'$\Delta C_{\ell}/C_{\ell}(\textrm{spec})$')
    plt.xlabel(r'Bin pairs')

    props = dict(boxstyle='square', lw=1.2,facecolor='white', alpha=1.)

    ax.text(0.82, 0.95, label, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

    plt.legend(loc='upper left',ncol=2, frameon=True,prop={'size':12})
    plt.savefig('plots/photoz/pz_xi_'+test+'.png',bbox_inches='tight')
    plt.close()

    if noparam:
      return

    param='sigma8'

    lines=[]
    for i,pz1 in enumerate(pz0):
      if boot:
        sig8_sig=pz.pz_spec_validation.calc_bootstrap_sig8(pz1,test,config.pztestdir,param,notomo=False)
        print sig8_sig
      else:
        sig8_sig=0.

      mean_notomo0 = Table.read(config.pztestdir+test+'/out/sim_data_notomo_'+pz1.name+'/means.txt', format='ascii')
      mean0 = Table.read(config.pztestdir+test+'/out/sim_data_'+pz1.name+'/means.txt', format='ascii')
      #print mean0
      for row in mean0:
        if param in row['parameter']:
          print 'true'
          mean=row['mean']
          low=row['std_dev']
          high=mean+row['std_dev']
          print pz1.name,(mean-1.)/low
      lines.append((pz1.name,mean,low,high,low,high,r'$\sigma_{8}$  Tomographic '+pz1.name))
      for row in mean_notomo0:
        if param in row['parameter']:
          mean_notomo=row['mean']
          low_notomo=row['std_dev']
          high_notomo=mean_notomo+row['std_dev']
      lines.append((pz1.name,mean_notomo,low_notomo,high_notomo,low_notomo,high_notomo,r'$\sigma_{8}$  2D '+pz1.name))
    #print lines
      
    #From StackOverflow  thanks, Paul Ivanov
    def rainbow_text(x,y,ls,lc,**kw):
        """
        Take a list of strings ``ls`` and colors ``lc`` and place them next to each
        other, with text ls[i] being shown in color lc[i].

        This example shows how to do both vertical and horizontal text, and will
        pass all keyword arguments to plt.text, so you can set the font size,
        family, etc.
        """
        from matplotlib.offsetbox import HPacker, TextArea, AnnotationBbox
        
        ax = pylab.gca()
        texts = [TextArea(s,textprops=dict(color=c, family='serif')) for (s,c) in zip(ls,lc)]
        txt = HPacker(children=texts, 
                align="baseline", 
                pad=0, sep=0) 

        def txt_offset(*kl):
          return ax.transData.transform_point((x, y))
        txt.set_offset(txt_offset)

        
        ax.add_artist(txt)

    bold_lines = [
    "Fiducial DES-SV cosmic shear",
    "CFHTLenS (H13) original conservative scales",
    "Planck (TT+LowP)",
    ]

    def plot_set(lines, color, start):
      col=['r','g','b','c','y','b']
      n=len(lines)
      yscale = 5.
      y = np.arange(n)
      x = [l[1] for l in lines]
      e=np.zeros((2,n))
      e95=np.zeros((2,n))
      for i,l in enumerate(lines):
        e[0,i] = l[1]-l[3]
        e[1,i] = l[2]-l[1]
        e95[0,i] = l[1]-l[5]
        e95[1,i] = l[4]-l[1]

      for i,l in enumerate(lines):
        #print i,l
        pylab.errorbar(l[1],(start-i)*yscale,xerr=l[2], fmt=col[i/2]+'.')

      for i,line in enumerate(lines):
        text = line[6]
        if text in bold_lines:
          weight = "bold"
          print text, " bold"
        else:
          weight = "normal"
        if text.startswith("Planck") and text.endswith("DES-SV"):
          text = [text[:-6], " "+text[-6:]]
          colors = [color, 'b']
          rainbow_text(0.95+0.02, (start-i)*yscale-0.5,text, colors, family='serif', weight=weight)
        else:
          pylab.text(1.15+0.005, (start-i)*yscale-0.5, text, color=col[i/2], family='serif', weight=weight)
        #pylab.plot([line[3]+0.01, 1.04], [(start-i)*yscale, (start-i)*yscale], ':', color='darkgray')


    pylab.figure(figsize=(12,12))

    pylab.axvspan(1.-sig8_sig, 1.+sig8_sig,alpha=0.2, color='gray')
    plot_set(lines, 'b', len(lines))
    pylab.ylim(-2, 52)
    pylab.xlim(0.85, 1.25)


    pylab.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',      # ticks along the bottom edge are off
        right='off',         # ticks along the top edge are off
        labelleft='off' # labels along the bottom edge are off
    )

    pylab.tick_params(
        axis='x',          # changes apply to the x-axis
        which='major',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        length=10,
        labelsize=20,
    )

    pylab.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        length=5,
    )
    a = pylab.gca()
    fontProperties = {'family':'serif', 'size':20}
    a.set_xticklabels(a.get_xticks(), fontProperties)


    pylab.xlabel("$\sigma_8$", fontsize=20)
    pylab.savefig('plots/photoz/pz_sig8_'+test+'.png')
    pylab.close()


    # print 'tomo'

    # np.set_printoptions(precision=2)
    # for i in xrange(4):
    #   for j in xrange(4):
    #     if j==0:
    #       print '\\'+SVA1PZPaper.get_pz_name(i+1,test).lower()+' & ',
    #     if i==j:
    #       print '- & ',
    #     else:
    #       print str(np.around(np.log(sum0a[i]/sum0a[j]),decimals=2))+' ('+str(np.around(np.log(sum0b[i]/sum0b[j]),decimals=2))+') & ',

    #   print ''

    #   print ''

    # np.set_printoptions(edgeitems=3,infstr='inf',linewidth=75, nanstr='nan', precision=8,suppress=False, threshold=1000, formatter=None)

    return 

  @staticmethod
  def plot_nofz(pz0,test,spec=True):

    col=['k','r','g','b','r','g','b']
    ax=plt.subplot(2,1,1)
    for i in xrange(len(pz0.pz)-1):
      plt.plot(pz0.bin,pz0.pz[i+1,:],color=col[i+1],linestyle='-',linewidth=1.,drawstyle='steps-mid',label='')
      plt.axvline(x=np.average(pz0.bin,weights=pz0.pz[i+1,:]), ymin=0., ymax = 1, linewidth=2, color=col[i+1])
      print i+1,np.average(pz0.bin,weights=pz0.pz[i+1,:])-np.average(pz0.bin,weights=pz0.spec[i+1,:])
      if spec:
        plt.plot(pz0.bin,pz0.spec[i+1,:],color=col[i+1],linestyle=':',linewidth=3.,drawstyle='steps-mid',label='')
        plt.axvline(x=np.average(pz0.bin,weights=pz0.spec[i+1,:]), ymin=0., ymax = 1, linewidth=1, color=col[i+1])
    ax.set_xticklabels([])
    plt.ylabel(r'$n(z)$')
    ax=plt.subplot(2,1,2)
    plt.plot(pz0.bin,pz0.pz[0,:],color=col[0],linestyle='-',linewidth=1.,drawstyle='steps-mid',label='')
    plt.axvline(x=np.average(pz0.bin,weights=pz0.pz[0,:]), ymin=0., ymax = 1, linewidth=2, color=col[0])
    print 0,np.average(pz0.bin,weights=pz0.pz[0,:])-np.average(pz0.bin,weights=pz0.spec[0,:])
    if spec:
      plt.plot(pz0.bin,pz0.spec[0,:],color=col[0],linestyle=':',linewidth=3.,drawstyle='steps-mid',label='')
      plt.axvline(x=np.average(pz0.bin,weights=pz0.spec[0,:]), ymin=0., ymax = 1, linewidth=1, color=col[0])
    plt.xlabel(r'$z$')
    plt.ylabel(r'$n(z)$')
    # plt.xscale('log')
    # plt.xlim((0,5.))
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig('plots/photoz/pz_nofz_'+test+'.png',bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_nofz_comp_pz(pzlist,label='',spec=True,notomo=False):

    col=['k','r','g','b','c','r','g']

    plt.figure(figsize=(8,16))
    for i in xrange(len(pzlist[0].pz)):
      if notomo&(i>0):
        continue
      ax=plt.subplot(len(pzlist[0].pz),1,i+1)
      for pz0i,pz0 in enumerate(pzlist):
        if spec:
          plt.plot(pz0.bin,pz0.pz[i,:],color=col[pz0i+1],linestyle='-',linewidth=1.,drawstyle='steps-mid',label=pz0.name)
          plt.axvline(x=np.average(pz0.bin,weights=pz0.pz[i,:]), ymin=0., ymax = 1, linewidth=2, color=col[pz0i+1])
          if pz0i==0:
            plt.plot(pzlist[0].bin,pzlist[0].spec[i,:],color=col[0],linestyle=':',linewidth=2.,drawstyle='steps-mid',label='')
            plt.axvline(x=np.average(pz0.bin,weights=pzlist[0].spec[i,:]), ymin=0., ymax = 1, linewidth=2, color=col[0])
        else:
          plt.plot(pz0.bin,pz0.pz[i,:],color=col[pz0i],linestyle='-',linewidth=1.,drawstyle='steps-mid',label=pz0.name)
          plt.axvline(x=np.average(pz0.bin,weights=pz0.pz[i,:]), ymin=0., ymax = 1, linewidth=2, color=col[pz0i])          
        print i,np.average(pz0.bin,weights=pz0.pz[i,:])-np.average(pzlist[0].bin,weights=pzlist[0].spec[0,:])
      props = dict(boxstyle='square', lw=1.2,facecolor='white', alpha=1.)
      if i==0:
        ax.text(0.73, 0.95, 'Non-tomographic', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
      else:
        ax.text(0.9, 0.95, 'Bin '+str(i), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
      plt.ylabel(r'$n(z)$')
      ax.minorticks_on()
      if i<len(pzlist[0].pz)-1:
        ax.set_xticklabels([])
      # plt.xscale('log')
      plt.xlim((0,1.5))
    plt.xlabel(r'$z$')
    plt.legend(loc='upper left',ncol=2, frameon=True,prop={'size':12})
    plt.subplots_adjust(hspace=0,wspace=0)
    if label!='':
      label+='_'
    plt.savefig('plots/photoz/pz_nofz_'+label+str(len(pzlist[0].pz)-1)+'_weight-'+str(pz0.wt)+'.png',bbox_inches='tight')
    plt.close()

    return


  @staticmethod
  def sig_crit_spec(pz0,bins,bootstrap,point,lensbins,edge):

    gs = gridspec.GridSpec(bins+1,1)
    plt.figure(figsize=(8,20))

    col=['k','r','g','b','c','r','g']

    for ibin in xrange(bins+1):
      print ibin
      ax1=plt.subplot(gs[ibin,0])        
      ax1.minorticks_on()
      for ipz,pz in enumerate(pz0):

        z_l=np.load('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_z.npy')
        if pz.wt:
          data=np.load('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_dat_weighted.npy')
        else:
          data=np.load('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_dat.npy')
        if bootstrap:
          if pz.wt:
            var=np.load('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_var_weighted.npy')
          else:
            var=np.load('text/spec_'+str(ibin)+'_'+pz.name+'_invsigcrit_var.npy')
        else:
          var=np.zeros(len(data))

        ls='o'
        if (ibin==0):
          ax1.errorbar(z_l[z_l<edge[ibin-1]]+1.1*ipz*(z_l[1]-z_l[0])/5.,data[z_l<edge[ibin-1]],yerr=var[z_l<edge[ibin-1]],linestyle='', marker=ls,color=col[ipz+1],label=pz.name)
        else:
          ax1.errorbar(z_l[z_l<edge[ibin-1]]+1.1*ipz*(z_l[1]-z_l[0])/5.,data[z_l<edge[ibin-1]],yerr=var[z_l<edge[ibin-1]],linestyle='', marker=ls,color=col[ipz+1],label='')
        ls='-'
        ax1.plot([-.05,1.3],[0,0],marker='',color='k',label='')

        if ibin==0:
          bintext=r'$'+str(np.around(np.min(pz.spec_full),2))+'<z_{\mathrm{source}}<'+str(np.around(np.max(pz.spec_full),2))+'$'
        else:
          bintext=r'$'+str(np.around(edge[ibin-1],2))+'<z_{\mathrm{source}}<'+str(np.around(edge[ibin],2))+'$'
        props = dict(boxstyle='square', lw=1.2,facecolor='white', alpha=1.)
        ax1.text(0.7, 0.95, bintext, transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)
        
        if (ibin==bins):
          plt.xlabel(r'$z_{\mathrm{lens}}$')
        else:
          ax1.set_xticklabels([])

        plt.ylabel(r'$\Delta\langle \Sigma_{\mathrm{crit}}^{-1}\rangle/\langle \Sigma_{\mathrm{crit}}^{-1}\rangle_{\mathrm{spec}}$')
        plt.ylim((-1.,1.))
        ax1.xaxis.set_major_locator(MultipleLocator(.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(.05))              
        plt.xlim((0.,1.1))
        plt.legend(loc='upper left',ncol=2, fancybox=True, shadow=True)

    gs.update(wspace=0.,hspace=0.)
    plt.savefig('plots/photoz/spec_point-'+str(point)+'_bins-'+str(bins)+'_bootstrap-'+str(bootstrap)+'_invsigcrit_weight-'+str(pz.wt)+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def sig_crit_spec2(pz0,point):

    gs = gridspec.GridSpec(len(pz0)+1,1)
    plt.figure(figsize=(8,16))

    col=['k','r','g','b','c','r','g']

    for ipz,pz1 in enumerate(pz0):
      ax1=plt.subplot(gs[ipz,0])
      ax1.minorticks_on()
      z_l=np.load('text/spec_'+pz1.name+'_invsigcrit_zl.npy')
      if pz1.wt:
        data=np.load('text/spec_'+pz1.name+'_invsigcrit_dat_weighted.npy')
      else:
        data=np.load('text/spec_'+pz1.name+'_invsigcrit_dat.npy')

      plt.imshow(data, interpolation='bilinear', origin='lower', extent=[np.min(z_l), np.max(z_l), np.min(z_l), np.max(z_l)])
      plt.xlabel(r'$z_{\mathrm{lens}}$')
      plt.ylabel(r'$z_{\mathrm{source}}$')
      if ipz==0:
        plt.title(r'$\Delta\langle \Sigma_{\mathrm{crit}}^{-1}\rangle/\langle \Sigma_{\mathrm{crit}}^{-1}\rangle_{\mathrm{spec}}$')

      cb = plt.colorbar(orientation='vertical')
      cb.set_label(pz1.name)
    gs.update(wspace=0.,hspace=0.)
    plt.savefig('plots/photoz/spec_point-'+str(point)+'_invsigcrit_weight-'+str(pz1.wt)+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def spec_loop_hist2d(pz0):


    for i,x in enumerate(['mc','mean','peak']):
      gs = gridspec.GridSpec(len(pz0)+1,1)
      plt.figure(figsize=(8,16))
      for ipz,pz in enumerate(pz0):
        ax1=plt.subplot(gs[ipz,0])
        ax1.minorticks_on()
        plt.hist2d(pz.spec_full,getattr(pz,'z_'+x+'_full'),weights=pz.w,bins=200,range=((0,1.5),(0,1.5)))
        plt.ylabel(r'$z_\mathrm{'+x+'}$')
        if (ipz==len(pz0)-1):
          plt.xlabel(r'$z_{\mathrm{spec}}$')
        else:
          ax1.set_xticklabels([])
        cb = plt.colorbar(orientation='vertical')
        cb.set_label(pz.name)
      gs.update(wspace=0.,hspace=0.)
      plt.savefig('plots/photoz/spec_hist_'+x+'_weighted.png', bbox_inches='tight')
      plt.close()

      gs = gridspec.GridSpec(len(pz0)+1,1)
      plt.figure(figsize=(8,16))
      for ipz,pz in enumerate(pz0):
        ax1=plt.subplot(gs[ipz,0])
        ax1.minorticks_on()
        plt.hist2d(pz.spec_full,getattr(pz,'z_'+x+'_full'),bins=200,range=((0,1.5),(0,1.5)))
        if (ipz==len(pz0)-1):
          plt.xlabel(r'$z_{\mathrm{spec}}$')
        else:
          ax1.set_xticklabels([])
        plt.ylabel(r'$z_\mathrm{'+x+'}$')
        cb = plt.colorbar(orientation='vertical')
        cb.set_label(pz.name)
      gs.update(wspace=0.,hspace=0.)
      plt.savefig('plots/photoz/spec_hist_'+x+'.png', bbox_inches='tight')
      plt.close()

    return


  @staticmethod
  def plot_pzrw(cat,pz,bins,w,label,edge):

    plt.figure(0,figsize=(5,10))
    ax=plt.subplot(1,2,1)

    col=['r','b','g']
    plt.hist(pz,bins=100,color='k',linestyle=('solid'),linewidth=1.,label='Full sample',histtype='step',normed=True)
    for i in range(cat.sbins):
      plt.hist(pz[bins==i],bins=100,color=col[i],linestyle=('solid'),linewidth=1.,label=r'$'+"{0:.2f}".format(edge[i])+'<$'+label+'$<'+"{0:.2f}".format(edge[i+1])+'$',histtype='step',weights=w[bins==i],normed=True)
      plt.hist(pz[bins==i],bins=100,color=col[i],linestyle=('dashed'),linewidth=1.,label='',histtype='step',normed=True)
    plt.legend(loc='upper right')
    #plt.ylim((0,2.5))
    plt.xlabel('z')
    plt.ylabel('n(z)')

    ax=plt.subplot(1,2,2)

    plt.axvline(x=1)
    for i in range(cat.sbins):
      plt.hist(w[bins==i],bins=50,alpha=.5,color=col[i],label=r'$'+"{0:.2f}".format(edge[i])+'<$'+label+'$<'+"{0:.2f}".format(edge[i+1])+'$',normed=True,range=[-1, 5])
    plt.legend(loc='upper right')
    #plt.xlim((-1,4))
    plt.xlabel('w')
    plt.ylabel('n(w)')
    plt.savefig('plots/split/pzrw_'+cat.name+'_'+label+'.png', bbox_inches='tight')
    plt.close()

    return


  @staticmethod
  def plot_IA(r,out,err,label):

    col=['r','b','g','c']
    name=['gp','gx','ee','xx']
    for i in xrange(len(out)):
      plt.errorbar(r[out[i]>0]*(1.+.2*i),out[i][out[i]>0],yerr=err[i][out[i]>0],color=col[i],linestyle='',marker='o',label=name[i])
      plt.errorbar(r[out[i]<0]*(1.+.2*i),-out[i][out[i]<0],yerr=err[i][out[i]<0],color=col[i],linestyle='',marker='s',label='')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.xlabel('R [Mpc/h]')
    plt.savefig('plots/IA_'+label+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def imshow_symlog(my_matrix, logthresh=3):
    # From Freddie Nfbnm on stack overflow  http://stackoverflow.com/questions/11138706/

    plt.imshow( my_matrix,interpolation='nearest',origin='lower',vmin=np.min(my_matrix), vmax=np.max(my_matrix),norm=matplotlib.colors.SymLogNorm(10**-logthresh) )

    maxlog=int(np.ceil( np.log10(np.max(my_matrix)) ))
    minlog=int(np.ceil( np.log10(-np.min(my_matrix)) ))

    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]+[0.0]+[(10**x) for x in xrange(-logthresh,maxlog+1)] )

    plt.colorbar(ticks=tick_locations)

    return
    