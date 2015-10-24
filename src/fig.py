import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
plt.style.use('/home/troxel/SVA1/SVA1StyleSheet.mplstyle')
from matplotlib.colors import LogNorm
import pylab

import catalog
import config
import pz

class plot_methods(object):
  """
  Plotting routines used by various modules.
  """

  @staticmethod
  def plot_hist(x1,bins=500,name='',label=''):

    plt.figure()
    plt.hist(x1,bins=bins)
    plt.ylabel(r'$n$')
    s=config.lbl.get(label,None)
    if config.log_val.get(label,None):
      s='log '+s
    plt.xlabel(s)
    plt.minorticks_on()
    plt.savefig('plots/hist/hist_'+name+'_'+label+'.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_2D_hist(x1,y1,bins=500,xname='',yname='',xlabel='',ylabel=''):

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
  def plot_hexbin(x1,cat,mask=None,bins=500,name='',label=''):

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    plt.figure()
    plt.hexbin(cat.ra[mask],cat.dec[mask],x1,gridsize=bins, cmap=plt.cm.afmhot)
    cb = plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.minorticks_on()
    s=config.lbl.get(label,None)
    if config.log_val.get(label,None):
      s='log '+s
    cb.set_label(s)
    plt.savefig('plots/footprint/hexbin_'+name+'_'+label+'.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.hexbin(cat.ra[mask],cat.dec[mask],x1,gridsize=bins,bins='log', cmap=plt.cm.afmhot)
    cb = plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.minorticks_on()
    s=config.lbl.get(label,None)
    if config.log_val.get(label,None):
      s='log '+s
    cb.set_label(s)
    plt.savefig('plots/footprint/hexbin_'+name+'_'+label+'_log.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_footprint(cat,mask=None,label='',bins=100):

    dra=np.max(cat.ra[mask])-np.min(cat.ra[mask])
    ddec=np.max(cat.dec[mask])-np.min(cat.dec[mask])

    plt.figure()
    a=plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=((np.min(cat.ra[mask])-.1*dra,np.max(cat.ra[mask])+.1*dra),(np.min(cat.dec[mask])-.1*ddec,np.max(cat.dec[mask])+.1*ddec)),normed=True, cmap=plt.cm.afmhot)
    cb = plt.colorbar()
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/footprint_'+cat.name+'_'+label+'.png', dpi=500, bbox_inches='tight')
    plt.close()    

    plt.figure()
    plt.hist2d(cat.ra[mask],cat.dec[mask],bins=(int(dra*bins),int(ddec*bins)),range=((np.min(cat.ra[mask])-.1*dra,np.max(cat.ra[mask])+.1*dra),(np.min(cat.dec[mask])-.1*ddec,np.max(cat.dec[mask])+.1*ddec)),normed=True,norm=LogNorm(), cmap=plt.cm.afmhot)
    cb = plt.colorbar()
    plt.gca().set_aspect('equal', 'box')
    plt.savefig('plots/footprint/footprint_'+cat.name+'_'+label+'_log.png', dpi=500, bbox_inches='tight')
    plt.close()    

    return

  @staticmethod
  def plot_whisker(x,y,e1,e2,name='',label='',scale=.01,key=''):

    plt.figure()
    Q = plt.quiver(x,y,e1,e2,units='width',pivot='middle',headwidth=0,width=.0005)
    plt.quiverkey(Q,0.2,0.2,scale,str(scale)+' '+key,labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
    plt.savefig('plots/footprint/whisker_'+name+'_'+label+'.png', dpi=500,bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_lin_split(x,e1,e2,e1err,e2err,m1,m2,b1,b2,cat,val,log=False):

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
    plt.savefig('plots/split/lin_split_'+cat.name+'_'+val+'_bs-'+str(cat.bs)+'.png', bbox_inches='tight')
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
  def fig_create_xi(cat,corr,theta,out,err,k,ga,gb):

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
    plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_bs-'+str(cat.bs)+'.png', bbox_inches='tight')
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
    plt.savefig('plots/xi/xi_'+corr+'_'+ga+'_'+gb+'_'+cat.name+'_bs-'+str(cat.bs)+'_log.png', bbox_inches='tight')
    plt.close()

    return

  @staticmethod
  def plot_pz_sig8(test,label='',boot=False,tomobins=3):

    from astropy.table import Table

    tomobins=np.sum(np.arange(tomobins+1))

    theta=np.loadtxt(config.pzrootdir+test+'/out/sim_data_notomo_spec_skynet/ell.txt')

    cov=np.loadtxt(config.pzrootdir+test+'/out/sim_data_notomo_spec_skynet/covmat.txt')
    cov=np.sqrt(cov[4:])/np.sqrt(7.)
    xi00=np.loadtxt(config.pzrootdir+test+'/out/sim_data_notomo_spec_skynet/data.txt')[2:]
    sig_notomo=np.mean(cov/xi00)

    cov=np.loadtxt(config.pzrootdir+test+'/out/sim_data_spec_skynet/covmat.txt')
    tmp=cov[(cov[:,2]==cov[:,0])&(cov[:,3]==cov[:,1])]
    cov=np.zeros((tomobins,7))
    for i in range(tomobins):
      cov[i,:]=np.sqrt(tmp[i,4:])/np.sqrt(7.)
    xi00=np.loadtxt(config.pzrootdir+test+'/out/sim_data_spec_skynet/data.txt')[:,2:]
    sig_tomo=np.mean(cov/xi00,axis=1)
    print sig_tomo

    ratio=np.zeros(tomobins+1)
    sig=np.zeros(tomobins+1)
    data=np.loadtxt(config.pzrootdir+test+'/out/sim_data_notomo_skynet/data.txt')[2:]
    data0=np.loadtxt(config.pzrootdir+test+'/out/sim_data_notomo_spec_skynet/data.txt')[2:]
    ratio[0]=np.mean((data[:]-data0[:])/data0[:])
    if boot:
      sig[0]=pz.pz_spec_validation.calc_bootstrap(test,config.pzrootdir,tomobins,notomo=True)
    else:
      sig[0]=0. 

    data=np.loadtxt(config.pzrootdir+test+'/out/sim_data_skynet/data.txt')[:,2:]
    data0=np.loadtxt(config.pzrootdir+test+'/out/sim_data_spec_skynet/data.txt')[:,2:]
    for bin in range(tomobins):
      ratio[bin+1]=np.mean((data[bin,:]-data0[bin,:])/data0[bin,:])
    if boot:
      sig[1:]=pz.pz_spec_validation.calc_bootstrap(test,config.pzrootdir,tomobins,notomo=False)
    else:
      sig[1:]=0.

    fig, ax = plt.subplots()

    plt.errorbar(np.arange(tomobins+1),ratio,yerr=sig,marker='o',linestyle='',color='k')

    plt.fill_between(0-.4+.8*np.arange(100)/100.,-sig_notomo*np.ones(100),sig_notomo*np.ones(100),interpolate=True,color='k',alpha=0.2)
    for i in range(tomobins):
      plt.fill_between(i+1-.4+.8*np.arange(100)/100.,-sig_tomo[i]*np.ones(100),sig_tomo[i]*np.ones(100),interpolate=True,color='k',alpha=0.2)
    plt.plot(np.arange(tomobins+3)-1,np.zeros((tomobins+3)), marker='', linestyle='-',color='k',label='')
    ax.xaxis.set_major_locator(MultipleLocator(1.))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    plt.xticks(np.arange(tomobins+3)-1,np.append([' ','2D','11','21','22','31','32','33','41','42','43','44','51','52','53','54','55','61','62','63','64','65','66'][:tomobins+2],[' ']))
    plt.ylim((-.3,.3))
    plt.xlim((-1,tomobins+1))
    plt.ylabel(r'$\Delta C_{\ell}/C_{\ell}(\textrm{spec})$')
    plt.xlabel(r'Bin pairs')

    props = dict(boxstyle='square', lw=1.2,facecolor='white', alpha=1.)

    ax.text(0.82, 0.95, label, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

    plt.savefig('plots/pz_xi_'+test+'.png')
    plt.close()


    param='sigma8'

    if boot:
      sig8_sig=pz.pz_spec_validation.calc_bootstrap_sig8(test,config.pzrootdir,param,notomo=False)
    else:
      sig8_sig=0.

    lines=[]
    mean_notomo0 = Table.read(config.pzrootdir+test+'/out/sim_data_notomo_skynet/means.txt', format='ascii')
    mean0 = Table.read(config.pzrootdir+test+'/out/sim_data_skynet/means.txt', format='ascii')
    low_notomo0 = Table.read(config.pzrootdir+test+'/out/sim_data_notomo_skynet/means.txt', format='ascii')
    low0 = Table.read(config.pzrootdir+test+'/out/sim_data_skynet/means.txt', format='ascii')
    high_notomo0 = Table.read(config.pzrootdir+test+'/out/sim_data_notomo_skynet/means.txt', format='ascii')
    high0 = Table.read(config.pzrootdir+test+'/out/sim_data_skynet/means.txt', format='ascii')
    for row in mean0:
      if param in row['parameter']:
        mean=row['mean']
        low=mean-row['std_dev']
        high=mean+row['std_dev']
    for row in low0:
        low95=mean-row['std_dev']
    for row in high0:
        high95=mean+row['std_dev']
    for row in mean_notomo0:
      if param in row['parameter']:
        mean_notomo=row['mean']
        low_notomo=mean_notomo-row['std_dev']
        high_notomo=mean_notomo+row['std_dev']
    for row in low_notomo0:
        low95_notomo=mean-row['std_dev']
    for row in high_notomo0:
        high95_notomo=mean+row['std_dev']
    lines.append(('skynet',mean,low,high,low95,high95,r'$\sigma_{8}$ - Tomographic'))
    lines.append(('skynet',mean_notomo,low_notomo,high_notomo,low95_notomo,high95_notomo,r'$\sigma_{8}$ - 2D'))
    
    #From StackOverflow - thanks, Paul Ivanov
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
      n=len(lines)
      yscale = 1.5
      y = np.arange(n)
      x = [l[1] for l in lines]
      e=np.zeros((2,n))
      e95=np.zeros((2,n))
      for i,l in enumerate(lines):
        e[0,i] = l[1]-l[3]
        e[1,i] = l[2]-l[1]
        e95[0,i] = l[1]-l[5]
        e95[1,i] = l[4]-l[1]

      pylab.errorbar(x,(start-y)*yscale,xerr=e, fmt=color+'.')

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
          pylab.text(1.05+0.005, (start-i)*yscale-0.5, text, color=color, family='serif', weight=weight)
        pylab.plot([line[3]+0.01, 1.04], [(start-i)*yscale, (start-i)*yscale], ':', color='darkgray')


    pylab.figure(figsize=(12,12))

    pylab.axvspan(1.-sig8_sig, 1.+sig8_sig,alpha=0.2, color='gray')
    plot_set(lines, 'b', len(lines))
    pylab.ylim(-2, 52)
    pylab.xlim(0.95, 1.1)


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
    pylab.savefig('plots/pz_sig8_'+test+'.png')
    pylab.close()


    # print 'tomo'

    # np.set_printoptions(precision=2)
    # for i in xrange(4):
    #   for j in xrange(4):
    #     if j==0:
    #       print '\\'+SVA1PZPaper.get_pz_name(i+1,test).lower()+' & ',
    #     if i==j:
    #       print '-- & ',
    #     else:
    #       print str(np.around(np.log(sum0a[i]/sum0a[j]),decimals=2))+' ('+str(np.around(np.log(sum0b[i]/sum0b[j]),decimals=2))+') & ',

    #   print ''

    #   print ''

    # np.set_printoptions(edgeitems=3,infstr='inf',linewidth=75, nanstr='nan', precision=8,suppress=False, threshold=1000, formatter=None)

    return 

  @staticmethod
  def plot_nofz(pz0,test):

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
  def imshow_symlog(my_matrix, logthresh=3):
    # From Freddie Nfbnm on stack overflow - http://stackoverflow.com/questions/11138706/

    plt.imshow( my_matrix,interpolation='nearest',origin='lower',vmin=np.min(my_matrix), vmax=np.max(my_matrix),norm=matplotlib.colors.SymLogNorm(10**-logthresh) )

    maxlog=int(np.ceil( np.log10(np.max(my_matrix)) ))
    minlog=int(np.ceil( np.log10(-np.min(my_matrix)) ))

    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]+[0.0]+[(10**x) for x in xrange(-logthresh,maxlog+1)] )

    plt.colorbar(ticks=tick_locations)

    return 
