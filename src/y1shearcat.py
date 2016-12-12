import os
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
dirname=os.path.split(__file__)[0]
style_file=os.path.join(dirname, "SVA1StyleSheet.mplstyle")
plt.style.use(style_file)
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import src.sys_split as sys_split
import src.config as config
import src.catalog as catalog

class y1(object):

    @staticmethod
    def load_data(i3file,mcalfile,goldfile):

        mcal = catalog.CatalogStore('matched_metacal',cutfunc=catalog.CatalogMethods.matched_metacal_cut(),cattype='ng',catfile=mcalfile,goldfile=goldfile)

        mask = mcal.size/mcal.psffwhm>0.5
        catalog.CatalogMethods.match_cat(mcal,mask)

        mcal.bs    = True
        mcal.wt    = True
        mcal.lbins = 20

        #np.save('mcal_coadds.npy',np.vstack((mcal.coadd,np.ones(len(mcal.coadd)),np.ones(len(mcal.coadd)),np.zeros(len(mcal.coadd)),np.zeros(len(mcal.coadd)),mcal.w)).T)

        i3 = catalog.CatalogStore('matched_i3',cutfunc=catalog.CatalogMethods.matched_i3_cut(),cattype='i3',catfile=i3file,goldfile=goldfile)

        i3.bs    = True
        i3.wt    = True
        i3.lbins = 20
        i3.w     = 1./(2*0.22**2+i3.cov11+i3.cov22)

        #np.save('i3_coadds.npy',np.vstack((i3.coadd,i3.m1,i3.m2,i3.c1,i3.c2,i3.w)).T)

        return i3,mcal

class y1_plots(object):

    @staticmethod
    def mean_e(cat1,cat2):

        txt.write_methods.heading('Linear Splits',cat,label='y1_paper',create=False)

        y1_plots.evssnr(cat1,cat2)
        y1_plots.evsradius(cat1,cat2)
        y1_plots.evspsf1(cat1,cat2)
        y1_plots.evspsf2(cat1,cat2)
        y1_plots.evspsfsize(cat1,cat2)

        return

    @staticmethod
    def mean_e_subplot(cat,n,val,fig):

        array=getattr(cat,val)
        tmp,tmp,arr1,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=split_gals_lin_along_base([cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name,log=config.log_val.get(val,False),plot=False)

        plt.figure(fig)
        ax=plt.subplot(2,1,n)
        plt.errorbar(arr1,e1,yerr=e1err,marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
        plt.errorbar(arr1,m1*arr1+b1,marker='',linestyle='-',color='r')
        plt.errorbar(arr1+(x[1]-x[0])/5.,e2,yerr=e2err,marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
        plt.errorbar(arr1,m2*arr1+b2,marker='',linestyle='-',color='b')
        ax.minorticks_on()
        plt.ylabel(r'$\langle e \rangle$')
        # plt.axhline(.004,color='k')
        # plt.axhline(-.004,color='k')
        plt.xlabel(config.lbl.get(val,val.replace('_','-')))
        plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})

        return

    @staticmethod
    def evssnr(cat1,cat2):

        plt.figure(1)

        y1_plots.mean_e_subplot(cat1,0,'snr',1)
        y1_plots.mean_e_subplot(cat2,1,'snr',1)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_snr.pdf', bbox_inches='tight')
        plt.close(1)

        return

    @staticmethod
    def evsradius(cat1,cat2):

        plt.figure(2)

        y1_plots.mean_e_subplot(cat1,0,'rgp',2)
        y1_plots.mean_e_subplot(cat2,1,'size',2)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_radius.pdf', bbox_inches='tight')
        plt.close(2)

        return

    @staticmethod
    def evspsf1(cat1,cat2):

        plt.figure(3)

        y1_plots.mean_e_subplot(cat1,0,'psf1',3)
        y1_plots.mean_e_subplot(cat2,1,'psf1',3)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_psf1.pdf', bbox_inches='tight')
        plt.close(3)

        return

    @staticmethod
    def evspsf2(cat1,cat2):

        plt.figure(4)

        y1_plots.mean_e_subplot(cat1,0,'psf2',4)
        y1_plots.mean_e_subplot(cat2,1,'psf2',4)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_psf2.pdf', bbox_inches='tight')
        plt.close(4)

        return

    @staticmethod
    def evspsfsize(cat1,cat2):

        plt.figure(5)

        y1_plots.mean_e_subplot(cat1,0,'psffwhm',5)
        y1_plots.mean_e_subplot(cat2,1,'psffwhm',5)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_psfsize.pdf', bbox_inches='tight')
        plt.close(5)

        return

    @staticmethod
    def whisker(psfdir,cat1,cat2):

        cols = ['ra','dec','e1','e2','ccd','col','row','psf_e1','psf_e2','mag']
        psf  = catalog.CatalogStore('psf',cutfunc=catalog.CatalogMethods.final_null_cuts_ra(),cols=cols,cattype='psf',catdir=psfdir)

        plt.figure(6)

        y1_plots.whiskerplot(cat1,'e',6)
        y1_plots.whiskerplot(cat2,'psf',6)
        y1_plots.whiskerplot(cat2,'dpsf',6)

        plt.subplots_adjust(hspace=0,wspace=0)
        plt.savefig('plots/y1/lin_split_psfsize.pdf', bbox_inches='tight')
        plt.close(6)

        return

#   @staticmethod
#   def plot_lin_split(x,e1,e2,e1err,e2err,m1,m2,b1,b2,name,val,log=False,label='',e=True,val2=None,trend=True):

#     plt.figure()
#     if e:
#       l1=r'$\langle e_1 \rangle$'
#       l2=r'$\langle e_2 \rangle$'
#       plt.errorbar(x,e1,yerr=e1err,marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
#     else:
#       l1=r'$\langle e_1 \rangle$'      
#       plt.errorbar(x,e1,yerr=e1err,marker='o',linestyle='',color='r',label='')
#     if trend:
#       plt.errorbar(x,m1*x+b1,marker='',linestyle='-',color='r')
#     if e:
#       plt.errorbar(x+(x[1]-x[0])/5.,e2,yerr=e2err,marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
#       if trend:
#         plt.errorbar(x,m2*x+b2,marker='',linestyle='-',color='b')
#       plt.ylabel(r'$\langle e \rangle$')
#     else:
#       plt.ylabel(r'$\langle $'+config.lbl.get(val2,val2)+r'$ \rangle$')
#     if e:
#       plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})
#       plt.axhline(.004,color='k')
#       plt.axhline(-.004,color='k')
#     if config.log_val.get(val,False):
#       plt.xlabel('log '+config.lbl.get(val,val.replace('_','-')))
#     else:
#       plt.xlabel(config.lbl.get(val,val.replace('_','-')))
#     y1=np.min(np.minimum(e1,e2))
#     if e:   
#       y2=np.max(np.maximum(e1,e2))
#     else:
#       y2=y1
#     # plt.ylim((np.min([y1-(y2-y1)/10.,-.005]),np.max([y2+(y2-y1)/10.,.005])))
#     plt.minorticks_on()
#     if val2 is not None:
#       val+='-'+val2
#     plt.savefig('plots/split/lin_split_'+name+'_'+val+'_'+label.replace('_','-')+'.png', bbox_inches='tight')
#     plt.close()

#     return




# class SVA1(object):

#   @staticmethod
#   def Fig_2ptPaper_amp(cat,xip,dxip,dxiperr):

#     chi2st=999999
#     amp=0.
#     for a in xrange(-200,200):
#       chi2=0.
#       for i in xrange(cat.tbins):
#         chi2=np.sum((dxip-a*xip/100.)**2./dxiperr**2.)
#       if chi2<chi2st:
#         chi2st=chi2
#         amp=a/100.

#     return amp

#   @staticmethod
#   def Fig_2ptPaper_subplot(cat,fig,r,c,n,array,label):

#     print ' '
#     print ' '
#     print '..........'+label+'..........'
#     print ' '
#     print ' '
#     theta,edge,edgemean,xip,xiperr,xim,ximerr,dxip0,dxip1,dxip2,xiperr0,xiperr1,xiperr2,chi2p0,chi2p1,chi2p2,dxim0,dxim1,dxim2,ximerr0,ximerr1,ximerr2,chi2m0,chi2m1,chi2m2=split_systematics.split_gals_2pt_along_(array,cat,label)
#     #edge=edge.astype(int)
#     A1b=SVA1.Fig_2ptPaper_amp(cat,xip,dxip1,xiperr1)
#     A2b=SVA1.Fig_2ptPaper_amp(cat,xip,dxip2,xiperr2)
#     A1bm=SVA1.Fig_2ptPaper_amp(cat,xim,dxim1,ximerr1)
#     A2bm=SVA1.Fig_2ptPaper_amp(cat,xim,dxim2,ximerr2)
    
#     cat.edgest=np.vstack((cat.edgest,edge))
#     cat.meanst=np.vstack((cat.meanst,edgemean))
#     cat.astp=np.vstack((cat.astp,np.array([A1a,A2a,A1b,A2b])))
#     cat.astm=np.vstack((cat.astm,np.array([A1am,A2am,A1bm,A2bm])))
#     cat.chi2stp=np.vstack((cat.chi2stp,np.array([chi2p0,chi2p1,chi2p2])))
#     cat.chi2stm=np.vstack((cat.chi2stm,np.array([chi2m0,chi2m1,chi2m2])))

#     print 'theta',theta
#     print 'xip',xip
#     print 'xiperr',xiperr
#     print 'dxip1',dxip1
#     print 'dxiperr1',xiperr1
#     print 'dxip2',dxip2
#     print 'dxiperr2',xiperr2
#     print 'chi',chi2p1,chi2p2,chi2p0
#     print 'Aa',A1a,A2a
#     print 'Ab',A1b,A2b
#     plt.figure(fig)
#     ax=plt.subplot(r,c,n)
#     ax.fill_between(theta,-xiperr/xip,xiperr/xip,facecolor='gray',alpha=0.4)
#     plt.errorbar(theta,np.zeros((len(theta))),marker='',linestyle='-',color='k')
#     plt.errorbar(theta,A1b*np.ones((len(theta))),marker='',linestyle='-',color='r')
#     plt.errorbar(theta*(1),dxip1/xip,yerr=xiperr1/xip,marker='v',linestyle='',color='r')
#     plt.errorbar(theta,A2b*np.ones((len(theta))),marker='',linestyle='-',color='b')
#     plt.errorbar(theta*(1.2),dxip2/xip,yerr=xiperr2/xip,marker='^',linestyle='',color='b')
#     t=ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=11,verticalalignment='top')
#     t.set_bbox(dict(color='white', alpha=0.5, edgecolor='white'))
#     #plt.yscale('log')
#     plt.xscale('log')
#     plt.xlim((1,500))
#     plt.ylim((-1.49,1.49))
#     ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#     #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
#     if n<5:
#       ax.set_xticklabels([])
#     else:
#       plt.xlabel(r'$\theta$ (arcmin)')
#     if n%2==0:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\Delta\xi_+/\xi_+$')
#     if n==1:
#       ax.set_title('im3shape7.2')
#     elif n==2:
#       ax.set_title('ngmix010')

#     plt.figure(fig+1)
#     ax=plt.subplot(r,c,n)
#     ax.fill_between(theta,-np.abs(ximerr/xim),np.abs(ximerr/xim),facecolor='gray',alpha=0.4)
#     plt.errorbar(theta,np.zeros((len(theta))),marker='',linestyle='-',color='k')
#     plt.errorbar(theta,A1bm*np.ones((len(theta))),marker='',linestyle='-',color='r')
#     plt.errorbar(theta*(1),dxim1/xim,yerr=ximerr1/xim,marker='v',linestyle='',color='r')
#     plt.errorbar(theta,A2bm*np.ones((len(theta))),marker='',linestyle='-',color='b')
#     plt.errorbar(theta*(1.2),dxim2/xim,yerr=ximerr2/xim,marker='^',linestyle='',color='b')
#     t=ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=11,verticalalignment='top')
#     t.set_bbox(dict(color='white', alpha=0.5, edgecolor='white'))    
#     #plt.yscale('log')
#     plt.xscale('log')
#     plt.xlim((1,500))
#     plt.ylim((-1.49,1.49))
#     ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#     #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
#     if n<5:
#       ax.set_xticklabels([])
#     else:
#       plt.xlabel(r'$\theta$ (arcmin)')
#     if n%2==0:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\Delta\xi_{-}/\xi_{-}$')
#     if n==1:
#       ax.set_title('im3shape7.2')
#     elif n==2:
#       ax.set_title('ngmix010')
#     return


#   @staticmethod
#   def Fig_2ptPaper_subplotold(cat,fig,r,c,n,array,label):

#     plt.figure(fig)
#     print ' '
#     print ' '
#     print '..........'+label+'..........'
#     print ' '
#     print ' '
#     theta,edge,xip,xiperr,dxip1,dxip2,xiperr1,xiperr2,chi2=split_systematics.split_gals_2pt_along_(array,cat,label)
#     edge=edge.astype(int)
#     A1=SVA1.Fig_2ptPaper_amp(cat,xip,dxip1,xiperr1)
#     A2=SVA1.Fig_2ptPaper_amp(cat,xip,dxip2,xiperr2)
#     print 'theta',theta
#     print 'xip',xip
#     print 'xiperr',xiperr
#     print 'dxip1',dxip1
#     print 'dxiperr1',xiperr1
#     print 'dxip2',dxip2
#     print 'dxiperr2',xiperr2
#     print 'chi',chi2
#     print 'A',A1,A2
#     ax=plt.subplot(r,c,n)
#     plt.errorbar(theta[xip>0],xip[xip>0]*theta[xip>0],yerr=xiperr[xip>0]*theta[xip>0],marker='.',linestyle='',color='k')
#     plt.errorbar(theta[xip<0],xip[xip<0]*theta[xip<0],yerr=xiperr[xip<0]*theta[xip<0],marker='.',linestyle='',color='k')
#     plt.errorbar(theta,np.abs(xip*A1*theta),marker='',linestyle=':',color='r')
#     plt.errorbar(theta[dxip1>0]*(1.2),dxip1[dxip1>0]*theta[dxip1>0],yerr=xiperr1[dxip1>0]*theta[dxip1>0],marker='v',linestyle='',color='r',label='snr=['+str(edge[0])+','+str(edge[1])+']')
#     plt.errorbar(theta[dxip1<0]*(1.2),-dxip1[dxip1<0]*theta[dxip1<0],yerr=xiperr1[dxip1<0]*theta[dxip1<0],marker='1',linestyle='',color='r')
#     plt.errorbar(theta,np.abs(xip*A2*theta),marker='',linestyle=':',color='b')
#     plt.errorbar(theta[dxip2>0]*(1.4),dxip2[dxip2>0]*theta[dxip2>0],yerr=xiperr2[dxip2>0]*theta[dxip2>0],marker='^',linestyle='',color='b',label='snr=['+str(edge[1])+','+str(edge[2])+']')
#     plt.errorbar(theta[dxip2<0]*(1.4),-dxip2[dxip2<0]*theta[dxip2<0],yerr=xiperr2[dxip2<0]*theta[dxip2<0],marker='2',linestyle='',color='b')
#     #plt.yscale('log')
#     plt.xscale('log')
#     plt.xlim((1,500))
#     plt.ylim((1e-6,9e-4))
#     #plt.legend(loc='upper center',ncol=2, frameon=False,prop={'size':10})
#     if n<5:
#       ax.set_xticklabels([])
#     else:
#       plt.xlabel(r'$\theta$ (arcmin)')
#     if n%2==0:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\theta\xi_+$')

#     return


#   @staticmethod
#   def Fig_ShearPaper_subplot(cat,fig,r,c,n,array,label):

#     plt.figure(fig)
#     print ' '
#     print ' '
#     print '..........'+label+'..........'
#     print ' '
#     print ' '
#     # lbinst=cat.lbins
#     # cat.lbins=50
#     # arr1,me1,me2,e1err,e2err,slp1,slp2,b1,b2=split_systematics.split_gals_lin_along_(array,cat,label)
#     # cat.lbins=1
#     # tmp,mm1,mm2,tmp,tmp,tmp,tmp,tmp,tmp=split_systematics.split_gals_lin_along_(array,cat,label)
#     # cat.lbins=lbinst
#     arr1,me1,me2,e1err,e2err,tmp,tmp,tmp,tmp=split_systematics.split_gals_lin_along_(array,cat,label)
#     print 'psf',arr1
#     print 'me',me1,me2
#     print 'err',e1err,e2err
#     print 'slp',slp1,slp2
#     print 'b',b1,b2
#     ax=plt.subplot(r,c,n)
#     if n==3:
#       plt.errorbar(arr1,me1,yerr=e1err,marker='.',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
#       plt.errorbar(arr1,me2,yerr=e2err,marker='.',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
#     else:
#       plt.errorbar(arr1,me1,yerr=e1err,marker='.',linestyle='',color='r')
#       plt.errorbar(arr1,me2,yerr=e2err,marker='.',linestyle='',color='b')
#     plt.errorbar(arr1,slp1*arr1+b1,marker='',linestyle='-',color='r')
#     plt.errorbar(arr1,slp2*arr1+b2,marker='',linestyle='-',color='b')
#     plt.errorbar(arr1,mm1*np.ones((len(arr1))),marker='',linestyle=':',color='r')
#     plt.errorbar(arr1,mm2*np.ones((len(arr1))),marker='',linestyle=':',color='b')
#     plt.legend(loc='upper right',ncol=2, frameon=False,prop={'size':12})
#     ax.minorticks_on()
#     if n<4:
#       ax.set_xticklabels([])
#       plt.ylim((-.00175,.00175))
#     else:
#       plt.ylim((-.0005,.00175))
#       if label=='psfe1':
#         plt.xlabel(r'PSF $e_1$')
#         plt.xlim((-.019,.0325))
#       if label=='psfe2':
#         plt.xlabel(r'PSF $e_2$')
#         plt.xlim((-.019,.029))
#       if label=='psffwhm':
#         plt.xlabel(r'PSF FWHM')
#         plt.xlim((.51,.7))
#     if (n==2)|(n==5):
#       ax.set_yticklabels([])
#       #plt.xlim((-.0075,.0225))
#     elif (n==3)|(n==6):
#       ax.set_yticklabels([])
#       #plt.xlim((-.0075,.0225))
#     else:
#       plt.ylabel(r'$\langle e \rangle$')
#       #plt.xlim((-.0025,.0175))

#     return


#   @staticmethod
#   def Fig_ShearPaper_subplot2(cat,ng,i3,fig,n,label,mask):

#     #if cat.cat==ng.cat:
#     #  mask1=np.in1d(ng.coadd,i3.coadd,assume_unique=True)
#     #  mask2=np.in1d(i3.coadd,ng.coadd,assume_unique=True)
#     #  tmp=i3.w[mask2]
#     #  sort=np.argsort(tmp)[np.argsort(np.argsort(ng.coadd[mask1]))]
#     #  mask=mask1
#     #  cat.w[mask1]=tmp[sort]
#     #else:
#     #  mask=np.ones((len(cat.coadd))).astype(bool)

#     print ' '
#     print ' '
#     print '..........psf-shear '+label+'..........'
#     print ' '
#     print ' '
#     theta,ggp,ggm,ggperr,ggmerr,gpp,gpm,gpperr,gpmerr,ppp,ppm,ppperr,ppmerr,alphap,alphaperr,leakagep,leakageperr,alpham,alphamerr,leakagem,leakagemerr=xi_2pt_shear_test_methods.xi_2pt_psf(cat,mask)
#     print 'ggp',ggp
#     print 'ggperr',ggperr
#     print 'gpp',gpp
#     print 'gpperr',gpperr
#     print 'ppp',ppp
#     print 'ppperr',ppperr
#     print 'alpha',alphap,alpham
#     print 'leakage',leakagep,leakagem

#     plt.figure(fig)

#     ax=plt.subplot(3,2,1+n)
#     if n==0:
#       ax.set_title(CatalogMethods.cat_name(cat.cat))
#     else:
#       ax.set_title(CatalogMethods.cat_name(cat.cat))
#     plt.errorbar(theta,ggp,yerr=ggperr,marker='.',linestyle='',color='k',label='gg')
#     plt.errorbar(theta,np.abs(gpp),yerr=gpperr,marker='.',linestyle='',color='r',label='gp')
#     plt.errorbar(theta,ppp,yerr=ppperr,marker='.',linestyle='',color='b',label='pp')
#     plt.legend(loc='upper right',ncol=3, frameon=False,prop={'size':12})
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlim((1,500))
#     plt.ylim(5e-7,6e-4)
#     ax.set_xticklabels([])
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\xi_+$')

#     ax=plt.subplot(3,2,3+n)
#     plt.errorbar(theta,alphap,marker='',linestyle=':',color='k')
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.5,.5)
#     ax.set_xticklabels([])
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\alpha$')

#     ax=plt.subplot(3,2,5+n)
#     plt.errorbar(theta,leakagep,yerr=leakageperr,marker='',linestyle=':',color='k')
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.5,.5)
#     plt.xlabel(r'$\theta$ (arcmin)')
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'Leakage into $\xi_+$')


#     cat.tmpsize=np.average(cat.radius[mask])
#     alphap=xi_2pt_shear_test_methods.calc_alpha(gpp,ppp,cat)
#     leakagep=xi_2pt_shear_test_methods.calc_psf_leakage(alphap,ppp,cat)/ggp
#     print alphap
#     print leakagep
#     print 'mean size',np.mean(cat.radius)
#     cat.tmpsize=1.


#     ax=plt.subplot(3,2,3+n)
#     plt.errorbar(theta,alphap,marker='',linestyle='-',color='k')
#     #ax.fill_between(theta,alphap-alphaperr,alphap+alphaperr,facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.5,.5)
#     ax.set_xticklabels([])
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\alpha$')

#     ax=plt.subplot(3,2,5+n)
#     plt.errorbar(theta,leakagep,yerr=leakageperr,marker='',linestyle='-',color='k')
#     #ax.fill_between(theta,leakagep-leakageperr,leakagep+leakageperr,facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.5,.5)
#     plt.xlabel(r'$\theta$ (arcmin)')
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'Leakage into $\xi_+$')

#     plt.figure(fig+1)

#     ax=plt.subplot(3,2,1+n)
#     plt.errorbar(theta,ggm,yerr=ggmerr,marker='.',linestyle='',color='k',label='gg')
#     plt.errorbar(theta,gpm,yerr=gpmerr,marker='.',linestyle='',color='r',label='gp')
#     plt.errorbar(theta,ppm,yerr=ppmerr,marker='.',linestyle='',color='b',label='pp')
#     plt.legend(loc='upper right',ncol=3, frameon=False,prop={'size':12})
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlim((1,500))
#     plt.ylim(5e-7,6e-4)
#     ax.set_xticklabels([])
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\xi_-$')

#     ax=plt.subplot(3,2,3+n)
#     plt.errorbar(theta,alpham,marker='',linestyle='-',color='k')
#     ax.fill_between(theta,alpham-alphamerr,alpham+alphamerr,facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.1,.49)
#     ax.set_xticklabels([])
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'$\alpha$')

#     ax=plt.subplot(3,2,5+n)
#     plt.errorbar(theta,leakagem,yerr=leakagemerr,marker='',linestyle='-',color='k')
#     ax.fill_between(theta,leakagem-leakagemerr,leakagem+leakagemerr,facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.xlim((1,500))
#     plt.ylim(-.1,.49)
#     plt.xlabel(r'$\theta$ (arcmin)')
#     if n==1:
#       ax.set_yticklabels([])
#     else:
#       plt.ylabel(r'Leakage into $\xi_-$')

#     return


#   @staticmethod
#   def Fig1_2ptPaper(i3,ng):


#     print '...........'
#     print 'Figure 1'
#     print '...........'

#     SVA1.Fig_2ptPaper_subplot(i3,1,3,2,1,i3.snr,'Signal-to-Noise')
#     SVA1.Fig_2ptPaper_subplot(ng,1,3,2,2,ng.snr,'Signal-to-Noise')
#     SVA1.Fig_2ptPaper_subplot(i3,1,3,2,3,i3.radius,'Size')
#     SVA1.Fig_2ptPaper_subplot(ng,1,3,2,4,ng.radius,'Size')
#     SVA1.Fig_2ptPaper_subplot(i3,1,3,2,5,i3.colour,'Colour')
#     SVA1.Fig_2ptPaper_subplot(ng,1,3,2,6,ng.colour,'Colour')
#     plt.figure(1)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig1p.png', bbox_inches='tight')
#     plt.close(1)   
#     plt.figure(2)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig1m.png', bbox_inches='tight')
#     plt.close(2)  

#     print 'edge final',i3.edgest,ng.edgest
#     print 'mean final',i3.meanst,ng.meanst
#     print 'Astp final',i3.astp,ng.astp
#     print 'Astm final',i3.astm,ng.astm
#     print 'chi2p final',i3.chi2stp,ng.chi2stp
#     print 'chi2m final',i3.chi2stm,ng.chi2stm

#     tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
#     tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
#     tmpi3=np.around(tmpi3, decimals=2)
#     tmpng=np.around(tmpng, decimals=2)
#     np.set_printoptions(suppress=True)
#     print tmpi3
#     print tmpng    


#     return

#   @staticmethod
#   def Fig2_2ptPaper(i3,ng):


#     print '...........'
#     print 'Figure 2'
#     print '...........'

    
#     SVA1.Fig_2ptPaper_subplot(i3,2,3,2,1,i3.ra,'RA')
#     SVA1.Fig_2ptPaper_subplot(ng,2,3,2,2,ng.ra,'RA')
#     SVA1.Fig_2ptPaper_subplot(i3,2,3,2,3,i3.dec,'Dec')
#     SVA1.Fig_2ptPaper_subplot(ng,2,3,2,4,ng.dec,'Dec')
#     SVA1.Fig_2ptPaper_subplot(i3,2,3,2,5,i3.ebv,'E(B-V)')
#     SVA1.Fig_2ptPaper_subplot(ng,2,3,2,6,ng.ebv,'E(B-V)')
#     plt.figure(2)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig2p.png', bbox_inches='tight')
#     plt.close(2)   
#     plt.figure(3)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig2m.png', bbox_inches='tight')
#     plt.close(3)   

#     print 'edge final',i3.edgest,ng.edgest
#     print 'mean final',i3.meanst,ng.meanst
#     print 'Astp final',i3.astp,ng.astp
#     print 'Astm final',i3.astm,ng.astm
#     print 'chi2p final',i3.chi2stp,ng.chi2stp
#     print 'chi2m final',i3.chi2stm,ng.chi2stm

#     tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
#     tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
#     tmpi3=np.around(tmpi3, decimals=2)
#     tmpng=np.around(tmpng, decimals=2)
#     np.set_printoptions(suppress=True)
#     print tmpi3
#     print tmpng    


#     return

#   @staticmethod
#   def Fig3_2ptPaper(i3,ng):


#     print '...........'
#     print 'Figure 3'
#     print '...........'

    
#     SVA1.Fig_2ptPaper_subplot(i3,3,3,2,1,i3.exptime,'Exposure Time')
#     SVA1.Fig_2ptPaper_subplot(ng,3,3,2,2,ng.exptime,'Exposure Time')
#     SVA1.Fig_2ptPaper_subplot(i3,3,3,2,3,i3.maglimit,'Mag Limit')
#     SVA1.Fig_2ptPaper_subplot(ng,3,3,2,4,ng.maglimit,'Mag Limit')
#     SVA1.Fig_2ptPaper_subplot(i3,3,3,2,5,i3.airmass,'Air Mass')
#     SVA1.Fig_2ptPaper_subplot(ng,3,3,2,6,ng.airmass,'Air Mass')
#     plt.figure(3)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig3p.png', bbox_inches='tight')
#     plt.close(3)   
#     plt.figure(4)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig3m.png', bbox_inches='tight')
#     plt.close(4)  

#     print 'edge final',i3.edgest,ng.edgest
#     print 'mean final',i3.meanst,ng.meanst
#     print 'Astp final',i3.astp,ng.astp
#     print 'Astm final',i3.astm,ng.astm
#     print 'chi2p final',i3.chi2stp,ng.chi2stp
#     print 'chi2m final',i3.chi2stm,ng.chi2stm

#     tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
#     tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
#     tmpi3=np.around(tmpi3, decimals=2)
#     tmpng=np.around(tmpng, decimals=2)
#     np.set_printoptions(suppress=True)
#     print tmpi3
#     print tmpng    

#     return


#   @staticmethod
#   def Fig4_2ptPaper(i3,ng):


#     print '...........'
#     print 'Figure 4'
#     print '...........'

  
#     SVA1.Fig_2ptPaper_subplot(i3,4,3,2,1,i3.skysigma,'Sky Sigma')
#     SVA1.Fig_2ptPaper_subplot(ng,4,3,2,2,ng.skysigma,'Sky Sigma')
#     SVA1.Fig_2ptPaper_subplot(i3,4,3,2,3,i3.skybrite,'Sky Brightness')
#     SVA1.Fig_2ptPaper_subplot(ng,4,3,2,4,ng.skybrite,'Sky Brightness')
#     SVA1.Fig_2ptPaper_subplot(i3,4,3,2,5,i3.fwhm,'FWHM')
#     SVA1.Fig_2ptPaper_subplot(ng,4,3,2,6,ng.fwhm,'FWHM')
#     plt.figure(4)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig4p.png', bbox_inches='tight')
#     plt.close(4)   
#     plt.figure(5)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig4m.png', bbox_inches='tight')
#     plt.close(5)   
        
#     print 'edge final',i3.edgest,ng.edgest
#     print 'mean final',i3.meanst,ng.meanst
#     print 'Astp final',i3.astp,ng.astp
#     print 'Astm final',i3.astm,ng.astm
#     print 'chi2p final',i3.chi2stp,ng.chi2stp
#     print 'chi2m final',i3.chi2stm,ng.chi2stm

#     tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
#     tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
#     tmpi3=np.around(tmpi3, decimals=2)
#     tmpng=np.around(tmpng, decimals=2)
#     np.set_printoptions(suppress=True)
#     print tmpi3
#     print tmpng    


#     return

#   @staticmethod
#   def Fig5_2ptPaper(i3,ng):


#     print '...........'
#     print 'Figure 5'
#     print '...........'

    
#     SVA1.Fig_2ptPaper_subplot(i3,5,3,2,1,i3.psf1,r'PSF $e_1$')
#     SVA1.Fig_2ptPaper_subplot(ng,5,3,2,2,ng.psf1,r'PSF $e_1$')
#     SVA1.Fig_2ptPaper_subplot(i3,5,3,2,3,i3.psf2,r'PSF $e_2$')
#     SVA1.Fig_2ptPaper_subplot(ng,5,3,2,4,ng.psf2,r'PSF $e_2$')
#     SVA1.Fig_2ptPaper_subplot(i3,5,3,2,5,i3.psffwhm,'PSF FWHM')
#     SVA1.Fig_2ptPaper_subplot(ng,5,3,2,6,ng.psffwhm,'PSF FWHM')
#     plt.figure(5)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig5p.png', bbox_inches='tight')
#     plt.close(5)   
#     plt.figure(6)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('2ptPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_pzreweight-'+str(i3.use_zrw)+'_ztype-'+i3.ztyp+'_Fig5m.png', bbox_inches='tight')
#     plt.close(6)   
        
#     print 'edge final',i3.edgest,ng.edgest
#     print 'mean final',i3.meanst,ng.meanst
#     print 'Astp final',i3.astp,ng.astp
#     print 'Astm final',i3.astm,ng.astm
#     print 'chi2p final',i3.chi2stp,ng.chi2stp
#     print 'chi2m final',i3.chi2stm,ng.chi2stm

#     tmpi3=np.hstack((i3.edgest,i3.meanst,i3.chi2stp))
#     tmpng=np.hstack((ng.edgest,ng.meanst,ng.chi2stp))
#     tmpi3=np.around(tmpi3, decimals=2)
#     tmpng=np.around(tmpng, decimals=2)
#     np.set_printoptions(suppress=True)
#     print tmpi3
#     print tmpng    


#     return

#   @staticmethod
#   def Fig1_ShearPaper(i3,ng):


#     print '...........'
#     print 'Figure 1'
#     print '...........'

#     SVA1.Fig_ShearPaper_subplot(i3,1,2,3,1,i3.psf1,'psfe1')    
#     SVA1.Fig_ShearPaper_subplot(i3,1,2,3,2,i3.psf2,'psfe2')    
#     SVA1.Fig_ShearPaper_subplot(i3,1,2,3,3,i3.psffwhm,'psffwhm')    
#     SVA1.Fig_ShearPaper_subplot(ng,1,2,3,4,ng.psf1,'psfe1')    
#     SVA1.Fig_ShearPaper_subplot(ng,1,2,3,5,ng.psf2,'psfe2')    
#     SVA1.Fig_ShearPaper_subplot(ng,1,2,3,6,ng.psffwhm,'psffwhm')    
#     plt.figure(1)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.lbins)+'_Fig1.png', bbox_inches='tight')
#     plt.close(1)
    
#     return

#   @staticmethod
#   def Fig2_ShearPaper(i3,ng,mask,mask2):


#     print '...........'
#     print 'Figure 2'
#     print '...........'

#     SVA1.Fig_ShearPaper_subplot2(i3,ng,i3,1,0,'im3shape',mask)
#     SVA1.Fig_ShearPaper_subplot2(ng,ng,i3,1,1,'ngmix',mask2)
#     plt.figure(1)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_bs-'+str(i3.bs)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_wt-'+str(i3.wt)+'_Fig2_p.png', bbox_inches='tight')
#     plt.close(1)
#     plt.figure(2)
#     plt.subplots_adjust(hspace=0,wspace=0)
#     plt.savefig('ShearPaper_jk-'+str(i3.use_jk)+'_nbins-'+str(i3.tbins)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_Fig2_m.png', bbox_inches='tight')
#     plt.close(2)
    
#     return

#   @staticmethod
#   def xi2cosebi(cat,theta,ggp,ggm):

#     if cat.sep[0] not in [1]:
#       print 'must use theta min == 1'
#       return np.zeros((cat.cbins)),np.zeros((cat.cbins)),np.zeros((cat.cbins))
#     if cat.sep[1] not in [30,60,100,200,400]:
#       print 'must use theta max in [30,60,100,200,400]'
#       return np.zeros((cat.cbins)),np.zeros((cat.cbins)),np.zeros((cat.cbins))
#     if cat.tbins not in [20,300]:
#       print 'must use theta bins in [20,300]'
#       return np.zeros((cat.cbins)),np.zeros((cat.cbins)),np.zeros((cat.cbins))

#     np.savetxt('cosebis/COSEBIstmpXi.out',np.column_stack((theta,ggp,ggm)))

#     cosebicall='./cr/EB_data '+str(cat.sep[1])+' '+str(cat.sep[0])+' '+str(cat.tbins)+' ./cosebis/COSEBIstmpXi.out ./cosebis/ '+str(cat.cbins)+' 0'
#     os.system(cosebicall)
#     #loadtxt(os.popen(cosebicall).read())
#     cosebi=np.loadtxt('COSEBIstmpXi.in')

#     return cosebi[:,0],cosebi[:,1],cosebi[:,2]

#   @staticmethod
#   def test_2pt(cat):

#     print '...........'
#     print 'test'
#     print '...........'

#     mask=np.ones((len(cat.coadd)))

#     theta,ggp,ggm,ggperr,ggmerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(cat,mask.astype(bool),mask)
#     plotting_methods.fig_open_xi_test_p(theta,ggp,ggperr,cat,'xi+')
#     plotting_methods.fig_open_xi_test_m(theta,ggm,ggmerr,cat,'xi-')
#     plotting_methods.fig_close_xi_test(cat,'test')
#     np.savetxt('2pcf/'+CatalogMethods.cat_name(cat.cat)+'_jk-'+str(cat.use_jk)+'_nbins-'+str(cat.tbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_t0-'+str(cat.sep[0])+'_t1-'+str(cat.sep[1])+'_pzreweight-'+str(cat.use_zrw)+'_ztype-'+cat.ztyp+'_xi.txt',np.column_stack((theta,ggp,ggperr,ggm,ggmerr)))

#     cbin,ce,cb=SVA1.xi2cosebi(cat,theta,ggp,ggm)

#     np.savetxt('cosebis/'+CatalogMethods.cat_name(cat.cat)+'_jk-'+str(cat.use_jk)+'_nbins-'+str(cat.tbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_t0-'+str(cat.sep[0])+'_t1-'+str(cat.sep[1])+'_pzreweight-'+str(cat.use_zrw)+'_ztype-'+cat.ztyp+'_cosebi.txt',np.column_stack((cbin,ce,cb)))

#     return

#   @staticmethod
#   def test_2pt_full(cat):


#     print '...........'
#     print 'test full'
#     print '...........'

#     #cat.slop=1.
    
#     cat.tbins=20
#     cat.sep=np.array([1,30])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,60])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,100])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,200])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,400])
#     SVA1.test_2pt(cat)

#     cat.tbins=300
#     cat.sep=np.array([1,30])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,60])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,100])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,200])
#     SVA1.test_2pt(cat)
#     cat.sep=np.array([1,400])
#     SVA1.test_2pt(cat)

#     return

#   @staticmethod
#   def ratiotest_2pt(i3,ng,label):


#     print '...........'
#     print 'ratiotest'
#     print '...........'

#     mask=np.ones((len(i3.coadd)))

#     theta,i3xip,i3xiperr,i3xim,i3ximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(i3,mask.astype(bool),mask)
#     theta,ngxip,ngxiperr,ngxim,ngximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ng,mask.astype(bool),mask)

#     slpp,bp=np.polyfit(np.log(theta), ngxip/i3xip, 1)
#     slpm,bm=np.polyfit(np.log(theta), ngxim/i3xim, 1)

#     print ngxim, i3xim,ngxim/i3xim

#     plt.figure(10)
#     plt.errorbar(theta,ngxip/i3xip, yerr=None, xerr=None,marker='o',linestyle='',color='green',label='xip')
#     plt.errorbar(theta*1.02,ngxim/i3xim, yerr=None, xerr=None,marker='o',linestyle='',color='blue',label='xim')
#     plt.errorbar(theta,slpp*np.log(theta)+bp, yerr=None, xerr=None,marker='',linestyle=':',color='green',label=str(np.around(np.mean(ngxip/i3xip),decimals=2))+' ('+str(np.around(np.mean(ngxip[0:len(theta)*8/10]/i3xip[0:len(theta)*8/10]),decimals=2))+')')
#     plt.errorbar(theta,slpm*np.log(theta)+bm, yerr=None, xerr=None,marker='',linestyle=':',color='blue',label=str(np.around(np.mean(ngxim/i3xim),decimals=2))+' ('+str(np.around(np.mean(ngxim[0:len(theta)*8/10]/i3xim[0:len(theta)*8/10]),decimals=2))+')')
#     plt.xscale('log')
#     plt.ylabel(r'$\xi(ngmix010)/\xi(im3shape7.2)$')
#     plt.xlabel(r'$\theta$ (arcmin)')
#     plt.xlim((1,500))
#     plt.ylim((-1,3))
#     plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=4, fancybox=True, shadow=True)
#     plt.savefig('nbins-'+str(i3.tbins)+'_biasorsens-'+str(i3.bs)+'_'+str(ng.bs)+'_weight-'+str(i3.wt)+'_'+str(ng.wt)+'_t0-'+str(i3.sep[0])+'_t1-'+str(i3.sep[1])+'_ztype-'+i3.ztyp+'_xi_ratio_'+label+'.png', bbox_inches='tight')
#     plt.close(10)


    
#     return


#   @staticmethod
#   def pz_pdf_stackW(cat):

#     from DES_pdf_stacker_v2 import return_stacked_pdf

#     def labelstr(bin):
#       if bin==0:
#         return '0.3<z<0.644W'
#       if bin==1:
#         return '0.644<z<0.901W'
#       if bin==2:
#         return '0.901<z<1.3W'
#       return ' '

#     binsarray=np.load('/home/zuntz/photoz/bins/bins.npy')
#     binsarray=binsarray[np.argsort(binsarray[:,0])]
#     mask=np.in1d(binsarray[:,0],cat.coadd,assume_unique=True)
#     binsarray=binsarray[mask]
#     mask=np.in1d(cat.coadd,binsarray[:,0],assume_unique=True)
#     tmp=cat.w[mask]
#     tmp2=cat.coadd[mask]

#     for bin in xrange(3):
#       print bin
#       bingals=binsarray[:,0][(binsarray[:,1]==bin+1)]
#       mask=np.in1d(tmp2,bingals,assume_unique=True)
#       tpz2=return_stacked_pdf(bingals,'TPZ2',tmp[mask])
#       bpz=return_stacked_pdf(bingals,'BPZ',tmp[mask])
#       tpz=return_stacked_pdf(bingals,'TPZ',tmp[mask])
#       zebra=return_stacked_pdf(bingals,'ZEBRA',tmp[mask])
#       skynet=return_stacked_pdf(bingals,'SKYNET',tmp[mask])
#       skynet2=return_stacked_pdf(bingals,'SKYNET_MAG_AUTO',tmp[mask])
#       annz=return_stacked_pdf(bingals,'ANNZ',tmp[mask])
#       np.save('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',tpz2)
#       np.save('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',bpz)
#       np.save('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',tpz)
#       np.save('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',skynet)
#       np.save('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',skynet2)
#       np.save('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',annz)
#       np.save('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',zebra)
#       if bin==0:
#         tpz2st=np.zeros((len(tpz2['z']),4))
#         bpzst=np.zeros((len(bpz['z']),4))
#         tpzst=np.zeros((len(tpz['z']),4))
#         zebrast=np.zeros((len(zebra['z']),4))
#         skynet2st=np.zeros((len(skynet2['z']),4))
#         skynetst=np.zeros((len(skynet['z']),4))
#         annzst=np.zeros((len(annz['z']),4))
#         tpz2st[:,0]=tpz2['z']
#         bpzst[:,0]=bpz['z']
#         tpzst[:,0]=tpz['z']
#         zebrast[:,0]=zebra['z']
#         skynet2st[:,0]=skynet2['z']
#         skynetst[:,0]=skynet['z']
#         annzst[:,0]=annz['z']
#       tpz2st[:,1+bin]=tpz2['pdf']
#       bpzst[:,1+bin]=bpz['pdf']
#       tpzst[:,1+bin]=tpz['pdf']
#       zebrast[:,1+bin]=zebra['pdf']
#       skynet2st[:,1+bin]=skynet2['pdf']
#       skynetst[:,1+bin]=skynet['pdf']
#       annzst[:,1+bin]=annz['pdf']

#     np.save('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',tpz2st)
#     np.save('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',bpzst)
#     np.save('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',tpzst)
#     np.save('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',zebrast)
#     np.save('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',skynet2st)
#     np.save('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',skynetst)
#     np.save('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.npy',annzst)

#     np.savetxt('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',tpz2st)
#     np.savetxt('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',bpzst)
#     np.savetxt('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',tpzst)
#     np.savetxt('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',zebrast)
#     np.savetxt('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',skynet2st)
#     np.savetxt('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',skynetst)
#     np.savetxt('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_fullW.txt',annzst)

#     return


#   @staticmethod
#   def pz_pdf_stack(cat):

#     from DES_pdf_stacker_v2 import return_stacked_pdf

#     def labelstr(bin):
#       if bin==0:
#         return '0.3<z<0.644'
#       if bin==1:
#         return '0.644<z<0.901'
#       if bin==2:
#         return '0.901<z<1.3'
#       return ' '

#     binsarray=np.load('/home/zuntz/photoz/bins/bins.npy')
#     binsarray=binsarray[np.argsort(binsarray[:,0])]
#     mask=np.in1d(binsarray[:,0],cat.coadd,assume_unique=True)
#     binsarray=binsarray[mask]
#     mask=np.in1d(cat.coadd,binsarray[:,0],assume_unique=True)
#     tmp=cat.w[mask]
#     tmp2=cat.coadd[mask]

#     for bin in xrange(3):
#       print bin
#       bingals=binsarray[:,0][(binsarray[:,1]==bin+1)]
#       mask=np.in1d(tmp2,bingals,assume_unique=True)
#       tpz2=return_stacked_pdf(bingals,'TPZ2',tmp[mask])
#       bpz=return_stacked_pdf(bingals,'BPZ',tmp[mask])
#       tpz=return_stacked_pdf(bingals,'TPZ',tmp[mask])
#       zebra=return_stacked_pdf(bingals,'ZEBRA',tmp[mask])
#       skynet=return_stacked_pdf(bingals,'SKYNET',tmp[mask])
#       skynet2=return_stacked_pdf(bingals,'SKYNET_MAG_AUTO',tmp[mask])
#       annz=return_stacked_pdf(bingals,'ANNZ',tmp[mask])
#       np.save('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',tpz2)
#       np.save('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',bpz)
#       np.save('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',tpz)
#       np.save('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',skynet)
#       np.save('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',skynet2)
#       np.save('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',annz)
#       np.save('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_'+labelstr(bin)+'.npy',zebra)
#       if bin==0:
#         tpz2st=np.zeros((len(tpz2['z']),4))
#         bpzst=np.zeros((len(bpz['z']),4))
#         tpzst=np.zeros((len(tpz['z']),4))
#         zebrast=np.zeros((len(zebra['z']),4))
#         skynet2st=np.zeros((len(skynet2['z']),4))
#         skynetst=np.zeros((len(skynet['z']),4))
#         annzst=np.zeros((len(annz['z']),4))
#         tpz2st[:,0]=tpz2['z']
#         bpzst[:,0]=bpz['z']
#         tpzst[:,0]=tpz['z']
#         zebrast[:,0]=zebra['z']
#         skynet2st[:,0]=skynet2['z']
#         skynetst[:,0]=skynet['z']
#         annzst[:,0]=annz['z']
#       tpz2st[:,1+bin]=tpz2['pdf']
#       bpzst[:,1+bin]=bpz['pdf']
#       tpzst[:,1+bin]=tpz['pdf']
#       zebrast[:,1+bin]=zebra['pdf']
#       skynet2st[:,1+bin]=skynet2['pdf']
#       skynetst[:,1+bin]=skynet['pdf']
#       annzst[:,1+bin]=annz['pdf']

#     np.save('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',tpz2st)
#     np.save('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',bpzst)
#     np.save('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',tpzst)
#     np.save('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',zebrast)
#     np.save('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',skynet2st)
#     np.save('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',skynetst)
#     np.save('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_full.npy',annzst)

#     np.savetxt('photoz/tpz2_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',tpz2st)
#     np.savetxt('photoz/bpz_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',bpzst)
#     np.savetxt('photoz/tpz_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',tpzst)
#     np.savetxt('photoz/zebra_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',zebrast)
#     np.savetxt('photoz/skynet2_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',skynet2st)
#     np.savetxt('photoz/skynet_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',skynetst)
#     np.savetxt('photoz/annz_'+CatalogMethods.cat_name(cat.cat)+'_full.txt',annzst)

#     return

#   @staticmethod
#   def sim_spread(cat,theta,xip0,xiperr,xim0,ximerr,xip,xim,ce,cb): 

#     for i in xrange(cat.num_patch):
#       if xip[i,0]==0:
#         continue
#       plt.figure(10)
#       plt.errorbar(theta,xip[i,:]*theta, yerr=None, xerr=None,label='')
#       plt.figure(11)
#       plt.errorbar(theta,xim[i,:]*theta, yerr=None, xerr=None,label='')
#       plt.figure(12)
#       plt.errorbar([1,2,3,4,5],ce[i,:], yerr=None, xerr=None,label='')
#       plt.figure(13)
#       plt.errorbar([1,2,3,4,5],cb[i,:], yerr=None, xerr=None,label='')

#     plt.figure(10)
#     plt.xscale('log')
#     plt.ylabel('xip*theta')
#     plt.xlabel('theta')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_xip.png', bbox_inches='tight')
#     plt.close(10)
#     plt.figure(11)
#     plt.xscale('log')
#     plt.ylabel('xim')
#     plt.xlabel('theta')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_xim.png', bbox_inches='tight')
#     plt.close(11)
#     plt.figure(12)
#     plt.xscale('log')
#     plt.ylabel('ce')
#     plt.xlabel('bin')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_ce.png', bbox_inches='tight')
#     plt.close(12)
#     plt.figure(13)
#     plt.xscale('log')
#     plt.ylabel('cb')
#     plt.xlabel('bin')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_cb.png', bbox_inches='tight')
#     plt.close(13)

#     plt.figure(14)
#     ax=plt.subplot(111)
#     plt.errorbar(theta,xip0*theta, yerr=xiperr*theta, xerr=None,label='')
#     ax.fill_between(theta,np.min(xip[xip[:,0]!=0],axis=0)*theta,np.max(xip[xip[:,0]!=0],axis=0)*theta,facecolor='gray',alpha=0.25)
#     ax.fill_between(theta,xip0*theta-np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0)/1.2)*theta,xip0*theta+np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0)/1.2)*theta,facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.ylabel('xip*theta')
#     plt.xlabel('theta')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_var_xip.png', bbox_inches='tight')
#     plt.close(14)

#     plt.figure(15)
#     ax=plt.subplot(111)
#     plt.errorbar(theta,xim0*theta, yerr=ximerr*theta, xerr=None,label='')
#     ax.fill_between(theta,np.min(xim[xim[:,0]!=0],axis=0)*theta,np.max(xim[xim[:,0]!=0],axis=0)*theta,facecolor='gray',alpha=0.25)
#     ax.fill_between(theta,xim0*theta-np.sqrt(np.var(xim[xim[:,0]!=0]*theta,axis=0)),xim0*theta+np.sqrt(np.var(xim[xim[:,0]!=0]*theta,axis=0)),facecolor='gray',alpha=0.5)
#     plt.xscale('log')
#     plt.ylabel('xim*theta')
#     plt.xlabel('theta')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_var_xim.png', bbox_inches='tight')
#     plt.close(15)

#     plt.figure(16)
#     ax=plt.subplot(111)
#     plt.errorbar(theta,xiperr/xip0, yerr=None, xerr=None,label='xiperr/xip')
#     plt.errorbar(theta,np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0))/xip0, yerr=None, xerr=None,label='sqrt(var)/xip')
#     plt.errorbar(theta,np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0))/xiperr, yerr=None, xerr=None,label='sqrt(var)/xiperr')
#     plt.errorbar(theta,np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0)/1.2)/xip0, yerr=None, xerr=None,label='sqrt(var)/xip')
#     plt.errorbar(theta,np.sqrt(np.var(xip[xip[:,0]!=0],ddof=1,axis=0)/1.2)/xiperr, yerr=None, xerr=None,label='sqrt(var)/xiperr')
#     plt.legend(loc='upper left')
#     plt.ylim((0,5))
#     plt.xscale('log')
#     plt.ylabel('')
#     plt.xlabel('theta')
#     plt.savefig(CatalogMethods.cat_name(cat.cat)+'_nbins-'+str(cat.lbins)+'_biasorsens-'+str(cat.bs)+'_weight-'+str(cat.wt)+'_jk-'+str(cat.use_jk)+'_sim_spread_var_xip_ratio.png', bbox_inches='tight')
#     plt.close(16)

#     return   

#   @staticmethod
#   def band_comp(i3v8,ng,ngr,ngi,ngz): 

#     maskr=np.in1d(ngr.coadd,i3v8.coadd,assume_unique=True)&np.in1d(ngr.coadd,ng.coadd,assume_unique=True)&np.in1d(ngr.coadd,ngi.coadd,assume_unique=True)&np.in1d(ngr.coadd,ngz.coadd,assume_unique=True)
#     maski3=np.in1d(i3v8.coadd,ngr.coadd[maskr],assume_unique=True)
#     maskng=np.in1d(ng.coadd,ngr.coadd[maskr],assume_unique=True)
#     maski=np.in1d(ngi.coadd,ngr.coadd[maskr],assume_unique=True)
#     maskz=np.in1d(ngz.coadd,ngr.coadd[maskr],assume_unique=True)

#     # i3v8.e1=i3v8.e1[maski3]
#     # i3v8.e2=i3v8.e2[maski3]
#     # ng.e1=ng.e1[maskng]
#     # ng.e2=ng.e2[maskng]
#     # ngr.e1=ngr.e1[maskr]
#     # ngr.e2=ngr.e2[maskr]
#     # ngi.e1=ngi.e1[maski]
#     # ngi.e2=ngi.e2[maski]
#     # ngz.e1=ngz.e1[maskz]
#     # ngz.e2=ngz.e2[maskz]

#     print 'number of gal',len(i3v8.e1)

#     sorti3=np.argsort(i3v8.coadd[maski3])[np.argsort(np.argsort(ngr.coadd[maskr]))]
#     sortng=np.argsort(ng.coadd[maskng])[np.argsort(np.argsort(ngr.coadd[maskr]))]
#     sorti=np.argsort(ngi.coadd[maski])[np.argsort(np.argsort(ngr.coadd[maskr]))]
#     sortz=np.argsort(ngz.coadd[maskz])[np.argsort(np.argsort(ngr.coadd[maskr]))]

#     print '|i3-ng|',np.mean(i3v8.e1[sorti3]-i3v8.c1[sorti3]-ng.e1[sortng]),np.mean(i3v8.e2[sorti3]-i3v8.c2[sorti3]-ng.e2[sortng])
#     print '|i3-ngr|',np.mean(i3v8.e1[sorti3]-i3v8.c1[sorti3]-ngr.e1),np.mean(i3v8.e2[sorti3]-i3v8.c2[sorti3]-ngr.e2)
#     print '|i3-ngi|',np.mean(i3v8.e1[sorti3]-i3v8.c1[sorti3]-ngi.e1[sorti]),np.mean(i3v8.e2[sorti3]-i3v8.c2[sorti3]-ngi.e2[sorti])
#     print '|i3-ngz|',np.mean(i3v8.e1[sorti3]-i3v8.c1[sorti3]-ngz.e1[sortz]),np.mean(i3v8.e2[sorti3]-i3v8.c2[sorti3]-ngz.e2[sortz])
#     print '|ngr-ng|',np.mean(ngr.e1-ng.e1[sortng]),np.mean(ngr.e2-ng.e2[sortng])
#     print '|ngi-ng|',np.mean(ngi.e1[sorti]-ng.e1[sortng]),np.mean(ngi.e2[sorti]-ng.e2[sortng])
#     print '|ngz-ng|',np.mean(ngz.e1[sortz]-ng.e1[sortng]),np.mean(ngz.e2[sortz]-ng.e2[sortng])

#     mask=np.ones((len(i3.coadd)))

#     theta,i3xip,i3xiperr,i3xim,i3ximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(i3,mask.astype(bool),mask)
#     theta,ngxip,ngxiperr,ngxim,ngximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ng,mask.astype(bool),mask)

#     theta,ngxip,ngxiperr,ngxim,ngximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ngr,mask.astype(bool),mask)
#     theta,ngxip,ngxiperr,ngxim,ngximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ngi,mask.astype(bool),mask)
#     theta,ngxip,ngxiperr,ngxim,ngximerr,chi2p,chi2m,tmp,tmp,tmp,tmp,tmp,tmp,tmp=xi_2pt_shear_test_methods.xi_2pt(ngz,mask.astype(bool),mask)

#     return

#   @staticmethod
#   def e_2d_hist(cat1,cat2,ngr,nbins=100):

#     mask1=np.in1d(cat1.coadd,cat2.coadd,assume_unique=True)&(np.in1d(cat1.coadd,ngr.coadd,assume_unique=True))
#     mask2=np.in1d(cat2.coadd,cat1.coadd,assume_unique=True)&(np.in1d(cat2.coadd,ngr.coadd,assume_unique=True))

#     e11=cat1.e1[mask1]
#     e12=cat1.e2[mask1]
#     e21=cat2.e1[mask2]
#     e22=cat2.e2[mask2]

#     sort1=np.argsort(cat1.coadd[mask1])[np.argsort(np.argsort(cat2.coadd[mask2]))]
#     sort2=np.ones(len(e22)).astype(bool)

#     plt.figure()
#     plt.hist2d(e11[sort1],e21[sort2],bins=nbins,norm=LogNorm())
#     plt.xlabel('e1 - '+CatalogMethods.cat_name(cat1.cat))
#     plt.ylabel('e1 - '+CatalogMethods.cat_name(cat2.cat))
#     plt.xlim((-1,1))
#     plt.ylim((-1,1))
#     plt.colorbar()
#     plt.savefig(CatalogMethods.cat_name(cat1.cat)+'_'+CatalogMethods.cat_name(cat2.cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(cat1.bs)+'_'+str(cat2.bs)+'_e11_2dhist.png', bbox_inches='tight')
#     plt.close()

#     plt.figure()
#     plt.hist2d(e12[sort1],e22[sort2],bins=nbins,norm=LogNorm())
#     plt.xlabel('e2 - '+CatalogMethods.cat_name(cat1.cat))
#     plt.ylabel('e2 - '+CatalogMethods.cat_name(cat2.cat))
#     plt.xlim((-1,1))
#     plt.ylim((-1,1))
#     plt.colorbar()
#     plt.savefig(CatalogMethods.cat_name(cat1.cat)+'_'+CatalogMethods.cat_name(cat2.cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(cat1.bs)+'_'+str(cat2.bs)+'_e22_2dhist.png', bbox_inches='tight')
#     plt.close()

#     plt.figure()
#     plt.hist2d(e11[sort1],e22[sort2],bins=nbins,norm=LogNorm())
#     plt.xlabel('e1 - '+CatalogMethods.cat_name(cat1.cat))
#     plt.ylabel('e2 - '+CatalogMethods.cat_name(cat2.cat))
#     plt.xlim((-1,1))
#     plt.ylim((-1,1))
#     plt.colorbar()    
#     plt.savefig(CatalogMethods.cat_name(cat1.cat)+'_'+CatalogMethods.cat_name(cat2.cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(cat1.bs)+'_'+str(cat2.bs)+'_e12_2dhist.png', bbox_inches='tight')
#     plt.close()

#     plt.figure()
#     plt.hist2d(e12[sort1],e21[sort2],bins=nbins,norm=LogNorm())
#     plt.xlabel('e2 - '+CatalogMethods.cat_name(cat1.cat))
#     plt.ylabel('e1 - '+CatalogMethods.cat_name(cat2.cat))
#     plt.xlim((-1,1))
#     plt.ylim((-1,1))
#     plt.colorbar()
#     plt.savefig(CatalogMethods.cat_name(cat1.cat)+'_'+CatalogMethods.cat_name(cat2.cat)+'_nbins-'+str(nbins)+'_biasorsens-'+str(cat1.bs)+'_'+str(cat2.bs)+'_e21_2dhist.png', bbox_inches='tight')
#     plt.close()

#     return
