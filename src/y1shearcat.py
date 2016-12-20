import os
import numpy as np
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
import src.txt as txt
import src.fig as fig0
import src.field as field

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

    @staticmethod
    def load_psf_data(psfdir):

        cols = ['e1','e2','ccd','col','row','psf1','psf2','mag']
        psf  = catalog.CatalogStore('psf',cutfunc=catalog.CatalogMethods.final_null_cuts_ra_flag(),cols=cols,cattype='psf',catdir=psfdir)

        return psf

class y1_plots(object):

    @staticmethod
    def mean_e(cat1,cat2):

        txt.write_methods.heading('Linear Splits',cat1,label='y1_paper',create=False)

        y1_plots.evssnr(cat1,cat2)
        y1_plots.evsradius(cat1,cat2)
        y1_plots.evspsf1(cat1,cat2)
        y1_plots.evspsf2(cat1,cat2)
        y1_plots.evspsfsize(cat1,cat2)

    @staticmethod
    def mean_e_epoch(cat1,cat2):

        txt.write_methods.heading('Linear Splits on epoch quantities',cat1,label='y1_paper',create=False)

        y1_plots.evsrow(cat1,cat2)



    @staticmethod
    def mean_e_subplot(cat,n,val,fig):

        array=getattr(cat,val)
        name=fig0.plot_methods.get_filename_str(cat)
        if isinstance(cat,catalog.CatalogStore):
            mask=catalog.CatalogMethods.check_mask(cat.coadd,None)
        tmp,tmp,arr2,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=sys_split.split_gals_lin_along_base([cat.cat,cat.bs,cat.wt,cat.e1,cat.e2,cat.m1,cat.m2,cat.c1,cat.c2,cat.w],val,array,mask,name,log=config.log_val.get(val,False),plot=False)

        if config.log_val.get(val,False):
            arr1=10**arr2
        else:
            arr1=arr2

        plt.figure(fig)
        ax=plt.subplot(2,1,n)
        plt.errorbar(arr1,e1,yerr=e1err,marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
        plt.errorbar(arr1,m1*arr2+b1,marker='',linestyle='-',color='r')
        plt.errorbar(arr1+(arr1[1]-arr1[0])/5.,e2,yerr=e2err,marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
        plt.errorbar(arr1,m2*arr2+b2,marker='',linestyle='-',color='b')
        ax.minorticks_on()
        plt.ylabel(r'$\langle e \rangle$')
        if config.log_val.get(val,False):
            plt.xscale('log')
        # plt.axhline(.004,color='k')
        # plt.axhline(-.004,color='k')
        plt.ylim((-0.002,0.002))
        plt.xlabel(config.lbl.get(val,val.replace('_','-')))
        if n==1:
            plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})


    @staticmethod
    def evssnr(cat1,cat2):

        plt.figure(1)

        y1_plots.mean_e_subplot(cat1,0,'snr',1)
        y1_plots.mean_e_subplot(cat2,1,'snr',1)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_snr.pdf', bbox_inches='tight')
        plt.close(1)

        return

    @staticmethod
    def evsradius(cat1,cat2):

        plt.figure(2)

        y1_plots.mean_e_subplot(cat1,0,'rgp',2)
        y1_plots.mean_e_subplot(cat2,1,'size',2)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_radius.pdf', bbox_inches='tight')
        plt.close(2)

        return

    @staticmethod
    def evspsf1(cat1,cat2):

        plt.figure(3)

        y1_plots.mean_e_subplot(cat1,0,'psf1',3)
        y1_plots.mean_e_subplot(cat2,1,'psf1',3)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf1.pdf', bbox_inches='tight')
        plt.close(3)

        return

    @staticmethod
    def evspsf2(cat1,cat2):

        plt.figure(4)

        y1_plots.mean_e_subplot(cat1,0,'psf2',4)
        y1_plots.mean_e_subplot(cat2,1,'psf2',4)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf2.pdf', bbox_inches='tight')
        plt.close(4)

        return

    @staticmethod
    def evspsfsize(cat1,cat2):

        plt.figure(5)

        y1_plots.mean_e_subplot(cat1,0,'psffwhm',5)
        y1_plots.mean_e_subplot(cat2,1,'psffwhm',5)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psfsize.pdf', bbox_inches='tight')
        plt.close(5)

        return


    def evsrow(i3epoch,mcepoch):

        plt.figure(6)

        y1_plots.mean_e_subplot(i3epoch,0,'row',6)
        y1_plots.mean_e_subplot(mcepoch,1,'row',6)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_ccdrow.pdf', bbox_inches='tight')
        plt.close(6)


    @staticmethod
    def psf_whisker(psf):

        psf.dpsf1 = psf.psf1-psf.e1
        psf.dpsf2 = psf.psf2-psf.e2

        y1_plots.whiskerplot(psf,'psf',6)
        y1_plots.whiskerplot(psf,'dpsf',7)

        return

    @staticmethod
    def e_whisker(cat1dir,cat1,cat2dir,cat2):

        cols=['coadd','expnum','ccd','row','col']
        epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
        y1_plots.whiskerplot(epoch1,'e',8)

        cols=['coadd','expnum','ccd','row','col']
        epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
        y1_plots.whiskerplot(epoch2,'e',9)

        return

    @staticmethod
    def whiskerplot(cat,col,fig):

        if col == 'psf':
            key = r'e_{PSF}'
        if col == 'e':
            key = r'e'
        if col == 'dpsf':
            key = r'\Delta e_{PSF}'

        scale=0.02

        y,x,mw,e1,e2,e=field.field.whisker_calc(cat,col=col)
        pos0=0.5*np.arctan2(e2/mw,e1/mw)
        e/=mw
        for i in range(len(x)):
            y[i,:,:],x[i,:,:]=field.field_methods.ccd_to_field(i,y[i,:,:]-2048,x[i,:,:]-1024)

        print 'y,x',y[i,:,:],x[i,:,:]

        plt.figure(fig)
        print np.shape(x),np.shape(y),np.shape(np.sin(pos0)*e),np.shape(np.cos(pos0)*e)
        Q = plt.quiver(np.ravel(y),np.ravel(x),np.ravel(np.sin(pos0)*e),np.ravel(np.cos(pos0)*e),units='width',pivot='middle',headwidth=0,width=.0005)
        plt.quiverkey(Q,0.2,0.2,scale,str(scale)+' '+key,labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
        plt.savefig('plots/y1/whisker_'+col+'.pdf', dpi=500, bbox_inches='tight')
        plt.close(fig)

        return

