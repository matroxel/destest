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
import src.corr as corr

class y1(object):

    @staticmethod
    def load_data(i3file,mcalfile,goldfile):

        goldcols = []

        mcal = catalog.CatalogStore('matched_metacal',cutfunc=catalog.CatalogMethods.matched_metacal_cut(),cutfunclive=catalog.CatalogMethods.matched_metacal_cut_live(),cattype='mcal',catfile=mcalfile,goldfile=goldfile,goldcols=goldcols)

        mcal.bs    = True
        mcal.wt    = False
        mcal.lbins = 20

        #np.save('mcal_coadds.npy',np.vstack((mcal.coadd,np.ones(len(mcal.coadd)),np.ones(len(mcal.coadd)),np.zeros(len(mcal.coadd)),np.zeros(len(mcal.coadd)),mcal.w)).T)

        i3 = catalog.CatalogStore('matched_i3',cutfunc=catalog.CatalogMethods.matched_i3_cut(),cattype='i3',catfile=i3file,goldfile=goldfile,goldcols=goldcols)

        i3.bs    = True
        i3.wt    = True
        i3.lbins = 20

        #np.save('i3_coadds.npy',np.vstack((i3.coadd,i3.m1,i3.m2,i3.c1,i3.c2,i3.w)).T)

        return i3,mcal

    @staticmethod
    def load_psf_data(psfdir):

        cols = ['e1','e2','ccd','col','row','psf1','psf2','mag']
        psf  = catalog.CatalogStore('psf',cutfunc=catalog.CatalogMethods.final_null_cuts_ra_flag(),cols=cols,cattype='psf',catdir=psfdir)

        return psf

class y1_plots(object):

    @staticmethod
    def save_obj(obj, name ):
        with open('obj/'+ name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load_obj(name ):
        with open('obj/' + name + '.pkl', 'rb') as f:
            return pickle.load(f)

    @staticmethod
    def mean_e(cat1,cat2,replace=False):

        txt.write_methods.heading('Linear Splits',cat1,label='y1_paper',create=False)

        y1_plots.evspsf1(cat1,cat2,replace=replace)
        y1_plots.evspsf2(cat1,cat2,replace=replace)
        y1_plots.evssnr(cat1,cat2,replace=replace)
        y1_plots.evsradius(cat1,cat2,replace=replace)
        y1_plots.evspsfsize(cat1,cat2,replace=replace)

    @staticmethod
    def mean_e_epoch(cat1,cat2):

        txt.write_methods.heading('Linear Splits on epoch quantities',cat1,label='y1_paper',create=False)

        y1_plots.evsrow(cat1,cat2)



    @staticmethod
    def mean_e_subplot(cat,n,val,fig,replace=False):

        name=fig0.plot_methods.get_filename_str(cat)
        name = '/text/lin_'+name+'_'+val+'.txt'

        if replace|(os.path.exists(name)):
            array=getattr(cat,val)
            if isinstance(cat,catalog.CatalogStore):
                mask=catalog.CatalogMethods.check_mask(cat.coadd,None)
            tmp,tmp,arr2,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=sys_split.split_gals_lin_along_base(cat,val,array,mask,name,log=config.log_val.get(val,False),plot=False)

            if config.log_val.get(val,False):
                arr1=10**arr2
            else:
                arr1=arr2

            d = {

            'arr1' : arr1,
            'arr2' : arr2,
            'e1' : e1,
            'e2' : e2,
            'e1err' : e1err,
            'e2err' : e2err,
            'm1' : m1,
            'm2' : m2,
            'b1' : b1,
            'b2' : b2
            }

            y1_plots.save_obj(d,name)
 
        else:

            d = y1_plots.load_obj(name)

        plt.figure(fig)
        ax=plt.subplot(2,1,n)
        plt.errorbar(d['arr1'],d['e1'],yerr=d['e1err'],marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
        plt.errorbar(d['arr1'],d['m1']*d['arr2']+d['b1'],marker='',linestyle='-',color='r')
        plt.errorbar(d['arr1']+(d['arr1'][1]-d['arr1'][0])/5.,d['e2'],yerr=d['e2err'],marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
        plt.errorbar(d['arr1'],d['m2']*d['arr2']+d['b2'],marker='',linestyle='-',color='b')
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
    def evssnr(cat1,cat2,replace=False):

        plt.figure(1)

        y1_plots.mean_e_subplot(cat1,0,'snr',1,replace=replace)
        y1_plots.mean_e_subplot(cat2,1,'snr',1,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_snr.pdf', bbox_inches='tight')
        plt.close(1)

        return

    @staticmethod
    def evsradius(cat1,cat2,replace=False):

        plt.figure(2)

        y1_plots.mean_e_subplot(cat1,0,'rgpp_rp',2,replace=replace)
        y1_plots.mean_e_subplot(cat2,1,'size',2,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_radius.pdf', bbox_inches='tight')
        plt.close(2)

        return

    @staticmethod
    def evspsf1(cat1,cat2,replace=False):

        plt.figure(3)

        y1_plots.mean_e_subplot(cat1,0,'psf1',3,replace=replace)
        y1_plots.mean_e_subplot(cat2,1,'psf1',3,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf1.pdf', bbox_inches='tight')
        plt.close(3)

        return

    @staticmethod
    def evspsf2(cat1,cat2,replace=False):

        plt.figure(4)

        y1_plots.mean_e_subplot(cat1,0,'psf2',4,replace=replace)
        y1_plots.mean_e_subplot(cat2,1,'psf2',4,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf2.pdf', bbox_inches='tight')
        plt.close(4)

        return

    @staticmethod
    def evspsfsize(cat1,cat2,replace=False):

        plt.figure(5)

        y1_plots.mean_e_subplot(cat1,0,'psffwhm',5,replace=replace)
        y1_plots.mean_e_subplot(cat2,1,'psffwhm',5,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psfsize.pdf', bbox_inches='tight')
        plt.close(5)

        return


    def evsrow(i3epoch,mcepoch,replace=False):

        plt.figure(6)

        y1_plots.mean_e_subplot(i3epoch,0,'row',6,replace=replace)
        y1_plots.mean_e_subplot(mcepoch,1,'row',6,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_ccdrow.pdf', bbox_inches='tight')
        plt.close(6)


    @staticmethod
    def psf_whisker(psf,replace=False):

        psf.dpsf1 = psf.psf1-psf.e1
        psf.dpsf2 = psf.psf2-psf.e2

        y1_plots.whiskerplot(psf,'psf',6,replace=replace)
        y1_plots.whiskerplot(psf,'dpsf',7,replace=replace)

        return

    @staticmethod
    def e_whisker(cat1dir,cat1,cat2dir,cat2,replace=False):

        cols=['coadd','expnum','ccd','row','col']
        epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
        y1_plots.whiskerplot(epoch1,'e',8)

        cols=['coadd','expnum','ccd','row','col']
        epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
        y1_plots.whiskerplot(epoch2,'e',9)

        return

    @staticmethod
    def whiskerplot(cat,col,fig,replace=False):

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

        

    @staticmethod 
    def tangential_shear_plot(i3, metacal, centers, centers_mask=None):

        mask = None

        i3_theta,i3_out,i3_err,i3_chi2 = corr.xi_2pt.xi_2pt(centers, i3, corr='NG', 
                                                            maska=centers_mask, maskb=mask, ran=True)
        i3_gammat = i3_out[0]
        i3_gammax = i3_out[2]
        i3_gammat_err = i3_err[0]
        i3_gammax_err = i3_err[2]
        
        mc_theta,mc_out,mc_err,mc_chi2 = corr.xi_2pt.xi_2pt(centers, metacal, corr='NG', 
                                                            maska=centers_mask, maskb=mask, ran=True)
        mc_gammat = mc_out[0]
        mc_gammax = mc_out[2]
        mc_gammat_err = mc_err[0]
        mc_gammax_err = mc_err[2]
        plt.figure()
        plt.errorbar(i3_theta, i3_gammat, i3_gammat_err, fmt='r.', label='Im3shape')
        plt.errorbar(mc_theta, mc_gammat, mc_gammat_err, fmt='b.', label='Metacal')
        plt.xscale('log')
#        plt.yscale('log', nonposy='clip')
        plt.savefig('plots/y1/special_gammat.pdf', dpi=500, bbox_inches='tight')
        plt.close()
        print "We think the imaginary gammat is gammax but not sure!"
        plt.figure()
        plt.errorbar(i3_theta, i3_gammax, i3_gammax_err, fmt='r.', label='Im3shape')
        plt.errorbar(mc_theta, mc_gammax, mc_gammax_err, fmt='b.', label='Metacal')
        plt.xscale('log')
#        plt.yscale('log', nonposy='clip')
        plt.savefig('plots/y1/special_gammax.pdf', dpi=500, bbox_inches='tight')
        plt.close()

    @staticmethod
    def b_mode_plot(i3, metacal):
        bp = corr.bandpowers()

        #We might want to move this as we may need it more than once
        i3_theta,i3_out,i3_err,i3_chi2 = corr.xi_2pt.xi_2pt(i3)
        i3_xip = i3_out[0]
        i3_xim = i3_out[1]
        i3_xiperr = i3_err[0]
        i3_ximerr = i3_err[1]
        i3_E,i3_B,i3_E_err,i3_B_err = bp.bandpowersEB(i3_xip,i3_xim,i3_xiperrm,i3_ximerr)

                #We might want to move this as we may need it more than once
        mc_theta,mc_out,mc_err,mc_chi2 = corr.xi_2pt.xi_2pt(metacal)
        mc_xip = mc_out[0]
        mc_xim = mc_out[1]
        mc_xiperr = mc_err[0]
        mc_ximerr = mc_err[1]
        mc_E,mc_B,mc_E_err,mc_B_err = bp.bandpowersEB(mc_xip,mc_xim,mc_xiperrm,mc_ximerr)

        ell = [bp.lm(i) for i in xrange(bp.nell)]

        plt.figure()
        plt.errorbar(ell, i3_B, i3_B_err, fmt='r.', label='Im3shape')
        plt.errorbar(ell, mc_B, mc_B_err, fmt='b.', label='Metacal')
        plt.savefig('plots/y1/b_modes.pdf', dpi=500, bbox_inches='tight')
        plt.close()






