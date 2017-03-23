import os
import numpy as np
import pickle
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
import fitsio as fio
import src.sys_split as sys_split
import src.config as config
import src.lin as lin
import src.catalog as catalog
import src.txt as txt
import src.fig as fig0
import src.field as field
import src.corr as corr

#pzbins = [0.2,0.43,0.63,0.9,1.3]

def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)

class y1(object):

    @staticmethod
    def load_data(i3file,mcalfile,goldfile,bpzfile,bpz0file,i3pickle,mcalpickle,load_pickle=True):


        if (load_pickle)&(os.path.exists(mcalpickle))&(os.path.exists(i3pickle)):

            mcal = load_obj(mcalpickle)
            i3   = load_obj(i3pickle)

        else:

            bpz = catalog.PZStore('bpz',setup=True,pztype='bpz',filetype='fits',file=bpz0file,sheared=True,nofzfile=bpzfile)

            goldcols = []

            mcal = catalog.CatalogStore('metacalibration',cutfunc=catalog.CatalogMethods.matched_metacal_cut(),cutfunclive=catalog.CatalogMethods.matched_metacal_cut_live(),cattype='mcal',catfile=mcalfile,goldfile=goldfile,goldcols=goldcols)

            # manually cut out S82
            catalog.CatalogMethods.match_cat(mcal,mcal.dec<-35)

            # manually cut out parts that don't contribute to selection correction to save memory
            catalog.CatalogMethods.match_cat(mcal,0==((mcal.flags_select&2**0)
                                                    |(mcal.flags_select&2**1)
                                                    |(mcal.flags_select&2**2)
                                                    |(mcal.flags_select&2**3)
                                                    |(mcal.flags_select&2**6)
                                                    |(mcal.flags_select&2**7)
                                                    |(mcal.flags_select&2**30)))
            # manually cut out parts that do contribute to selection correction but are never selected to save memory
            catalog.CatalogMethods.match_cat(mcal,(~((mcal.flags_select&2**4!=0)
                                                    &(mcal.flags_select_1p&2**4!=0)
                                                    &(mcal.flags_select_1m&2**4!=0)
                                                    &(mcal.flags_select_2p&2**4!=0)
                                                    &(mcal.flags_select_2m&2**4!=0)))
                                                    &(~((mcal.flags_select&2**5!=0)
                                                    &(mcal.flags_select_1p&2**5!=0)
                                                    &(mcal.flags_select_1m&2**5!=0)
                                                    &(mcal.flags_select_2p&2**5!=0)
                                                    &(mcal.flags_select_2m&2**5!=0)))
                                                    &(~((mcal.flags_select&2**8!=0)
                                                    &(mcal.flags_select_1p&2**8!=0)
                                                    &(mcal.flags_select_1m&2**8!=0)
                                                    &(mcal.flags_select_2p&2**8!=0)
                                                    &(mcal.flags_select_2m&2**8!=0))))
            # 'dec_cut':2**0,
            # 'flags_badregion':2**1,
            # 'flags_gold':2**2,
            # 'mcal_flags':2**3,
            # 'size_cut':2**4,     # used in corrections
            # 'snr_cut':2**5,      # used in corrections
            # 'mof_flags':2**6,
            # 'mof_lowflux':2**7,
            # 'mcal_lowflux':2**8, # used in corrections
            # 'not_measured':2**30,

            mcal.bs    = True
            mcal.wt    = False
            mcal.lbins = 20
            mcal.add_pz(bpz,sheared=True,bounds=[0.2,1.3])

            #np.save('mcal_coadds.npy',np.vstack((mcal.coadd,np.ones(len(mcal.coadd)),np.ones(len(mcal.coadd)),np.zeros(len(mcal.coadd)),np.zeros(len(mcal.coadd)),mcal.w)).T)

            i3 = catalog.CatalogStore('im3shape',cutfunc=catalog.CatalogMethods.matched_i3_cut(),cattype='i3',catfile=i3file,goldfile=goldfile,goldcols=goldcols)

            # manually cut out S82
            catalog.CatalogMethods.match_cat(i3,i3.dec<-35)

            i3.bs    = True
            i3.wt    = True
            i3.lbins = 20
            i3.add_pz(bpz,sheared=False,bounds=[0.2,1.3])

            #np.save('i3_coadds.npy',np.vstack((i3.coadd,i3.m1,i3.m2,i3.c1,i3.c2,i3.w)).T)

            bpz = None

            save_obj(mcal,mcalpickle)
            save_obj(i3,i3pickle)

        return i3,mcal

    @staticmethod
    def load_psf_data(psfdir):

        # MAX_CENTROID_SHIFT = 1.0

        # NOT_USED = 1
        # MEAS_BAD_MEASUREMENT = 2
        # MEAS_CENTROID_SHIFT = 4
        # PSFEX_BAD_MEASUREMENT = 8
        # PSFEX_CENTROID_SHIFT = 16
        # PSFEX_FAILURE = 32
        # RESERVED = 64
        # #DESDM_BAD_MEASUREMENT = 64  # Prior to y1a1-v13, this was the meaning of 64-256
        # #DESDM_CENTROID_SHIFT = 128
        # #DESDM_FAILURE = 256
        # #DESDM_FLAG_FACTOR = DESDM_BAD_MEASUREMENT / PSFEX_BAD_MEASUREMENT
        # BLACK_FLAG_FACTOR = 512 # blacklist flags are this times the original exposure blacklist flag
        #                         # blacklist flags go up to 64, so this uses up to 1<<15

        cols = ['e1','e2','ccd','col','row','psf1','psf2','mag','size','psf_size','flag']
        psf  = catalog.CatalogStore('psf',cutfunc=catalog.CatalogMethods.final_null_cuts_ra_flag(),cols=cols,cattype='psf',catdir=psfdir)

        return psf

class y1_plots(object):

    x_lims = {
    'snr' : (8,400),
    'psf1': (-0.024,0.044),
    'psf2': (-0.024,0.044),
    'psfsize': None,
    'psffwhm': (3.,4.8),
    'size': None,
    'rgp': (1.,2.2)
    }

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
        name = 'text/lin_'+name+'_'+val+'.pkl'

        if replace|(not os.path.exists(name)):

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

            save_obj(d,name)
 
        else:

            d = load_obj(name)

        if val in ['size','psfsize']:
            d['arr1']=np.sqrt(d['arr1'])

        plt.figure(fig)
        ax=plt.subplot(2,1,n)
        plt.errorbar(d['arr1'],d['e1'],yerr=d['e1err'],marker='o',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
        plt.errorbar(d['arr1'],d['m1']*d['arr2']+d['b1'],marker='',linestyle='-',color='r')
        plt.errorbar(d['arr1'],d['e2'],yerr=d['e2err'],marker='o',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
        plt.errorbar(d['arr1'],d['m2']*d['arr2']+d['b2'],marker='',linestyle='-',color='b')
        ax.minorticks_on()
        plt.ylabel(r'$\langle e \rangle$')
        if config.log_val.get(val,False):
            plt.xscale('log')
        plt.axhline(0.,color='k')
        plt.ylim((-0.0015,0.0015))
        if y1_plots.x_lims[val] is not None:
            plt.xlim(y1_plots.x_lims[val])
        plt.xlabel(config.lbl.get(val,val.replace('_','-')))
        plt.title(cat.name)
        if n==1:
            plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})


    @staticmethod
    def evssnr(cat1,cat2,replace=False):

        plt.figure(1)

        y1_plots.mean_e_subplot(cat1,1,'snr',1,replace=replace)
        y1_plots.mean_e_subplot(cat2,2,'snr',1,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_snr.pdf', bbox_inches='tight')
        plt.close(1)

        return

    @staticmethod
    def evsradius(cat1,cat2,replace=False):

        plt.figure(2)

        y1_plots.mean_e_subplot(cat1,1,'rgp',2,replace=replace)
        y1_plots.mean_e_subplot(cat2,2,'size',2,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_radius.pdf', bbox_inches='tight')
        plt.close(2)

        return

    @staticmethod
    def evspsf1(cat1,cat2,replace=False):

        plt.figure(3)

        y1_plots.mean_e_subplot(cat1,1,'psf1',3,replace=replace)
        y1_plots.mean_e_subplot(cat2,2,'psf1',3,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf1.pdf', bbox_inches='tight')
        plt.close(3)

        return

    @staticmethod
    def evspsf2(cat1,cat2,replace=False):

        plt.figure(4)

        y1_plots.mean_e_subplot(cat1,1,'psf2',4,replace=replace)
        y1_plots.mean_e_subplot(cat2,2,'psf2',4,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf2.pdf', bbox_inches='tight')
        plt.close(4)

        return

    @staticmethod
    def evspsfsize(cat1,cat2,replace=False):

        plt.figure(5)

        y1_plots.mean_e_subplot(cat1,1,'psffwhm',5,replace=replace)
        y1_plots.mean_e_subplot(cat2,2,'psfsize',5,replace=replace)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psfsize.pdf', bbox_inches='tight')
        plt.close(5)

        return


    def evsrow(i3epoch,mcepoch,replace=False):

        plt.figure(6)

        y1_plots.mean_e_subplot(i3epoch,1,'row',6,replace=replace)
        y1_plots.mean_e_subplot(mcepoch,2,'row',6,replace=replace)

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

    @staticmethod
    def footprint_sub(ra,dec,rasep,decsep,nside,fig):

        import skymapper as skm

        bc, ra, dec, vertices = skm.getCountAtLocations(ra, dec, nside=nside, return_vertices=True)

        # setup figure
        import matplotlib.cm as cm
        cmap = cm.YlOrRd
        ax = fig.add_subplot(111, aspect='equal')

        # setup map: define AEA map optimal for given RA/Dec
        proj = skm.createConicMap(ax, ra, dec, proj_class=skm.AlbersEqualAreaProjection)
        # add lines and labels for meridians/parallels (separation 5 deg)
        sep = 5
        meridians = np.arange(-90, 90+decsep, decsep)
        parallels = np.arange(0, 360+rasep, rasep)
        skm.setMeridianPatches(ax, proj, meridians, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
        skm.setParallelPatches(ax, proj, parallels, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
        skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
        skm.setParallelLabels(ax, proj, parallels, loc="top")

        # add vertices as polygons
        vmin, vmax = np.percentile(bc,[10,90])
        poly = skm.addPolygons(vertices, proj, ax, color=bc, vmin=vmin, vmax=vmax, cmap=cmap, zorder=3, rasterized=True)

        # add colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.0)
        cb = fig.colorbar(poly, cax=cax)
        cb.set_label('$n_g$ [arcmin$^{-2}$]')
        cb.solids.set_edgecolor("face")

        return

    @staticmethod
    def footprint_plot(cat):

        mask1 = mask[np.in1d(mask,np.where((cat.ra>-70)&(cat.ra<100))[0])]
        fig = plt.figure(figsize=(6.5,6))
        y1_plots.footprint_sub(cat.ra[mask1], cat.dec[mask1],10,10,1024,fig)
        plt.savefig('plots/y1/footprint.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def psf_star_dist(cat):

        mask = cat.flag==0

        unique_ccd_ids = np.max(cat.ccd[mask])*cat.filenum[mask]+cat.ccd[mask]
        u,c=np.unique(unique_ccd_ids,return_counts=True)
        plt.hist(c,bins=100,range=(0,700),histtype='stepfilled',color='k')
        plt.ylabel('Number of CCDs')
        plt.xlabel('Number of stars per CCD')
        plt.savefig('plots/y1/psf_star_dist.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def bin_mean_new(x,y,edge):
        # I need to fix main destest version to work like this, much faster

        ind0 = 0
        mean = np.zeros(len(edge)-1)
        err = np.zeros(len(edge)-1)
        for j in range(len(edge)-1):
            ind = np.searchsorted(x, edge[j+1])
            if ind == 0:
                ind = -1
            mean[j] = np.mean(y[ind0:ind])
            err[j]  = np.std(y[ind0:ind])/np.sqrt(ind-ind0)
            ind0 = ind

        return mean,err


    @staticmethod
    def psf_mag_res_plot(cat,replace=False,bins=60):

        name = 'text/psf_resid.pkl'

        if replace|(not os.path.exists(name)):

            edge=lin.linear_methods.find_bin_edges(cat.mag,bins)
            i = np.argsort(cat.mag)

            arr1, tmp    = y1_plots.bin_mean_new(cat.mag[i],cat.mag[i],edge)
            dT,   dTerr  = y1_plots.bin_mean_new(cat.mag[i],(2.*cat.psf_size**2-2.*cat.size**2)[i],edge)
            de1,  de1err = y1_plots.bin_mean_new(cat.mag[i],(cat.psf1-cat.e1)[i],edge)
            de2,  de2err = y1_plots.bin_mean_new(cat.mag[i],(cat.psf2-cat.e2)[i],edge)

            d = {

            'arr1'   : arr1,
            'dT'     : dT,
            'de1'    : de1,
            'de2'    : de2,
            'dTerr'  : dTerr,
            'de1err' : de1err,
            'de2err' : de2err,
            }

            save_obj(d,name)
 
        else:

            d = load_obj(name)

        ax=plt.subplot(2,1,1)
        plt.errorbar(d['arr1'],d['dT'],yerr=d['dTerr'],marker='o',linestyle='',color='k')
        ax.minorticks_on()
        plt.ylabel(r'$T_{\mathrm{PSF}}-T_{\mathrm{model}}~(\mathrm{arcsec}^{2})$')
        plt.axvline(cat.mag[cat.flag==0].min(),color='k')
        plt.axhline(0.,color='k')
        ax.set_xticklabels([])

        ax=plt.subplot(2,1,2)
        plt.errorbar(d['arr1'],d['de1'],yerr=d['de1err'],marker='o',linestyle='',color='r',label=r'$e_1$')
        plt.errorbar(d['arr1'],d['de2'],yerr=d['de2err'],marker='o',linestyle='',color='b',label=r'$e_2$')
        ax.minorticks_on()
        plt.ylabel(r'$e_{\mathrm{PSF}}-e_{\mathrm{model}}$')
        plt.axvline(cat.mag[cat.flag==0].min(),color='k')
        plt.axhline(0.,color='k')
        plt.xlabel('Magnitude')

        plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})
        plt.tight_layout()
        plt.savefig('plots/y1/psf_res_mag.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def bin_median(x,y):
        # I need to fix main destest version to work like this, much faster

        ind0 = 0
        jmax = len(np.unique(x))
        median = np.zeros(jmax)
        d = np.where(np.diff(x)>0)[0]
        for j in range(jmax):
            if j==jmax-1:
                ind = -1
            else:
                ind = d[j]
            median[j] = np.median(y[ind0:ind])
            ind0 = ind

        return median

    @staticmethod
    def psf_star_fwhm_subplot(size,psf_exp,exp,band,color):

        flist = exp['exp'][exp['filter']==band]
        flist = np.char.strip(flist,'DECam_').astype(int)
        flist = np.unique(flist)
        fmask = np.in1d(psf_exp,flist,assume_unique=False)
        print psf_exp,size[fmask],np.sum(np.isnan(size[fmask]))
        median = y1_plots.bin_median(psf_exp[fmask],size[fmask])
        print median,np.max(median),np.min(median)
        print 'median seeing for '+band+' band',np.median(median)
        plt.hist(median,bins=50,histtype='stepfilled',color=color,alpha=0.2)

        return median

    @staticmethod
    def psf_star_fwhm_dist(cat,expfile):

        mask = cat.flag==0
        exp = fio.FITS(expfile)[-1].read(columns=['filter','exp'])

        psf_exp = np.char.strip(cat.filename[mask],'_psf')
        psf_exp = np.char.strip(psf_exp,'DECam_').astype(int)
        i = np.argsort(psf_exp)
        psf_exp = psf_exp[i]

        rmedian = y1_plots.psf_star_fwhm_subplot(np.sqrt(cat.size[mask][i].min()/2)*2.355,psf_exp,exp,'r','r')
        imedian = y1_plots.psf_star_fwhm_subplot(np.sqrt(cat.size[mask][i].min()/2)*2.355,psf_exp,exp,'i','b')
        zmedian = y1_plots.psf_star_fwhm_subplot(np.sqrt(cat.size[mask][i].min()/2)*2.355,psf_exp,exp,'z','k')
        print 'median seeing for riz bands',np.median(np.append(np.append(rmedian,imedian),zmedian))
        plt.ylabel('Number of exposures')
        plt.xlabel('Seeing FWHM (arcsec)')
        plt.savefig('plots/y1/psf_fwhm_dist.pdf', bbox_inches='tight')
        plt.close()

        return