import os
import numpy as np
import pickle
import glob
import healpy as hp
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
    print "Pickling:", name
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    print "Unpickling:", name
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
            mcal.add_pz(bpz,sheared=True,bounds=[0.2,1.3])

            #np.save('mcal_coadds.npy',np.vstack((mcal.coadd,np.ones(len(mcal.coadd)),np.ones(len(mcal.coadd)),np.zeros(len(mcal.coadd)),np.zeros(len(mcal.coadd)),mcal.w)).T)

            i3 = catalog.CatalogStore('im3shape',cutfunc=catalog.CatalogMethods.matched_i3_cut(),cattype='i3',catfile=i3file,goldfile=goldfile,goldcols=goldcols)

            # manually cut out S82
            catalog.CatalogMethods.match_cat(i3,i3.dec<-35)

            i3.bs    = True
            i3.wt    = True
            i3.add_pz(bpz,sheared=False,bounds=[0.2,1.3])

            #np.save('i3_coadds.npy',np.vstack((i3.coadd,i3.m1,i3.m2,i3.c1,i3.c2,i3.w)).T)

            bpz = None

            save_obj(mcal,mcalpickle)
            save_obj(i3,i3pickle)

        return i3,mcal

    @staticmethod
    def load_extra(goldfile,i3pickle,mcalpickle,fwhmfile,airmassfile,skybritefile,maglimfile):

        def load(cat,goldfile,fwhmfile,airmassfile,skybritefile,maglimfile):
            gold = fio.FITS(goldfile)[-1].read(columns=['coadd_objects_id','ebv','mof_flux_r','mof_flux_i'])
            s1,s2 = catalog.CatalogMethods.sort2(cat.coadd,gold['coadd_objects_id'])
            catalog.CatalogMethods.match_cat(cat,s1)
            cat.ebv = gold['ebv'][s2]
            cat.colour = gold['mof_flux_r'][s2]-gold['mof_flux_i'][s2]

            pix = hp.ang2pix(4096, np.pi/2.-np.radians(cat.dec),np.radians(cat.ra), nest=False)
            tmp=fio.FITS(fwhmfile)[-1].read()
            map = np.zeros(12*4096**2)
            map[tmp['PIXEL']] = tmp['SIGNAL']
            cat.fwhm = map[pix]

            tmp=fio.FITS(airmassfile)[-1].read()
            map = np.zeros(12*4096**2)
            map[tmp['PIXEL']] = tmp['SIGNAL']
            cat.airmass = map[pix]

            tmp=fio.FITS(skybritefile)[-1].read()
            map = np.zeros(12*4096**2)
            map[tmp['PIXEL']] = tmp['SIGNAL']
            cat.skybrite = map[pix]

            tmp=fio.FITS(maglimfile)[-1].read()
            map = np.zeros(12*4096**2)
            map[tmp['PIXEL']] = tmp['SIGNAL']
            cat.maglim = map[pix]

            return cat

        mcal = load_obj(mcalpickle)
        mcal = load(mcal,goldfile,fwhmfile,airmassfile,skybritefile,maglimfile)
        save_obj(mcal,mcalpickle+'_extra')
        mcal = None

        i3   = load_obj(i3pickle)
        i3   = load(i3,goldfile,fwhmfile,airmassfile,skybritefile,maglimfile)
        save_obj(i3,i3pickle+'_extra')
        i3   = None

        return

    @staticmethod
    def load_field_centers(special_points_file, rm_maskfile, rm_randomsfile, special_pickle, load_pickle=True):
        if os.path.exists(special_pickle) and load_pickle:
            special = load_obj(special_pickle)
        else:
            special=catalog.CatalogStore("special", catfile=special_points_file, cattype='gal', cols=-1,
                cutfunc=catalog.CatalogMethods.final_null_cuts_ra(), ranfile=rm_randomsfile)
            special.add_pixel(4096, nest=False)

            rm_mask_hpix = fio.FITS(rm_maskfile)[1]['HPIX'][:]
            inmask=np.in1d(special.pix, rm_mask_hpix)
            catalog.CatalogMethods.match_cat(special,inmask)
            save_obj(special, special_pickle)
        return special


    @staticmethod
    def load_psf_data(psfdir,psfpickle,load_pickle=True):

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

        if (load_pickle)&(os.path.exists(psfpickle)):

            #pass
            psf = load_obj(psfpickle)

        else:

            cols = ['e1','e2','ccd','col','row','psf1','psf2','mag','size','psf_size','flag']
            psf  = catalog.CatalogStore('psf',cutfunc=catalog.CatalogMethods.final_null_cuts_ra(),cols=cols,cattype='psf',catdir=psfdir)
            psf.ccd-=1
            catalog.CatalogMethods.match_cat(psf,psf.flag<2)

            #save_obj(psf,psfpickle)

        return psf

    @staticmethod
    def epoch_data_dump(epochdir,cols,hdu):

        lenst=0
        # Loop over file(s) [in directory]
        for ifile,file in enumerate(glob.glob(epochdir+'*')):
          # File format may not be readable
          try:
            fits=fio.FITS(file)
          except IOError:
            print 'error loading fits file: ',file

          # Dump the requested columns into memory if everything is there
          try:
            tmparray=fits[hdu].read(columns=cols)
          except ValueError:
            print 'error loading fits file: ',file
          fits.close()

          # If first file, generate storage array to speed up reading in of data
          if lenst==0:
            array=np.empty((500000000), dtype=tmparray.dtype.descr)
            
          array[lenst:lenst+len(tmparray)]=tmparray

          lenst+=len(tmparray)
          print ifile,len(tmparray),lenst,file

        return array

    @staticmethod
    def get_nonepoch(i3pickle,mcalpickle):

        mcal = load_obj(mcalpickle)
        mask = catalog.CatalogMethods.get_cuts_mask(mcal,full=False)
        e1,e2,w,m1,m2 = lin.linear_methods.get_lin_e_w_ms(mcal,mask=mask)
        mcal = np.vstack((np.copy(mcal.coadd[mask]),e1/m1,e2/m2,w)).T
        mask=np.argsort(mcal[:,0])
        mcal = mcal[mask]

        i3   = load_obj(i3pickle)
        e1,e2,w,m1,m2 = lin.linear_methods.get_lin_e_w_ms(i3)
        i3 = np.vstack((np.copy(i3.coadd),e1/m1,e2/m2,w)).T
        mask=np.argsort(i3[:,0])
        i3 = i3[mask]

        return mcal,i3

    @staticmethod
    def load_epoch_data(i3dir,mcaldir,imagefile,i3pickle,mcalpickle,i3epochpickle,mcalepochpickle,load_pickle=True):

        if (load_pickle)&(os.path.exists(mcalepochpickle))&(os.path.exists(i3epochpickle)):

            mcalepoch = load_obj(mcalepochpickle)
            i3epoch   = load_obj(i3epochpickle)

        else:

            mcal,i3 = y1.get_nonepoch(i3pickle,mcalpickle)

            mcalepoch0 = y1.epoch_data_dump(mcaldir,['id','orig_row','orig_col','image_id'],-2)
            mask = np.where(np.in1d(mcal[:,0],mcalepoch0['id'],assume_unique=False))[0]
            mcal = mcal[mask]
            mask = np.where(np.in1d(mcalepoch0['id'],mcal[:,0],assume_unique=False))[0]
            mcalepoch0 = mcalepoch0[mask]

            mcalepoch=catalog.CatalogStore('mcal',setup=False,cattype='epoch',release='y1')
            mcalepoch.coadd = mcalepoch0['id']
            mcalepoch.row = mcalepoch0['orig_row']
            mcalepoch.col = mcalepoch0['orig_col']
            mcalepoch.imid = mcalepoch0['image_id']
            mcalepoch0 = None

            mask=np.argsort(mcalepoch.coadd)
            catalog.CatalogMethods.match_cat(mcalepoch,mask)
            diff=np.diff(mcalepoch.coadd)
            diff=np.where(diff!=0)[0]+1
            diff=np.append([0],diff)
            diff=np.append(diff,[None])

            mcalepoch.e1=np.zeros(len(mcalepoch.coadd))
            mcalepoch.e2=np.zeros(len(mcalepoch.coadd))
            mcalepoch.w=np.zeros(len(mcalepoch.coadd))
            for i in xrange(len(diff)-1):
                mcalepoch.e1[diff[i]:diff[i+1]]=mcal[i,1]
                mcalepoch.e2[diff[i]:diff[i+1]]=mcal[i,2]
                mcalepoch.w[diff[i]:diff[i+1]]=mcal[i,3]

            im = fio.FITS(imagefile)[-1].read()
            im = im[np.argsort(im['id'])]
            mask=np.argsort(mcalepoch.imid)
            catalog.CatalogMethods.match_cat(mcalepoch,mask)
            im = im[np.in1d(im['id'],mcalepoch.imid,assume_unique=False)]
            catalog.CatalogMethods.match_cat(mcalepoch,np.in1d(mcalepoch.imid,im['id'],assume_unique=False))
            diff=np.diff(mcalepoch.imid)
            diff=np.where(diff!=0)[0]+1
            diff=np.append([0],diff)
            diff=np.append(diff,[None])

            mcalepoch.ccd=np.zeros(len(mcalepoch.coadd)).astype(int)
            for i in xrange(len(diff)-1):
                mcalepoch.ccd[diff[i]:diff[i+1]]=im['ccd'][i]-1

            save_obj(mcalepoch,mcalepochpickle)

            i3epoch0 = y1.epoch_data_dump(i3dir,['coadd_objects_id','ccd','orig_row','orig_col'],-1)
            mask = np.where(np.in1d(i3[:,0],i3epoch0['coadd_objects_id'],assume_unique=False))[0]
            i3 = i3[mask]
            mask = np.where(np.in1d(i3epoch0['coadd_objects_id'],i3[:,0],assume_unique=False))[0]
            i3epoch0 = i3epoch0[mask]

            i3epoch=catalog.CatalogStore('i3',setup=False,cattype='epoch',release='y1')
            i3epoch.coadd = i3epoch0['coadd_objects_id']
            i3epoch.ccd = i3epoch0['ccd']-1
            i3epoch.row = i3epoch0['orig_row']
            i3epoch.col = i3epoch0['orig_col']
            i3epoch0 = None

            mask=np.argsort(i3epoch.coadd)
            catalog.CatalogMethods.match_cat(i3epoch,mask)
            diff=np.diff(i3epoch.coadd)
            diff=np.where(diff!=0)[0]+1
            diff=np.append([0],diff)
            diff=np.append(diff,[None])

            i3epoch.e1=np.zeros(len(i3epoch.coadd))
            i3epoch.e2=np.zeros(len(i3epoch.coadd))
            i3epoch.w=np.zeros(len(i3epoch.coadd))
            for i in xrange(len(diff)-1):
                i3epoch.e1[diff[i]:diff[i+1]]=i3[i,1]
                i3epoch.e2[diff[i]:diff[i+1]]=i3[i,2]
                i3epoch.w[diff[i]:diff[i+1]]=i3[i,3]

            save_obj(i3epoch,i3epochpickle)

        return mcalepoch, i3epoch

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
    def mean_e(cat1,cat2,psfcat,replace=False):

        txt.write_methods.heading('Linear Splits',cat1,label='y1_paper',create=False)

        y1_plots.evspsf1(cat1,cat2,psfcat,replace=replace)
        y1_plots.evspsf2(cat1,cat2,psfcat,replace=replace)
        y1_plots.evssnr(cat1,cat2,replace=replace)
        y1_plots.evsradius(cat1,cat2,replace=replace)
        y1_plots.evspsfsize(cat1,cat2,replace=replace)

    # @staticmethod
    # def mean_e_epoch(cat1,cat2):

    #     txt.write_methods.heading('Linear Splits on epoch quantities',cat1,label='y1_paper',create=False)

    #     y1_plots.evsrow(cat1,cat2)

    @staticmethod
    def mean_e_subplot(cat,n,val,fig,replace=False,psf=None):

        name=fig0.plot_methods.get_filename_str(cat)
        name = 'text/lin_'+name+'_'+val+'.pkl'

        if replace|(not os.path.exists(name)):

            array=getattr(cat,val)
            if isinstance(cat,catalog.CatalogStore):
                mask=catalog.CatalogMethods.check_mask(cat.coadd,None)
            tmp,tmp,arr2,arr1err,e1,e2,e1err,e2err,m1,m2,b1,b2,m1err,m2err,b1err,b2err=sys_split.split_gals_lin_along_base(cat,val,array,mask,name,log=config.log_val.get(val,False),plot=False)

            if val in ['psf1','psf2']:
                edge=lin.linear_methods.find_bin_edges(getattr(psf,val),config.cfg.get('lbins',10))
                i=np.argsort(getattr(psf,val))
                catalog.CatalogMethods.match_cat(psf,i)
                if val == 'psf1':
                    mean,err=y1_plots.bin_mean_new(psf.psf1,psf.e1-psf.psf1,edge)
                    xmean,xerr=y1_plots.bin_mean_new(psf.psf1,psf.psf1,edge)
                else:
                    mean,err=y1_plots.bin_mean_new(psf.psf2,psf.e2-psf.psf2,edge)
                    xmean,xerr=y1_plots.bin_mean_new(psf.psf2,psf.psf2,edge)
                m,b,merr,berr = lin.fitting.lin_fit(xmean,mean,err)
            else:
                m=b=merr=berr=None

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
            'm1err' : m1err,
            'm2err' : m2err,
            'b1' : b1,
            'b2' : b2,
            'b1err' : b1err,
            'b2err' : b2err,
            'm' : m,
            'b' : b,
            'merr' : merr,
            'berr' : berr
            }

            save_obj(d,name)
 
        else:

            d = load_obj(name)


        xmin,xmax = y1_plots.x_lims[val]
        if val in ['size','psfsize']:
            d['arr1']=np.sqrt(d['arr1'])
            xmin=np.sqrt(xmin)
            xmax=np.sqrt(xmax)
        if config.log_val.get(val,False):
            xmin=np.log10(xmin)
            xmax=np.log10(xmax)
        y1min = np.minimum(np.minimum((d['m1']+d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']-d['b1err'],(d['m1']-d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']-d['b1err']),np.minimum((d['m1']-d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']+d['b1err'],(d['m1']+d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']+d['b1err']))
        y1max = np.maximum(np.maximum((d['m1']+d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']-d['b1err'],(d['m1']-d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']-d['b1err']),np.maximum((d['m1']-d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']+d['b1err'],(d['m1']+d['m1err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1']+d['b1err']))
        y2min = np.minimum(np.minimum((d['m2']+d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']-d['b2err'],(d['m2']-d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']-d['b2err']),np.minimum((d['m2']-d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']+d['b2err'],(d['m2']+d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']+d['b2err']))
        y2max = np.maximum(np.maximum((d['m2']+d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']-d['b2err'],(d['m2']-d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']-d['b2err']),np.maximum((d['m2']-d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']+d['b2err'],(d['m2']+d['m2err'])*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2']+d['b2err']))

        plt.figure(fig)
        ax=plt.subplot(2,1,n)
        plt.errorbar(d['arr1'],d['e1'],yerr=d['e1err'],marker='.',linestyle='',color='r',label=r'$\langle e_1 \rangle$')
        plt.errorbar(np.concatenate(([xmin],d['arr1'],[xmax])),d['m1']*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1'],marker='',linestyle='-',color='r')
        plt.fill_between(np.concatenate(([xmin],d['arr2'],[xmax])),y1min,y1max,interpolate=True,color='r',alpha=0.2)
        if val == 'psf1':
            plt.errorbar(np.concatenate(([xmin],d['arr1'],[xmax])),d['m']*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b1'],marker='',linestyle='--',color='r',lw=2)
        plt.errorbar(d['arr1'],d['e2'],yerr=d['e2err'],marker='.',linestyle='',color='b',label=r'$\langle e_2 \rangle$')
        plt.errorbar(np.concatenate(([xmin],d['arr1'],[xmax])),d['m2']*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2'],marker='',linestyle='-',color='b')
        plt.fill_between(np.concatenate(([xmin],d['arr2'],[xmax])),y2min,y2max,interpolate=True,color='b',alpha=0.2)
        if val == 'psf2':
            plt.errorbar(np.concatenate(([xmin],d['arr1'],[xmax])),d['m']*np.concatenate(([xmin],d['arr2'],[xmax]))+d['b2'],marker='',linestyle='--',color='b',lw=2)
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
    def evspsf1(cat1,cat2,psfcat,replace=False):

        plt.figure(3)

        y1_plots.mean_e_subplot(cat1,1,'psf1',3,replace=replace,psf=psfcat)
        y1_plots.mean_e_subplot(cat2,2,'psf1',3,replace=replace,psf=psfcat)

        plt.tight_layout()
        plt.savefig('plots/y1/lin_split_psf1.pdf', bbox_inches='tight')
        plt.close(3)

        return

    @staticmethod
    def evspsf2(cat1,cat2,psfcat,replace=False):

        plt.figure(4)

        y1_plots.mean_e_subplot(cat1,1,'psf2',4,replace=replace,psf=psfcat)
        y1_plots.mean_e_subplot(cat2,2,'psf2',4,replace=replace,psf=psfcat)

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


    # def evsrow(i3epoch,mcepoch,replace=False):

    #     plt.figure(6)

    #     y1_plots.mean_e_subplot(i3epoch,1,'row',6,replace=replace)
    #     y1_plots.mean_e_subplot(mcepoch,2,'row',6,replace=replace)

    #     plt.tight_layout()
    #     plt.savefig('plots/y1/lin_split_ccdrow.pdf', bbox_inches='tight')
    #     plt.close(6)


    # @staticmethod
    # def psf_whisker(psf,replace=False):

    #     psf.dpsf1 = psf.psf1-psf.e1
    #     psf.dpsf2 = psf.psf2-psf.e2

    #     y1_plots.whiskerplot(psf,'psf',6,replace=replace)
    #     y1_plots.whiskerplot(psf,'dpsf',7,replace=replace)

    #     return

    # @staticmethod
    # def e_whisker(cat1dir,cat1,cat2dir,cat2,replace=False):

    #     cols=['coadd','expnum','ccd','row','col']
    #     epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
    #     y1_plots.whiskerplot(epoch1,'e',8)

    #     cols=['coadd','expnum','ccd','row','col']
    #     epoch1=catalog.CatalogStore('epoch1',cutfunc=None,cattype='i3',cols=cols,catdir=cat1dir,release='y1')
    #     y1_plots.whiskerplot(epoch2,'e',9)

    #     return

    # @staticmethod
    # def whiskerplot(cat,col,fig,replace=False):

    #     if col == 'psf':
    #         key = r'e_{PSF}'
    #     if col == 'e':
    #         key = r'e'
    #     if col == 'dpsf':
    #         key = r'\Delta e_{PSF}'

    #     scale=0.02

    #     y,x,mw,e1,e2,e=field.field.whisker_calc(cat,col=col)
    #     pos0=0.5*np.arctan2(e2/mw,e1/mw)
    #     e/=mw
    #     for i in range(len(x)):
    #         y[i,:,:],x[i,:,:]=field.field_methods.ccd_to_field(i,y[i,:,:]-2048,x[i,:,:]-1024)

    #     print 'y,x',y[i,:,:],x[i,:,:]

    #     plt.figure(fig)
    #     print np.shape(x),np.shape(y),np.shape(np.sin(pos0)*e),np.shape(np.cos(pos0)*e)
    #     Q = plt.quiver(np.ravel(y),np.ravel(x),np.ravel(np.sin(pos0)*e),np.ravel(np.cos(pos0)*e),units='width',pivot='middle',headwidth=0,width=.0005)
    #     plt.quiverkey(Q,0.2,0.2,scale,str(scale)+' '+key,labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
    #     plt.savefig('plots/y1/whisker_'+col+'.pdf', dpi=500, bbox_inches='tight')
    #     plt.close(fig)

    #     return


    @staticmethod 
    def tangential_psf_data(psf, centers):
        name = "text/special_psf_gamma.pkl"
        if os.path.exists(name):
            return
        theta,out,err,chi2 = corr.xi_2pt.xi_2pt(centers, psf, corr='NG', 
                                                maska=None, maskb=None, ran=True)
    
        data = theta,out,err,chi2
        pickle.dump(data, open(name,"w"))
        

    @staticmethod 
    def tangential_shear_plot(i3, metacal, centers, ran, centers_mask=None, load=True):
        name = "text/special_gamma_ran_{}.pkl".format(ran)

        if os.path.exists(name) and load:
            data = pickle.load(open(name))
            i3_data, mc_data = data
            i3_theta,i3_out,i3_err,i3_chi2 = i3_data
            mc_theta,mc_out,mc_err,mc_chi2 = mc_data
        else:

            mask = None

            i3_theta,i3_out,i3_err,i3_chi2 = corr.xi_2pt.xi_2pt(centers, i3, corr='NG', 
                                                                maska=centers_mask, maskb=mask, ran=ran)
        
            mc_theta,mc_out,mc_err,mc_chi2 = corr.xi_2pt.xi_2pt(centers, metacal, corr='NG', 
                                                                maska=centers_mask, maskb=mask, ran=ran)

            i3_data = i3_theta,i3_out,i3_err,i3_chi2
            mc_data = mc_theta,mc_out,mc_err,mc_chi2
            data = (i3_data,mc_data)
            pickle.dump(data, open(name,"w"))

        i3_gammat = i3_out[0]
        i3_gammax = i3_out[2]
        i3_gammat_err = i3_err[0]
        i3_gammax_err = i3_err[2]
        mc_gammat = mc_out[0]
        mc_gammax = mc_out[2]
        mc_gammat_err = mc_err[0]
        mc_gammax_err = mc_err[2]

        chi2_i3_t = (i3_gammat**2/i3_gammat_err**2).sum()
        chi2_i3_x = (i3_gammax**2/i3_gammax_err**2).sum()
        chi2_mc_t = (mc_gammat**2/mc_gammat_err**2).sum()
        chi2_mc_x = (mc_gammax**2/mc_gammax_err**2).sum()

        print i3_gammat/i3_gammat_err
        print i3_gammax/i3_gammax_err
        print mc_gammat/mc_gammat_err
        print mc_gammax/mc_gammax_err

        print i3_gammat

        print "chi2 i3_t = ", chi2_i3_t
        print "chi2 i3_x = ", chi2_i3_x
        print "chi2 mc_t = ", chi2_mc_t
        print "chi2 mc_x = ", chi2_mc_x

        plt.figure()
        plt.errorbar(i3_theta, i3_gammat, i3_gammat_err, fmt='r.', label='Im3shape')
        plt.errorbar(mc_theta, mc_gammat, mc_gammat_err, fmt='b.', label='Metacal')
        plt.xscale('log')
        plt.axhline()
        plt.xlabel(r"$\theta$ / arcmin") 
        plt.ylabel(r"$\gamma_t$") 
        plt.legend()
        plt.savefig('plots/y1/special_gammat_ran_{}.pdf'.format(ran), dpi=500, bbox_inches='tight')
        plt.close()
        print "We think the imaginary gammat is gammax but not sure!"
        plt.figure()
        plt.errorbar(i3_theta, i3_gammax, i3_gammax_err, fmt='r.', label='Im3shape')
        plt.errorbar(mc_theta, mc_gammax, mc_gammax_err, fmt='b.', label='Metacal')
        plt.xscale('log')
        plt.axhline()
        plt.xlabel(r"$\theta$ / arcmin") 
        plt.ylabel(r"$\gamma_t$") 
        plt.legend()
        plt.savefig('plots/y1/special_gammax_ran_{}.pdf'.format(ran), dpi=500, bbox_inches='tight')
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

        skm.addFootprint('DES', proj, ax, zorder=10, edgecolor='#2222B2', facecolor='None', lw=2)
        
        def add_label(ra_label,dec_label,label):
            x,y=proj(ra_label,dec_label)
            ax.text(x,y,label,fontsize=14)
        
        add_label(45,-64,"SPT")
        add_label(350,4.0, "Stripe 82")
        add_label(50,-26, "D04")
        add_label(38.5,-1.5, "D04")

        #Fiddle with axis limit so we can see the whole range
        xmin,_=proj(100,-30)
        _,xmax=ax.get_xlim()
        ax.set_xlim(xmin,xmax)
        ymin,ymax=ax.get_ylim()
        r=ymax-ymin
        ymin-=r/10.
        ymax+=r/5.
        ax.set_ylim(ymin,ymax)
        return

    @staticmethod
    def footprint_plot(cat):

        mask = catalog.CatalogMethods.get_cuts_mask(cat,full=False)
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
        plt.hist(c,bins=100,range=(0,700),histtype='stepfilled',color='k', alpha=0.5)
        plt.ylabel('Number of CCD exposures')
        plt.xlabel('Stars per CCD exposure')
        plt.savefig('plots/y1/psf_star_dist.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def bin_mean_new(x,y,edge):
        # I need to fix main destest version to work like this, much faster

        ind0 = np.searchsorted(x, edge[0])
        mean = np.zeros(len(edge)-1)
        err = np.zeros(len(edge)-1)
        for j in range(len(edge)-1):
            if j%100==0:
                print j
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
            dT,   dTerr  = y1_plots.bin_mean_new(cat.mag[i],-(2.*cat.psf_size**2-2.*cat.size**2)[i],edge)
            de1,  de1err = y1_plots.bin_mean_new(cat.mag[i],-(cat.psf1-cat.e1)[i],edge)
            de2,  de2err = y1_plots.bin_mean_new(cat.mag[i],-(cat.psf2-cat.e2)[i],edge)

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

        print "Adding magnitude offset 5.3"
        mag_offset = 5.3

        ax=plt.subplot(2,1,1)
        plt.errorbar(mag_offset+d['arr1'],d['dT'],yerr=d['dTerr'],marker='.',linestyle='',color='k')
        ax.minorticks_on()
        plt.ylabel(r'$T_{\mathrm{PSF}}-T_{\mathrm{model}}~(\mathrm{arcsec}^{2})$')
        plt.axvline(mag_offset+cat.mag[cat.flag==0].min(),color='k')
        plt.axhline(0.,color='k')
        plt.axvspan(mag_offset+9, mag_offset+cat.mag[cat.flag==0].min(), facecolor='k', alpha=0.2)
        plt.xlim(15,23)

        ax.set_xticklabels([])

        ax=plt.subplot(2,1,2)
        plt.errorbar(mag_offset+d['arr1'],d['de1'],yerr=d['de1err'],marker='.',linestyle='',color='r',label=r'$e_1$')
        plt.errorbar(mag_offset+d['arr1'],d['de2'],yerr=d['de2err'],marker='.',linestyle='',color='b',label=r'$e_2$')
        ax.minorticks_on()
        plt.ylabel(r'$e_{\mathrm{PSF}}-e_{\mathrm{model}}$')
        plt.axvline(mag_offset+cat.mag[cat.flag==0].min(),color='k')
        plt.axhline(0.,color='k')
        plt.axvspan(mag_offset+9, mag_offset+cat.mag[cat.flag==0].min(), facecolor='k', alpha=0.2)
        plt.xlabel('Magnitude')
        plt.xlim(15,23)

        plt.legend(loc='upper right',ncol=1, frameon=True,prop={'size':12})
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
        plt.hist(median,bins=50,histtype='stepfilled',color=color,alpha=0.5,label=band)

        return median

    @staticmethod
    def psf_star_fwhm_dist(cat,expfile):

        mask = cat.flag==0
        exp = fio.FITS(expfile)[-1].read(columns=['filter','exp'])

        psf_exp = np.char.strip(cat.filename[mask],'_psf')
        psf_exp = np.char.strip(psf_exp,'DECam_').astype(int)
        i = np.argsort(psf_exp)
        psf_exp = psf_exp[i]

        ax = plt.subplot()
        fwhm = cat.size[mask][i]*2.355
        zmedian = y1_plots.psf_star_fwhm_subplot(fwhm,psf_exp,exp,'z','k')
        imedian = y1_plots.psf_star_fwhm_subplot(fwhm,psf_exp,exp,'i','b')
        rmedian = y1_plots.psf_star_fwhm_subplot(fwhm,psf_exp,exp,'r','r')
        print 'median seeing for riz bands',np.median(np.append(np.append(rmedian,imedian),zmedian))
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')
        plt.ylabel('Number of exposures')
        plt.xlabel('Seeing FWHM (arcsec)')
        plt.savefig('plots/y1/psf_fwhm_dist.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def psf_e_fov(cat,replace=False):
        import scipy.stats as stats

        name = 'text/psf_focal.pkl'

        if replace|(not os.path.exists(name)):

            mask = cat.flag==0
            x,y=field.field_methods.get_field_pos(cat)
            x=x[mask]
            y=y[mask]

            psf1,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,cat.psf1[mask],bins=1000)
            psf2,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,cat.psf2[mask],bins=1000)
            dpsf1,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,cat.psf1[mask]-cat.e1[mask],bins=1000)
            dpsf2,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,cat.psf2[mask]-cat.e2[mask],bins=1000)

            d = {
            'psf1'  : psf1,
            'psf2'  : psf2,
            'dpsf1' : dpsf1,
            'dpsf2' : dpsf2
            }

            save_obj(d,name)
 
        else:

            d = load_obj(name)

        plt.figure(figsize=(20,30))
        fig, ax = plt.subplots(nrows=2, ncols=2)
        im = ax[0,0].imshow(d['psf1'].T,vmin=-0.06, vmax=0.06, cmap = plt.get_cmap('PuOr'))
        im = ax[1,0].imshow(d['psf2'].T,vmin=-0.06, vmax=0.06, cmap = plt.get_cmap('PuOr'))
        im = ax[0,1].imshow(d['dpsf1'].T*10,vmin=-0.06, vmax=0.06, cmap = plt.get_cmap('PuOr'))
        im = ax[1,1].imshow(d['dpsf2'].T*10,vmin=-0.06, vmax=0.06, cmap = plt.get_cmap('PuOr'))
        fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.5)
        for ax_ in ax.flat:
            ax_.axis('off')
            ax_.set_aspect(1)
        ax[0,0].text(-100,250, r'$e_{1}^{\mathrm{PSF}}$', verticalalignment='center', horizontalalignment='center', rotation='vertical', fontdict={"size":14})
        ax[0,0].text(250,575, 'Mean', verticalalignment='center', horizontalalignment='center', fontdict={"size":12})
        ax[1,0].text(-100,250, r'$e_{2}^{\mathrm{PSF}}$', verticalalignment='center', horizontalalignment='center', rotation='vertical', fontdict={"size":14})
        ax[0,1].text(250,575, 'Mean residual ()', verticalalignment='center', horizontalalignment='center', fontdict={"size":12})
        plt.subplots_adjust(wspace=0.1, hspace=0.1, right=0.7, bottom=0.1, top=0.99)
        plt.savefig('plots/y1/psf_focal.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def mean_e_row_plot(cat,replace=False,bins=500):

        name = 'text/mean_e_row.pkl'

        if replace|(not os.path.exists(name)):

            edge=lin.linear_methods.find_bin_edges(cat.col,bins)
            i = np.argsort(cat.col)

            arr1, tmp  = y1_plots.bin_mean_new(cat.col[i],cat.col[i],edge)
            e1,  e1err = y1_plots.bin_mean_new(cat.col[i],cat.e1[i],edge)
            e2,  e2err = y1_plots.bin_mean_new(cat.col[i],cat.e2[i],edge)

            d = {

            'arr1'  : arr1,
            'e1'    : e1,
            'e2'    : e2,
            'e1err' : e1err,
            'e2err' : e2err,
            }

            save_obj(d,name)
 
        else:

            d = load_obj(name)

        # plt.errorbar(d['arr1'],d['e1'],yerr=d['e1err'],marker='.',linestyle='',color='r',label=r'$e_1$')
        # plt.errorbar(d['arr1'],d['e2'],yerr=d['e2err'],marker='.',linestyle='',color='b',label=r'$e_2$')
        plt.errorbar(d['arr1'],d['e1'],marker='.',linestyle='',color='r',label=r'$e_1$')
        plt.errorbar(d['arr1'],d['e2'],marker='.',linestyle='',color='b',label=r'$e_2$')
        plt.minorticks_on()
        plt.ylabel(r'$\langle e \rangle$')
        plt.axhline(0.,color='k')
        plt.fill_between([0.,30.],-0.004*np.ones(2),0.004*np.ones(2),interpolate=True,color='k',alpha=0.2)
        plt.fill_between([2017.,2048.],-0.004*np.ones(2),0.004*np.ones(2),interpolate=True,color='k',alpha=0.2)
        plt.xlim((0.,2048.))
        plt.ylim((-0.0015,0.0015))
        plt.xlabel('CCD pixel column')

        plt.legend(loc='lower right',ncol=1, frameon=True,prop={'size':12})
        plt.savefig('plots/y1/mean_e_row.pdf', bbox_inches='tight')
        plt.close()

        return


    @staticmethod
    def mean_e_fov(mcal,i3,replace=False):
        import scipy.stats as stats

        name = 'text/mean_e_focal.pkl'

        if replace|(not os.path.exists(name)):

            x,y=field.field_methods.get_field_pos(mcal)
            cate1,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,mcal.e1,bins=500)
            cate2,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,mcal.e2,bins=500)

            x,y=field.field_methods.get_field_pos(i3)
            cat2e1,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,i3.e1,bins=500)
            cat2e2,tmp,tmp,tmp = stats.binned_statistic_2d(x,y,i3.e2,bins=500)

            d = {
            'cate1'  : cate1,
            'cate2'  : cate2,
            'cat2e1'  : cat2e1,
            'cat2e2'  : cat2e2
            }

            save_obj(d,name)
 
        else:

            d = load_obj(name)

        plt.figure(figsize=(20,30))
        fig, ax = plt.subplots(nrows=2, ncols=2)
        im = ax[0,0].imshow(d['cate1'].T,vmin=-0.02, vmax=0.02, cmap = plt.get_cmap('PuOr'))
        im = ax[1,0].imshow(d['cate2'].T,vmin=-0.02, vmax=0.02, cmap = plt.get_cmap('PuOr'))
        im = ax[0,1].imshow(d['cat2e1'].T,vmin=-0.02, vmax=0.02, cmap = plt.get_cmap('PuOr'))
        im = ax[1,1].imshow(d['cat2e2'].T,vmin=-0.02, vmax=0.02, cmap = plt.get_cmap('PuOr'))
        fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.5)
        for ax_ in ax.flat:
            ax_.axis('off')
            ax_.set_aspect(1)
        ax[0,0].text(-100,250, r'$e_{1}$', verticalalignment='center', horizontalalignment='center', rotation='vertical', fontdict={"size":14})
        ax[0,0].text(250,575, r'$\textsc{metacal}$', verticalalignment='center', horizontalalignment='center', fontdict={"size":12})
        ax[1,0].text(-100,250, r'$e_{2}$', verticalalignment='center', horizontalalignment='center', rotation='vertical', fontdict={"size":14})
        ax[0,1].text(250,575, r'$\textsc{im3shape}$', verticalalignment='center', horizontalalignment='center', fontdict={"size":12})
        plt.subplots_adjust(wspace=0.1, hspace=0.1, right=0.7, bottom=0.1, top=0.99)
        plt.savefig('plots/y1/mean_e_focal.pdf', bbox_inches='tight')
        plt.close()

        return

    @staticmethod
    def mean_e_mock():

        g = np.zeros((100,8,4,2))
        for seed in range(100):
          for rlsn in range(8):
            for zbin in range(4):
                f = '/global/cscratch1/sd/seccolf/y1_patch/seed'+str(seed+1)+'/kgg-s'+str(seed+1)+'-f2z'+str(zbin+1)+'_c'+str(rlsn+1)+'.fits'
                print f
                fmap = fio.FITS(f)[-1].read(columns=['Q_STOKES','U_STOKES'])
                g[seed,rlsn,zbin,0]=np.mean(fmap['Q_STOKES'])
                g[seed,rlsn,zbin,1]=np.mean(fmap['U_STOKES'])

        cnt = 0
        g0 = np.zeros((800,4,2))
        for seed in range(100):
          for rlsn in range(8):
            for zbin in range(4):
                g0[cnt,zbin,0]=g[seed,rlsn,zbin,0]
                g0[cnt,zbin,1]=g[seed,rlsn,zbin,1]
            cnt+=1
        cov1=np.sum((g0-np.mean(g0,axis=0))*(g0-np.mean(g0,axis=0)),axis=0)*(len(g0)-1.)/len(g0)*(len(g0)-1-1)/(len(g0)-1)
        np.save('text/mean_e_flask.npy',g0)

        return
