import numpy as np
import os
import fitsio as fio
import healpy as hp
import src.catalog as catalog
import src.config as config
import src.fig as fig
import src.lin as lin
import src.sys_split as sys_split
import src.corr as corr
import src.field as field
import src.pz as pz
import src.cosmo as cosmo
import src.y1shearcat as y1
import src.txt as txt

#400 seconds to read mcal
#250 seconds to read i3



goldfile  = '/global/cscratch1/sd/tvarga/WLCAT/release/Y1A1_GOLD_1_0_3_wide_match_3.fits'
i3file    = '/global/cscratch1/sd/tvarga/WLCAT/release/Y1A1_GOLD_1_0_3_im3shape_3_psfex_2_nbc_2_match_2.fits'
mcalfile  = '/global/cscratch1/sd/tvarga/WLCAT/release/Y1A1_GOLD_1_0_3_metacalibration_2_psfex_2_match_2.fits'
i3epochdir= '/project/projectdirs/des/wl/desdata/wlpipe/im3shape_y1a1_v3/bord/epoch/'
mcalepoch = '/global/cscratch1/sd/tvarga/WLCAT/release/Y1A1_GOLD_1_0_3_metacalibration_2_psfex_2_match_2.fits'
psfdir    = '/global/cscratch1/sd/troxel/psf_cats/'

i3,mcal  = y1.y1.load_data(i3file,mcalfile,goldfile)
mcal.m1=np.ones(len(mcal.coadd))
mcal.m2=np.ones(len(mcal.coadd))


i3epoch = catalog.CatalogStore("i3epoch", 
    cattype="i3epoch", catdir=i3epochdir, 
    maxiter=10, cols=["coadd", "row", "col",  "e1", "e2"])


def add_i3_cal_to_epoch(epoch, i3):
    x=np.in1d(epoch.coadd,i3.coadd,assume_unique=False)
    catalog.CatalogMethods.match_cat(epoch,x)
    mask=np.argsort(epoch.coadd)
    catalog.CatalogMethods.match_cat(epoch,mask)

    diff=np.diff(epoch.coadd)
    diff=np.where(diff!=0)[0]+1
    # diff = np.concatenate([[0], diff, [None]])
    diff=np.append([0],diff)
    diff=np.append(diff,[None])

    epoch.m1=np.zeros(len(epoch.coadd))
    epoch.m2=np.zeros(len(epoch.coadd))
    epoch.c1=np.zeros(len(epoch.coadd))
    epoch.c2=np.zeros(len(epoch.coadd))
    epoch.w=np.zeros(len(epoch.coadd))

    for i in xrange(len(diff)-1):
        epoch.m1[diff[i]:diff[i+1]]=i3.m1[i]
        epoch.m2[diff[i]:diff[i+1]]=i3.m2[i]
        epoch.c1[diff[i]:diff[i+1]]=i3.c1[i]
        epoch.c2[diff[i]:diff[i+1]]=i3.c2[i]
        epoch.w[diff[i]:diff[i+1]]=i3.w[i]
    

add_i3_cal_to_epoch(i3epoch, i3)



txt.write_methods.heading('---------------',mcal,label='y1_paper',create=True)

y1.y1_plots.mean_e(i3,mcal)
# psf  = y1.y1.load_psf_data(psfdir)
# y1.y1_plots.psf_whisker(psf)
# y1.y1_plots.e_whisker(i3epoch,i3,mcalepoch,mcal)
print "Warning: no metacal epoch file yet!"
y1.y1_plots.mean_e_epoch(i3epoch, i3epoch)

