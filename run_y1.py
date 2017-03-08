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

goldfile  = '/global/cscratch1/sd/troxel/y1a1-gold-mof-badregion.fits'
i3file    = '/global/cscratch1/sd/troxel/y1a1-im3shape_v4-nbc_v2_matched_v3.fits'
mcalfile  = '/global/cscratch1/sd/troxel/mcal-y1a1-combined-riz-blind-v2-matched.fits'
bpzfile   = '/global/cscratch1/sd/troxel/BPZ_ngmix_mof_slr_HiZ_combined_matched.fits'
i3epochdir= '/project/projectdirs/des/wl/desdata/wlpipe/im3shape_y1a1_v3/bord/epoch/'
mcalepoch = '/global/cscratch1/sd/tvarga/WLCAT/release/Y1A1_GOLD_1_0_3_metacalibration_2_psfex_2_match_2.fits'
rmfile    = '/global/cscratch1/sd/troxel/redmagicv6.4.11/y1a1_gold_1.0.2c-full_redmapper_v6.4.11_redmagic_combined_troxel.fit'
psfdir    = '/global/cscratch1/sd/troxel/psf_cats/'
special_points_file = '/global/cscratch1/sd/zuntz/y1a1_special_field_points.fits'

i3,mcal  = y1.y1.load_data(i3file,mcalfile,goldfile)
# rm = catalog.CatalogStore('rm',cutfunc=catalog.CatalogMethods.final_null_cuts_ra(),cattype='gal',catfile=rmfile,cols=['coadd','ra','dec'])
# special=catalog.CatalogStore("special", catfile=special_points_file, cattype='gal', cols=-1, cutfunc=catalog.CatalogMethods.final_null_cuts_ra())


# i3epoch = catalog.CatalogStore("i3epoch", 
#     cattype="i3epoch", catdir=i3epochdir, 
#     maxiter=10, cols=["coadd", "row", "col",  "e1", "e2"])


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
    

def make_randoms(name, cat, rans_per_object):
    if os.path.exists(name+"_random.fits.gz"):
        print "Using existing random", name+"_random.fits.gz"
    else:
        print "Making new random (one off)"
        catalog.CatalogMethods.create_random_cat_from_cat(cat, nran=len(cat.e1)*rans_per_object, label=name+"_")
    f = fitsio.FITS(name+"_random.fits.gz")
    data = f[1].read_columns(['ra', 'dec'])
    cat.ran_ra = data['ra']
    cat.ran_dec = data['dec']
    f.close()



# add_i3_cal_to_epoch(i3epoch, i3)

# add_randoms(i3, 1)
# add_randoms(metacal, 1)


txt.write_methods.heading('---------------',mcal,label='y1_paper',create=True)

y1.y1_plots.mean_e(i3,mcal)

# psf  = y1.y1.load_psf_data(psfdir)
# y1.y1_plots.psf_whisker(psf)
# y1.y1_plots.e_whisker(i3epoch,i3,mcalepoch,mcal)
# print "Warning: no metacal epoch file yet!"
# y1.y1_plots.mean_e_epoch(i3epoch, i3epoch)

